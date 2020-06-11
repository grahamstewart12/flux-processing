### ============================================================================
# Biomet data gapfilling =======================================================
### ============================================================================

# Purpose: 

# References:

# Franz, D., Koebsch, F., Larmanou, E., Augustin, J., & Sachs, T. (2016). High 
# net CO2 and CH4 release at a eutrophic shallow lake on a formerly drained fen. 
# Biogeosciences, 13(10), 3051–3070. https://doi.org/10.5194/bg-13-3051-2016

# Mauder, M., Cuntz, M., Drüe, C., Graf, A., Rebmann, C., Schmid, H. P., et al. 
# (2013). A strategy for quality and uncertainty assessment of long-term 
# eddy-covariance measurements. Agricultural and Forest Meteorology, 169, 
# 122–135. https://doi.org/10.1016/j.agrformet.2012.09.006

# Vuichard, N., & Papale, D. (2015). Filling the gaps in meteorological 
# continuous data measured at FLUXNET sites with ERA-Interim reanalysis. Earth 
# System Science Data, 7(2), 157–171. https://doi.org/10.5194/essd-7-157-2015


# Input(s):

# Output(s):

start_time <- Sys.time()

# Load the required packages
#suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
library(REddyProc)
library(lubridate)
library(tidyverse)
# Packages needed but not loaded: solartime, bigleaf, broom

# Load reference files
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")
source("~/Desktop/DATA/Flux/tools/reference/var_attributes.R")
source("~/Desktop/DATA/Flux/tools/reference/gf_control.R")

# Load functions
source("~/Desktop/DATA/Flux/tools/engine/functions/clean.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/dates_and_times.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/flag.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/latest_version.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/utilities.R")


### Helper functions ===========================================================

get_best_lag <- function(y, x, lag.max = 12) {
  
  dy <- y - dplyr::lag(y)
  dx <- x - dplyr::lag(x)
  
  ccf <- broom::tidy(
    ccf(dx, dy, lag.max = lag.max, na.action = na.pass, plot = FALSE)
  )
  
  ccf %>%
    dplyr::arrange(dplyr::desc(acf)) %>%
    purrr::pluck("lag", 1) %>%
    magrittr::multiply_by(-1)
}

apply_lag <- function(x, lag) {
  
  if (lag >= 0) {
    out <- dplyr::lag(x, lag)
  } else {
    out <- dplyr::lead(x, -lag)
  }
  
  out
}

select_clean <- function(data, vars) {
  
  vars <- purrr::map_chr(vars, rlang::as_string)
  df <- dplyr::select(data, timestamp, dplyr::all_of(vars))
  
  get_qc_var <- function(x) {
    rlang::as_string(x) %>%
      stringr::str_c("qc_", .) %>%
      rlang::sym()
  }
  
  for (i in 2:length(df)) {
    name <- rlang::sym(names(df)[i])
    qc_name <- get_qc_var(name)
    df <- dplyr::mutate(
      df, !!name := clean(!!name, dplyr::pull(data, !!qc_name))
    )
  }
  
  df
}

roll_mean_real <- function(x, n = 1L) {
  
  m <- RcppRoll::roll_mean(x, n, fill = NA, na.rm = TRUE)
  
  # Only real x values get a roll_mean value
  na_ind <- which(is.na(x))
  m[na_ind] <- NA
  
  # Fill with original values where roll_mean is cut off at the ends
  dplyr::coalesce(m, x)
}

plot_harmonies <- function(data, aux_data, aux_suffix) {
  
  sfx <- aux_suffix
  sfx_ <- stringr::str_c("_", sfx)
  
  aux_data %>%
    dplyr::rename_with(~ stringr::str_c(., sfx_), -timestamp) %>%
    dplyr::left_join(data, by = "timestamp") %>%
    tidyr::pivot_longer(-timestamp, names_to = "var", values_to = "value") %>%
    dplyr::mutate(
      site = dplyr::if_else(stringr::str_detect(var, sfx_), sfx, "site"),
      var = stringr::str_remove(var, sfx_)
    ) %>%
    tidyr::pivot_wider(names_from = "site", values_from = "value") %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(!!rlang::sym(sfx), site)) +
    ggplot2::facet_wrap(~ var, scales = "free") +
    ggplot2::geom_hex(
      ggplot2::aes(alpha = log(..count..)), fill = "steelblue", bins = 20
    ) +
    ggplot2::geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE, size = 0.75, linetype = 2, 
      color = "black"
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.4, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

debias_init <- function(data, aux_data, var, diff = FALSE, 
                        type = c("lm", "lm0", "ratio"), lag = 0, ctrl) {
  
  # Handle vars with different names between datasets 
  if (length(var) > 1) {
    var1 <- rlang::enquo(var[1])
    var2 <- rlang::enquo(var[2])
  } else {
    var1 <- rlang::enquo(var)
    var2 <- rlang::enquo(var)
  }
  
  # Get information from control list if provided
  if (!missing(ctrl)) {
    diff <- purrr::pluck(ctrl, var[1], "db_diff")
    type <- purrr::pluck(ctrl, var[1], "db_type")
    lag <- purrr::pluck(ctrl, var[1], "db_lag")
  } else {
    type <- rlang::arg_match(type)
  }
  
  tbl <- tibble::tibble(
    y = dplyr::pull(data, !!var1),
    x = dplyr::pull(aux_data, !!var2)
  )
  
  # Difference if requested (useful for highly autocorrelated variables)
  if (diff) tbl <- dplyr::mutate_all(tbl, ~ . - dplyr::lag(., 1))
  
  # Apply lag if indicated
  if (lag < 0) {
    tbl <- dplyr::mutate(tbl, x = dplyr::lead(x, -lag))
  } else {
    tbl <- dplyr::mutate(tbl, x = dplyr::lag(x, lag))
  } 
  
  # Only include existing pairs of values 
  tbl <- dplyr::mutate(
    tbl,
    y = dplyr::if_else(is.na(x), NA_real_, y),
    x = dplyr::if_else(is.na(y), NA_real_, x)
  )
  tbl <- tidyr::drop_na(tbl)
  
  # Fit debiasing model
  if (type == "lm") {
    fit <- lm(y ~ x, data = tbl)
    names(fit$coefficients) <- c("b", "m")
  } else if (type == "lm0") {
    fit <- lm(y ~ 0 + x, data = tbl)
    fit$coefficients <- c(b = 0, m = unname(coef(fit)))
  }  else if (type == "ratio") {
    # Manually calculate coefficient for ratios
    beta <- tbl %>%
      dplyr::filter(sign(y) == sign(x) | sign(y) == 0 | sign(x) == 0) %>%
      dplyr::summarize(
        beta = sum(abs(y)) / sum(abs(x)), .groups = "drop_last"
      ) %>%
      purrr::pluck(1, 1)
    # Assign coefficient as offset in lm
    fit <- lm(y ~ 0 + x, data = tbl)
    fit$coefficients <- c(b = 0, m = beta)
  }
  
  # Add lag attribute for fit on differenced variables
  if (diff) attr(fit, "lag") <- lag
  
  fit
}

debias <- function(data, aux_data, var, fit, diff = FALSE, lag = 0, ctrl) {
  
  # Function to help find QC variable
  get_qc_var <- function(x) {
    rlang::as_string(x) %>% stringr::str_c("qc_", .) %>% rlang::sym()
  }
  
  # Handle vars with different names between datasets 
  if (length(var) > 1) {
    var1 <- rlang::enquo(var[1])
    var2 <- rlang::enquo(var[2])
  } else {
    var1 <- rlang::enquo(var)
    var2 <- rlang::enquo(var)
  }
  
  # Get information from control list if provided
  if (!missing(ctrl)) {
    diff <- purrr::pluck(ctrl, var[1], "db_diff")
    lag <- purrr::pluck(ctrl, var[1], "db_lag")
  }
  
  tbl <- tibble::tibble(
    y = dplyr::pull(data, !!var1),
    x = dplyr::pull(aux_data, !!var2)
  )
  
  # Apply lag if indicated
  if (lag < 0) {
    tbl <- dplyr::mutate(tbl, x = dplyr::lead(x, -lag))
  } else {
    tbl <- dplyr::mutate(tbl, x = dplyr::lag(x, lag))
  } 
  
  # If differences, reconstruct data using time series fill
  if (diff) {
    out <- fill_along(
      tbl$y, tbl$x, coef(fit)[2], lag = 0, align = TRUE, prioritize = "fore"
    )
    #out <- dplyr::if_else(is.na(tbl$x - dplyr::lag(tbl$x, 1)), NA_real_, out)
  } else {
    out <- coef(fit)[1] + coef(fit)[2] * tbl$x
  }
  
  out
}

fill_linear <- function(x, n_max = NULL) {
  
  # Option for safely returning original vector
  if (!is.null(n_max)) {
    if (is.na(n_max) | n_max == 0) return(x)
  }
  
  y <- seq_along(x)
  out <- approx(y, x, y)$y
  
  if (!is.null(n_max)) {
    x_tmp <- tidyr::replace_na(x, 9999)
    rl <- rle(c(x_tmp))
    n <- rep(rl$lengths, times = rl$lengths)
    out <- dplyr::if_else(n > n_max & is.na(x), NA_real_, out)
  }
  
  out
}

fill_along <- function(y, x, coef, lag = 0, align = FALSE, prioritize) {
  
  # Essentially forecasting/hindcasting, using a secondary variable as reference
  
  prioritize <- rlang::arg_match(prioritize, c("fore", "hind", "both"))
  
  x <- dplyr::lag(x, lag)
  dx <- x - dplyr::lag(x)
  y_lag <- dplyr::lag(y, 1)
  y_f <- y
  
  # Set NA indices and initialize loop
  na_i <- which(is.na(y) & !is.na(y_lag))
  n_na <- length(which(is.na(y)))
  iter <- 0
  
  # Fill gaps moving forward in time (i.e. forecast)
  while (n_na > 0) {
    y_f[na_i] <- y_lag[na_i] + coef * dx[na_i]
    
    y_lag <- dplyr::lag(y_f, 1)
    
    na_i <- which(is.na(y_f) & !is.na(y_lag))
    if (length(which(is.na(y_f))) == n_na) break
    n_na <- length(which(is.na(y_f)))
  }
  
  dx <- dplyr::lead(x, 1) - x
  y_lead <- dplyr::lead(y, 1)
  y_b <- y
  
  # Reset NA indices and re-initialize loop
  na_i <- which(is.na(y) & !is.na(y_lag))
  n_na <- length(which(is.na(y)))
  iter <- 0
  
  # Fill gaps moving backward in time (i.e. hindcast)
  while (n_na > 0) {
    y_b[na_i] <- y_lead[na_i] - coef * dx[na_i]
    
    y_lead <- dplyr::lead(y_b, 1)
    
    na_i <- which(is.na(y_b) & !is.na(y_lead))
    if (length(which(is.na(y_b))) == n_na) break
    n_na <- length(which(is.na(y_b)))
  }
  
  # Interpolation is forecast, hindcast, or average of both
  if (prioritize == "fore") {
    out <- dplyr::coalesce(y_f, y_b)
  } else if (prioritize == "hind") {
    out <- dplyr::coalesce(y_b, y_f)
  } else if (prioritize == "both") {
    out <- purrr::pmap_dbl(list(y_f, y_b), purrr::lift_vd(mean, na.rm = TRUE))
  }
  
  if (align) {
    gaps <- dplyr::if_else(is.na(y), 1L, 0L)
    weights <- gaps %>%
      tidy_rle() %>%
      dplyr::mutate(yf = y_f, d = dplyr::lead(yf) - yf) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(
        n = seq_along(values),
        di = dplyr::last(d) - dplyr::nth(d, -2),
        pf = di * (n / lengths),
        pf = dplyr::if_else(values == 1, pf, 0),
        pf = dplyr::if_else(lengths == 1, 0, pf),
        pf = tidyr::replace_na(pf, 0)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::pull(pf)
    out <- out + weights
  }
  
  out
}

plot_fill_vars <- function(data, fill_name) {
  
  data %>%
    tidyr::pivot_longer(-timestamp, values_to = fill_name) %>%
    ggplot2::ggplot(ggplot2::aes(timestamp, !!rlang::sym(fill_name))) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::scale_x_datetime(date_labels = "%b") +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_bw()
}

blend_vars <- function(data, aux_data, x, var1, var2, var3 = NULL, 
                       diff = FALSE) {
  
  # In wetlands the "depth" of soil temperature changes with the water level
  # - could be more accurate to allow "depth" to vary over time 
  # - this function creates a blend two time series, where the weight given to 
  #   each varies daily based on a regression with the dependent variable
  
  x <- rlang::enquo(x)
  var1 <- rlang::enquo(var1)
  var2 <- rlang::enquo(var2)
  var3 <- rlang::enquo(var3)
  
  # THIS IS HORRIBLE but allowing n vars as [...] arg too complicated for now
  if (rlang::quo_is_null(var3)) {
    blend <- data %>%
      dplyr::mutate(date = lubridate::date(timestamp), y = !!x) %>% 
      dplyr::left_join(
        dplyr::select(aux_data, timestamp, x1 = !!var1, x2 = !!var2), 
        by = "timestamp"
      ) %>%
      dplyr::select(date, timestamp, y, x1, x2)
    
    if (diff) blend <- dplyr::mutate_if(blend, is.numeric, ~ . - dplyr::lag(.))
    
    blend <- blend %>%
      tidyr::drop_na() %>%
      tidyr::nest(data = c(timestamp, y, x1, x2)) %>%
      dplyr::mutate(
        lm = purrr::map(data, ~ lm(y ~ 0 + x1 + x2, data = .)), 
        coefs = purrr::map(lm, coef), 
        p1 = purrr::map_dbl(coefs, pluck, 1) %>% 
          RcppRoll::roll_mean(7, na.rm = TRUE, fill = NA), 
        p2 = purrr::map_dbl(coefs, pluck, 2) %>%
          RcppRoll::roll_mean(7, na.rm = TRUE, fill = NA), 
        tot = purrr::pmap_dbl(list(p1, p2), sum)
      ) %>% 
      dplyr::mutate(dplyr::across(c(p1, p2), ~ pmax(pmin(. / tot, 1), 0))) %>%
      tidyr::unnest(data) %>% 
      # Smooth daily boundaries
      dplyr::transmute(
        timestamp = timestamp,
        p1 = RcppRoll::roll_mean(p1, 21, na.rm = TRUE, fill = NA), 
        p2 = RcppRoll::roll_mean(p2, 21, na.rm = TRUE, fill = NA)
      )
    
    blend <- aux_data %>%
      dplyr::select(timestamp, x1 = !!var1, x2 = !!var2) %>%
      dplyr::left_join(blend, by = "timestamp") %>%
      dplyr::mutate(dplyr::across(
        c(p1, p2), ~ dplyr::if_else(is.nan(.), NA_real_, .)
      )) %>%
      tidyr::fill(p1, p2, .direction = "downup") %>%
      dplyr::mutate(b = x1 * p1 + x2 * p2)
  } else {
    
    blend <- data %>%
      dplyr::transmute(
        timestamp = timestamp, date = lubridate::date(timestamp), y = !!x
      ) %>% 
      dplyr::left_join(dplyr::select(
        aux_data, timestamp, x1 = !!var1, x2 = !!var2, x3 = !!var3
      ), by = "timestamp")
    
    if (diff) blend <- dplyr::mutate_if(blend, is.numeric, ~ . - dplyr::lag(.))
    
    blend <- blend %>%
      tidyr::drop_na() %>%
      tidyr::nest(data = c(timestamp, y, x1, x2, x3)) %>%
      dplyr::mutate(
        lm = purrr::map(data, ~ lm(y ~ 0 + x1 + x2 + x3, data = .)), 
        coefs = purrr::map(lm, coef), 
        p1 = purrr::map_dbl(coefs, pluck, 1) %>% 
          RcppRoll::roll_mean(7, na.rm = TRUE, fill = NA), 
        p2 = purrr::map_dbl(coefs, pluck, 2) %>%
          RcppRoll::roll_mean(7, na.rm = TRUE, fill = NA), 
        p3 = purrr::map_dbl(coefs, pluck, 3) %>%
          RcppRoll::roll_mean(7, na.rm = TRUE, fill = NA), 
        tot = purrr::pmap_dbl(list(p1, p2, p3), sum)
      ) %>% 
      dplyr::mutate(dplyr::across(
        c(p1, p2, p3), ~ pmax(pmin(. / tot, 1), 0)
      )) %>%
      tidyr::unnest(data) %>% 
      # Smooth daily boundaries
      dplyr::transmute(
        timestamp = timestamp,
        p1 = RcppRoll::roll_mean(p1, 21, na.rm = TRUE, fill = NA), 
        p2 = RcppRoll::roll_mean(p2, 21, na.rm = TRUE, fill = NA),
        p3 = RcppRoll::roll_mean(p3, 21, na.rm = TRUE, fill = NA)
      )
    
    blend <- aux_data %>%
      dplyr::select(timestamp, x1 = !!var1, x2 = !!var2, x3 = !!var3) %>%
      dplyr::left_join(blend, by = "timestamp") %>%
      dplyr::mutate(dplyr::across(
        c(p1, p2, p3), ~ dplyr::if_else(is.nan(.), NA_real_, .)
      )) %>%
      tidyr::fill(p1, p2, p3, .direction = "downup") %>%
      dplyr::mutate(b = x1 * p1 + x2 * p2 + x3 * p3)
  }
  
  dplyr::pull(blend, b)
}

flatten_period <- function(x, n, min_step = 1e-4, max_iter = 50) {
  range <- range(x, na.rm = TRUE)
  
  out <- x
  offset <- mean(out[1:n], na.rm = TRUE)
  step <- min_step + 1
  iter <- 0
  
  while (step > min_step & iter < max_iter) {
    offset <- mean(out[1:n], na.rm = TRUE)
    out <- 1 - abs(offset - out)
    out <- scales::rescale(out, range)
    step <- mean(abs(offset - out[1:n]), na.rm = TRUE)
    iter <- iter + 1
  }
  
  out
}

gather_fill_data <- function(data, var, backup, order, max_i, ctrl) {
  
  var_name <- rlang::as_string(var)
  var <- rlang::ensym(var)
  
  # Get information from control list if provided
  if (!missing(ctrl)) {
    backup <- purrr::pluck(ctrl, var_name, "gf_backup")
    order <- purrr::pluck(ctrl, var_name, "gf_order")
    max_i <- purrr::pluck(ctrl, var_name, "gf_max_i")
  }
  
  list <- list(x = dplyr::pull(data, !!var))
  
  # Gap filling method modules
  
  # "Backup" variable
  if (!is.na(backup)) {
    b_var <- rlang::ensym(backup)
    # Simple de-bias in case sensors are slightly off
    coef <- data %>% 
      dplyr::select(!!var, !!b_var) %>% 
      tidyr::drop_na() %>% 
      dplyr::summarize(sum(!!var) / sum(!!b_var), .groups = "drop_last") %>% 
      purrr::pluck(1)
    list <- append(list, list(b = dplyr::pull(data, !!b_var) * coef))
  }
  
  # Potential radiation
  if (stringr::str_detect(order, "p")) {
    pot <- dplyr::pull(data, sw_in_pot)
    
    # - purpose is mostly to correct negative nighttime values
    pot <- dplyr::if_else(pot == 0, pot, NA_real_)
    list <- append(list, list(p = pot))
  }
  
  # Interpolation
  if (stringr::str_detect(order, "i")) {
    i <- fill_linear(dplyr::pull(data, !!var), max_i)
    # Remove included original values
    i <- dplyr::if_else(!is.na(dplyr::pull(data, !!var)), NA_real_, i)
    list <- append(list, list(i = i))
  }
  
  # MDC
  if (stringr::str_detect(order, "m")) {
    m_var <- var_name %>% stringr::str_c(., "_f") %>% rlang::sym()
    m_qc_var <- var_name %>% stringr::str_c(., "_fqc") %>% rlang::sym()
    m <- clean(
      dplyr::pull(data, !!m_var), dplyr::pull(data, !!m_qc_var), 
      value = 3, na.as = NA
    )
    # Remove included original values
    m <- dplyr::if_else(!is.na(dplyr::pull(data, !!var)), NA_real_, m)
    list <- append(list, list(m = m))
  }
  
  # De-biased auxilliary data
  if (stringr::str_detect(order, "a")) {
    a_var <- var_name %>% stringr::str_c(., "_d") %>% rlang::sym()
    list <- append(list, list(a = dplyr::pull(data, !!a_var)))
  }
  
  # De-biased ERA data
  if (stringr::str_detect(order, "e")) {
    e_var <- var_name %>% stringr::str_c(., "_df") %>% rlang::sym()
    list <- append(list, list(e = dplyr::pull(data, !!e_var)))
  }
  
  list
}

plan_fill <- function(list, backup, order, ctrl, var) {
  #browser()
  # Get information from control list if provided
  if (!missing(ctrl)) {
    backup <- purrr::pluck(ctrl, var, "gf_backup")
    order <- purrr::pluck(ctrl, var, "gf_order")
  }
  
  # Add backup as first method if backup var is detected
  if (!is.na(backup)) order <- stringr::str_c("b", order)
  
  # Add "x" (original data) to front of list
  order <- stringr::str_c("x", order)
  
  names <- order %>% stringr::str_split("") %>% purrr::pluck(1)
  
  # Gap-focused approach (Plan A)
  gaps <- list %>%
    # Replace non-missing values with the gapfill variable name
    tibble::as_tibble() %>%
    # Make sure fill vars are in the correct order
    dplyr::select(dplyr::all_of(names)) %>%
    dplyr::mutate(
      gap = dplyr::if_else(is.na(x), 1L, 0L),
      gap_id = dplyr::pull(tidy_rle(gap), id),
      gap_len = dplyr::pull(tidy_rle(gap), lengths)
    ) %>%
    dplyr::group_by(gap_id) %>%
    dplyr::summarize(
      dplyr::across(-gap, ~ length(na.omit(.x))), .groups = "drop_last"
    ) %>%
    dplyr::ungroup()
  
  cover <- gaps %>%
    dplyr::mutate(
      # Do gap-fill vars cover at least part of a gap?
      any = purrr::imap_dfr(
        ., ~ dplyr::if_else(.x != 0, .y, NA_character_)
      ) %>% 
        dplyr::select(-dplyr::starts_with("gap")) %>% 
        purrr::pmap_chr(dplyr::coalesce), 
      # Do gap-fill vars cover entire gaps?
      all = purrr::imap_dfr(
        ., ~ dplyr::if_else(.x == gap_len, .y, NA_character_)
      ) %>% 
        dplyr::select(-dplyr::starts_with("gap")) %>% 
        purrr::pmap_chr(dplyr::coalesce)
    ) %>% 
    dplyr::mutate(
      plan = purrr::pmap_chr(dplyr::select(., all, any), dplyr::coalesce)
    )
  
  plan_a <- rep(cover$all, times = cover$gap_len)
  
  # Helper function for replacing non-NA values
  replace_real <- function(x, replace) {
    #x <- vctrs::vec_cast(x, class(replace))
    dplyr::if_else(!is.na(x), replace, as.character(x))
  }
  
  # Point-focused approach (Plan B)
  # Replace non-missing values with the gapfill variable name
  tbl <- list %>% 
    purrr::imap(~ replace_real(.x, .y)) %>%
    tibble::as_tibble()
  
  # Make sure fill vars are in the correct order
  inorder <- dplyr::select(tbl, dplyr::all_of(names))
  
  # Coalesce to form "Plan B"
  plan_b <- dplyr::coalesce(!!!as.list(inorder))
  
  plan <- dplyr::coalesce(plan_a, plan_b)
  
  plan
}

gapfill_biomet <- function(fill_data, fmeth) {
  
  # This method is MUCH slower than coalescing vectors
  
  out <- rep(NA_real_, length(fmeth))
  for (i in seq_along(fmeth)) {
    if (is.na(fmeth[i])) next
    out[i] <- purrr::pluck(fill_data, fmeth[i], i)
  }
  
  out
}

qc_biomet_fmeth <- function(fmeth, mdc_fqc) {
  
  if (missing(mdc_fqc)) mdc_fqc <- rep(NA_integer_, length(fmeth))
  
  dplyr::case_when(
    fmeth %in% c("x", "b") ~ 0L,
    fmeth %in% c("p", "i") ~ 1L,
    fmeth == "m" ~ as.integer(mdc_fqc),
    fmeth == "a" ~ 2L,
    fmeth == "e" ~ 3L,
    TRUE ~ NA_integer_
  )
}

plot_filled <- function(data, var, type = c("fmeth", "fqc")) {
  
  type <- rlang::arg_match(type)
  var_name <- rlang::as_string(var)
  var <- rlang::ensym(var)
  
  var_f <- var_name %>% stringr::str_c(., "_f") %>% rlang::sym()
  
  var_fmeth <- var_name %>% stringr::str_c(., "_fmeth") %>% rlang::sym()
  
  if (type == "fmeth") {
    color <- var_fmeth
  } else if (type == "fqc") {
    color <- var_name %>% stringr::str_c(., "_fqc") %>% rlang::sym()
    data <- dplyr::mutate(data, !!color := factor(!!color))
  }
  
  type <- rlang::sym(type)
  
  filled <- dplyr::mutate(
    data, 
    !!var_f := dplyr::if_else(!is.na(!!var_fmeth), !!var_f, NA_real_),
    !!type := !!color
  )
  
  filled <- data %>%
    dplyr::filter(!is.na(!!var_fmeth)) %>%
    dplyr::mutate(!!type := !!color)
  
  data %>%
    ggplot2::ggplot(ggplot2::aes(timestamp, !!var_f)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::geom_point(
      data = filled, ggplot2::aes(color = !!type), size = 0.75, na.rm = TRUE
    ) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_bw()
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)
path_in <- file.path(wd, "processing", "03_biomet_qc", "data")

# Input file - biomet output with QC flags
data_input <- latest_version(path_in)
# Input file - biomet data from nearby site
aux_input <- latest_version(
  stringr::str_replace(path_in, settings$site, md$closest_site[1])
)
# Input file - processed ERA data for site location
era_input <- latest_version(
  file.path("~/Desktop", "DATA", "Flux", "JLL", "all", "era"), "era_proc"
)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing", "04_biomet_gapfill")


### Import and initialize data =================================================

cat("Importing data files...")

# Load the data
data <- readr::read_csv(
  data_input, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
)

# Set timestamp to local time zone to allow alignment with ERA data
data <- data %>%
  dplyr::mutate(
    timestamp = lubridate::force_tz(timestamp, tzone = md$tz_name)
  ) %>%
  # Add timestamp components
  add_time_comps()

# Load the auxilliary data
aux <- readr::read_csv(
  aux_input, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
)

# Set timestamp to local time zone
aux <- dplyr::mutate(
  aux,
  timestamp = lubridate::force_tz(timestamp, tzone = md$tz_name)
)

# Load the ERA data
era <- readr::read_csv(
  era_input, col_types = readr::cols(.default = readr::col_guess()), 
  progress = FALSE
)

# Create empty list to hold saved plots
plots <- list()

cat("done.\n")
cat("De-biasing auxilliary data...")

# Add lags to control list 
aux_ctrl <- purrr::list_modify(
  control,
  g = list(db_lag = get_best_lag(data$g, aux$g)),
  swc = list(db_lag = get_best_lag(data$swc, aux$swc)),
  ts = list(db_lag = get_best_lag(data$ts, aux$ts))
)
#purrr::map_dbl(aux_ctrl, pluck, "db_lag") # check lags

# Select and clean variables to be filled 
data_c <- data %>% select_clean(aux_vars) %>% tibble::as_tibble()
aux_c <- aux %>% select_clean(aux_vars) %>% tibble::as_tibble()

# Check harmonization
plots$aux_vars <- plot_harmonies(data_c, aux_c, "aux")

# Fit biomet and aux variables to generate de-biasing coefficients
aux_fits <- purrr::map(
  aux_vars, ~ debias_init(data_c, aux_c, ., ctrl = aux_ctrl)
)

# De-bias the auxilliary data
aux_d <- aux_vars %>% 
  purrr::map2(aux_fits, ~ debias(data_c, aux_c, .x, .y, ctrl = aux_ctrl)) %>%
  # Interpolate small gaps
  purrr::imap(~ fill_linear(.x, purrr::pluck(aux_ctrl, .y, "gf_max_i"))) %>%
  rlang::set_names(stringr::str_c(aux_vars, "_d")) %>%
  tibble::as_tibble() %>%
  dplyr::bind_cols(dplyr::select(data_c, timestamp), .)

# Check results
plots$aux_debias <- plot_fill_vars(aux_d, "aux_d")

cat("done.\n")
cat("Writing de-biased auxilliary data...")

# Save de-biased auxilliary data with documentation
aux_d_out <- file.path(
  path_out, "aux", paste0("biomet_aux_db_", tag_out, ".csv")
)
readr::write_csv(aux_d, aux_d_out)

# Create documentation for debiased auxilliary output
aux_d_docu <- purrr::prepend(
  settings, 
  list(
    files = c(data_input, aux_input),
    aux_site = purrr::pluck(md, "closest_site", 1),
    control = aux_ctrl,
    d_coefs = purrr::map(aux_fits, purrr::pluck, "coefficients")
  ),
  before = length(settings)
)
aux_d_docu_out <- stringr::str_replace(aux_d_out, ".csv", ".txt")
# Save documentation
sink(aux_d_docu_out)
print(aux_d_docu)
sink()

cat("done.\n")

### Data from external sources =================================================

cat("Harmonizing ERA data...")

# Prepare ERA data to be used for gap-filling
era_c <- era %>%
  tibble::as_tibble() %>%
  # Make copies for vars with multiple versions 
  dplyr::mutate(ta_ep = ta, ta_bm = ta, vpd_ep = vpd, vpd_bm = vpd) %>%
  # Note: ts0 = skin, ts1 = 0-7cm, ts2 = 7-28cm, ts3 = 28-100cm
  dplyr::select(
    timestamp, ta_ep, ta_bm, vpd_ep, vpd_bm, pa, sw_in, sw_out, lw_in, lw_out, 
    p_rain, swc1, swc2, swc3, ts1, ts2, ts3, ws
  ) %>%
  # Subset current year
  dplyr::right_join(dplyr::select(data, timestamp), by = "timestamp") %>% 
  tidyr::drop_na()

# Harmonize temporal resolution with ERA data
data_h <- era_c %>%
  dplyr::select(timestamp) %>%
  dplyr::left_join(data, by = "timestamp") %>%
  # Smooth differenced variables 
  # - balances inflated differences due to sensor noise (doesn't exist in ERA)
  dplyr::mutate(
    swc = roll_mean_real(swc, 7), 
    ts = roll_mean_real(ts, 11)
  ) %>%
  # p_rain is a sum, so need to add half-hours before joining
  dplyr::select(-p_rain) %>%
  dplyr::left_join(
    data %>%
      dplyr::mutate(timestamp = lubridate::round_date(timestamp, "hour")) %>%
      dplyr::group_by(timestamp) %>%
      dplyr::summarize(p_rain = sum(p_rain), .groups = "drop_last"),
    by = "timestamp"
  )

# Combine ERA soil var levels to account for changing WTD throughout year
era_c <- era_c %>%
  dplyr::bind_cols(dplyr::select(data_h, swc, ts)) %>%
  dplyr::mutate(dplyr::across(
    c(swc1, swc2, swc3), ~ apply_lag(.x, get_best_lag(swc, .x))
  )) %>%
  dplyr::mutate(dplyr::across(
    c(ts1, ts2, ts3), ~ apply_lag(.x, get_best_lag(ts, .x))
  )) %>%
  dplyr::mutate(
    swc = blend_vars(data_h, ., swc, swc2, swc3),
    ts = blend_vars(data_h, ., ts, ts2, ts3)
  )

# Set primary var if options are available (ta, vpd)
era_ctrl <- purrr::list_modify(
  control,
  swc = list(db_lag = get_best_lag(data_h$swc, era_c$swc)),
  ts = list(db_lag = get_best_lag(data_h$ts, era_c$ts))
)

# Harmonize cleaned data
data_hc <- era_c %>%
  dplyr::select(timestamp) %>%
  dplyr::left_join(data_c, by = "timestamp") %>%
  # Smooth differenced variables 
  # - balances inflated differences due to sensor noise (not the case in ERA)
  dplyr::mutate(
    swc = roll_mean_real(swc, 7),  
    ts = roll_mean_real(ts, 11)
  ) %>%
  # p_rain is a sum, so need to add half-hours before joining
  dplyr::select(-p_rain) %>%
  dplyr::left_join(
    data_c %>%
      dplyr::bind_cols(dplyr::select(data, sol_ang)) %>%
      dplyr::mutate(timestamp = lubridate::round_date(timestamp, "hour")) %>%
      dplyr::group_by(timestamp) %>%
      # Solar angle is used to reconstruct time series
      dplyr::summarize(
        p_rain = sum(p_rain), 
        sol_ang_h = sum(dplyr::if_else(sol_ang < 0, 0, sol_ang)) / 2, 
        .groups = "drop_last"
      ),
    by = "timestamp"
  )

# Plot harmonization
plots$era_vars <- data_hc %>% 
  dplyr::select(-ppfd_in, -sol_ang_h) %>% 
  plot_harmonies(era_c, "era")

# Calculate average distribution (in half-hrs) of hourly p_rain
# - needed to re-distribute p_rain when re-constructing original time series
p_rain_dist <- data_c %>%
  dplyr::mutate(
    timestamp = lubridate::round_date(timestamp, "hour"),
    p_rain_ind = dplyr::if_else(p_rain > 0, 1, 0)
  ) %>%
  dplyr::group_by(timestamp) %>%
  dplyr::summarize(
    p_rain_hrs = dplyr::na_if(sum(p_rain_ind, na.rm = TRUE), 0), 
    .groups = "drop_last"
  ) %>% 
  dplyr::summarize(mean(p_rain_hrs, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::pull()

cat("done.\n")
cat("De-biasing ERA data...")

# Fit biomet and aux variables to generate de-biasing coefficients
era_fits <- purrr::map(
  era_vars, ~ debias_init(data_hc, era_c, ., ctrl = era_ctrl)
)

# Create list with de-biasing coefficients
(era_coef <- purrr::map(era_fits, pluck, "coefficients"))

# De-bias the era data
era_d <- era_vars %>% 
  purrr::map2(era_fits, ~ debias(data_hc, era_c, .x, .y, ctrl = era_ctrl)) %>%
  rlang::set_names(stringr::str_c(era_vars, "_d")) %>%
  tibble::as_tibble() %>%
  dplyr::bind_cols(dplyr::select(data_hc, timestamp), .)

# Plot results
plots$era_debias <- plot_fill_vars(era_d, "era_d")

# Get coefficient of ppfd_in/sw_in relationship
frac_ppfd_era <- data_hc %>% 
  dplyr::select(ppfd_in) %>% 
  dplyr::bind_cols(dplyr::select(era_d, sw_in_d)) %>%
  tidyr::drop_na() %>% 
  dplyr::summarize(sum(sw_in_d) / sum(ppfd_in), .groups = "drop_last") %>% 
  purrr::pluck(1)
frac_ppfd_era

# Reconstruct original time series
era_df <- era_d %>%
  dplyr::rename_all(str_remove, "_d$") %>%
  dplyr::bind_cols(dplyr::select(data_hc, sol_ang_h)) %>%
  dplyr::right_join(
    dplyr::select(data, timestamp, sol_ang), by = "timestamp"
  ) %>%
  dplyr::arrange(timestamp) %>%
  dplyr::mutate(
    # Fill ERA solar angle ahead
    sol_ang_h = dplyr::if_else(is.na(sol_ang_h), lead(sol_ang_h), sol_ang_h),
    # Ratio of solar angles determines how to distribute sw_in to half-hours
    sol_ang_ratio = tidyr::replace_na(
      dplyr::if_else(sol_ang < 0, 0, sol_ang) / sol_ang_h, 0
    ),
    # p_rain is reconstructed according to typical hourly distribution
    p_rain = dplyr::case_when(
      round(p_rain_dist) == 2 & !is.na(p_rain) ~ p_rain / 2,
      round(p_rain_dist) == 2 & is.na(p_rain) ~ lead(p_rain) / 2,
      round(p_rain_dist) == 1 & !is.na(p_rain) ~ p_rain,
      round(p_rain_dist) == 1 & is.na(p_rain) & !is.na(lead(p_rain)) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  # Fill every other point with linear interpolation
  dplyr::mutate(dplyr::across(c(-timestamp, -p_rain), fill_linear)) %>%
  # Estimate first point using next available difference
  dplyr::mutate(dplyr::across(
    c(-timestamp, -p_rain), 
    ~ dplyr::if_else(
      dplyr::row_number() == 1, 
      dplyr::lead(.x, 1) - (dplyr::lead(.x, 2) - dplyr::lead(.x, 1)), .x
    )
  )) %>%
  dplyr::mutate(dplyr::across(c(sw_in, sw_out), ~ . * sol_ang_ratio)) %>%
  dplyr::mutate(
    # Adhere to local time zone
    timestamp = lubridate::with_tz(timestamp, md$tz_name),
    # Calculate ppfd_in equivalent of sw_in
    ppfd_in = sw_in / frac_ppfd_era
  ) %>%
  dplyr::select(-dplyr::contains("sol_ang")) %>%
  dplyr::rename_with(~ stringr::str_c(., "_df"), -timestamp)

# Plot results
plots$era_reconst <- era_df %>% 
  dplyr::select(-ppfd_in_df) %>% 
  plot_fill_vars("era_df")

cat("done.\n")
cat("Writing de-biased ERA data...")

# Save de-biased external data with documentation
era_df_out <- file.path(
  path_out, "era", paste0("era_db_reconst", tag_out, ".csv")
)
readr::write_csv(era_df, era_df_out)

# Create documentation for processed output
era_df_docu <- purrr::prepend(
  settings, 
  list(
    files = c(data_input, era_input),
    control = era_ctrl,
    d_coefs = purrr::map(era_fits, purrr::pluck, "coefficients"),
    p_rain_dist = p_rain_dist,
    frac_ppfd = frac_ppfd_era
  ),
  before = length(settings)
)
era_df_docu_out <- stringr::str_replace(era_df_out, ".csv", ".txt")

# Save documentation
sink(era_df_docu_out)
print(era_df_docu)
sink()

cat("done.\n")

### Prepare for gap-filling ====================================================

cat("Generating filled data with MDC algorithm...")

# Gather vars that will run MDC
mdc_names <- gf_vars %>% 
  tibble::enframe(value = "var") %>% 
  dplyr::mutate(var = purrr::simplify(var)) %>%
  dplyr::left_join(tibble::enframe(control, name = "var"), by = "var") %>%
  tidyr::unnest_wider(value) %>%
  dplyr::filter(stringr::str_detect(gf_order, "m")) %>%
  dplyr::pull(var)

# Initialize REddyProc class (necessary for MDC gap-filling)
biomet_mdc <- suppressMessages(suppressWarnings(sEddyProc$new(
  md$site_code, as.data.frame(data_c), mdc_names, "timestamp", 
  LatDeg = md$lat, LongDeg = md$lon, TimeZoneHour = md$tz
)))

# Run MDC for each variable (no qc var because vars are already cleaned)
# - see ?sEddyProc_sFillInit for description of output
for (i in seq_along(mdc_names)) {
  suppressMessages(suppressWarnings(biomet_mdc$sMDSGapFill(
    mdc_names[i], V1 = "none", FillAll = FALSE, isVerbose = FALSE
  )))
}
mdc_results <- biomet_mdc$sExportResults() %>% tibble::as_tibble()

cat("done.\n")
cat("Gathering all fill variables...")

# Initialize gap-filling data frame
biomet_f <- data_c %>%
  # Add backup variables
  dplyr::mutate(
    ppfd_in_sw = clean(data$ppfd_in_sw, data$qc_sw_in),
    sw_in_ppfd = clean(data$sw_in_ppfd, data$qc_ppfd_in)
  ) %>%
  dplyr::mutate(
    dplyr::across(c(ppfd_in_sw, sw_in_ppfd), clean, data$qc_biomet_all)
  ) %>%
  # Add debiased aux variables
  dplyr::left_join(aux_d, by = "timestamp") %>%
  # Add debiased era variables
  dplyr::left_join(era_df, by = "timestamp") %>%
  # Add MDC variables
  dplyr::bind_cols(mdc_results)

cat("done.\n")

### Fill gaps in biomet variables ==============================================

cat("Customizing gap filling based on available data...")

# Organize gapfill variables for each biomet variable
fill_data <- gf_vars %>% 
  purrr::map(~ gather_fill_data(biomet_f, ., ctrl = era_ctrl)) %>%
  rlang::set_names(gf_vars)

# Create reference vector indicating which gapfill variable to use
gf_meths <- purrr::imap(fill_data, ~ plan_fill(.x, ctrl = era_ctrl, var = .y))

# Create a gapfill QC vector based on method
gf_qcs <- gf_meths %>% 
  purrr::map(qc_biomet_fmeth) %>% 
  # MDC QC values may vary
  purrr::list_modify(!!!purrr::imap(
    magrittr::extract(gf_meths, mdc_names), 
    ~ qc_biomet_fmeth(
      .x, dplyr::pull(biomet_f, !!rlang::sym(stringr::str_c(.y, "_fqc")))
    )
  ))

cat("done.\n")
cat("Filling gaps...")

# Coalesce gapfill vars into one final vector
filled_vars <- purrr::map2(fill_data, gf_meths, gapfill_biomet)

# Add fill vars, methods, and QC variables to main data frame
data <- dplyr::bind_cols(
  data,
  filled_vars %>% 
    tibble::as_tibble() %>% 
    dplyr::rename_with(~ stringr::str_c(., "_f")),
  gf_meths %>% 
    tibble::as_tibble() %>% 
    dplyr::rename_with(~ stringr::str_c(., "_fmeth")) %>% 
    dplyr::mutate(dplyr::across(.fns = ~ dplyr::na_if(.x, "x"))),
  gf_qcs %>% 
    tibble::as_tibble() %>% 
    dplyr::rename_with(~ stringr::str_c(., "_fqc"))
)

cat("done.\n")
cat("Cleaning up results...")

# Set nighttime shortwave radiation to 0 and recode fmeth/fqc
# Used only when potential radiation is 0 (i.e. nighttime) AND sun is not
#   actively rising or setting, OR when var < 0
data <- data %>% 
  dplyr::mutate(dark = magrittr::and(
    sun_time == "night", lag(sun_time) != "set" & lead(sun_time) != "rise"
  )) %>%
  # Only the ends will be NA, and they are at midnight
  dplyr::mutate(dark = tidyr::replace_na(dark, TRUE)) %>%
  dplyr::mutate(dplyr::across(
    c(ppfd_in_f, sw_in_f, sw_out_f), ~ dplyr::if_else(dark | .x < 0, 0, .x)
  )) %>%
  dplyr::mutate(dplyr::across(
    c(ppfd_in_fmeth, sw_in_fmeth, sw_out_fmeth), ~ dplyr::if_else(dark, "p", .x)
  )) %>%
  dplyr::mutate(dplyr::across(
    c(ppfd_in_fqc, sw_in_fqc, sw_out_fqc), ~ dplyr::if_else(dark, 0L, .x)
  )) %>%
  dplyr::select(-dark)

# Set SWC to its maximum at saturation
data <- dplyr::mutate(
  data, 
  swc_f = pmin(swc_f, md$swc_sat), 
  swc_fmeth = dplyr::if_else(swc_f > md$swc_sat, "s", swc_fmeth),
  swc_fqc = dplyr::if_else(swc_f > md$swc_sat, 0L, swc_fqc)
)

# Recalculate net radiation using filled components
data <- dplyr::mutate(
  data, 
  netrad_f = sw_in_f + lw_in_f - (sw_out_f + lw_out_f)
)

# Fill RH using filled/converted VPD
data <- dplyr::mutate(
  data,
  rh_ep_f = dplyr::coalesce(
    clean(rh_ep, qc_rh_ep), bigleaf::VPD.to.rH(vpd_ep_f / 10, ta_ep_f) * 100
  ),
  rh_bm_f = dplyr::coalesce(
    clean(rh_bm, qc_rh_bm), bigleaf::VPD.to.rH(vpd_bm_f / 10, ta_ep_f) * 100
  )
)

# Additional variable calculations using gap-filled vars
data <- data %>%
  dplyr::mutate(s = 0) %>%
  dplyr::mutate(
    # Nighttime (Mauder et al. 2013 definition; overwrites original)
    night = sw_in_f <= 20,
    # Albedo
    albedo = sw_out_f / sw_in_f,
    # Water suface temperature (Stefan-Boltzmann law; Franz et al. 2016)
    tw_surf = (lw_out_f / (0.960 * 5.67e-8))^(1 / 4) - 273.15,
    # Potential evapotranspiration (Priestley & Taylor 1972)
    et_pot = as.data.frame(.) %>% 
      bigleaf::potential.ET("ta_ep_f", "pa_f", "netrad_f", "g_f", "s") %>%
      dplyr::pull(ET_pot)
    # Likely fog indicator
    #fog = detect_fog(rh_f, lw_in_f, lw_out_f, sw_in_f, p_rain_f)
  ) %>%
  dplyr::select(-s)

cat("done.\n")
cat("Plotting results...")

# Plot results
plots$gf_results <- data %>% 
  dplyr::select(timestamp, dplyr::ends_with("_f")) %>% 
  plot_fill_vars("var_f")

plots_fmeth <- gf_vars %>% 
  purrr::map(~ plot_filled(data, ., "fmeth")) %>% 
  rlang::set_names(stringr::str_c(gf_vars, "_fmeth"))
plots <- append(plots, plots_fmeth)

plots_fqc <- gf_vars %>% 
  purrr::map(~ plot_filled(data, ., "fqc")) %>% 
  rlang::set_names(stringr::str_c(gf_vars, "_fqc"))
plots <- append(plots, plots_fqc)

cat("done.\n")


### Save gap-filled data and all accessory files ===============================

cat("Writing filled data and diagnostic plots...")

# Set "official" vars
data <- dplyr::mutate(
  data,
  # Air temperature (from CO2 gas analyzer)
  ta = ta_ep,
  ta_f = ta_ep_f,
  # Relative humidity (from Biomet)
  rh = rh_bm,
  rh_f = rh_bm_f,
  # Vapor pressure deficit (calculated using "official" TA & RH)
  vpd = vpd_bm,
  vpd_f = vpd_bm_f
)

# Save full gap-filled dataset
data_out <- file.path(path_out, "data", paste0("biomet_gf_", tag_out, ".csv"))
data %>%
  remove_time_comps(-timestamp) %>%
  # Force timestamp back to UTC for data storage
  dplyr::mutate(timestamp = lubridate::force_tz(timestamp, "UTC")) %>%
  readr::write_csv(data_out)

# Create documentation for processed data output
data_docu <- purrr::prepend(
  settings, 
  list(
    files = c(data_input, aux_input, era_input),
    aux_site = md$closest_site[1],
    vars = gf_vars,
    control = control
  ),
  before = length(settings)
)
data_docu_out <- stringr::str_replace(data_out, ".csv", ".txt")
# Save documentation
sink(data_docu_out)
print(data_docu)
sink()

# Save Biomet gapfill details dataset
biomet_f_out <- file.path(
  path_out, "details", paste0("biomet_gf_details_", tag_out, ".csv")
)
biomet_f %>%
  # Force timestamp back to UTC for data storage
  dplyr::mutate(timestamp = lubridate::force_tz(timestamp, "UTC")) %>%
  readr::write_csv(biomet_f_out)

# Save one pdf document with all diagnostic/summary plots
plot_path <- file.path(
  path_out, "plots", paste0("biomet_gf_plots_", tag_out, ".pdf")
)
pdf(plot_path)
purrr::map(plots, print)
dev.off()

# Save models and uncertainties
aux_fit_df <- aux_fits %>% 
  purrr::map(broom::glance) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(var = names(aux_fits), data = "aux")
era_fit_df <- era_fits %>% 
  purrr::map(broom::glance) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(var = names(era_fits), data = "era") 
fit_df <- dplyr::bind_rows(aux_fit_df, era_fit_df) %>% 
  dplyr::select(var, data, dplyr::everything())

fit_df_out <- file.path(
  path_out, "fits", paste0("biomet_db_fits_", tag_out, ".csv")
)
readr::write_csv(fit_df, fit_df_out)

end_time <- Sys.time()
elapsed_time <- round(unclass(end_time - start_time), 1)

cat("done.\n")
cat(
  "Finished processing in ", elapsed_time, " ", attr(elapsed_time, "units"), 
  ".\n", sep = ""
)
cat("Output located in", path_out, "\n")

# Uncertainty something like this
#biomet_c %>% 
#  mutate(ta_gf = gapfill_biomet(
#    fill_data[[1]], plan_fill_all(fill_data[[1]], "ta", "imae")
#  )) %>% 
#  do(lm = lm(ta_ep ~ ta_gf, data = .)) %>% 
#  pluck(1, 1) %>% summary()

# Finished
