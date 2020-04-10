### ============================================================================
# Biomet data gapfilling =======================================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

# Load the required packages
suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
library(openeddy)
library(REddyProc)
library(lubridate)
library(tidyverse)
# Packages needed but not loaded: solartime, bigleaf, broom

source("~/Desktop/DATA/Flux/tools/engine/biomet_gapfill.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")
source("~/Desktop/DATA/Flux/tools/reference/var_attributes.R")
source("~/Desktop/DATA/Flux/tools/reference/gf_control.R")


### Helper functions ===========================================================

acf_1 <- function(x) {
  acf <- acf(x, plot = FALSE, na.action = na.pass)$acf
  acf %>% array_branch() %>% flatten_dbl() %>% pluck(2)
}

get_best_lag <- function(y, x, lag.max = 12) {
  
  dy <- y - dplyr::lag(y)
  dx <- x - dplyr::lag(x)
  
  ccf <- ccf(dx, dy, lag.max = lag.max, na.action = na.pass, plot = FALSE)
  
  tbl <- tibble::tibble(
    lag = ccf %>% 
      purrr::pluck("lag") %>%
      purrr::array_branch() %>%
      purrr::simplify(),
    acf = ccf %>% 
      purrr::pluck("acf") %>%
      purrr::array_branch() %>%
      purrr::simplify()
  )
  
  lag <- tbl %>%
    dplyr::arrange(desc(acf)) %>%
    dplyr::slice(1) %>%
    purrr::pluck("lag")
  
  (-lag)
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
  
  vars <- map_chr(vars, rlang::as_string)
  df <- dplyr::select_at(data, vars(timestamp, all_of(vars)))
  
  get_qc_var <- function(x) {
    rlang::as_string(x) %>%
      stringr::str_c("qc_", .) %>%
      rlang::sym()
  }
  
  for (i in 2:length(df)) {
    name <- names(df)[i]
    qc_name <- get_qc_var(name)
    df[, name] <- clean(df[, name], dplyr::pull(data, !!qc_name))
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

plot_harmonies <- function(data, data_aux, aux_suffix) {
  
  sfx <- aux_suffix
  sfx_ <- stringr::str_c("_", sfx)
  
  data_aux %>%
    dplyr::rename_at(vars(-timestamp), ~stringr::str_c(., sfx_)) %>%
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
  #browser()
  
  tbl <- tibble::tibble(
    y = dplyr::pull(data, !!var1),
    x = dplyr::pull(aux_data, !!var2)
  )
  
  # Difference if requested (useful for highly autocorrelated variables)
  if (diff) tbl <- dplyr::mutate_all(tbl, ~ . - dplyr::lag(., 1))
  
  # Apply lag if indicated
  if (lag < 0) {
    tbl <- mutate(tbl, x = dplyr::lead(x, -lag))
  } else {
    tbl <- mutate(tbl, x = dplyr::lag(x, lag))
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
      dplyr::summarize(beta = sum(abs(y)) / sum(abs(x))) %>%
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
    tbl <- mutate(tbl, x = dplyr::lead(x, -lag))
  } else {
    tbl <- mutate(tbl, x = dplyr::lag(x, lag))
  } 
  
  # If differences, reconstruct data using time series fill
  if (diff) {
    out <- fill_along(tbl$y, tbl$x, coef(fit)[2], lag = 0, align = TRUE)
    #out <- dplyr::if_else(is.na(tbl$x - dplyr::lag(tbl$x, 1)), NA_real_, out)
  } else {
    out <- coef(fit)[1] + coef(fit)[2] * tbl$x
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
      dplyr::mutate_at(dplyr::vars(p1, p2), ~ pmax(pmin(. / tot, 1), 0)) %>%
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
      dplyr::mutate_at(
        dplyr::vars(p1, p2), ~ dplyr::if_else(is.nan(.), NA_real_, .)
      ) %>%
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
      dplyr::mutate_at(dplyr::vars(p1, p2, p3), ~ pmax(pmin(. / tot, 1), 0)) %>%
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
      dplyr::mutate_at(
        dplyr::vars(p1, p2, p3), ~ dplyr::if_else(is.nan(.), NA_real_, .)
      ) %>%
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
      dplyr::summarize(sum(!!var) / sum(!!b_var)) %>% 
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
    dplyr::select_at(dplyr::vars(tidyselect::all_of(names))) %>%
    dplyr::mutate(
      gap = dplyr::if_else(is.na(x), 1L, 0L),
      gap_id = dplyr::pull(tidy_rle(gap), id),
      gap_len = dplyr::pull(tidy_rle(gap), lengths)
    ) %>%
    dplyr::group_by(gap_id) %>%
    dplyr::summarize_at(
      dplyr::vars(-gap, -dplyr::group_cols()), ~ length(na.omit(.))
    ) %>%
    dplyr::ungroup()
  
  cover <- gaps %>%
    dplyr::mutate(
      # Do gap-fill vars cover at least part of a gap?
      any = purrr::imap_dfr(
        ., ~ dplyr::if_else(.x != 0, .y, NA_character_)
      ) %>% 
        dplyr::select_at(dplyr::vars(-dplyr::starts_with("gap"))) %>% 
        purrr::pmap_chr(dplyr::coalesce), 
      # Do gap-fill vars cover entire gaps?
      all = purrr::imap_dfr(
        ., ~ dplyr::if_else(.x == gap_len, .y, NA_character_)
      ) %>% 
        dplyr::select_at(dplyr::vars(-dplyr::starts_with("gap"))) %>% 
        purrr::pmap_chr(dplyr::coalesce)
    ) %>% 
    dplyr::mutate(
      plan = purrr::pmap_chr(dplyr::select(., all, any), dplyr::coalesce)
    )
  
  plan_a <- rep(cover$all, times = cover$gap_len)
  
  # Point-focused approach (Plan B)
  # Replace non-missing values with the gapfill variable name
  tbl <- list %>% 
    purrr::imap(~ replace_real(.x, .y)) %>%
    tibble::as_tibble()
  
  # Make sure fill vars are in the correct order
  inorder <- dplyr::select_at(tbl, dplyr::vars(tidyselect::all_of(names)))
  
  # Coalesce to form "Plan B"
  plan_b <- dplyr::coalesce(!!!as.list(inorder))
  
  plan <- dplyr::coalesce(plan_a, plan_b)
  
  plan
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
path_in <- file.path(wd, "processing_data", "04_biomet_qc", "output")

# Input file - biomet output with QC flags
biomet_input <- latest_version(path_in, "biomet_qc")
# Input file - biomet data from nearby sites
aux_input <- latest_version(
  stringr::str_replace(path_in, settings$site, md$closest_site[1]), 
  "biomet_qc"
)
# Input file - processed ERA data for site location
era_input <- latest_version(
  file.path("~/Desktop", "DATA", "Flux", "JLL", "all", "output"), "era_proc"
)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing_data", "05_biomet_gapfill", "output")


### Perform Biomet data initialization =========================================

# Load the Biomet file
cat("Importing Biomet file.\n")
biomet <- read.csv(biomet_input, stringsAsFactors = FALSE)

# Add timestamp components, some initialization
biomet <- biomet %>%
  mutate(
    timestamp = ymd_hms(timestamp, tz = md$tz_name),
    # Calculate vpd (easier to harmonize than RH)
    # - can I use ta_ep here? using ta for now since from the same instrument
    vpd = bigleaf::rH.to.VPD(rh / 100, ta) * 10,
    qc_vpd = combine_flags(qc_rh, qc_ta),
    qc_vpd_ep = combine_flags(qc_rh, qc_ta_ep)
  ) %>%
  add_time_comps()

# Create empty list to hold saved plots
plots <- list()


### Import & de-bias auxilliary data to be used in gap-filling =================

cat("Importing auxilliary file.\n")
# Load the auxilliary file
aux <- read.csv(aux_input, stringsAsFactors = FALSE)
# Prepare to be used for gap-filling
aux <- mutate(
  aux,
  timestamp = ymd_hms(timestamp, tz = md$tz_name),
  # Can I use ta_ep here? using ta for now since from the same instrument
  vpd = bigleaf::rH.to.VPD(rh / 100, ta) * 10,
  qc_vpd = combine_flags(qc_rh, qc_ta),
  # Match appropriate TA vars 
  !!sym(md$ta_var) := !!sym(pluck(site_metadata, md$closest_site[1], "ta_var"))
)

cat("De-biasing auxilliary data.\n")
# Set primary var if options are available (ta, vpd)
aux_vars <- list_modify(aux_vars, ta = md$ta_var)

# Add lags to control list 
aux_ctrl <- list_modify(
  control,
  g = list(db_lag = get_best_lag(biomet$g, aux$g)),
  swc = list(db_lag = get_best_lag(biomet$swc, aux$swc)),
  ts = list(db_lag = get_best_lag(biomet$ts, aux$ts))
)
map_dbl(aux_ctrl, pluck, "db_lag") # check lags

# Select and clean variables to be filled 
biomet_c <- biomet %>% select_clean(aux_vars) %>% as_tibble()
aux_c <- aux %>% select_clean(aux_vars) %>% as_tibble()

# Check harmonization
plots$aux_vars <- plot_harmonies(biomet_c, aux_c, "aux")
#plots$aux_vars

# Fit biomet and aux variables to generate de-biasing coefficients
aux_fits <- map(aux_vars, ~ debias_init(biomet_c, aux_c, ., ctrl = aux_ctrl))

# Check de-biasing coefficients
#map(aux_fits, pluck, "coefficients") %>% enframe() %>% unnest_wider(value)

# De-bias the auxilliary data
aux_d <- aux_vars %>% 
  map2(aux_fits, ~ debias(biomet_c, aux_c, .x, .y, ctrl = aux_ctrl)) %>%
  # Interpolate small gaps
  imap(~ fill_linear(.x, pluck(aux_ctrl, .y, "gf_max_i"))) %>%
  set_names(str_c(aux_vars, "_d")) %>%
  as_tibble() %>%
  bind_cols(select(biomet_c, timestamp), .)

# Check results
plots$aux_debias <- plot_fill_vars(aux_d, "aux_d")
#plots$aux_debias

cat("Writing de-biased auxilliary data.\n")
# Save de-biased auxilliary data with documentation
aux_d_out <- file.path(path_out, paste0("biomet_aux_db_", tag_out, ".csv"))
write.csv(aux_d, aux_d_out, row.names = FALSE)

# Create documentation for processed Biomet output
aux_d_docu <- append(
  settings, list(
    files = c(biomet_input, aux_input),
    aux_site = md$closest_site[1],
    control = aux_ctrl,
    d_coefs = map(aux_fits, pluck, "coefficients")
  )
)
aux_d_docu_out <- str_replace(aux_d_out, ".csv", ".txt")
# Save documentation
sink(aux_d_docu_out)
aux_d_docu
sink()


### Data from external sources =================================================

cat("Importing ERA file.\n")
# ERA
era <- read.csv(era_input, stringsAsFactors = FALSE)
era <- mutate(era, timestamp = ymd_hms(timestamp, tz = "UTC"))

cat("Harmonizing ERA data.\n")
# Prepare to be used for gap-filling
era_c <- era %>%
  as_tibble() %>%
  # Note: ts0 = skin, ts1 = 0-7cm, ts2 = 7-28cm, ts3 = 28-100cm
  select(
    timestamp, ta, pa_ep, sw_in, sw_out, lw_in, lw_out, vpd, p_rain, 
    swc1, swc2, swc3, ts0, ts1, ts2, ts3
  ) %>%
  mutate(
    # Estimate g at different levels 
    g1 = ts0 - ts2, g2 = ts1 - ts2, g3 = ts2 - ts3,
    # Match appropriate TA vars 
    !!sym(md$ta_var) := ta
  ) %>%
  # Update g using swc as proxy for wtd
  mutate_at(vars(g1, g2, g3), ~ . * scales::rescale(1 - swc3)) %>%
  # Adjust swc for offset when wetland is likely inundated
  mutate_at(
    vars(swc2, swc3), flatten_period, 
    n = pluck(md, "drydown_start", as.character(settings$year)) * 24
  ) %>%
  # Subset current year
  right_join(select(biomet, timestamp), by = "timestamp") %>% 
  drop_na()

# Harmonize biomet temporal resolution with ERA data
biomet_h <- era_c %>%
  select(timestamp) %>%
  left_join(biomet, by = "timestamp") %>%
  # Smooth differenced variables 
  # - balances inflated differences due to sensor noise (doesn't exist in ERA)
  mutate(
    g = roll_mean_real(g, 3), swc = roll_mean_real(swc, 7), 
    ts = roll_mean_real(ts, 11)
  ) %>%
  # p_rain is a sum, so need to add half-hours before joining
  select(-p_rain) %>%
  left_join(
    biomet %>%
      mutate(timestamp = round_date(timestamp, "hour")) %>%
      group_by(timestamp) %>%
      summarize(p_rain = sum(p_rain)),
    by = "timestamp"
  )

# Combine ERA soil var levels to account for changing WTD throughout year
era_c <- era_c %>%
  bind_cols(select(biomet_h, g, swc, ts)) %>%
  mutate_at(vars(g1, g2, g3), ~ apply_lag(., get_best_lag(g, .))) %>%
  mutate_at(vars(swc1, swc2, swc3), ~ apply_lag(., get_best_lag(swc, .))) %>%
  mutate_at(vars(ts1, ts2, ts3), ~ apply_lag(., get_best_lag(ts, .))) %>%
  mutate(
    g = blend_vars(biomet_h, ., g, g1, g2, g3, diff = TRUE),
    swc = blend_vars(biomet_h, ., swc, swc2, swc3),
    ts = blend_vars(biomet_h, ., ts, ts2, ts3)
)

# Set primary var if options are available (ta, vpd)
era_vars <- list_modify(era_vars, ta = md$ta_var)
era_ctrl <- list_modify(
  control,
  g = list(db_lag = get_best_lag(biomet_h$g, era_c$g)),
  swc = list(db_lag = get_best_lag(biomet_h$swc, era_c$swc)),
  ts = list(db_lag = get_best_lag(biomet_h$ts, era_c$ts))
)
#map_dbl(era_ctrl, pluck, "db_lag") # check lags

# Harmonize cleaned biomet
biomet_hc <- era_c %>%
  select(timestamp) %>%
  left_join(biomet_c, by = "timestamp") %>%
  # Smooth differenced variables 
  # - balances inflated differences due to sensor noise (not the case in ERA)
  mutate(
    g = roll_mean_real(g, 3), swc = roll_mean_real(swc, 7),  
    ts = roll_mean_real(ts, 11)
  ) %>%
  # p_rain is a sum, so need to add half-hours before joining
  select(-p_rain) %>%
  left_join(
    biomet_c %>%
      bind_cols(select(biomet, sol_ang)) %>%
      mutate(timestamp = round_date(timestamp, "hour")) %>%
      group_by(timestamp) %>%
      # Solar angle is used to reconstruct time series
      summarize(
        p_rain = sum(p_rain), 
        sol_ang_h = sum(if_else(sol_ang < 0, 0, sol_ang)) / 2
      ),
    by = "timestamp"
  )

# Check harmonization
plots$era_vars <- biomet_hc %>% 
  select(-ppfd_in, -sol_ang_h) %>% 
  plot_harmonies(era_c, "era")
#plots$era_vars

# Calculate average distribution (in half-hrs) of hourly p_rain
# - needed to re-distribute p_rain when re-constructing original time series
p_rain_dist <- biomet_c %>%
  mutate(
    timestamp = round_date(timestamp, "hour"),
    p_rain_ind = if_else(p_rain > 0, 1, 0)
  ) %>%
  group_by(timestamp) %>%
  summarize(p_rain_hrs = na_if(sum(p_rain_ind, na.rm = TRUE), 0)) %>% 
  summarize(mean(p_rain_hrs, na.rm = TRUE)) %>%
  pull()

cat("De-biasing ERA data.\n")
# Fit biomet and aux variables to generate de-biasing coefficients
era_fits <- map(era_vars, ~ debias_init(biomet_hc, era_c, ., ctrl = era_ctrl))

# Create list with de-biasing coefficients
(era_coef <- map(era_fits, pluck, "coefficients"))

# De-bias the era data
era_d <- era_vars %>% 
  map2(era_fits, ~ debias(biomet_hc, era_c, .x, .y, ctrl = era_ctrl)) %>%
  set_names(str_c(era_vars, "_d")) %>%
  as_tibble() %>%
  bind_cols(select(biomet_hc, timestamp), .)

# Check results
plots$era_debias <- plot_fill_vars(era_d, "era_d")
#plots$era_debias

# Get coefficient of ppfd_in/sw_in relationship
frac_ppfd_era <- biomet_hc %>% 
  select(ppfd_in) %>% 
  bind_cols(select(era_d, sw_in_d)) %>%
  drop_na() %>% 
  summarize(sum(sw_in_d) / sum(ppfd_in)) %>% 
  pluck(1)
frac_ppfd_era

# Reconstruct original time series
era_df <- era_d %>%
  rename_all(str_remove, "_d$") %>%
  bind_cols(select(biomet_hc, sol_ang_h)) %>%
  right_join(select(biomet, timestamp, sol_ang), by = "timestamp") %>%
  mutate(
    # Fill ERA solar angle ahead
    sol_ang_h = if_else(is.na(sol_ang_h), lead(sol_ang_h), sol_ang_h),
    # Ratio of solar angles determines how to distribute sw_in to half-hours
    sol_ang_ratio = replace_na(if_else(sol_ang < 0, 0, sol_ang) / sol_ang_h, 0),
    # p_rain is reconstructed according to typical hourly distribution
    p_rain = dplyr::case_when(
      round(p_rain_dist) == 2 & !is.na(p_rain) ~ p_rain / 2,
      round(p_rain_dist) == 2 & is.na(p_rain) ~ lead(p_rain) / 2,
      round(p_rain_dist) == 1 & !is.na(p_rain) ~ p_rain,
      round(p_rain_dist) == 1 & is.na(p_rain) & !is.na(lead(p_rain)) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate_at(vars(-timestamp, -p_rain), fill_linear) %>%
  mutate_at(vars(sw_in, sw_out), ~ . * sol_ang_ratio) %>%
  mutate(
    # Adhere to local time zone
    timestamp = with_tz(timestamp, md$tz_name),
    # Calculate ppfd_in equivalent of sw_in
    ppfd_in = sw_in / frac_ppfd_era
  ) %>%
  select_at(vars(-contains("sol_ang"))) %>%
  rename_at(vars(-timestamp), ~str_c(., "_df"))

# Check results
plots$era_reconst <- era_df %>% select(-ppfd_in_df) %>% plot_fill_vars("era_df")
#plots$era_reconst

# Check the final dataset
#era_df %>% summarize(first(timestamp), last(timestamp)) # start/end?

cat("Writing de-biased ERA data.\n")
# Save de-biased external data with documentation
era_df_out <- file.path(path_out, paste0("era_db_reconst", tag_out, ".csv"))
write.csv(era_df, era_df_out, row.names = FALSE)

# Create documentation for processed output
era_df_docu <- append(
  settings, 
  list(
    files = c(biomet_input, era_input),
    control = era_ctrl,
    d_coefs = map(era_fits, pluck, "coefficients"),
    p_rain_dist = p_rain_dist,
    frac_ppfd = frac_ppfd_era
  )
)
era_df_docu_out <- str_replace(era_df_out, ".csv", ".txt")

# Save documentation
sink(era_df_docu_out)
era_df_docu
sink()


### Prepare for gap-filling ====================================================

# Set primary var if options are available (ta, vpd)
gf_vars <- list_modify(gf_vars, ta = md$ta_var)

cat("Generating filled data with MDC algorithm.\n")
# Gather vars that will run MDC
mdc_names <- gf_vars %>% 
  enframe(value = "var") %>% 
  mutate(var = simplify(var)) %>%
  left_join(enframe(control, name = "var"), by = "var") %>%
  unnest_wider(value) %>%
  filter(str_detect(gf_order, "m")) %>%
  pull(var)

# Initialize REddyProc class (necessary for MDC gap-filling)
biomet_mdc <- suppressMessages(suppressWarnings(sEddyProc$new(
  md$site_code, as.data.frame(biomet_c), mdc_names, "timestamp", 
  LatDeg = md$lat, LongDeg = md$lon, TimeZoneHour = md$tz
)))

# Run MDC for each variable (no qc var because vars are already cleaned)
# - see ?sEddyProc_sFillInit for description of output
for (i in seq_along(mdc_names)) {
  suppressMessages(suppressWarnings(biomet_mdc$sMDSGapFill(
    mdc_names[i], V1 = "none", FillAll = FALSE, isVerbose = FALSE
  )))
}
mdc_results <- biomet_mdc$sExportResults() %>% as_tibble()

cat("Gathering all fill variables.\n")
# Initialize gap-filling data frame
ta_alt <- pluck(control, md$ta_var, "gf_backup")
biomet_f <- biomet_c %>%
  # Add backup variables
  mutate(
    !!sym(ta_alt) := clean(
      pull(biomet, !!sym(ta_alt)), pull(biomet, !!sym(str_c("qc_", ta_alt)))
    ),
    vpd_ep = clean(biomet$vpd_ep, biomet$qc_vpd_ep),
    ppfd_in_sw = clean(biomet$ppfd_in_sw, biomet$qc_sw_in),
    sw_in_ppfd = clean(biomet$sw_in_ppfd, biomet$qc_ppfd_in)
  ) %>%
  mutate_at(vars(ppfd_in_sw, sw_in_ppfd), clean, biomet$qc_biomet_all) %>%
  # Add debiased aux variables
  left_join(aux_d, by = "timestamp") %>%
  # Add debiased era variables
  left_join(era_df, by = "timestamp") %>%
  # Add MDC variables
  bind_cols(mdc_results)


### Fill gaps in biomet variables ==============================================

# Organize gapfill variables for each biomet variable
fill_data <- gf_vars %>% 
  map(~ gather_fill_data(biomet_f, ., ctrl = era_ctrl)) %>%
  set_names(gf_vars)

cat("Customizing gap filling based on available data.\n")
# Create reference vector indicating which gapfill variable to use
gf_meths <- imap(fill_data, ~ plan_fill(.x, ctrl = era_ctrl, var = .y))

# Create a gapfill QC vector based on method
gf_qcs <- gf_meths %>% 
  map(qc_biomet_fmeth) %>% 
  # MDC QC values may vary
  list_modify(!!!imap(
    magrittr::extract(gf_meths, mdc_names), 
    ~ qc_biomet_fmeth(.x, pull(biomet_f, !!sym(str_c(.y, "_fqc"))))
  ))

cat("Filling gaps.\n")
# Coalesce gapfill vars into one final vector
filled_vars <- map2(fill_data, gf_meths, gapfill_biomet)

# Add fill vars, methods, and QC variables to main data frame
biomet <- bind_cols(
  biomet,
  filled_vars %>% as_tibble() %>% rename_all(~ str_c(., "_f")),
  gf_meths %>% 
    as_tibble() %>% 
    rename_all(~ str_c(., "_fmeth")) %>% 
    mutate_all(na_if, "x"),
  gf_qcs %>% as_tibble() %>% rename_all(~ str_c(., "_fqc"))
)

cat("Cleaning up results.\n")
# Set nighttime shortwave radiation to 0 and recode fmeth/fqc
# Used only when potential radiation is 0 (i.e. nighttime) AND sun is not
#   actively rising or setting, OR when var < 0
biomet <- biomet %>% 
  mutate(night2 = magrittr::and(
    sun_time == "night", lag(sun_time) != "set" & lead(sun_time) != "rise"
  )) %>%
  # Only the ends will be NA, and they are at midnight
  mutate(night2 = replace_na(night2, TRUE)) %>%
  mutate_at(
    vars(ppfd_in_f, sw_in_f, sw_out_f), ~ if_else(night2 | . < 0, 0, .)
  ) %>%
  mutate_at(
    vars(ppfd_in_fmeth, sw_in_fmeth, sw_out_fmeth), ~ if_else(night2, "p", .)
  ) %>%
  mutate_at(
    vars(ppfd_in_fqc, sw_in_fqc, sw_out_fqc), ~ if_else(night2, 0L, .)
  ) %>%
  select(-night2)

# Set SWC to its maximum at saturation
biomet <- mutate(
  biomet, 
  swc_f = pmin(swc_f, md$swc_sat), 
  swc_fmeth = if_else(swc_f > md$swc_sat, "s", swc_fmeth),
  swc_fqc = if_else(swc_f > md$swc_sat, 0L, swc_fqc)
)

# Recalculate net radiation using filled components
biomet <- mutate(
  biomet, 
  netrad_f = sw_in_f + lw_in_f - (sw_out_f + lw_out_f), 
  netrad_fqc = if_else(swc_f > md$swc_sat, 0L, swc_fqc)
)

# Fill RH using filled/converted VPD
ta_gf_var <- str_c(md$ta_var, "_f")
# Make sure that both TA vars are filled
ta_db_fit <- biomet %>%
  mutate(!!sym(ta_alt) := pull(biomet, !!sym(ta_alt))) %>%
  do(fit = lm(
    !!sym(ta_alt) ~ !!sym(ta_gf_var), na.action = na.exclude, data = .
  )) %>%
  pluck(1, 1)
biomet <- biomet %>% 
  mutate(
    # Fill alternate TA var
    !!sym(str_c(ta_alt, "_f")) := predict(
      ta_db_fit, select(biomet, !!sym(ta_gf_var))
    ),
    # Reconstruct RH
    rh_f = coalesce(
      clean(rh, qc_rh), bigleaf::VPD.to.rH(vpd_f / 10, ta_f) * 100
    )
  )

# Additional variable calculations using gap-filled vars
biomet <- biomet %>% 
  mutate(
    # Albedo
    albedo = sw_out_f / sw_in_f,
    # Water suface temperature
    tw_surf = (lw_out / (0.960 * 5.67e-8))^(1 / 4) - 273.15,
    # Potential evapotranspiration
    et_pot = pull(
      bigleaf::potential.ET(., "ta_ep_f", "pa_ep_f", "netrad_f", "g_f"), ET_pot
    ),
    # Likely fog indicator
    fog = detect_fog(rh_f, lw_in_f, lw_out_f, sw_in_f, p_rain_f)
  )

cat("Plotting results.\n")
# Plot results
plots$gf_results <- biomet %>% 
  select_at(vars(timestamp, ends_with("_f"))) %>% 
  plot_fill_vars("var_f")

plots_fmeth <- gf_vars %>% 
  map(~plot_filled(biomet, ., "fmeth")) %>% 
  set_names(str_c(gf_vars, "_fmeth"))
plots <- append(plots, plots_fmeth)

plots_fqc <- gf_vars %>% 
  map(~plot_filled(biomet, ., "fqc")) %>% 
  set_names(str_c(gf_vars, "_fqc"))
plots <- append(plots, plots_fqc)


### Save gap-filled data and all accessory files ===============================

cat("Writing filled data and diagnostic plots.\n")
# Save full gap-filled Biomet dataset
biomet_out <- file.path(path_out, paste0("biomet_gf_", tag_out, ".csv"))
write.csv(biomet, biomet_out, row.names = FALSE)

# Create documentation for processed Biomet output
biomet_docu <- append(
  settings, list(
    files = c(biomet_input, aux_input, era_input),
    aux_site = md$closest_site[1],
    vars = gf_vars,
    control = control
  )
)
biomet_docu_out <- str_replace(biomet_out, ".csv", ".txt")
# Save documentation
sink(biomet_docu_out)
biomet_docu
sink()

# Save Biomet gapfill details dataset
biomet_f_out <- str_replace(biomet_out, "_gf_", "_gf_details_")
write.csv(biomet_f, biomet_f_out, row.names = FALSE)

# Save one pdf document with all diagnostic/summary plots
plot_path <- file.path(path_out, paste0("biomet_gf_plots_", tag_out, ".pdf"))
pdf(plot_path)
map(plots, ~ .)
dev.off()

# Save models and uncertainties
aux_fit_df <- aux_fits %>% 
  map(broom::glance) %>% 
  bind_rows() %>% 
  mutate(var = names(aux_fits), data = "aux")
era_fit_df <- era_fits %>% 
  map(broom::glance) %>% 
  bind_rows() %>% 
  mutate(var = names(era_fits), data = "era") 
fit_df <- bind_rows(aux_fit_df, era_fit_df) %>% select(var, data, everything())

fit_df_out <- file.path(path_out, paste0("biomet_db_fits_", tag_out, ".csv"))
write.csv(fit_df, fit_df_out, row.names = FALSE)

# Uncertainty something like this
#biomet_c %>% 
#  mutate(ta_gf = gapfill_biomet(
#    fill_data[[1]], plan_fill_all(fill_data[[1]], "ta", "imae")
#  )) %>% 
#  do(lm = lm(ta_ep ~ ta_gf, data = .)) %>% 
#  pluck(1, 1) %>% summary()

# Finished
