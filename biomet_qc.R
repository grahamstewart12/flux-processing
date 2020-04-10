### ============================================================================
# Biomet data quality control ==================================================
### ============================================================================

# Purpose: 

# References:

# Estévez, J., Gavilán, P., Giráldez, J.V., 2011. Guidelines on validation 
# procedures for meteorological data from automatic weather stations. Journal of 
# Hydrology 402, 144–154. https://doi.org/10.1016/j.jhydrol.2011.02.031

# Olson, R.J., Holladay, S.K., Cook, R.B., Falge, E., Baldocchi, D., Gu, L., 
# 2004. FLUXNET: Database of fluxes, site characteristics, and flux-community 
# information. Oak Ridge National Laboratory, Oak Ridge, TN. 
# https://doi.org/10.2172/1184413

# ONEFlux Processing Pipeline
# https://github.com/AmeriFlux/ONEFlux/blob/master/oneflux_steps/qc_auto/src

# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/Desktop/RESEARCH/fluxtools", quiet = TRUE)
library(openeddy)
library(lubridate)
library(tidyverse)

#source("~/Desktop/DATA/Flux/tools/engine/biomet_qc_v2.R")
source("~/Desktop/DATA/Flux/tools/reference/var_attributes.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")


### Helper functions ===========================================================

acf_1 <- function(x) {
  acf <- acf(x, plot = FALSE, na.action = na.pass)$acf
  acf %>% 
    purrr::array_branch() %>% 
    purrr::flatten_dbl() %>% 
    purrr::pluck(2)
}

validate_flags <- function(data, var1, var2) {
  
  var1_name <- rlang::ensym(var1) %>% rlang::as_string()
  var2_name <- rlang::ensym(var2) %>% rlang::as_string()
  
  # Save the original flag names
  var1_names <- data %>% purrr::pluck(var1_name) %>% names()
  var2_names <- data %>% purrr::pluck(var2_name) %>% names()
  
  # Extract flags for each var
  var1_flags <- data %>%
    purrr::pluck(var1_name) %>%
    dplyr::rename_all(
      stringr::str_remove_all, stringr::str_c("qc_", var1_name, "_")
    )
  var2_flags <- data %>% 
    purrr::pluck(var2_name) %>%
    dplyr::rename_all(
      stringr::str_remove_all, stringr::str_c("qc_", var2_name, "_")
    )
  
  # Subset flags used in both vars
  flags <- stringr::str_subset(names(var1_flags), names(var2_flags))
  
  # Check by FLAG, not by var
  for (i in seq_along(flags)) {
    
    # Not applied to null or gap tests
    if (flags[i] %in% c("null", "gap")) next
    
    # Find points that are flagged for both vars 
    val_ind <- which(var1_flags[, flags[i]] == 2 & var2_flags[, flags[i]] == 2)
    
    # Downgrade flags for these points
    var1_flags[val_ind, flags[i]] <- 1L
    var2_flags[val_ind, flags[i]] <- 1L
  }
  
  # Give back original names and return modified list
  out <- data
  names(var1_flags) <- var1_names
  purrr::pluck(out, var1_name) <- var1_flags
  names(var2_flags) <- var2_names
  purrr::pluck(out, var2_name) <- var2_flags
  out
}

biomet_auto_flags <- function(data, var, p_rain) {
  
  # Might be better to break this up into units (i.e. individual tests)
  # - check biomet_qc_v3 for some of these
  
  var_name <- rlang::ensym(var) %>% rlang::as_string()
  var <- rlang::enquo(var)
  p_rain <- rlang::enquo(p_rain)
  
  cat(var_name, "...", sep = "")
  
  tbl <- data %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      timestamp = timestamp,
      year = lubridate::year(timestamp - 900),
      month = lubridate::month(timestamp - 900),
      date = lubridate::date(timestamp - 900),
      x = !!var,
      dx = abs(x - dplyr::lag(x)),
      p = !!p_rain
    )
  
  # Tests from Taylor & Loescher 2013 (summary of NEON's QC strategy)
  
  # Notes
  # - order in paper: 1) range, 2) sigma, 3) delta, 4) step, 5) null
  # - general principle is thresholds are determined by mean +/- 2sd
  # - however this implementation is slightly different in that it aims to
  # -  detect aberrant sensor behavior, not determine empirical threshold
  # - thus thresholds are computed and set at short (day-month) time windows 
  # - 
  # - they recommend visual inspection after flagging to remove false positives
  # - instead, I added secondary checks on tests prone to false positives
  # - double-difference checks step & range tests
  # - acf checks sigma test
  
  # Operates on points within days: range, delta, step
  pd <- tbl %>% 
    dplyr::select(x, dx, date, timestamp) %>% 
    tidyr::nest(data = c(timestamp, x, dx)) %>% 
    # Compute statistics
    dplyr::mutate(
      # Range test
      x = purrr::map(data, dplyr::pull, x), 
      x_mean = purrr::map_dbl(x, base::mean, na.rm = TRUE), 
      x_sigma = purrr::map_dbl(x, stats::sd, na.rm = TRUE), 
      x_lo = x_mean - 2 * x_sigma, 
      x_hi = x_mean + 2 * x_sigma,
      
      # Delta & step test
      dx = purrr::map(data, dplyr::pull, dx),
      dx_mean = purrr::map_dbl(dx, base::mean, na.rm = TRUE), 
      dx_sigma = purrr::map_dbl(dx, stats::sd, na.rm = TRUE), 
      dx_lo = dx_mean - 2 * dx_sigma, 
      dx_hi = dx_mean + 2 * dx_sigma,
      
      abs_mean = purrr::map_dbl(x, ~ mean(abs(.), na.rm = TRUE))
    ) %>% 
    dplyr::select(-x, -dx) %>% 
    tidyr::unnest(data) %>% 
    dplyr::arrange(timestamp) %>%
    # Flag point values
    dplyr::mutate(
      range = dplyr::if_else(x < x_lo | x > x_hi, 2L, 0L),
      # Secondary check on range test - not worth flagging
      spread = abs(x_hi - x_lo) / abs_mean,
      # - downgrade flag if range of limits is less than 0.5x mean
      # - some vars generally don't vary much relative to mean
      range = dplyr::if_else(spread < 0.5, pmax(range - 1L, 0L), range),
      
      # Secondary check on delta/step tests - not worth flagging
      # - set absolute limits on low/high based on mean +/- 1sd of all data
      dx_lo = pmin(dx_lo, mean(dx, na.rm = TRUE) - sd(dx, na.rm = TRUE)),
      dx_hi = pmax(dx_hi, mean(dx, na.rm = TRUE) + sd(dx, na.rm = TRUE)),
      
      delta = dplyr::if_else(dx < dx_lo, 2L, 0L),
      step = dplyr::if_else(dx > dx_hi, 2L, 0L),
      
      # Secondary check on step test - double difference (diff2)
      dx2 = diff2(x),
      # - if point is "in line" with adjacent points then downgrade flag
      # - use dx high for dx2 - marginally different from calculating dx2 lims
      step = dplyr::if_else(dx2 < dx_hi, pmax(step - 1L, 0L), step)
    ) %>%
    dplyr::select(range, delta, step)
  
  # Cross-check SWC & TS step flags against P_RAIN
  # - large jumps are reasonable if it rained
  if (var_name %in% c("swc", "ts")) {
    pd <- pd %>%
      dplyr::bind_cols(tbl) %>%
      dplyr::mutate_at(
        dplyr::vars(range, step), 
        ~ dplyr::if_else(around(p, ~ .x > 0), pmax(. - 1L, 0L), .)
      ) %>%
      dplyr::select(range, delta, step)
  }
  
  # Operates on days within months: sigma, spike
  dm <- tbl %>% 
    dplyr::select(x, dx, month, date, timestamp) %>% 
    tidyr::nest(data = c(timestamp, x, dx)) %>% 
    # Compute statistics
    dplyr::mutate(
      x = purrr::map(data, dplyr::pull, x),
      x_var = purrr::map_dbl(x, stats::sd, na.rm = TRUE),
      
      dx = purrr::map(data, dplyr::pull, dx),
      dx_var = purrr::map_dbl(dx, stats::sd, na.rm = TRUE),
      
      acf = purrr::map_dbl(x, acf_1)
    ) %>%
    dplyr::select(-x, -dx) %>%
    tidyr::unnest(data) %>%
    tidyr::nest(data = c(date, timestamp, x, dx, x_var, dx_var, acf)) %>%
    dplyr::mutate(
      # Sigma test
      x_var = purrr::map(data, dplyr::pull, x_var), 
      x_mean = purrr::map_dbl(x_var, base::mean, na.rm = TRUE), 
      x_sigma = purrr::map_dbl(x_var, stats::sd, na.rm = TRUE), 
      v_lo = x_mean - 2 * x_sigma, 
      v_hi = x_mean + 2 * x_sigma,
      
      # Spike test
      dx_var = purrr::map(data, dplyr::pull, dx_var), 
      dx_mean = purrr::map_dbl(dx_var, base::mean, na.rm = TRUE), 
      dx_sigma = purrr::map_dbl(dx_var, stats::sd, na.rm = TRUE), 
      k_lo = dx_mean - 2 * dx_sigma, 
      k_hi = dx_mean + 2 * dx_sigma
    ) %>% 
    dplyr::select(-x_var, -dx_var) %>% 
    tidyr::unnest(data) %>% 
    dplyr::mutate(
      # Secondary check on sigma/spike tests - not worth flagging
      # - set absolute limit on low based on mean - 1sd of all data
      v_lo = pmin(v_lo, mean(x_var, na.rm = TRUE) - sd(x_var, na.rm = TRUE)),
      v_hi = pmax(v_hi, mean(x_var, na.rm = TRUE) - sd(x_var, na.rm = TRUE)),
      k_lo = pmin(k_lo, mean(dx_var, na.rm = TRUE) - sd(dx_var, na.rm = TRUE)),
      k_hi = pmax(k_hi, mean(dx_var, na.rm = TRUE) - sd(dx_var, na.rm = TRUE)),
      
      v_flag = dplyr::if_else(x_var < v_lo | x_var > v_hi, 2L, 0L),
      
      k_flag = dplyr::if_else(dx_var < k_lo | dx_var > k_hi, 2L, 0L),
      
      # Secondary check on sigma/spike tests - autocorrelation
      # - if ACF is high enough then downgrade flag
      v_flag = dplyr::if_else(acf >= 0.90, pmax(v_flag - 1L, 0L), v_flag),
      k_flag = dplyr::if_else(acf >= 0.90, pmax(k_flag - 1L, 0L), k_flag)
    ) %>%
    dplyr::select(sigma = v_flag, spike = k_flag)
  
  # Cross-check SWC & TS step flags against P_RAIN
  # - large jumps are reasonable if it rained
  if (var_name %in% c("swc", "ts")) {
    dm <- dm %>%
      dplyr::bind_cols(dplyr::select(tbl, date, p)) %>%
      dplyr::group_by(date) %>%
      dplyr::mutate(p = sum(p, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_at(
        dplyr::vars(sigma, spike), 
        ~ dplyr::if_else(p > 0.1, pmax(. - 1L, 0L), .)
      ) %>%
      dplyr::select(sigma, spike)
  }
  
  
  # Operates on days within years: null
  # - no secondary checks since this already indicates underlying issues
  dy <- tbl %>% 
    dplyr::select(x, year, date, timestamp) %>% 
    tidyr::nest(data = c(timestamp, x)) %>% 
    dplyr::mutate(
      x = purrr::map(data, dplyr::pull, x),
      # Number of missing values
      n = purrr::map_dbl(x, ~dplyr::na_if(length(na.omit(.)), 0))
    ) %>%
    dplyr::select(-x) %>%
    tidyr::unnest(data) %>%
    tidyr::nest(data = c(date, timestamp, x, n)) %>%
    dplyr::mutate(
      # Null test
      n = purrr::map(data, dplyr::pull, n), 
      n_mean = purrr::map_dbl(n, base::mean, na.rm = TRUE), 
      n_sigma = purrr::map_dbl(n, stats::sd, na.rm = TRUE), 
      n_lo = n_mean - 2 * n_sigma
    ) %>% 
    dplyr::select(-n) %>% 
    tidyr::unnest(data) %>% 
    dplyr::mutate(null = dplyr::if_else(n < n_lo, 2L, 0L)) %>%
    dplyr::select(null)
  
  cat("done. ", sep = "")
  
  # Output tbl of all flags
  dplyr::bind_cols(pd, dm, dy) %>% dplyr::mutate_all(tidyr::replace_na, 0L)
}

flag_var_diffs <- function(data, var1, var2, qc_data, alpha = 0.001) {
  
  var1_name <- rlang::ensym(var1) %>% rlang::as_string()
  var1 <- rlang::enquo(var1)
  var2_name <- rlang::ensym(var2) %>% rlang::as_string()
  var2 <- rlang::enquo(var2)
  
  # Gather var1 flags
  var1_flags <- qc_data %>% 
    dplyr::select_at(
      dplyr::vars(tidyselect::contains(stringr::str_c("_", var1_name, "_")))
    ) %>%
    dplyr::mutate_all(~dplyr::if_else(. == 1, 0L, .))
  if (!stringr::str_detect(var1_name, "_ep")) {
    var1_flags <- dplyr::select_at(
      var1_flags, dplyr::vars(-tidyselect::contains("_ep_"))
    )
  }
  qc_var1 <- combine_flags(var1_flags)
  
  # Gather var2 flags
  var2_flags <- qc_data %>% 
    dplyr::select_at(
      dplyr::vars(tidyselect::contains(stringr::str_c("_", var2_name, "_")))
    ) %>%
    dplyr::mutate_all(~dplyr::if_else(. == 1, 0L, .))
  if (!stringr::str_detect(var2_name, "_ep")) {
    var2_flags <- dplyr::select_at(
      var2_flags, dplyr::vars(-tidyselect::contains("_ep_"))
    )
  }
  qc_var2 <- combine_flags(var2_flags)
  
  x <- dplyr::pull(data, !!var1)
  y <- dplyr::pull(data, !!var2)
  
  # Flag differences between vars
  flag <- flag_mahalanobis(x, y, alpha = alpha)
  
  # Var is not given diff flag if other var is already flagged
  # - indicates problem is elsewhere
  tbl <- tibble::tibble(
    flag1 = dplyr::if_else(qc_var2 == 2, 0L, as.integer(flag)),
    flag2 = dplyr::if_else(qc_var1 == 2, 0L, as.integer(flag))
  )
  
  names(tbl) <- stringr::str_c("qc_", c(var1_name, var2_name), "_vardiff")
  tbl
}

flag_rad <- function(x, pot) {
  
  # Via OneFlux:
  # https://github.com/LI-COR/ONEFlux/blob/master/oneflux_steps/qc_auto/src
  
  out <- rep(0L, length(x))
  value <- x - pot
  
  out <- dplyr::case_when(
    value > 0 & pot == 0 & x > 50 ~ 2L,
    value > 50 & pot > 200 ~ 2L,
    pot > 0 & (value / pot > 0.15) ~ 2L,
    TRUE ~ 0L
  )
  
  out
}

plot_flags <- function(data, var, qc_data, geom = c("point", "line")) {
  
  geom <- rlang::arg_match(geom)
  
  var_name <- rlang::ensym(var) %>% 
    rlang::as_string() %>% 
    stringr::str_c("_", ., "_")
  var <- rlang::enquo(var)
  
  # If qc_data is not provided, assume it is a column in data named "qc_var"
  if (!missing(qc_data)) {
    qc_data <- rlang::enquo(qc_data)
  } else {
    qc_var <- var_name %>% 
      stringr::str_c("qc", .) %>% 
      stringr::str_remove("_$") %>% rlang::sym()
    qc_data <- tibble::tibble(!!qc_var := dplyr::pull(data, !!qc_var)) %>%
      dplyr::rename_all(~ stringr::str_c(., "_all"))
  }
  
  
  # If qc_data is not a data frame, assume it is a column in data
  if (!inherits(rlang::eval_tidy(qc_data), "data.frame")) {
    qc_data <- tibble::tibble(!!qc_data := dplyr::pull(data, !!qc_data)) %>%
      dplyr::rename_all(~ stringr::str_c(., "_all"))
  }
  
  plot_data <- data %>%
    dplyr::select(timestamp, !!var) %>%
    dplyr::bind_cols(dplyr::select_at(
      rlang::eval_tidy(qc_data), dplyr::vars(dplyr::contains(var_name))
    ))
  
  if (!stringr::str_detect(var_name, "_ep")) {
    plot_data <- dplyr::select_at(
      plot_data, dplyr::vars(-dplyr::contains("_ep"))
    )
  }
  
  plot_data <- plot_data %>%
    coalesce_flags(-timestamp, -!!var) %>%
    dplyr::filter(flag != "none") %>%
    # Nicer flag names
    dplyr::mutate(flag = stringr::str_remove_all(flag, "^in_|^out_|^ep_"))
  
  plot <- data %>% ggplot2::ggplot(ggplot2::aes(timestamp, !!var))
  
  if (geom == "point") {
    plot <- plot + ggplot2::geom_point(alpha = 0.8, stroke = 0, na.rm = TRUE)
  } else if (geom == "line") {
    plot <- plot + ggplot2::geom_line(alpha = 0.8, na.rm = TRUE) 
  }
  
  plot +
    ggplot2::geom_point(
      data = plot_data, ggplot2::aes(color = flag), size = 0.75, na.rm = TRUE
    )  +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_bw()
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)
path_in <- file.path(wd, "processing_data", "03_combine_biomet", "output")

# Input file - biomet output with QC flags
biomet_input <- latest_version(path_in, "biomet_combined")
# Input file - biomet data from nearby sites
aux_input <- latest_version(
  stringr::str_replace(path_in, settings$site, md$closest_site[1]),
  "biomet_combined"
)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing_data", "04_biomet_qc", "output")

# Set variables to be flagged
vars <- rlang::exprs(
  ta, ta_ep, pa_ep, rh, ppfd_in, sw_in, sw_out, lw_in, lw_out, g, swc, ts
)

# List of vars with reps
rep_vars <- list(
  g = expr(list(g_1_1_1, g_2_1_1, g_3_1_1)),
  swc = expr(list(swc_1_1_1, swc_2_1_1, swc_3_1_1)),
  ts = expr(list(ts_1_1_1, ts_2_1_1, ts_3_1_1))
)


### Biomet data initialization =================================================

# Load the Biomet file
biomet <- read.csv(biomet_input, stringsAsFactors = FALSE)

# Add timestamp components
biomet <- biomet %>% 
  mutate(timestamp = ymd_hms(timestamp, tz = md$tz_name)) %>%
  add_time_comps()

# Load Biomet file from closest site
biomet_aux <- read.csv(aux_input, stringsAsFactors = FALSE)
biomet_aux <- mutate(
  biomet_aux, 
  # Parse timestamp
  timestamp = ymd_hms(timestamp, tz = md$tz_name),
  # Remove records during system failure
  p_rain = clean(p_rain, flag_biomet_system(biomet_aux))
)

# Initialize QC data frame
qc_biomet <- select(biomet, timestamp)

# Warn if both ta or rh vars are the same
#biomet %>% filter(ta == ta_ep) %>% summarize(n())
#biomet %>% filter(rh == rh_ep & rh != 100) %>% summarize(n())


### Automatic statistical flags ================================================

# Compute all auto flags
auto_flags <- vars %>% 
  map(~ biomet_auto_flags(biomet, !!., p_rain)) %>%
  set_names(map_chr(vars, rlang::as_string))

# Validate flags for analagous sensors
# - can do this for sw_in vs sw_out because same unit NOT same sensor
auto_flags <- auto_flags %>%
  list_modify(sw_in2 = pluck(auto_flags, "sw_in")) %>%
  validate_flags(ta, ta_ep) %>%
  validate_flags(ppfd_in, sw_in) %>%
  validate_flags(sw_out, sw_in2) %>%
  list_modify(sw_in2 = NULL)

# Give flags QC names
for (i in seq_along(vars)) {
  auto_flags[[i]] <- rename_all(
    auto_flags[[i]], ~ str_c("qc_", names(auto_flags)[i], "_", .)
  )
}

# Plot all flags for each var
auto_plots <- map2(vars, auto_flags, ~ plot_flags(biomet, !!.x, .y))
#map(auto_plots, ~ .)

# Add combined flag to QC dataset
for (i in seq_along(vars)) {
  var_qc_name <- str_c("qc_", rlang::as_string(vars[[i]]), "_auto")
  qc_biomet[, var_qc_name] <- combine_flags(auto_flags[[i]], clean_value = 1)
}


### Flag theoretically implausible values ======================================

# Add p_rain to vars list
all_vars <- append(vars, list(expr(p_rain)))

for (i in seq_along(all_vars)) {
  var_qc_name <- str_c("qc_", rlang::as_string(all_vars[[i]]), "_plaus")
  qc_biomet[, var_qc_name] <- drop_attributes(apply_thr(
    pull(biomet, !!all_vars[[i]]),
    pluck(var_attrs, rlang::as_string(all_vars[[i]]), "limits"),
    flag = "outside"
  ))
}


### Special potential incoming radiation flags =================================

qc_biomet <- mutate(
  qc_biomet,
  # Check SW_IN against potential radiation
  qc_sw_in_pot = flag_rad(biomet$sw_in, biomet$sw_in_pot),
  # Check converted PPFD_IN against potential radiation
  qc_ppfd_in_pot = flag_rad(biomet$sw_in_ppfd, biomet$sw_in_pot)
)


### Special precipitation flags ================================================

p_rain_rh_lim <- biomet %>% 
  filter(p_rain > 0, !is.na(rh)) %>% 
  summarize(mean(rh) - 2 * sd(rh)) %>%
  pluck(1, 1) %>% signif(1)
p_rain_kt_lim <- biomet %>% 
  filter(p_rain > 0, !is.na(kt)) %>% 
  summarize(mean(kt) + 2 * sd(kt)) %>%
  pluck(1, 1) %>% signif(1)

qc_biomet <- mutate(
  qc_biomet,
  # Unlikely precipitation - soft flags
  # - non-rain events are not flagged because no plausible explanation for that
  # - Estevez et al. 2011 thresholds: rh = 80, kt = 0.5
  # - allowing empirical thresholds here due to site differences
  qc_p_rain_rh = if_else(biomet$p_rain > 0 & biomet$rh < p_rain_rh_lim, 1L, 0L),
  qc_p_rain_kt = if_else(biomet$p_rain > 0 & biomet$kt > p_rain_kt_lim, 1L, 0L),
  qc_p_rain_lik = combine_flags(qc_p_rain_rh, qc_p_rain_kt),
  # Did nearby site experience similar conditions? (+/- 1 lag)  
  qc_p_rain_aux = if_else(magrittr::and(
    biomet$p_rain > 0, around(biomet_aux$p_rain, ~ .x == 0, .and = TRUE)
  ), 1L, 0L),
  # Hard flag if unlikely precipitation AND closest site had no rain
  qc_p_rain_err = if_else(qc_p_rain_lik + qc_p_rain_aux > 1, 2L, 0L) %>%
    replace_na(0L)
)


### Difference between analagous measurement flags =============================

# TA, RH, SW_IN
qc_biomet <- bind_cols(
  qc_biomet,
  flag_var_diffs(biomet, ta, ta_ep, qc_biomet),
  # Only consider daytime differences for incoming rad
  flag_var_diffs(
    mutate_at(biomet, vars(sw_in, ppfd_in), ~ if_else(night_pot, NA_real_, .)), 
    sw_in, ppfd_in, qc_biomet
  )
)


### Combine flags and check results ============================================

for (i in seq_along(all_vars)) {
  var_name <- rlang::as_string(all_vars[[i]])
  var_flags <- qc_biomet %>% 
    select_at(vars(contains(str_c("_", var_name, "_")))) %>%
    mutate_all(~ if_else(. == 1, 0L, .))
  if (!str_detect(var_name, "_ep")) {
    var_flags <- select_at(var_flags, vars(-contains("_ep_")))
  }
  
  # Combine flags, flag isolated points
  biomet[, str_c("qc_", var_name)] <- flag_around(combine_flags(var_flags), 1)
}

# Plot of all flags for each var
flag_plots <- all_vars %>%
  map(~ plot_flags(biomet, !!., qc_biomet, geom = "line")) %>%
  set_names(str_c(map_chr(all_vars, rlang::as_label), "_flags"))
#map(flag_plots, ~ .)

# How many flagged?
biomet %>%
  select_at(vars(all_of(str_c("qc_", map_chr(all_vars, rlang::as_string))))) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  filter(value == 2) %>%
  summarize(n_flag = n())

# Plot overall flag for each var
qc_plots <- all_vars %>%
  map(~ plot_flags(biomet, !!., geom = "line")) %>%
  set_names(str_c(map_chr(all_vars, rlang::as_label), "_qc"))
#map(qc_plots, ~ .)

# Combine all plot lists
plots <- append(flag_plots, qc_plots)


### Save output ================================================================

# Save processed Biomet dataset
biomet_out <- file.path(path_out, str_c("biomet_qc_", tag_out, ".csv"))
write.csv(biomet, biomet_out, row.names = FALSE)

# Create documentation for processed Biomet output
biomet_docu <- append(
  settings, 
  list(
    files = c(biomet_input, aux_input),
    aux_site = md$closest_site[1],
    plaus_lims = modify(var_attrs, "limits"),
    p_rain_lik = list(kt_lim = p_rain_kt_lim, rh_lim = p_rain_rh_lim)
  )
)
biomet_docu_out <- str_replace(biomet_out, ".csv", ".txt")
# Save documentation
sink(biomet_docu_out)
biomet_docu
sink()

# Save Biomet QC flags datasets
auto_biomet_out <- str_replace(biomet_out, "_qc_", "_qc_auto_flags_")
write.csv(bind_cols(auto_flags), auto_biomet_out, row.names = FALSE)
qc_biomet_out <- str_replace(biomet_out, "_qc_", "_qc_flags_")
write.csv(qc_biomet, qc_biomet_out, row.names = FALSE)

# Save one pdf document with all diagnostic/summary plots
plot_path <- file.path(path_out, paste0("biomet_qc_plots_", tag_out, ".pdf"))
pdf(plot_path)
map(plots, ~ .)
dev.off()

# Finished
