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

start_time <- Sys.time()

# Load the required packages
#suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
library(lubridate, warn.conflicts = FALSE)
library(tidyverse)

# Load reference files
source("~/Desktop/DATA/Flux/tools/reference/var_attributes.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")

# Load functions
source("~/Desktop/DATA/Flux/tools/engine/functions/dates_and_times.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/flag.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/latest_version.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/utilities.R")


### Helper functions ===========================================================

acf_1 <- function(x) {
  acf <- acf(x, plot = FALSE, na.action = na.pass)$acf
  acf %>% 
    purrr::array_branch() %>% 
    purrr::flatten_dbl() %>% 
    purrr::pluck(2)
}

tidy_acf <- purrr::compose(
  broom::tidy, ~ acf(.x, plot = FALSE, na.action = na.pass)
)

diff2 <- function(x, na.rm = FALSE, replace = 0) {
  
  d1 <- x - dplyr::lag(x, 1)
  d2 <- dplyr::lead(x, 1) - x
  
  # Don't remove NA if both d1 and d2 are NA
  if (na.rm) {
    # First look one further
    d1_out <- dplyr::if_else(is.na(d1) & !is.na(d2), x - dplyr::lag(x, 1), d1)
    d2_out <- dplyr::if_else(is.na(d2) & !is.na(d1), dplyr::lead(x, 2) - x, d2)
    
    # If still NA, fill with 'replace' value
    d1_out <- tidyr::replace_na(d1_out, replace[1])
    d2_out <- tidyr::replace_na(d2_out, replace[1])
  } else {
    d1_out <- d1
    d2_out <- d2
  }
  
  d1_out - d2_out
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
    dplyr::rename_with(
      ~ stringr::str_remove_all(.x, stringr::str_c("qc_", var1_name, "_"))
    )
  var2_flags <- data %>% 
    purrr::pluck(var2_name) %>%
    dplyr::rename_with(
      ~ stringr::str_remove_all(.x, stringr::str_c("qc_", var2_name, "_"))
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
  
  dt <- data %>%
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
  pd <- dt %>% 
    dplyr::select(x, dx, date, timestamp) %>%
    dplyr::group_by(date) %>%
    # Compute statistics
    dplyr::mutate(
      # Range test
      x_mean = mean(x, na.rm = TRUE), 
      x_sigma = sd(x, na.rm = TRUE),
      x_lo = x_mean - 2 * x_sigma, 
      x_hi = x_mean + 2 * x_sigma,
      
      # Delta & step test
      dx_mean = mean(dx, na.rm = TRUE), 
      dx_sigma = sd(dx, na.rm = TRUE),
      dx_lo = dx_mean - 2 * dx_sigma, 
      dx_hi = dx_mean + 2 * dx_sigma,
      
      abs_mean = mean(abs(x), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
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
  # - large jumps are reasonable if it rained (+/- 1 lag)
  if (var_name %in% c("swc", "ts")) {
    pd <- pd %>%
      dplyr::bind_cols(dt) %>%
      dplyr::mutate(dplyr::across(
        c(range, step), 
        ~ dplyr::if_else(around(p, ~ .x > 0), pmax(.x - 1L, 0L), .x)
      )) %>%
      dplyr::select(range, delta, step)
  }
  
  # Operates on days within months: sigma, spike
  dm <- dt %>% 
    dplyr::select(x, dx, month, date, timestamp) %>%
    dplyr::group_by(month, date) %>%
    # Compute statistics
    dplyr::summarize(
      x_var = sd(x, na.rm = TRUE),
      
      dx_var = sd(dx, na.rm = TRUE),
      
      acf = x %>% tidy_acf() %>% purrr::pluck("acf", 2)
    ) %>%
    dplyr::mutate(
      # Sigma test
      x_mean = mean(x_var, na.rm = TRUE), 
      x_sigma = sd(x_var, na.rm = TRUE), 
      v_lo = x_mean - 2 * x_sigma, 
      v_hi = x_mean + 2 * x_sigma,
      
      # Spike test
      dx_mean = mean(dx_var, na.rm = TRUE), 
      dx_sigma = sd(dx_var, na.rm = TRUE), 
      k_lo = dx_mean - 2 * dx_sigma, 
      k_hi = dx_mean + 2 * dx_sigma
    ) %>%
    dplyr::ungroup() %>%
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
    # Expand daily values to full time series
    dplyr::right_join(dplyr::select(dt, date), by = "date") %>%
    dplyr::select(sigma = v_flag, spike = k_flag)
  
  # Cross-check SWC & TS sigma, spike flags against P_RAIN
  # - large jumps are reasonable if it rained that day
  if (var_name %in% c("swc", "ts")) {
    dm <- dm %>%
      dplyr::bind_cols(dplyr::select(dt, date, p)) %>%
      dplyr::group_by(date) %>%
      dplyr::mutate(p = sum(p, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(
        c(sigma, spike), ~ dplyr::if_else(p > 0.1, pmax(.x - 1L, 0L), .x)
      )) %>%
      dplyr::select(sigma, spike)
  }
  
  
  # Operates on days within years: null
  # - no empirical checks since this already indicates underlying issues
  # - but disregard flag if missing data is <15%
  dy <- dt %>% 
    dplyr::select(x, year, date, timestamp) %>%
    dplyr::group_by(date) %>%
    # Number of non-missing values (exclude days that are entirely missing)
    dplyr::summarize(n = dplyr::na_if(length(na.omit(x)), 0)) %>%
    dplyr::mutate(
      # Null test
      n_mean = mean(n, na.rm = TRUE), 
      n_sigma = sd(n, na.rm = TRUE), 
      n_lo = n_mean - 2 * n_sigma,
      
      null = dplyr::if_else(n < n_lo, 2L, 0L),
      # Don't bother flagging if n missing is less than 8
      null = dplyr::if_else(n > 40, pmax(null - 1L, 0L), null)
    ) %>%
    # Expand daily values to full time series
    dplyr::right_join(dplyr::select(dt, date), by = "date") %>%
    dplyr::select(null)
  
  cat("done. ", sep = "")
  
  # Output tbl of all flags
  dplyr::bind_cols(pd, dm, dy) %>% 
    dplyr::mutate(dplyr::across(.fns = ~ tidyr::replace_na(.x, 0L)))
}

flag_var_diffs <- function(data, var1, var2, qc_data, alpha = 0.001) {
  
  var1_name <- rlang::ensym(var1) %>% rlang::as_string()
  var1 <- rlang::enquo(var1)
  var2_name <- rlang::ensym(var2) %>% rlang::as_string()
  var2 <- rlang::enquo(var2)
  
  # Gather var1 flags
  var1_flags <- qc_data %>% 
    dplyr::select(dplyr::contains(stringr::str_c("_", var1_name, "_"))) %>%
    dplyr::mutate(dplyr::across(.fns = ~ dplyr::if_else(.x == 1, 0L, .x)))
  if (!stringr::str_detect(var1_name, "_ep")) {
    var1_flags <- dplyr::select(var1_flags, -dplyr::contains("_ep_"))
  }
  qc_var1 <- combine_flags(var1_flags)
  
  # Gather var2 flags
  var2_flags <- qc_data %>% 
    dplyr::select(dplyr::contains(stringr::str_c("_", var2_name, "_"))) %>%
    dplyr::mutate(dplyr::across(.fns = ~ dplyr::if_else(.x == 1, 0L, .x)))
  if (!stringr::str_detect(var2_name, "_ep")) {
    var2_flags <- dplyr::select(var2_flags, -dplyr::contains("_ep_"))
  }
  qc_var2 <- combine_flags(var2_flags)
  
  x <- dplyr::pull(data, !!var1)
  y <- dplyr::pull(data, !!var2)
  
  # Flag differences between vars
  flag <- flag_distance(x, y, alpha = alpha)
  
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
      dplyr::rename_with(~ stringr::str_c(., "_all"))
  }
  
  
  # If qc_data is not a data frame, assume it is a column in data
  if (!inherits(rlang::eval_tidy(qc_data), "data.frame")) {
    qc_data <- tibble::tibble(!!qc_data := dplyr::pull(data, !!qc_data)) %>%
      dplyr::rename_with(~ stringr::str_c(., "_all"))
  }
  
  plot_data <- data %>%
    dplyr::select(timestamp, !!var) %>%
    dplyr::bind_cols(
      dplyr::select(rlang::eval_tidy(qc_data), dplyr::contains(var_name))
    )
  
  if (!stringr::str_detect(var_name, "_ep")) {
    plot_data <- dplyr::select(plot_data, -dplyr::contains("_ep"))
  }
  
  # Helper function for organizing flags
  coalesce_flags <- function(data, ..., prefix = "qc_[[:alnum:]]+\\_") {
    
    data %>%
      tibble::as_tibble() %>%
      tidyr::gather("flag", "value", ...) %>%
      dplyr::mutate(
        flag = stringr::str_remove(flag, prefix),
        value = dplyr::if_else(value > 1, flag, NA_character_)
      ) %>%
      tidyr::spread(flag, value) %>%
      dplyr::mutate(
        flag = dplyr::coalesce(!!! dplyr::select(., ...)),
        flag = tidyr::replace_na(flag, "none")
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
    plot <- plot + ggplot2::geom_line(alpha = 0.7, size = 0.5, na.rm = TRUE) 
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
path_in <- file.path(wd, "processing", "02_correct_eddypro", "data")

# Input file
data_input <- latest_version(path_in)
# Input file - data from nearby sites
aux_input <- latest_version(
  stringr::str_replace(path_in, settings$site, md$closest_site[1])
)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing", "03_biomet_qc")

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

cat("Importing data files...")

# Load the data
data <- readr::read_csv(
  data_input, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Add timestamp components
data <- add_time_comps(data)

# Load Biomet file from closest site
aux <- readr::read_csv(
  aux_input, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Initialize QC data frame
qc_biomet <- dplyr::select(data, timestamp)

# Warn if both ta or rh vars are the same
#data %>% filter(ta == ta_ep) %>% summarize(n())
#data %>% filter(rh == rh_ep & rh != 100) %>% summarize(n())

cat("done.\n")


### Automatic statistical flags ================================================

cat("Computing automatic sensor flags.\n")

# Compute all auto flags
auto_flags <- vars %>% 
  purrr::map(~ biomet_auto_flags(data, !!., p_rain)) %>%
  rlang::set_names(purrr::map_chr(vars, rlang::as_string))

# Validate flags for analagous sensors
# - can do this for sw_in vs sw_out because same unit NOT same sensor
auto_flags <- auto_flags %>%
  purrr::list_modify(sw_in2 = purrr::pluck(auto_flags, "sw_in")) %>%
  validate_flags(ta, ta_ep) %>%
  validate_flags(ppfd_in, sw_in) %>%
  validate_flags(sw_out, sw_in2) %>%
  purrr::list_modify(sw_in2 = NULL)

# Give flags QC names
for (i in seq_along(vars)) {
  auto_flags[[i]] <- dplyr::rename_with(
    auto_flags[[i]], ~ stringr::str_c("qc_", names(auto_flags)[i], "_", .)
  )
}

# Plot all flags for each var
auto_plots <- purrr::map2(vars, auto_flags, ~ plot_flags(data, !!.x, .y))

# Add combined flag to QC dataset
for (i in seq_along(vars)) {
  var_qc_name <- stringr::str_c("qc_", rlang::as_string(vars[[i]]), "_auto")
  qc_biomet[, var_qc_name] <- combine_flags(auto_flags[[i]], clean_value = 1)
}


### Flag theoretically implausible values ======================================

cat("\nFlagging unlikely conditions...")

# Add p_rain to vars list
all_vars <- append(vars, list(expr(p_rain)))

for (i in seq_along(all_vars)) {
  var_qc_name <- stringr::str_c(
    "qc_", rlang::as_string(all_vars[[i]]), "_plaus"
  )
  qc_biomet[, var_qc_name] <- flag_thr(
    dplyr::pull(data, !!all_vars[[i]]),
    purrr::pluck(var_attrs, rlang::as_string(all_vars[[i]]), "limits"),
    rule = "outside"
  )
}


### Special potential incoming radiation flags =================================

qc_biomet <- dplyr::mutate(
  qc_biomet,
  # Check SW_IN against potential radiation
  qc_sw_in_pot = flag_rad(data$sw_in, data$sw_in_pot),
  # Check converted PPFD_IN against potential radiation
  qc_ppfd_in_pot = flag_rad(data$sw_in_ppfd, data$sw_in_pot)
)


### Special precipitation flags ================================================

p_rain_rh_lim <- data %>% 
  dplyr::filter(p_rain > 0, !is.na(rh)) %>% 
  dplyr::summarize(mean(rh) - 2 * sd(rh)) %>%
  purrr::pluck(1, 1) %>% signif(1)
p_rain_kt_lim <- data %>% 
  dplyr::filter(p_rain > 0, !is.na(kt)) %>% 
  dplyr::summarize(mean(kt) + 2 * sd(kt)) %>%
  purrr::pluck(1, 1) %>% signif(1)

qc_biomet <- dplyr::mutate(
  qc_biomet,
  # Unlikely precipitation - soft flags
  # - non-rain events are not flagged because no plausible explanation for that
  # - Estevez et al. 2011 thresholds: rh = 80, kt = 0.5
  # - allowing empirical thresholds here due to site differences
  qc_p_rain_rh = dplyr::if_else(
    data$p_rain > 0 & data$rh < p_rain_rh_lim, 1L, 0L
  ),
  qc_p_rain_kt = dplyr::if_else(
    data$p_rain > 0 & data$kt > p_rain_kt_lim, 1L, 0L
  ),
  qc_p_rain_lik = combine_flags(qc_p_rain_rh, qc_p_rain_kt),
  # Did nearby site experience similar conditions? (+/- 2 lags)  
  qc_p_rain_aux = dplyr::if_else(magrittr::and(
    data$p_rain > 0, around(aux$p_rain, ~ .x == 0, .n = 2, .and = TRUE)
  ), 1L, 0L),
  # Hard flag if unlikely precipitation AND closest site had no rain
  qc_p_rain_err = dplyr::if_else(qc_p_rain_lik + qc_p_rain_aux > 1, 2L, 0L) %>%
    tidyr::replace_na(0L)
)

cat("done.\n")


### Difference between analagous measurement flags =============================

# TA, RH, SW_IN
qc_biomet <- dplyr::bind_cols(
  qc_biomet,
  flag_var_diffs(data, ta, ta_ep, qc_biomet),
  # Only consider daytime differences for incoming rad
  flag_var_diffs(
    dplyr::mutate(data, dplyr::across(
      c(sw_in, ppfd_in), ~ dplyr::if_else(night_pot, NA_real_, .x)
    )), 
    sw_in, ppfd_in, qc_biomet
  )
)


### Combine flags and check results ============================================

cat("Combining flags and plotting results...")

for (i in seq_along(all_vars)) {
  var_name <- rlang::as_string(all_vars[[i]])
  var_flags <- qc_biomet %>% 
    dplyr::select(dplyr::contains(stringr::str_c("_", var_name, "_"))) %>%
    dplyr::mutate(dplyr::across(.fns = ~ dplyr::if_else(.x == 1, 0L, .x)))
  if (!stringr::str_detect(var_name, "_ep")) {
    var_flags <- dplyr::select(var_flags, -dplyr::contains("_ep_"))
  }
  
  # Combine flags, flag isolated points
  qc_name <- stringr::str_c("qc_", var_name)
  data[, qc_name] <- flag_around(combine_flags(var_flags), 1)
}

# Plot of all flags for each var
flag_plots <- all_vars %>%
  purrr::map(~ plot_flags(data, !!., qc_biomet, geom = "line")) %>%
  rlang::set_names(
    stringr::str_c(purrr::map_chr(all_vars, rlang::as_label), "_flags")
  )

# How many flagged?
data %>%
  dplyr::select(dplyr::all_of(
    stringr::str_c("qc_", purrr::map_chr(all_vars, rlang::as_string))
  )) %>%
  tidyr::pivot_longer(dplyr::everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(value == 2) %>%
  dplyr::summarize(n_flag = dplyr::n())

# Plot overall flag for each var
qc_plots <- all_vars %>%
  purrr::map(~ plot_flags(data, !!., geom = "line")) %>%
  rlang::set_names(
    stringr::str_c(purrr::map_chr(all_vars, rlang::as_label), "_qc")
  )

# Combine all plot lists
plots <- append(flag_plots, qc_plots)

cat("done.\n")

### Save output ================================================================

cat("Writing output...")

# Save processed dataset
data_out <- file.path(
  path_out, "data", stringr::str_c("biomet_qc_", tag_out, ".csv")
)
readr::write_csv(remove_time_comps(data, -timestamp), data_out)

# Create documentation for processed output
data_docu <- append(
  settings, 
  list(
    files = c(data_input, aux_input),
    aux_site = purrr::pluck(md, "closest_site", 1),
    plaus_lims = purrr::modify(var_attrs, "limits"),
    p_rain_lik = list(kt_lim = p_rain_kt_lim, rh_lim = p_rain_rh_lim)
  )
)
data_docu_out <- stringr::str_replace(data_out, ".csv", ".txt")
# Save documentation
sink(data_docu_out)
data_docu
sink()

# Save Biomet QC flags datasets
auto_biomet_out <- file.path(
  path_out, "flags", stringr::str_c("biomet_qc_auto_flags_", tag_out, ".csv")
)
readr::write_csv(dplyr::bind_cols(auto_flags), auto_biomet_out)

qc_biomet_out <- file.path(
  path_out, "flags", stringr::str_c("biomet_qc_flags_", tag_out, ".csv")
)
readr::write_csv(qc_biomet, qc_biomet_out)

# Save one pdf document with all diagnostic/summary plots
plot_path <- file.path(
  path_out, "plots", paste0("biomet_qc_plots_", tag_out, ".pdf")
)
pdf(plot_path)
purrr::map(plots, ~ .)
dev.off()

end_time <- Sys.time()
elapsed_time <- round(unclass(end_time - start_time), 1)

cat("done.\n")
cat(
  "Finished processing in ", elapsed_time, " ", attr(elapsed_time, "units"), 
  ".\n", sep = ""
)
cat("Output located in", path_out, "\n")

# Finished
