### ============================================================================
# Daily averaging of quality-controlled data ===================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/Desktop/RESEARCH/fluxtools")
library(lubridate)
library(tidyverse)

source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")
source("~/Desktop/DATA/Flux/tools/engine/daily_averages.R")

script_init(settings, site_metadata)


### Load required input data ===================================================

# EddyPro
site <- read.csv(site_input, stringsAsFactors = FALSE)

# Add timestamp components
site <- site %>% 
  dplyr::mutate(timestamp = lubridate::ymd_hms(timestamp, tz = md$tz_name)) %>%
  add_time_comps()


### Variables that must be first computed sub-daily ============================

# Indicator for thermal stratification/convection (re: Franz2016)
# - complicated at this site because PT is below sediment interface
site <- site %>%
  #filter(!is.na(tw), !is.na(tw_surf), tw > 0, tw_surf > 0) %>%
  transmute(
    timestamp = timestamp,
    # Air-free water density
    paf = map_dbl(tw + 273.15, possibly(IAPWS95::DTp, 0.101325), NA_real_),
    paf_surf = map_dbl(
      tw_surf + 273.15, possibly(IAPWS95::DTp, 0.101325), NA_real_
    ),
    # Air-saturated water density
    pas = paf - 0.004612 + 0.000106 * tw,
    pas_surf = paf_surf - 0.004612 + 0.000106 * tw_surf,
    # Water density gradient
    # - positive = stratification, negative = convective mixing
    wtd_abs = wtd_f - 0.2590146,
    wdg = (pas - pas_surf) / -wtd_abs,
    # Indicator for stratification - consider ignoring if wtd_f < 0.3
    stratified = if_else(wdg > 0, 1, 0)
  ) %>%
  select(timestamp, wdg, stratified) %>%
  left_join(site)

### Daily aggregation ==========================================================

# Re-create "un-gapfilled" variables - probably not necessary
site %>% select_at(vars(ends_with("_f"), ends_with("_fmeth"))) %>% names(.)
clean <- site %>%
  dplyr::mutate(
    sw_in = if_else(is.na(sw_in_fmeth), sw_in_f, NA_real_),
    ta = if_else(is.na(ta_fmeth), ta_f, NA_real_),
    vpd = if_else(is.na(vpd_fmeth), vpd_f, NA_real_),
    p_rain = if_else(is.na(p_rain_fmeth), p_rain_f, NA_real_),
    ts = if_else(is.na(ts_fmeth), ts_f, NA_real_),
    wtd = if_else(is.na(wtd_fmeth), wtd_f, NA_real_)
  ) %>%
  dplyr::select(-dplyr::ends_with("_fmeth"), -dplyr::ends_with("_fqc"))
  
# Clean non-gapfilled variables that have available flags
clean %>% select_at(vars(starts_with("qc_"))) %>% names(.)
clean <- clean %>% 
  clean_all() %>%
  mutate(
    night = if_else(sw_in_f < 20, 1, 0),
    fc = clean(fc, qc_fc_final),
    # Create an fc version that can be directly compared w/ fch4
    fc_paired = if_else(is.na(fch4), NA_real_, fc),
    fch4 = clean(fch4, qc_fch4_final)
  ) %>%
  dplyr::select(-dplyr::starts_with("qc_"))

# Generate reference list for how variables should be aggregated
agg_ref <- list(
  # Vars that should be summed daily (all others will be averaged)
  sum = c(
    "night", "sw_in", "ppfd_in", "lw_in", "lw_out", "netrad", "p_rain",
    "stratified"
  ),
  # Vars that should be split by day/night 
  diurnal = c(
    "\\bfc\\b", "fc_paired", "\\ble\\b", "\\bh\\b", "fch4", "ustar", "\\bwd\\b", 
    "x_peak_H00", "x_70perc_H00", "phi", "wdg",
    "\\bta", "\\brh\\b", "\\btw", "\\bdo", "\\bts"
  ),
  # Averaged vars for which statistics should be calculated
  stats = c("\\bfc\\b", "fc_paired", "\\ble\\b", "\\bh\\b", "fch4")
)

# Perform aggregation
# 1. Averages
daily <- clean %>%
  dplyr::select(
    -timestamp, -dplyr::matches(stringr::str_c(agg_ref$sum, collapse = "|"))
  ) %>%
  dplyr::group_by(date) %>%
  summarize_all(list(avg = mean), na.rm = TRUE)

# 2. Sums
daily <- clean %>%
  dplyr::select(
    date, dplyr::matches(stringr::str_c(agg_ref$sum, collapse = "|"))
  ) %>%
  dplyr::group_by(date) %>%
  # For sums it is good to know how many missing records
  summarize_all(list(sum = sum, n = ~ length(na.omit(.))), na.rm = TRUE) %>%
  # "Clean" night_sum now since it will be used for diurnal vars
  dplyr::mutate(night_sum = dplyr::na_if(night_sum, 48) %>% dplyr::na_if(0)) %>%
  dplyr::left_join(daily)

# 3. Day/night averages
daily <- clean %>%
  dplyr::select(
    date, night, 
    dplyr::matches(stringr::str_c(agg_ref$diurnal, collapse = "|")), 
    -ta_ep, -tstar
  ) %>%
  dplyr::group_by(date, night) %>%
  summarize_all(list(avg = mean), na.rm = TRUE) %>%
  dplyr::mutate(night = dplyr::recode(night, `0` = "day", `1` = "night")) %>%
  tidyr::pivot_wider(
    names_from = night, 
    values_from = c(dplyr::ends_with("sum"), dplyr::ends_with("avg"))
  ) %>%
  dplyr::select(-dplyr::ends_with("_NA")) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(daily)

# 4. Detailed stats
daily <- clean %>%
  dplyr::select(
    date, dplyr::matches(stringr::str_c(agg_ref$stats, collapse = "|"))
  ) %>%
  dplyr::group_by(date) %>%
  summarize_all(list(
    min = min, max = max, sd = sd, med = median, n = ~ length(na.omit(.))
  ), na.rm = TRUE) %>%
  dplyr::left_join(daily)

# 5. Detailed stats - day/night
daily <- clean %>%
  dplyr::select(
    date, night, dplyr::matches(stringr::str_c(agg_ref$stats, collapse = "|"))
  ) %>%
  dplyr::group_by(date, night) %>%
  summarize_all(list(
    min = min, max = max, sd = sd, med = median, n = ~ length(na.omit(.))
  ), na.rm = TRUE) %>%
  dplyr::mutate(night = dplyr::recode(night, `0` = "day", `1` = "night")) %>%
  tidyr::pivot_wider(
    names_from = night, values_from = c(-"date", -"night")
  ) %>%
  dplyr::select(-dplyr::ends_with("_NA")) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(daily)

# Convert to day-friendly units
daily <- daily %>%
  # C fluxes to mg m-2 h-1
  dplyr::mutate(dplyr::across(
      c(
        dplyr::starts_with("fch4_"), dplyr::starts_with("fc_"), 
        -dplyr::contains("_n_"), -dplyr::ends_with("_n")
      ), 
      ~ . * 12.011 * 3600 / 1000
  )) %>%
  # C fluxes to mg m-2 d-1 (extrapolate averages)
  dplyr::mutate(
    fch4 = fch4_avg * 24,
    fc_day = fc_avg_day * 24 * (1 - night_sum / night_n),
    fc_night = fc_avg_night * 24 * (night_sum / night_n),
    fc = fc_day + fc_night,
    fc_paired_day = fc_paired_avg_day * 24 * (1 - night_sum / night_n), 
    fc_paired_night = fc_paired_avg_night * 24 * (night_sum / night_n),
    fc_paired = fc_paired_day + fc_paired_night
  )

# Clean up final data set
daily <- dplyr::mutate(daily, dplyr::across(
  is.numeric, list(~ dplyr::ifelse(!is.finite(.), NA_real_, .))
))

# Check time series of vars
daily %>%
  ggplot(aes(date, fc)) +
  geom_point(aes(color = fc_n < 8))
daily %>%
  ggplot(aes(date, fch4)) +
  geom_point(aes(color = fch4_n < 8))

# Save dataset
daily_out <- paste0(
  wd, "/09_daily_averages/output/", "daily_", tag_out, ".csv"
)
write.csv(daily, daily_out, row.names = FALSE)

# Finished
