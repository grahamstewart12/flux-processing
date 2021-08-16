### ============================================================================
# Correct EddyPro output data ==================================================
### ============================================================================

# Purpose: Combine multiple eddypro output and/or biomet files, write the result
# to file. EddyPro is not always able to produce a single output file for a
# given time period due to changes in system configuration, etc. This script
# enables the user to resolve this issue by merging multiple output files (of
# potentially varying formats) with appropriate documentation.

# Input(s): Biomet file(s) from EddyPro, any additional data from sensors not
# associated with flux tower (e.g. PTs)

# Output(s): Combined Biomet file with any additional data added, documentation

start_time <- Sys.time()

# Load the required packages
devtools::load_all("~/Projects/flux/dscalr")
library(lubridate)
library(tidyverse)

# Load reference files
# source("reference/site_metadata.R")


### Helper functions ===========================================================

# - mostly wrappers to simplify computation

potential_radiation <- function(timestamp, md) {
  
  yday <- lubridate::yday(timestamp)
  hour <- decimal_hour(timestamp)
  
  bigleaf::potential.radiation(
    yday, hour, md$tower$coords$lat, md$tower$coords$lon, md$tz$utc_offset
  )
}

sun_position <- function(timestamp, md) {
  
  # solartime function gets tz from timestamp, so need to encode attribute
  timestamp <- lubridate::force_tz(timestamp, md$tz$name)
  
  sun_pos <- solartime::computeSunPosition(
    timestamp, md$tower$coords$lat, md$tower$coords$lon
  )
  
  sun_pos %>% 
    tibble::as_tibble() %>% 
    dplyr::pull(elevation)
}

apply_offset <- function(x, offset) {
  
  # Do this safely - don't accidentally apply offsets twice
  if (!is.null(attr(x, "offset"))) {
    stop("Offset has already been applied.")
  } 
  
  out <- x + offset
  attr(out, "offset") <- offset
  out
}

sun_time <- function(sol_ang, shoulder_n) {
  
  dplyr::case_when(
    magrittr::and(
      sol_ang > 0, before(sol_ang, ~ .x < 0, .n = shoulder_n)
    ) ~ "rise",
    magrittr::and(
      sol_ang < 0, before(sol_ang, ~ .x > 0, .n = shoulder_n)
    ) ~ "set",
    sol_ang > 0 ~ "day",
    sol_ang < 0 ~ "night"
  )
}

clearness_index <- function(sw_in, sw_in_pot, night_pot) {
  
  kt <- pmin(pmax(sw_in, 0) / sw_in_pot, 1)
  dplyr::if_else(night_pot, NA_real_, kt)
}

flag_biomet_system <- function(data, vars = c(ta_bm, lw_in, lw_out)) {
  
  # Suggested vars chosen using the criteria:
  # 1) not directly dependent, 2) unlikely to have natural runs, 3) not repped
  vars <- rlang::enquo(vars)
  
  tbl <- dplyr::select(data, !!vars)
  
  tbl %>%
    dplyr::mutate(dplyr::across(.fns = flag_repeats)) %>%
    dplyr::mutate(dplyr::across(.fns = ~ dplyr::na_if(.x, 0))) %>%
    as.list() %>%
    purrr::pmap_int(~ sum(.)) %>%
    tidyr::replace_na(0L)
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- yaml::read_yaml(file.path("data", settings$site, "metadata.yml"))

# Set tag for creating output file names
tag <- make_tag(settings$site, settings$year)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("data", settings$site, settings$year)
path_in <- file.path(wd, "01_combine_eddypro")

# Input files
file_data <- file.path(
  wd, "01_combine_eddypro", "eddypro", paste0("eddypro_combined_", tag, ".csv")
)
# file_data <- latest_version(file.path(path_in, "eddypro"))
file_settings <- file.path(
  wd, "01_combine_eddypro", "settings", 
  paste0("processing_combined_", tag, ".csv")
)
# file_settings <- latest_version(file.path(path_in, "settings"))

# Set path for output files
path_out <- file.path(wd, "02_correct_eddypro")

# List of vars housed in the biomet system
biomet_vars <- rlang::exprs(
  ta_bm, rh_bm, ppfd_in, sw_in, sw_out, lw_in, lw_out, p_rain,
  matches(c("^ts_|^swc_|^g_"))
)


### Load & clean up input data =================================================

cat("Importing data files...")

# Read in data
data <- read_csv(
  file_data, guess_max = 7000, col_types = cols(), progress = FALSE
)

# Fix names
data <- data %>%
  # Simplify non-replicated variable names (remove "_1_1_1")
  rename_with(
    ~ str_replace(., "_\\d_1_1", ""),
    c(matches("_\\d_1_1"), -matches(c("^ts_|^swc_|^g_"))),
  ) %>%
  rename(
    # Only one air pressure var, so remove suffix
    pa = pa_ep,
    # Give biomet ta & rh vars a suffix 
    ta_bm = ta,
    rh_bm = rh
  )

# Flag/clean system errors indicated by simultaneous runs in biomet vars
data <- data %>% 
  mutate(qc_biomet_all = flag_biomet_system(.)) %>%
  mutate(across(c(!!!biomet_vars), ~ clean(.x, qc_biomet_all)))

cat("done.\n")


### Correct names & units ======================================================

cat("Applying names & units corrections...")

# TODO get rid of this once I make sure everything uses the fluxnet output
# - don't need to correct CH4 units since those vars are taken from full_output

# Check p_rain before converting
if (max(data$p_rain, na.rm = TRUE) < 1) {
  data <- mutate(data, p_rain = p_rain * 1000)
}

# Automatic correction of temperatures and percentages
data <- data %>%
  mutate(across(starts_with("swc_"), ~ if_else(.x > 1, .x / 100, .x))) %>%
  mutate(across(
    c(starts_with("ts_"), ta_ep, ta_bm), ~ if_else(.x > 150, .x - 273.15, .x)
  ))

# CH4 variables from nmol to umol
if (median(data$fch4, na.rm = TRUE) > 1) {
  data <- mutate(data, across(
    c(
      fch4, fch4_randunc_hf, sch4_single, ch4_molar_density, ch4_mixing_ratio, 
      ch4, ch4_meas_sigma
    ), ~ .x / 1000
  ))
}

cat("done.\n")


### Correct wind direction for magnetic declination ============================

mag_decl <- NA

# Read in EddyPro settings
ep_settings <- read_csv(file_settings, col_types = cols())

# Check if geographic north was already used before correcting
if (sum(ep_settings$use_geo_north) == 0) {
  
  # Calculate magnetic declination given location & time
  mag_decl <- pluck(
    oce::magneticField(
      md$tower$coords$lon, md$tower$coords$lat, settings$year + 0.5
    ), 
    "declination"
  )
  
  # Reassign magnetic wind direction, apply declination
  data <- mutate(
    data,
    wd_mag = wd,
    wd = (wd_mag - mag_decl) %% 360
  )
}


### Combined replicated variables ==============================================

cat("Combining replicated variables...")

data <- data %>%
  rowwise() %>%
  mutate(
    g = mean(c(g_1_1_1, g_2_1_1, g_3_1_1)),
    swc = mean(c(swc_1_1_1, swc_2_1_1, swc_3_1_1)),
    ts = mean(c(ts_1_1_1, ts_2_1_1, ts_3_1_1))
  ) %>%
  ungroup()

cat("done.\n")


### Correct radiation values for zero-offset ===================================

cat("Correcting radiation for zero-offset...")

# Radiation offset
offset <- data %>%
  # Calculate sw_in_pot and sol_ang at end of averaging period
  mutate(
    sw_in_pot = potential_radiation(timestamp, md),
    sol_ang = sun_position(timestamp, md),
    # Stricter sun time to exclude any possible period with light
    sun_time = sun_time(sol_ang, shoulder_n = 3)
  ) %>%
  filter(sw_in_pot == 0, sun_time == "night") %>%
  summarize(across(
    c(ppfd_in, sw_in, sw_out), mean, na.rm = TRUE)
  ) %>%
  as.list()

data <- mutate(
  data, 
  ppfd_in = apply_offset(ppfd_in, -offset$ppfd_in),
  sw_in = apply_offset(sw_in, -offset$sw_in),
  sw_out = apply_offset(sw_out, -offset$sw_out)
)

cat("done.\n")


### Add auxilliary data to main data frame =====================================

cat("Calculating additional variables...")

# Other variables that can be calculated using core data

# Get coefficient of ppfd_in/sw_in relationship
frac_ppfd <- data %>% 
  select(sw_in, ppfd_in) %>% 
  drop_na() %>% 
  summarize(sum(sw_in) / sum(ppfd_in)) %>% 
  pluck(1)

data <- mutate(
  data,
  # Calculate sw_in equivalent of ppfd_in
  sw_in_ppfd = ppfd_in * frac_ppfd,
  # Calculate ppfd_in equivalent of sw_in
  ppfd_in_sw = sw_in / frac_ppfd,
  # Re-calculate potential radiation
  sw_in_pot = potential_radiation(timestamp - 900, md),
  night_pot = sw_in_pot < 20,
  # Solar angle and sun position
  sol_ang = sun_position(timestamp - 900, md),
  sun_time = sun_time(sol_ang, shoulder_n = 1),
  # Clearness index
  kt = clearness_index(sw_in, sw_in_pot, night_pot)
)

cat("done.\n")


### Save combined file with documentation ======================================

cat("Writing output...")

# Save the corrected data as a .csv file
data_out <- file.path(
  path_out, "data", paste0("eddypro_corrected_", tag, ".csv")
)
write_csv(data, data_out)

# Create documentation for data output
data_docu <- prepend(
  settings, list(
    files = file_data, mag_decl = mag_decl, frac_ppfd = frac_ppfd, 
    offsets = offset
  ), before = length(settings)
)
data_docu_out <- str_replace(data_out, ".csv", ".txt")

# Save documentation
sink(data_docu_out)
print(data_docu)
sink()

end_time <- Sys.time()
elapsed_time <- round(unclass(end_time - start_time), 1)

cat("done.\n")
cat(
  "Finished processing in ", elapsed_time, " ", attr(elapsed_time, "units"), 
  ".\n", sep = ""
)
cat("Output located in", path_out, "\n")

# Finished
