### ============================================================================
# Combine Biomet files =========================================================
### ============================================================================

# Purpose: Combine multiple eddypro output and/or biomet files, write the result
# to file. EddyPro is not always able to produce a single output file for a
# given time period due to changes in system configuration, etc. This script
# enables the user to resolve this issue by merging multiple output files (of
# potentially varying formats) with appropriate documentation.

# Input(s): Biomet file(s) from EddyPro, any additional data from sensors not
# associated with flux tower (e.g. PTs)

# Output(s): Combined Biomet file with any additional data added, documentation

# Load the required packages
suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
library(lubridate)
library(tidyverse)

#source("~/Desktop/DATA/Flux/tools/engine/combine_biomet.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")


### Helper functions ===========================================================

# - mostly wrappers to simplify computation

potential_radiation <- function(timestamp, md) {
  
  yday <- lubridate::yday(timestamp)
  hour <- decimal_hour(timestamp)
  
  bigleaf::potential.radiation(yday, hour, md$lat, md$lon, md$tz)
}

sun_position <- function(timestamp, md) {
  
  sun_pos <- solartime::computeSunPosition(timestamp, md$lat, md$lon)
  
  sun_pos %>% 
    tibble::as_tibble() %>% 
    dplyr::pull(elevation)
}

apply_offset <- function(x, offset) {
  
  # Do this safely - don't accidentally apply offsets twice
  if (!is.null(attr(x, "offset"))) stop("Offset has already been applied.")
  
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

flag_biomet_system <- function(data, vars = c(ta, lw_in, lw_out)) {
  
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
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)
path_in <- file.path(wd, "processing", "01_combine_eddypro")

# Input file - biomet output with QC flags
biomet_input <- latest_version(path_in, "eddypro_combined")

# Input file - eddypro settings
ep_settings_input <- latest_version(path_in, "processing")
# [179] "use_geo_north=0"

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing", "02_correct_eddypro")

# List of vars housed in the biomet system
biomet_vars <- rlang::exprs(
  ta, rh, ppfd_in, sw_in, sw_out, lw_in, lw_out, p_rain,
  dplyr::matches(c("^ts_|^swc_|^g_"))
)


### Load & clean up the Biomet file ============================================

# Read in Biomet file
biomet <- readr::read_csv(
  biomet_input, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Parse POSIX timestamp
#biomet <- dplyr::mutate(
#  biomet, timestamp = lubridate::ymd_hms(timestamp, tz = md$tz_name)
#)

# Simplify non-replicated variable names (remove "_1_1_1")
biomet <- biomet %>%
  dplyr::rename_with(
    ~ stringr::str_replace(., "_\\d_1_1", ""),
    c(dplyr::matches("_\\d_1_1"), -dplyr::matches(c("^ts_|^swc_|^g_"))),
  )

# Flag/clean system errors indicated by simultaneous runs in biomet vars
biomet <- biomet %>% 
  dplyr::mutate(qc_biomet_all = flag_biomet_system(.)) %>%
  dplyr::mutate(dplyr::across(c(!!!biomet_vars), ~ clean(.x, qc_biomet_all)))


### Correct units ==============================================================

# Check p_rain before converting
if (max(biomet$p_rain, na.rm = TRUE) < 1) {
  biomet <- dplyr::mutate(biomet, p_rain = p_rain * 1000)
}

# Automatic correction of temperatures and percentages
biomet <- biomet %>%
  dplyr::mutate(dplyr::across(
    dplyr::starts_with("swc_"), ~ dplyr::if_else(.x > 1, .x / 100, .x)
  )) %>%
  dplyr::mutate(dplyr::across(
    c(dplyr::starts_with("ts_"), ta, ta_ep), 
    ~ dplyr::if_else(.x > 150, .x - 273.15, .x)
  ))

# CH4 variables from nmol to umol
if (median(biomet$fch4, na.rm = TRUE) > 1) {
  biomet <- dplyr::mutate(biomet, dplyr::across(
    c(
      fch4, fch4_randunc_hf, sch4_single, ch4_molar_density, ch4_mixing_ratio, 
      ch4, ch4_meas_sigma
    ), ~ .x / 1000
  ))
}


### Correct wind direction for magnetic declination ============================

mag_decl <- NA

# Read in EddyPro settings
ep_settings <- readr::read_csv(
  ep_settings_input, col_types = readr::cols(.default = readr::col_guess())
)

# Check if geographic north was already used before correcting
if (sum(ep_settings$use_geo_north) == 0) {
  
  # Calculate magnetic declination given location & time
  mag_decl <- purrr::pluck(
    oce::magneticField(md$lon, md$lat, settings$year + 0.5), "declination"
  )
  
  # Reassign magnetic wind direction, apply declination
  biomet <- dplyr::mutate(
    biomet,
    wd_mag = wd,
    wd = (wd_mag - mag_decl) %% 360
  )
}


### Combined replicated variables ==============================================

biomet <- biomet %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    g = mean(c(g_1_1_1, g_2_1_1, g_3_1_1), na.rm = TRUE),
    swc = mean(c(swc_1_1_1, swc_2_1_1, swc_3_1_1), na.rm = TRUE),
    ts = mean(c(ts_1_1_1, ts_2_1_1, ts_3_1_1), na.rm = TRUE)
  ) %>%
  dplyr::ungroup()


### Correct radiation values for zero-offset ===================================

# Radiation offset
offset <- biomet %>%
  # Calculate sw_in_pot and sol_ang at end of averaging period
  dplyr::mutate(
    sw_in_pot = potential_radiation(timestamp, md),
    sol_ang = sun_position(timestamp, md),
    # Stricter sun time to exclude any possible period with light
    sun_time = sun_time(sol_ang, shoulder_n = 3)
  ) %>%
  dplyr::filter(sw_in_pot == 0, sun_time == "night") %>%
  dplyr::summarize(dplyr::across(
    c(ppfd_in, sw_in, sw_out), mean, na.rm = TRUE)
  ) %>%
  as.list()

biomet <- dplyr::mutate(
  biomet, 
  ppfd_in = apply_offset(ppfd_in, -offset$ppfd_in),
  sw_in = apply_offset(sw_in, -offset$sw_in),
  sw_out = apply_offset(sw_out, -offset$sw_out)
)


### Add auxilliary data to Biomet data frame ===================================

# Other variables that can be calculated using core data

# Get coefficient of ppfd_in/sw_in relationship
frac_ppfd <- biomet %>% 
  dplyr::select(sw_in, ppfd_in) %>% 
  tidyr::drop_na() %>% 
  dplyr::summarize(sum(sw_in) / sum(ppfd_in)) %>% 
  purrr::pluck(1)

biomet <- dplyr::mutate(
  biomet,
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


### Save combined file with documentation ======================================

# Save the combined Biomet inputs as .csv file
biomet_out <- file.path(path_out, paste0("biomet_combined_", tag_out, ".csv"))
readr::write_csv(biomet, biomet_out)

# Create documentation for Biomet output
biomet_docu <- append(
  settings, list(
    files = biomet_input, mag_decl = mag_decl, frac_ppfd = frac_ppfd, 
    offsets = offset
  )
)
biomet_docu_out <- stringr::str_replace(biomet_out, ".csv", ".txt")

# Save documentation
sink(biomet_docu_out)
biomet_docu
sink()

# Finished
