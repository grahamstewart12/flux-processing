### ============================================================================
# Combine EddyPro outputs ======================================================
### ============================================================================

# Purpose: Combine multiple eddypro output and/or biomet files, write the result
# to file. EddyPro is not always able to produce a single output file for a
# given time period due to changes in system configuration, etc. This script
# enables the user to resolve this issue by merging multiple output files (of
# potentially varying formats) with appropriate documentation.

# Inputs:

# Load the required packages
devtools::load_all("~/Desktop/RESEARCH/fluxtools")
library(openeddy)
library(lubridate)
library(tidyverse)

# Load the source files
#source("~/Desktop/DATA/Flux/tools/engine/combine_eddypro.R")
source("~/Desktop/DATA/Flux/tools/reference/combine_eddypro_control.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")


### Initialize script settings & documentation =================================

# Set directory where outputs are stored
# This assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)

# Set the folders containing the EddyPro outputs to be combined
ep_dirs <- pluck(control, settings$site, settings$year, "ep_dirs")
ep_paths <- file.path(wd, "processing_data", "00_eddypro_output", ep_dirs)

# Load metadata file
md <- pluck(site_metadata, settings$site)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing_data", "02_combine_eddypro", "output")


### Load, combine, and save the EddyPro files ==================================

lists <- read_eddypro_files(ep_paths)

# SPECIAL CASE: subset data based on correct calibrations
dir_start <- pluck(control, settings$site, settings$year, "dir_start")
dir_end <- pluck(control, settings$site, settings$year, "dir_end")
lists <- lists %>%
  map(
    map_at, ep_dirs[1], filter, 
    !between(timestamp, as_datetime(dir_start[1]), as_datetime(dir_start[2]))
  ) %>%
  map(
    map_at, ep_dirs[2], filter, 
    between(timestamp, as_datetime(dir_start[1]), as_datetime(dir_start[2]))
  )

# Combine and save the data as .csv files
combined_lists <- combine_eddypro(lists, path_out, tag_out, settings$year)


### Save a comprehensive dataset for further analysis ==========================

# Prepare output dataset
fluxnet <- lists %>% 
  pluck("fluxnet") %>%
  bind_rows() %>%
  select(timestamp, everything()) %>%
  arrange(timestamp) %>%
  # Convert all names to lowercase (better for R)
  rename_all(tolower) %>%
  # Fix errors caused by duplicated column names
  tibble::as_tibble(.name_repair = "universal") %>%
  # Remove columns that contain nothing (all NA)
  select_if(~!all(is.na(.)))

# Make sure timestamp forms a regular sequence for the entire year
fluxnet <- data.frame(
  timestamp = create_timesteps(settings$year, 48, shift_by = 30)
) %>%
  left_join(fluxnet, by = "timestamp")
fluxnet %>% summarize(first(timestamp), last(timestamp)) # Start/end?

# Save the combined fluxnet files as .csv file
fn_out <- file.path(path_out, paste0("eddypro_combined_", tag_out, ".csv"))
write.csv(fluxnet, fn_out, row.names = FALSE)

# Create documentation
fn_docu <- settings
fn_docu$files <- ep_paths %>% 
  map(list.files, pattern = ".csv", full.names = TRUE, recursive = TRUE) %>% 
  simplify() %>% 
  str_subset("fluxnet")
fn_docu_out <- str_replace(fn_out, ".csv", ".txt")
# Save documentation
sink(fn_docu_out)
fn_docu
sink()

