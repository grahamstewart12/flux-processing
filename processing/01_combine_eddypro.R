### ============================================================================
# Combine EddyPro outputs ======================================================
### ============================================================================

# Purpose: Combine multiple eddypro output and/or biomet files, write the result
# to file. EddyPro is not always able to produce a single output file for a
# given time period due to changes in system configuration, etc. This script
# enables the user to resolve this issue by merging multiple output files (of
# potentially varying formats) with appropriate documentation.

# Inputs:

start_time <- Sys.time()

# Load the required packages
devtools::load_all("~/Projects/flux/dscalr")
library(lubridate, warn.conflicts = FALSE)
library(tidyverse, warn.conflicts = FALSE)

# Load reference files
# source("reference/combine_eddypro_control.R")
# source("reference/site_metadata.R")


### Helper functions ===========================================================

# Helper function for locating files
str_subset_eddypro <- function(files) {
  names <- basename(files)
  keep <- names %>% 
    stringr::str_subset("^eddypro_") %>%
    stringr::str_subset("\\.csv$") %>%
    # Remove spectral analysis files
    stringr::str_subset("spectra", negate = TRUE) %>%
    # Risky subset to remove "user stats" files from consideration
    stringr::str_subset("user", negate = TRUE)
  inds <- which(names %in% keep)
  files[inds]
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- yaml::read_yaml(file.path("data", settings$site, "metadata.yml"))

# Set directory where outputs are stored
# This assumes that the subdirectory structure is already present
wd <- file.path("data", settings$site, settings$year)

# Set the folders containing the EddyPro outputs to be combined
# ep_dirs <- pluck(
#   control, settings$site, as.character(settings$year), "ep_dirs"
# )
ep_paths <- list.dirs(
  file.path(wd, "00_raw_output", "eddypro"), recursive = FALSE
)
#ep_paths <- file.path(wd, "00_raw_output", "eddypro", ep_dirs)
ep_dirs <- basename(ep_paths)

# Set tag for creating output file names
tag <- make_tag(settings$site, settings$year)

# Set path for output files
path_out <- file.path(wd, "01_combine_eddypro")


### Load, combine, and save the EddyPro files ==================================

# Read EddyPro settings
ep_settings <- ep_paths %>%
  map(list.files, pattern = ".eddypro", full.names = TRUE) %>%
  map(read_eddypro_settings) %>%
  map(as_tibble) %>%
  map2(ep_dirs, ~ mutate(.x, dir = .y)) %>%
  map(relocate, dir)
# all.equal(ep_settings[[1]], ep_settings[[2]])

# Find all EddyPro data files
data_files <- ep_paths %>% 
  map(list.files, pattern = ".csv", full.names = TRUE, recursive = TRUE) %>% 
  map(str_subset_eddypro)
file_types <- map(data_files, get_eddypro_file_type)

# Read EddyPro data (files may be large--be patient)
# data_ep <- data_files %>% 
#   map2(file_types, ~ map2(.x, .y, ~ read_eddypro(.x, type = .y))) %>% 
#   map(rlang::set_names, file_types)

# Set the lists to hold different types of output files
data_ep <- ep_dirs %>% 
  map2(data_files, ~ vctrs::vec_init(list(), length(.y))) %>%
  rlang::set_names(ep_dirs) %>%
  map2(file_types, ~ rlang::set_names(.x, .y))

# Read files into the lists
for (i in seq_along(ep_dirs)) {
  
  cat("Reading files from folder ", basename(ep_paths[i]), ": ", sep = "")
  
  for (j in seq_along(data_files[[i]])) {
    cat(paste0(file_types[[i]][j], "..."))
    temp <- quietly(read_eddypro)(
      data_files[[i]][j], type = file_types[[i]][j], guess_max = 5000
    )
    pluck(data_ep, i, file_types[[i]][j]) <- temp$result
    cat("done. ")
  }
  
  cat("\n")
}

data_ep_t <- data_ep %>%
  map(discard, is.null) %>%
  transpose() %>%
  map(discard, is.null)

# Subset data based on correct calibrations if applicable
recal <- map(md$recal_periods, as_datetime)
# recal <- control %>% 
#   pluck(settings$site, as.character(settings$year), "recal_periods") %>%
#   map(as_datetime)

if (!all(is.na(recal))) {
  
  cat("Removing uncalibrated periods...")
  
  # Recal dir should be labeled with "_recal"
  # - only supports one recal dir (for now)
  recal_dir <- str_subset(ep_dirs, "recal")
  
  # Get non-recal dirs
  other_dirs <- str_subset(ep_dirs, "recal", negate = TRUE)
  
  # Build logical expressions for subsetting periods
  recal_logic <- 1:length(recal) %>% 
    map(~ str_glue("between2(timestamp, recal[[", .x, "]])")) %>% 
    glue::glue_collapse(sep = " | ") %>% 
    str_glue("(", ., ")")
  other_logic <- str_glue("!", recal_logic)
  
  # Subset the data
  data_ep_t <- data_ep_t %>%
    # Remove non-recal periods from recal data
    map(map_at, recal_dir, filter, !!rlang::parse_expr(recal_logic)) %>%
    # Remove recal periods from non-recal data
    map(map_at, other_dirs, filter, !!rlang::parse_expr(other_logic))
  
  cat("done.\n")
}

cat("Merging data...")

# Combine and save the data as .csv files
data_merged <- data_ep_t %>% 
  map(merge_eddypro) %>%
  # Subset records by year
  map(filter, year(timestamp - 900) %in% settings$year) 

# Combine the data
# combined_lists <- lists_sub %>%
#   map(bind_rows) %>%
#   map(arrange, timestamp) %>%
#   # Subset year
#   map(filter, year(timestamp - 900) %in% settings$year) 

cat("done.\n")

# Add units and varnames back to data
# combined_lists <- combined_lists %>%
#   map(dplyr::select, -timestamp) %>%
#   map2(varnames, ~ add_attr(.x, "varnames", .y)) %>%
#   map2(units, ~ add_attr(.x, "units", .y))

# Write combined data as single .csv files

# Bind settings together and write as a table
ep_settings_out <- file.path(
  path_out, "settings", paste0("processing_combined_", tag, ".csv")
)
write_csv(bind_rows(ep_settings), ep_settings_out)

files_out <- file.path(
  path_out, 
  names(data_merged), paste0(names(data_merged), "_combined_", tag, ".csv")
)

cat("\nWriting datasets: ")
for (i in seq_along(data_merged)) {
  
  cat(names(data_merged)[i], "...", sep = "")
  
  # Each EddyPro output type gets its own folder
  if (!dir.exists(dirname(files_out[i]))) {
    dir.create(dirname(files_out[i]))
  }
  
  write_eddypro(data_merged[[i]], files_out[i])
  
  cat("done. ")
}


### Save a comprehensive dataset for further analysis ==========================

cat("\nPreparing FLUXNET data for further analysis...")

# Prepare output dataset
fluxnet <- data_merged %>%
  pluck("fluxnet") %>%
  drop_all_attrs() %>%
  # mutate(timestamp = as_datetime(timestamp)) %>%
  # Convert all names to lowercase (better for R)
  rename_with(tolower) %>%
  # Fix errors caused by duplicated column names
  as_tibble(.name_repair = "unique") %>%
  # Remove columns that contain nothing (all NA)
  select(-where(~ all(is.na(.x))))

# Subset necessary data
# fluxnet <- fluxnet %>%
#   # Remove metadata and uncorrected data
#   select(
#     -starts_with("badm"), -starts_with("inst"), -starts_with("custom"), 
#     -ends_with("method"), -ends_with("uncorr"), -matches("stage[[:digit:]]")
#   ) %>%
#   # Remove stats & QC details (this is ultimately pulled from full_output)
#   select(
#     -ends_with("kid"), -ends_with("zcd"), -ends_with("corrdiff"), 
#     -ends_with("nsr"), -ends_with("ss"), -ends_with("itc"), 
#     -ends_with("test"), -ends_with("median"), -ends_with("p25"), 
#     -ends_with("p75"), -ends_with("spikes"), -ends_with("nrex"), 
#     -ends_with("skw"), -ends_with("kur")
#   )

# Make sure timestamp forms a regular sequence for the entire year
fluxnet <- tibble(
  timestamp = create_timesteps(settings$year, 48, shift_by = 30)
) %>%
  left_join(fluxnet, by = "timestamp")

cat("done.\n")
cat("Writing FLUXNET data...")

# Save the combined fluxnet files as .csv file
fn_out <- file.path(
  path_out, "eddypro", paste0("eddypro_combined_", tag, ".csv")
)
write_csv(fluxnet, fn_out)

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

end_time <- Sys.time()
elapsed_time <- round(unclass(end_time - start_time), 1)

cat("done.\n")
cat(
  "Finished processing in ", elapsed_time, " ", attr(elapsed_time, "units"), 
  ".\n", sep = ""
)
cat("Output located in", path_out, "\n")
