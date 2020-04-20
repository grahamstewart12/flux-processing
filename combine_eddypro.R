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
suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
library(openeddy, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(tidyverse, warn.conflicts = FALSE)

# Load the source files
source("~/Desktop/DATA/Flux/tools/reference/combine_eddypro_control.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")


### Initialize script settings & documentation =================================

# Set directory where outputs are stored
# This assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)

# Set the folders containing the EddyPro outputs to be combined
ep_dirs <- purrr::pluck(control, settings$site, settings$year, "ep_dirs")
ep_paths <- file.path(wd, "processing_data", "00_eddypro_output", ep_dirs)

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing_data", "02_combine_eddypro", "output")


### Load, combine, and save the EddyPro files ==================================

# Helper function for locating files
str_subset_eddypro <- function(files) {
  names <- basename(files)
  keep <- names %>% 
    stringr::str_subset("^eddypro_") %>%
    stringr::str_subset("\\.csv$") %>%
    # Risky subset to remove "user stats" files from consideration
    stringr::str_subset("user", negate = TRUE)
  inds <- which(names %in% keep)
  files[inds]
}

# Set the lists to hold different types of output files
lists <- list(
  full_output = vector("list", length = length(ep_paths)),
  metadata = vector("list", length = length(ep_paths)),
  biomet = vector("list", length = length(ep_paths)),
  fluxnet = vector("list", length = length(ep_paths)),
  qc_details = vector("list", length = length(ep_paths)),
  st7 = vector("list", length = length(ep_paths))
)
lists <- 1:6 %>%
  purrr::map(~ vector("list", length = length(ep_paths))) %>%
  rlang::set_names(
    c("full_output", "metadata", "biomet", "fluxnet", "qc_details", "st7")
  )

data_files <- vector("list", length = length(ep_paths))

read_eddypro <- function(file, timestamp = paste(date, time), ...) {
  
  timestamp <- rlang::enexpr(timestamp)
  
  # Default: skip = 0, units = TRUE
  data <- openeddy::read_eddy(file, ...)
  
  data %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(timestamp = lubridate::ymd_hm(!!timestamp)) %>%
    dplyr::filter(!is.na(timestamp))
}

# Read files into the lists
for (i in seq_along(ep_paths)) {
  
  cat("Reading files from folder ", basename(ep_paths[i]), ": ", sep = "")
  
  data_files[[i]] <- ep_paths[i] %>%
    list.files(pattern = ".csv", full.names = TRUE, recursive = TRUE) %>% 
    str_subset_eddypro()
  files <- data_files[[i]]
  files_col <- stringr::str_c(files, collapse = "|")
  
  # FULL OUTPUT FILE
  if (stringr::str_detect(files_col, "full_output")) {
    cat("full_output...")
    file <- stringr::str_subset(files, "full_output")
    temp <- read_eddypro(file, skip = 1)
    purrr::pluck(lists, "full_output", i) <- temp
    cat("done. ")
  } 
  
  # METADATA FILE
  if (stringr::str_detect(files_col, "metadata")) {
    cat("metadata...")
    file <- stringr::str_subset(files, "metadata")
    temp <- read_eddypro(file, units = FALSE)
    purrr::pluck(lists, "metadata", i) <- temp
    cat("done. ")
  }
  
  # BIOMET FILE
  if (stringr::str_detect(files_col, "biomet")) {
    cat("biomet...")
    file <- stringr::str_subset(files, "biomet")
    temp <- read_eddypro(file)
    purrr::pluck(lists, "biomet", i) <- temp
    cat("done. ")
  }
  
  # FLUXNET FILE
  if (stringr::str_detect(files_col, "fluxnet")) {
    cat("fluxnet...")
    file <- stringr::str_subset(files, "fluxnet")
    temp <- read_eddypro(file, TIMESTAMP_END, units = FALSE)
    purrr::pluck(lists, "fluxnet", i) <- temp
    cat("done. ")
  }
  
  # QC DETAILS FILE
  if (stringr::str_detect(files_col, "qc_details")) {
    cat("qc_details...")
    file <- stringr::str_subset(files, "qc_details")
    temp <- read_eddypro(file, skip = 1)
    purrr::pluck(lists, "qc_details", i) <- temp
    cat("done. ")
  }
  
  # STATS FILE
  if (stringr::str_detect(files_col, "st7")) {
    cat("stats...")
    file <- stringr::str_subset(files, "st7")
    temp <- read_eddypro(file, skip = 1, units = FALSE)
    purrr::pluck(lists, "st7", i) <- temp
    cat("done. ")
  }
  
  cat("\n")
}

lists_sub <- lists %>%
  purrr::map(rlang::set_names, basename(ep_paths)) %>%
  purrr::map(purrr::discard, is.null)

# Subset data based on correct calibrations if applicable
recal <- control %>% 
  purrr::pluck(settings$site, settings$year, "recal_periods") %>%
  purrr::map(lubridate::as_datetime)

if (!all(is.na(recal))) {
  
  cat("Removing uncalibrated periods...")
  
  # Recal dir should be labeled with "_recal"
  # - only supports one recal dir (for now)
  recal_dir <- stringr::str_subset(ep_dirs, "recal")
  
  # Get non-recal dirs
  other_dirs <- stringr::str_subset(ep_dirs, "recal", negate = TRUE)
  
  # Build logical expressions for subsetting periods
  recal_logic <- 1:length(recal) %>% 
    purrr::map(
      ~ stringr::str_glue("between2(timestamp, recal[[", .x, "]])")
    ) %>% 
    glue::glue_collapse(sep = " | ") %>% 
    stringr::str_glue("(", ., ")")
  other_logic <- stringr::str_glue("!", recal_logic)
  
  # Subset the data
  lists_sub <- lists_sub %>%
    # Remove non-recal periods from recal data
    purrr::map(
      purrr::map_at, recal_dir, dplyr::filter, !!rlang::parse_expr(recal_logic)
    ) %>%
    # Remove recal periods from non-recal data
    purrr::map(
      purrr::map_at, other_dirs, dplyr::filter, !!rlang::parse_expr(other_logic)
    )
}

cat("done.\n")
cat("Coalescing data...")

# Combine and save the data as .csv files

# Save original varnames and units
varnames <- lists_sub %>%
  purrr::map_depth(2, dplyr::select, -timestamp) %>%
  purrr::map_depth(2, ~ openeddy::varnames(., names = TRUE)) %>%
  purrr::map_depth(2, dplyr::bind_rows) %>%
  purrr::map(dplyr::bind_rows) %>%
  purrr::map(dplyr::slice, 1) %>%
  purrr::map(unlist, use.names = TRUE) %>%
  purrr::map(dplyr::na_if, "-") %>%
  purrr::map(~ dplyr::coalesce(., names(.))) %>%
  purrr::map(tidyr::replace_na, "-") %>%
  purrr::map(unname)
units <- lists_sub %>%
  purrr::map_depth(2, dplyr::select, -timestamp) %>%
  purrr::map_depth(2, ~ openeddy::units(., names = TRUE)) %>%
  purrr::map_depth(2, dplyr::bind_rows) %>%
  purrr::map(dplyr::bind_rows) %>%
  purrr::map(dplyr::slice, 1) %>%
  purrr::map(unlist, use.names = FALSE)

# Combine the data
combined_lists <- lists_sub %>%
  purrr::map(dplyr::bind_rows) %>%
  purrr::map(dplyr::arrange, timestamp) %>%
  # Subset year
  purrr::map(
    dplyr::filter, lubridate::year(timestamp - 900) %in% settings$year
  ) 

cat("done.\n")
cat("Resolving duplicate timestamps: ")

# Handle duplicated observations
for (i in seq_along(combined_lists)) {
  
  name <- names(combined_lists)[i]
  cat(name, "...", sep = "")
  
  combined_lists[[i]] <- combined_lists[[i]] %>% 
    dplyr::rowwise() %>% 
    # Number of NA values in each row
    dplyr::mutate(n_na = sum(is.na(dplyr::c_across(is.numeric)))) %>% 
    dplyr::group_by(timestamp) %>%
    # Sort by least number of NAs per timestamp
    dplyr::arrange(n_na, .by_group = TRUE) %>%
    # Select the first (i.e. most complete) observation
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n_na)
  
  cat("done. ")
}

# Add units and varnames back to data
combined_lists <- combined_lists %>%
  purrr::map(dplyr::select, -timestamp) %>%
  purrr::map2(varnames, ~ add_attr(.x, "varnames", .y)) %>%
  purrr::map2(units, ~ add_attr(.x, "units", .y))

# Write combined data as single .csv files

cat("\nWriting datasets: ")

for (i in seq_along(combined_lists)) {
  
  name <- names(combined_lists)[i]
  cat(name, "...", sep = "")
  
  if (name %in% c("fluxnet", "st7")) {
    # FLUXNET and stats outputs don't have units
    write.csv(
      combined_lists[[i]],
      file.path(path_out, paste0(name, "_combined_", tag_out, ".csv")),
      na = "-9999", row.names = FALSE
    )
  } else {
    openeddy::write_eddy(
      combined_lists[[i]],
      file.path(path_out, paste0(name, "_combined_", tag_out, ".csv")),
      col.names = varnames[[i]]
    )
  }
  
  cat("done. ")
}


### Save a comprehensive dataset for further analysis ==========================

cat("\nPreparing FLUXNET data for further analysis...")

# Prepare output dataset
fluxnet <- combined_lists %>%
  purrr::pluck("fluxnet") %>%
  drop_all_attributes() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(timestamp = lubridate::ymd_hm(TIMESTAMP_END)) %>%
  dplyr::select(timestamp, dplyr::everything()) %>%
  # Convert all names to lowercase (better for R)
  dplyr::rename_all(tolower) %>%
  # Fix errors caused by duplicated column names
  tibble::as_tibble(.name_repair = "universal") %>%
  # Remove columns that contain nothing (all NA)
  dplyr::select_if(~ !all(is.na(.)))

# Make sure timestamp forms a regular sequence for the entire year
fluxnet <- data.frame(
  timestamp = create_timesteps(settings$year, 48, shift_by = 30)
) %>%
  dplyr::left_join(fluxnet, by = "timestamp")

cat("done.\n")
cat("Writing FLUXNET data...")

# Save the combined fluxnet files as .csv file
fn_out <- file.path(path_out, paste0("eddypro_combined_", tag_out, ".csv"))
write.csv(fluxnet, fn_out, row.names = FALSE)

# Create documentation
fn_docu <- settings
fn_docu$files <- ep_paths %>% 
  purrr::map(
    list.files, pattern = ".csv", full.names = TRUE, recursive = TRUE
  ) %>% 
  purrr::simplify() %>% 
  stringr::str_subset("fluxnet")
fn_docu_out <- stringr::str_replace(fn_out, ".csv", ".txt")
# Save documentation
sink(fn_docu_out)
fn_docu
sink()

end_time <- Sys.time()
elapsed_time <- round(unclass(end_time - start_time), 3)

cat("done.\n")
cat(
  "Finished processing in ", elapsed_time, " ", attr(elapsed_time, "units"), 
  ".\n", sep = ""
)
cat("Output located in", path_out, "\n")
