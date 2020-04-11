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
    lists[["full_output"]][[i]] <- openeddy::read_eddy(file, skip = 1) %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::mutate(timestamp = lubridate::ymd_hm(paste(date, time))) %>%
      dplyr::filter(!is.na(timestamp))
    cat("done. ")
  } 
  # METADATA FILE
  if (stringr::str_detect(files_col, "metadata")) {
    cat("metadata...")
    file <- stringr::str_subset(files, "metadata")
    lists[["metadata"]][[i]] <- openeddy::read_eddy(file, units = FALSE) %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::mutate(timestamp = lubridate::ymd_hm(paste(date, time))) %>%
      dplyr::filter(!is.na(timestamp))
    cat("done. ")
  }
  # BIOMET FILE
  if (stringr::str_detect(files_col, "biomet")) {
    cat("biomet...")
    file <- stringr::str_subset(files, "biomet")
    lists[["biomet"]][[i]] <- openeddy::read_eddy(file) %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::mutate(timestamp = lubridate::ymd_hm(paste(date, time))) %>%
      dplyr::filter(!is.na(timestamp))
    cat("done. ")
  }
  # FLUXNET FILE
  if (stringr::str_detect(files_col, "fluxnet")) {
    cat("fluxnet...")
    file <- stringr::str_subset(files, "fluxnet")
    lists[["fluxnet"]][[i]] <- file %>% 
      openeddy::read_eddy(units = FALSE) %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::mutate(timestamp = lubridate::ymd_hm(TIMESTAMP_END)) %>%
      dplyr::filter(!is.na(timestamp))
    cat("done. ")
  }
  # QC DETAILS FILE
  if (stringr::str_detect(files_col, "qc_details")) {
    cat("qc_details...")
    file <- stringr::str_subset(files, "qc_details")
    lists[["qc_details"]][[i]] <- file %>%
      openeddy::read_eddy(skip = 1) %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::mutate(timestamp = lubridate::ymd_hm(paste(date, time))) %>%
      dplyr::filter(!is.na(timestamp))
    cat("done. ")
  }
  # STATS FILE
  if (stringr::str_detect(files_col, "st7")) {
    cat("stats...")
    file <- stringr::str_subset(files, "st7")
    lists[["st7"]][[i]] <- file %>%
      openeddy::read_eddy(skip = 1, units = FALSE) %>%
      tibble::as_tibble(.name_repair = "unique") %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::mutate(timestamp = lubridate::ymd_hm(paste(date, time))) %>%
      dplyr::filter(!is.na(timestamp))
    cat("done. ")
  }
  cat("\n")
}

#ep_path_nms <- basename(ep_paths)

lists <- lists %>%
  purrr::map(rlang::set_names, basename(ep_paths)) %>%
  purrr::map(purrr::discard, is.null)
#lists <- read_eddypro_files(ep_paths)

# Subset data based on correct calibrations
# - I HATE this but it works
dir_start <- purrr::pluck(control, settings$site, settings$year, "dir_start")
dir_end <- purrr::pluck(control, settings$site, settings$year, "dir_end")
lists <- lists %>%
  purrr::map(
    purrr::map_at, ep_dirs[1], dplyr::filter, 
    !dplyr::between(
      timestamp, 
      lubridate::as_datetime(dir_start[1]), lubridate::as_datetime(dir_start[2])
    )
  ) %>%
  purrr::map(
    purrr::map_at, ep_dirs[2], dplyr::filter, 
    dplyr::between(
      timestamp, 
      lubridate::as_datetime(dir_start[1]), lubridate::as_datetime(dir_start[2])
    )
  )

# Combine and save the data as .csv files

cat("Coalescing data.\n")

# Save original varnames and units
varnames <- lists %>%
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
units <- lists %>%
  purrr::map_depth(2, dplyr::select, -timestamp) %>%
  purrr::map_depth(2, ~ openeddy::units(., names = TRUE)) %>%
  purrr::map_depth(2, dplyr::bind_rows) %>%
  purrr::map(dplyr::bind_rows) %>%
  purrr::map(dplyr::slice, 1) %>%
  purrr::map(unlist, use.names = FALSE)

# Combine the data
combined_lists <- lists %>%
  purrr::map(dplyr::bind_rows) %>%
  purrr::map(dplyr::arrange, timestamp) %>%
  # Subset year
  purrr::map(
    dplyr::filter, lubridate::year(timestamp - 900) %in% settings$year
  ) %>%
  purrr::map(dplyr::select, -timestamp) %>%
  purrr::map2(varnames, ~ add_attr(.x, "varnames", .y)) %>%
  purrr::map2(units, ~ add_attr(.x, "units", .y))

cat("Writing datasets: ")

for (i in seq_along(combined_lists)) {
  name <- names(combined_lists)[i]
  cat(name, "...", sep = "")
  
  if (name %in% c("fluxnet", "st7")) {
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

#combined_lists <- combine_eddypro(lists, path_out, tag_out, settings$year)


### Save a comprehensive dataset for further analysis ==========================

# Prepare output dataset
fluxnet <- lists %>% 
  purrr::pluck("fluxnet") %>%
  dplyr::bind_rows() %>%
  dplyr::select(timestamp, dplyr::everything()) %>%
  dplyr::arrange(timestamp) %>%
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
#fluxnet %>% summarize(first(timestamp), last(timestamp)) # Start/end?

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

