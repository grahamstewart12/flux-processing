### ============================================================================
# Spatial analysis of flux area ================================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

# Load the required packages
suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
devtools::load_all("~/R Projects/footprints")
library(progress)
library(lubridate)
library(tidyverse)

source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")


### Helper functions ===========================================================

enlist_names <- function(x) {
  
  x %>%
    rlang::set_names() %>%
    purrr::map(~ NULL)
}

progress_info <- function(len) {
  
  progress_bar$new(
    total = len, 
    format = paste0(
      "[:spin] Completed: :current (:percent)  ", 
      "Elapsed: :elapsed  Remaining: :eta"
    )
  )
}

cross_keys <- function(...) {
  
  dots <- rlang::list2(...)
  
  values <- dots %>%
    tidyr::crossing() %>%
    dplyr::summarize(value = prod(dplyr::c_across()))
}

cross_grids <- function(...) {
  
  dots <- rlang::list2(...)
  
  # Assume that all matrices have same dimensions
  dims <- dots %>% purrr::map(dim) %>% purrr::pluck(1)
  
  crossed <- dots %>%
    purrr::map(as.vector) %>%
    # Drop elements that only contain missing values
    purrr::discard(~ all(is.na(.x))) %>%
    purrr::lift_dl(stringr::str_c)(sep = "_")
  
  matrix(crossed, nrow = dims[1], ncol = dims[2])
}

### Initialize script settings & documentation =================================

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)
path_in <- file.path(wd, "processing_data", "06_footprint", "output")

# Set the file names
# Half-hourly water level
wtd_path <- file.path(
  "~/Desktop", "DATA", "PT", "output", settings$site, 
  paste0(settings$site, "_", settings$year, ".csv")
)

# Cover type image
img_path <- file.path(
  dirname(wd), "imagery", "image_classification", "output", "classified.tif"
)

# DEM
dem_path <- file.path(
  dirname(wd), "raster_processing", paste0(settings$site, "_DEM_proc.tif")
)

# WTD-corrected DEM
demw_path <- file.path(
  dirname(wd), "raster_processing", paste0(settings$site, "_DEM_wtdcorr.tif")
)

# TODO import other imagery into a big list for footprint weights 
# - e.g. NDVI, NDWI, texture stuff

aoi_path <- file.path(dirname(wd), paste0(settings$site, "_area_proc"))

fp_path <- latest_version(path_in, "half_hourly_ffp", "")

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)

# Set path for output files
path_out <- file.path(wd, "processing_data", "13_footprint_cover", "output")


### Load required input data ===================================================

# Import water level data
site_wtd <- readr::read_csv(
  wtd_path, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Import footprint grid (including AOI grid)
grid <- read_grid(file.path(fp_path, "grid"), names = c("x", "y", "aoi"))

# Import site image
img <- raster::raster(img_path)
#img_key <- c("wtr" = 1, "veg" = 2, "frs" = 3)

# Import DEM
dem <- raster::raster(dem_path)

demw <- raster::raster(demw_path)
#flood_key <- c("dry" = -1, "wet" = 1)


### Prepare input data =======================================

# Get measurement elevation (from NON-CORRECTED DEM)
#meas_elev <- raster::extract(dem, data.frame(x = md$x_utm, y = md$y_utm))

# Map imagery onto same grid as footprint
img_grid <- snap_to_grid(img, grid, c(md$x_utm, md$y_utm))
# use bilinear interpolation for DEM?
dem_grid <- snap_to_grid(dem, grid, c(md$x_utm, md$y_utm))
demw_grid <- snap_to_grid(demw, grid, c(md$x_utm, md$y_utm))

# Reclassify cover types with labels
#img_grid <- with_matrix(
#  img_grid, ~ dplyr::recode(.x, `1` = "wtr", `2` = "veg", `3` = "frs")
#)

# Create key for combined cover types & flooding
key <- img_key %>% 
  purrr::map(~ .x * flood_key) %>% 
  purrr::imap(~ rlang::set_names(.x, stringr::str_c(.y, "_", names(.)))) %>%
  purrr::flatten() %>%
  purrr::simplify()


### Retrieve footprints, calculate cover =======================================

# Get list of footprint files
fp_files <- list.files(file.path(fp_path, "footprints"), full.names = TRUE)

# Initialize loop
# should I read all matrices first and then perform calculations??

# Empty list with names from footprint files (i.e. the timestamps)
fp_timestamps <- fp_files %>%
  basename() %>%
  stringr::str_sub(1, -5) %>%
  lubridate::ymd_hms()
cover <- vector("list", length = length(fp_files))
wtd <- site_wtd %>%
  dplyr::right_join(
    tibble::enframe(fp_timestamps, name = NULL, value = "timestamp")
  ) %>%
  dplyr::pull(wtd_f)
p <- progress_info(length(fp_files))

for (i in 1:length(fp_files)) {
  browser()
  # Read and scale footprint matrix
  fp_temp <- read_matrix(fp_files[i]) # 6.61 ms

  # Mask footprint to 85% to account for rapid expansion upon approaching 100%
  # - recommended: between 80% and 90% (Kljun et al. 2015)
  # - this means data should be first filtered for phi > 0.85
  fp_mask <- mask_source_area(fp_temp, p = 0.85, mask_value = 0) # 4.17 ms
  fp_mask <- fp_mask * grid$aoi # 0.192 ms
  
  # Normalize (???)
  #fp_temp <- raster::calc(fp_temp, function(x) x / sum(x, na.rm = TRUE))
  
  # Calculate water level over AOI
  #aoi_wtd <- md$elev_tower + wtd[i] - demw_grid # 0.131 ms
  # Give wtd a placeholder value if it is NA, so that other covers can be summed
  aoi_wtd <- md$elev_tower + tidyr::replace_na(wtd[i], 0.1) - demw_grid
  # Classify flooded area
  aoi_flooded <- with_matrix(aoi_wtd, ~ dplyr::if_else(.x >= 0, 1, 0)) 
  # 3.15 ms
  #aoi_flooded <- ifelse(aoi_wtd >= 0, "wet", "dry") # 7.80 ms
  
  # Can use cross_grids() if want to combine more than 2 matrices
  # - this is kind of awkward but flexible
  # - automatically detects if one matrix is all NA and drops it
  # - ensures that all cells are classified properly
  #aoi_cover <- cross_grids(img_grid, aoi_flooded) # 5.52 ms
  aoi_cover <- img_grid * sign(aoi_wtd)
  
  # Calculate cover-type weights
  cover_wt <- summarize_cover(fp_mask, aoi_cover, levels = key) # 3.79 ms
  # CHECK: sum of these should equal sum(fp_temp)
  
  # Calculate footprint-weighted WTD & flooded area
  wtd_wt <- summarize_cover(fp_mask, aoi_wtd, type = "numeric") # 0.605 ms
  flood_wt <- summarize_cover(fp_mask, aoi_flooded, type = "numeric") # 0.605 ms
  
  # Add all weights to list 
  cover[[i]] <- c(cover_wt, "wtd_wt" = wtd_wt, "flood_wt" = flood_wt)

  p$tick()
}

# Gather everything into one data frame
fp_cover <- cover %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(timestamp = fp_timestamps) %>%
  dplyr::relocate(timestamp)

 # Sum individual cover types and wet/dry
fp_cover <- fp_cover %>%
  dplyr::mutate(
    # If no WTD data, only cover types were calculated so don't overwrite these
    wtr = dplyr::if_else(is.na(wtr), wtr_wet + wtr_dry, .),
    veg = dplyr::if_else(is.na(veg), veg_wet + veg_dry, .),
    frs = dplyr::if_else(is.na(frs), frs_wet + frs_dry, .),
    wet = wtr_wet + veg_wet + frs_wet,
    dry = wtr_dry + veg_dry + frs_dry
  )

# Write output to file
covers_out <- file.path(path_out, paste0("footprint_cover_", tag_out, ".csv"))
readr::write_csv(all_covers, covers_out)


