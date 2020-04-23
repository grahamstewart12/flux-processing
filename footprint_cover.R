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

# Performs a vectorized function on x and returns x as a matrix
with_matrix <- function(x, .f) {
  
  dims <- dims(x)
  mat_fun <- rlang::as_function(.f)
  
  vec <- mat_fun(as.vector(x))
  
  matrix(vec, nrow = dims[1], ncol = dims[2])
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


### Load required input data ===================================================

# Import water level data
wtd <- readr::read_csv(wtd_path)

# Import footprint grid (including AOI grid)
grid <- read_grid(file.path(fp_path, "grid"), names = c("x", "y", "aoi"))

# Import site image
img <- raster::raster(img_path)
#img_key <- c("wtr" = 1, "veg" = 2, "frs" = 3)
#img_key <- tibble::tibble(name = c("wtr", "veg", "frs"), value = c(1, 2, 3))
#plot_matrix(img_grid)
#img <- raster::setValues(img, factor(raster::values(img)))

# Import DEM
dem <- raster::raster(dem_path)
#plot_matrix(dem_grid)

demw <- raster::raster(demw_path)
#flood_key <- c("dry" = -1, "wet" = 1)
#flood_key <- tibble::tibble(name = c("wet", "dry"), value = c(1, -1))
#plot_matrix(demw_grid)
browser()

### Prepare input data =======================================

# Get measurement elevation (from NON-CORRECTED DEM)
#meas_elev <- raster::extract(dem, data.frame(x = md$x_utm, y = md$y_utm))

# Map imagery onto same grid as footprint
img_grid <- snap_to_grid(img, grid, c(md$x_utm, md$y_utm))
# use bilinear interpolation for DEM?
dem_grid <- snap_to_grid(dem, grid, c(md$x_utm, md$y_utm))
demw_grid <- snap_to_grid(demw, grid, c(md$x_utm, md$y_utm))

# Reclassify cover types with labels
img_grid <- with_matrix(
  img_grid, ~ dplyr::recode(.x, `1` = "wtr", `2` = "veg", `3` = "frs")
)

# Create key for combined cover types & flooding
#key <- img_key %>% 
#  purrr::map(~ .x * flood_key) %>%
#  purrr::flatten_dbl() %>%
#  rlang::set_names(
#    names(img_key) %>% 
#      purrr::map(~ stringr::str_c(.x, names(flood_key), sep = "_")) %>%
#      purrr::flatten_chr()
#  )


### Retrieve footprints, calculate cover =======================================

# Get list of footprint files
fp_files <- list.files(file.path(fp_path, "footprints"), full.names = TRUE)

# Initialize loop
# should I read all matrices first and then perform calculations??

# Empty list with names from footprint files (i.e. the timestamps)
fp_timestamps <- fp_files %>%
  basename() %>%
  stringr::str_sub(1, -5) %>%
  lubridate::ymd_hms(tz = md$tz_name)
cover <- vector("list", length = length(fp_files))
wtd <- site %>%
  dplyr::filter(timestamp %in% fp_timestamps) %>%
  dplyr::pull(wtd_f)
p <- progress_info(length(fp_files))

for (i in 1:length(fp_files)) {
  
  # Read and scale footprint matrix
  fp_temp <- read_matrix(fp_files[i])
  
  # Mask footprint to 85% to account for rapid expansion upon approaching 100%
  # - recommended: between 80% and 90% (Kljun et al. 2015)
  # - this means data should be first filtered for phi > 0.85
  fp_mask <- mask_source_area(fp_temp, p = 0.85)
  #fp_temp <- fp_percent(fp_temp, 0.85, fill = 0)
  
  # Normalize (???)
  #fp_temp <- raster::calc(fp_temp, function(x) x / sum(x, na.rm = TRUE))
  
  # Calculate water level over AOI
  aoi_wtd <- md$elev_well + wtd[i] - demw_re
  # Classify flooded area
  aoi_flooded <- ifelse(aoi_wtd >= 0, "wet", "dry")
  
  # Can use cross_grids() if want to combine more than 2 matrices
  # - this is kind of awkward but flexible
  # - automatically detects if one matrix is all NA and drops it
  # - ensures that all cells are classified properly
  aoi_cover <- cross_grids(img_grid, aoi_flooded)
  
  # Calculate cover-type weights
  cover_wt <- summarize_cover(fp_temp, aoi_cover)
  
  # Calculate footprint-weighted WTD
  wtd_wt <- summarize_cover(fp_temp, aoi_wtd, type = "numeric")
  
  # Add all weights to list 
  cover[[i]] <- c(cover_wt, "wtd_wt" = wtd_wt)
  
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

# Check patch coverage distributions
covers %>%
  select(-wtd, -wtd_wt, -wtr_dry, -veg_dry) %>%
  gather("type", "wt", -timestamp) %>%
  ggplot(aes(wt, fill = type, color = type)) +
  #facet_wrap(~ type, scales = "free") +
  geom_density(alpha = 0.2)

# Write output to file
covers_out <- file.path(
  wd, "processing_data", "13_footprint_cover", "output", 
  paste0("footprint_cover_", tag_out, ".csv")
)
write_eddy(all_covers, covers_out)


