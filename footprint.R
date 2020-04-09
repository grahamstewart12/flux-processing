### ============================================================================
# Footprint modeling ===========================================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

# Load the required packages
suppressWarnings(devtools::load_all("~/Desktop/RESEARCH/fluxtools"))
library(openeddy)
library(lubridate)
library(tidyverse)

source("~/Desktop/DATA/Flux/tools/engine/footprint.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")

# Initialize script settings & documentation
script_init(settings, site_metadata)


### Load required input data ===================================================

cat("Importing data files.\n")
site <- read.csv(site_input, stringsAsFactors = FALSE)

# Add timestamp components
site <- site %>%
  mutate(timestamp = ymd_hms(timestamp, tz = md$tz_name)) %>%
  add_time_comps()

# Import site AOI polygon
aoi <- rgdal::readOGR(aoi_input, verbose = FALSE)

# Import ERA data, get variables
era <- read.csv(era_input, stringsAsFactors = FALSE)
era <- era %>%
  mutate(timestamp = ymd_hms(timestamp) %>% with_tz(md$tz_name)) %>%
  select(timestamp, blh)

# Add to main data frame
site <- site %>%
  left_join(era, by = "timestamp") %>%
  # Interpolate hourly -> half-hourly
  mutate(blh = approx(seq_along(blh), blh, seq_along(blh))$y)


### Create a vectorized area boundary ==========================================
# - used for fetch filtering in QC

# Set the desired resolution (width of wind_dir bins)
# - smoothing not necessary, but removes some bias if area deliminated manually
#boundary <- roi_boundary(
#  aoi, res = 5, c(md$x_utm, md$y_utm), smooth = TRUE, smooth_method = "densify"
#)

# Write to .csv file
#boundary_out <- file.path(path_out, paste0("boundary_vec_", tag_out, ".csv"))
#write.csv(boundary, boundary_out, row.names = FALSE)

# Create documentation for phi output
#boundary_docu <- settings
#boundary_docu$files <- c(aoi_input)
#boundary_docu$settings <- list(res = 5, smooth = TRUE)
#boundary_docu_out <- str_replace(boundary_out, ".csv", ".txt")
# Save documentation
#sink(boundary_docu_out)
#boundary_docu 
#sink()


### Calculate 1-D footprints (fetch lengths) ===================================

if (control$fetch) {
  cat("Calculating half-hourly fetch lengths.\n")
  # Calculate half-hourly fetch lengths
  fetch <- fp_fetch(
    site, ws = "ws", ustar = "ustar", zeta = "zl", mo_length = "mo_length",
    z = md$tower_height, zd = md$displacement, zo = md$roughness_length, 
    percent = c(10, 30, 50, 70, 90, "peak"), method = control$fetch_model
  )
  
  # Rename with method ID tag
  fetch <- rename_all(
    fetch, list(~ stringr::str_c(., tolower(control$fetch_model)))
  )
  
  fetch <- bind_cols(select(site, timestamp), fetch)
  
  fetch_out <- file.path(
    path_out, paste0("fetch_", tag_out, ".csv")
  )
  write.csv(fetch, fetch_out, row.names = FALSE)
  
  # Create documentation for fetch output
  fetch_docu <- settings
  fetch_docu$files <- c(site_input)
  fetch_docu$method <- control$fetch_model
  fetch_docu_out <- str_replace(fetch_out, ".csv", ".txt")
  # Save documentation
  sink(fetch_docu_out)
  fetch_docu 
  sink()
}


### Calculate 2-D footprints and coverage ======================================

if (control$fp) {
  # Create new folder for output files
  fp_out <- file.path(path_out, paste0("half_hourly_", tag_out))
  dir.create(fp_out)
  
  # Select only necessary variables, remove missing data
  site_fp <- site %>% 
    select(timestamp, wd, ustar, mo_length, v_sigma) %>%
    drop_na()
  
  # Set up grid
  grid <- fp_grid(fetch = (md$tower_height - md$displacement) * 100)
  
  
  # Calculate footprints, create raster images, write to file
  # WARNING: this takes ~1 hr to run
  
  # Initialize loop
  len <- nrow(site_fp)
  phi_hh <- vector("list", length = len)
  progress <- progress_estimated(len)
  
  # Calculate footprint matrices
  for (i in 1:len) {
    
    # Calculate footprint, add to list
    fp_temp <- fp_calculate(
      wd = site_fp$wd[i],
      ustar = site_fp$ustar[i],
      mo_length = site_fp$mo_length[i],
      v_sigma = site_fp$v_sigma[i],
      z = md$tower_height,
      zd = md$displacement,
      zo = md$roughness_length,
      grid = grid,
      model = "H00"
    )
    
    # Set output name as the timestamp
    fp_temp_name <- format(site_fp$timestamp[i], "%Y-%m-%d-%H%M%S")
    
    # Rasterize
    fp_temp_rst <- fp_rasterize(
      fp_temp, grid = grid, coords = c(md$x_utm, md$y_utm), 
      crs = raster::crs(aoi)
    )
    
    # Clip to AOI
    fp_temp_rst <- raster::mask(fp_temp_rst, aoi) %>% raster::trim()
    
    # Calculate AOI coverage
    phi_hh[i] <- raster::cellStats(fp_temp_rst, sum)
    names(phi_hh)[i] <- fp_temp_name
    
    # Scale for more efficient storage as integer values
    fp_temp_rst <- fp_temp_rst * 1e9
    
    # Write to file
    raster::writeRaster(
      fp_temp_rst, file.path(fp_out, fp_temp_name), 
      format = "GTiff", overwrite = TRUE, datatype = "INT4S"
    )
    
    progress$tick()$print()
    
  }
  
  
  # Write file with just footprint data so this doesn't have to be run again
  fp_basics <- phi_hh %>% 
    tibble::enframe(name = "timestamp", value = "phi") %>% 
    mutate(timestamp = ymd_hms(timestamp, tz = md$tz_name)) %>%
    unnest(phi) %>%
    right_join(select(site, timestamp))
  
  fp_basics_out <- file.path(
    path_out, paste0("footprint_basics_", tag_out, ".csv")
  )
  write.csv(fp_basics, fp_basics_out, row.names = FALSE)
  
  # Create documentation for phi output
  fp_basics_docu <- settings
  fp_basics_docu$files <- c(site_input, aoi_input)
  fp_basics_docu$model_params <- list(
    method = "H00", fetch = attr(grid, "fetch"), res = attr(grid, "res")
  )
  fp_basics_docu_out <- str_replace(fp_basics_out, ".csv", ".txt")
  # Save documentation
  sink(fp_basics_docu_out)
  fp_basics_docu 
  sink()
}


### Calculate 2-D footprints - FFP method ======================================

if (control$ffp) {
  # See Kljun et al. 2015
  #source("~/Desktop/DATA/Flux/tools/FFP_R/calc_footprint_FFP.R")
  #source("~/Desktop/DATA/Flux/tools/FFP_R/calc_footprint_FFP_climatology.R")
  
  # Create new folder for output files
  ffp_out <- file.path(path_out, paste0("half_hourly_ffp_", tag_out))
  dir.create(ffp_out)
  
  # Select only necessary variables, remove missing data
  site_ffp <- site %>% 
    select(timestamp, wd, ustar, mo_length, v_sigma, blh) %>%
    drop_na()
  # Pre-filter known invalid cases (not using - just calculate all fp then see)
  # - NOTE: for ustar Kljun only indicate > 0.1, but adding upper limit of 2.0
  # - Olson et al. 2004 (FLUXNET) suggest 6.0, but here there are none >1 & <2
  # - this suggests spurious values above 2
  #filter(
  #  between(ustar, 0.1, 2), 
  #  ((md$tower_height - md$displacement) / mo_length) >= -15.5
  #)
  
  # Set up grid
  grid <- fp_grid(fetch = (md$tower_height - md$displacement) * 100)
  
  # Calculate footprints, create raster images, write to file
  # WARNING: this takes ~1 hr to run
  
  # Initialize loop
  len <- nrow(site_ffp)
  ffp_phi_hh <- vector("list", length = len)
  progress <- progress_estimated(len)
  
  # Calculate footprint matrices
  for (i in 1:len) {
    
    # Calculate footprint, add to list
    ffp_temp <- fp_calculate(
      wd = site_ffp$wd[i],
      ustar = site_ffp$ustar[i],
      mo_length = site_ffp$mo_length[i],
      v_sigma = site_ffp$v_sigma[i],
      blh = site_ffp$blh[i],
      z = md$tower_height,
      zd = md$displacement,
      zo = md$roughness_length,
      grid = grid,
      model = "K15"
    )
    
    # Set output name as the timestamp
    ffp_temp_name <- format(site_ffp$timestamp[i], "%Y-%m-%d-%H%M%S")
    
    # Rasterize
    ffp_temp_rst <- fp_rasterize(
      ffp_temp, grid = grid, coords = c(md$x_utm, md$y_utm), 
      crs = raster::crs(aoi)
    )
    
    # Clip to AOI
    ffp_temp_rst <- raster::mask(ffp_temp_rst, aoi) %>% raster::trim()
    
    # Calculate AOI coverage
    ffp_phi_hh[i] <- raster::cellStats(ffp_temp_rst, sum)
    names(ffp_phi_hh)[i] <- ffp_temp_name
    
    # Scale for more efficient storage as integer values
    ffp_temp_rst <- ffp_temp_rst * 1e9
    
    # Write to file
    raster::writeRaster(
      ffp_temp_rst, file.path(ffp_out, ffp_temp_name), 
      format = "GTiff", overwrite = TRUE, datatype = "INT4S"
    )
    
    progress$tick()$print()
  }
  
  # Write file with just footprint data so this doesn't have to be run again
  ffp_basics <- ffp_phi_hh %>% 
    tibble::enframe(name = "timestamp", value = "phi_ffp") %>% 
    mutate(timestamp = ymd_hms(timestamp, tz = md$tz_name)) %>%
    unnest(phi_ffp) %>%
    right_join(select(site, timestamp, ustar, mo_length)) %>%
    # Classify phi according to Gockede et al. 2008
    mutate(
      phi_ffp_class = case_when(
        # Homogenous measurements
        between(phi_ffp, 0.95, 1.00) ~ 0L,
        # Representative measurements
        between(phi_ffp, 0.80, 0.95) ~ 1L,
        # Acceptable measurements
        between(phi_ffp, 0.50, 0.80) ~ 2L,
        # Disturbed measurements
        between(phi_ffp, 0.00, 0.50) ~ 3L
      ),
      zeta = ((md$tower_height - md$displacement) / mo_length),
      ffp_inval = case_when(
        !between(ustar, 0.1, 2) ~ 1L,
        zeta < -15.5 ~ 1L,
        TRUE ~ 0L
      )
    ) %>% 
    select(-ustar, -mo_length, -zeta)
  
  ffp_basics_out <- file.path(
    path_out, paste0("footprint_basics_ffp_", tag_out, ".csv")
  )
  write.csv(ffp_basics, ffp_basics_out, row.names = FALSE)
  
  # Create documentation for phi output
  ffp_basics_docu <- settings
  ffp_basics_docu$files <- c(site_input, aoi_input)
  ffp_basics_docu$model_params <- list(
    method = "K15", fetch = attr(grid, "fetch"), res = attr(grid, "res")
  )
  ffp_basics_docu_out <- str_replace(ffp_basics_out, ".csv", ".txt")
  # Save documentation
  sink(ffp_basics_docu_out)
  ffp_basics_docu 
  sink()
}




