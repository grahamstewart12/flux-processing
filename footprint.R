### ============================================================================
# Footprint modeling ===========================================================
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

progress_info <- function(len) {
  
  progress_bar$new(
    total = len, 
    format = paste0(
      "[:spin] Completed: :current (:percent)  ", 
      "Elapsed: :elapsed  Remaining: :eta"
    )
  )
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site, settings$year)
path_in <- file.path(wd, "processing_data", "05_biomet_gapfill", "output")
  
# Input file - biomet output with QC flags
site_input <- latest_version(path_in, "biomet_gf")
  
# Input file - processed ERA data for site location
era_input <- latest_version(
  file.path("~/Desktop", "DATA", "Flux", "JLL", "all", "output"), "era_proc"
)
  
# AOI
aoi_input <- file.path(dirname(wd), paste0(settings$site, "_area_proc"))
  
# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)
  
# Set path for output files
path_out <- file.path(wd, "processing_data", "06_footprint", "output")


### Load required input data ===================================================

cat("Importing data files.\n")
site <- read.csv(site_input, stringsAsFactors = FALSE)

# Add timestamp components
site <- site %>%
  dplyr::mutate(timestamp = lubridate::ymd_hms(timestamp, tz = md$tz_name)) %>%
  add_time_comps()

# Import site AOI polygon
aoi <- rgdal::readOGR(aoi_input, verbose = FALSE)

# Import ERA data, get variables
era <- read.csv(era_input, stringsAsFactors = FALSE)
era <- era %>%
  dplyr::mutate(
    timestamp = lubridate::ymd_hms(timestamp) %>% lubridate::with_tz(md$tz_name)
  ) %>%
  dplyr::select(timestamp, blh)

# Add to main data frame
site <- site %>%
  dplyr::left_join(era, by = "timestamp") %>%
  # Interpolate hourly -> half-hourly
  dplyr::mutate(blh = approx(seq_along(blh), blh, seq_along(blh))$y)


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
  fetch <- dplyr::rename_with(
    fetch, ~ stringr::str_c(., tolower(control$fetch_model))
  )
  
  fetch <- dplyr::bind_cols(dplyr::select(site, timestamp), fetch)
  
  fetch_out <- file.path(
    path_out, paste0("fetch_", tag_out, ".csv")
  )
  write.csv(fetch, fetch_out, row.names = FALSE)
  
  # Create documentation for fetch output
  fetch_docu <- settings
  fetch_docu$files <- c(site_input)
  fetch_docu$method <- control$fetch_model
  fetch_docu_out <- stringr::str_replace(fetch_out, ".csv", ".txt")
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
    dplyr::select(timestamp, wd, ustar, mo_length, v_sigma) %>%
    tidyr::drop_na()
  
  # Set up grid
  grid <- grid_init(fetch = (md$tower_height - md$displacement) * 100)
  
  # Convert AOI to grid
  aoi_grid <- aoi_to_grid(aoi, grid, c(md$x_utm, md$y_utm))
  
  # Calculate footprints, create raster images, write to file
  # WARNING: this takes ~1 hr to run
  
  # Initialize loop
  phi_hh <- vector("list", length = nrow(site_fp))
  p <- progress_info(nrow(site_fp))
  
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
    
    p$tick()
    
  }
  
  
  # Write file with just footprint data so this doesn't have to be run again
  fp_basics <- phi_hh %>% 
    tibble::enframe(name = "timestamp", value = "phi") %>% 
    dplyr::mutate(
      timestamp = lubridate::ymd_hms(timestamp, tz = md$tz_name)
    ) %>%
    tidyr::unnest(phi) %>%
    dplyr::right_join(dplyr::select(site, timestamp))
  
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
  fp_basics_docu_out <- stringr::str_replace(fp_basics_out, ".csv", ".txt")
  # Save documentation
  sink(fp_basics_docu_out)
  fp_basics_docu 
  sink()
}


### Calculate 2-D footprints - FFP method ======================================

if (control$ffp) {
  # See Kljun et al. 2015
  
  # Create new folder for output files
  ffp_out <- file.path(path_out, paste0("half_hourly_ffp_", tag_out))
  dir.create(ffp_out, showWarnings = FALSE)
  
  # Create sub-folders for grid and footprint matrices
  dir.create(file.path(ffp_out, "grid"), showWarnings = FALSE)
  dir.create(file.path(ffp_out, "footprints"), showWarnings = FALSE)
  
  # Select only necessary variables, remove missing data
  site_ffp <- site %>% 
    dplyr::select(timestamp, wd, ustar, mo_length, v_sigma, blh) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      zd = md$displacement, 
      zo = md$roughness_length
    )
  # Pre-filter known invalid cases (not using - just calculate all fp then see)
  # - NOTE: for ustar Kljun only indicate > 0.1, but adding upper limit of 2.0
  # - Olson et al. 2004 (FLUXNET) suggest 6.0, but here there are none >1 & <2
  # - this suggests spurious values above 2
  #filter(
  #  between(ustar, 0.1, 2), 
  #  ((md$tower_height - md$displacement) / mo_length) >= -15.5
  #)
  
  # Set up grid
  grid <- grid_init(fetch = (md$tower_height - md$displacement) * 100)
  
  # Convert AOI to grid
  aoi_grid <- aoi_to_grid(aoi, grid, c(md$x_utm, md$y_utm))
  
  # Extent of AOI for trimming footprint grid
  extent_trim <- get_trim_extent(aoi_grid)
  
  # Trim grid to AOI area
  # - no need to calculate footprint weights outside of AOI 
  # - these are implicit in phi
  grid <- purrr::map(grid, trim_matrix, extent_trim)
  
  # Trim AOI
  aoi_grid <- trim_matrix(aoi_grid)
  
  # Write grid to file
  write_matrix(grid$x, file.path(ffp_out, "grid", "x"), trunc = NA)
  write_matrix(grid$y, file.path(ffp_out, "grid", "y"), trunc = NA)
  write_matrix(aoi_grid, file.path(ffp_out, "grid", "aoi"), trunc = NA)
  
  # Calculate footprints, create raster images, write to file
  
  # Initialize loop
  ffp_phi_hh <- vector("list", length = nrow(site_ffp))
  p <- progress_info(nrow(site_ffp))
  
  # Calculate footprint matrices
  for (i in 1:nrow(site_ffp)) {
    
    # Calculate footprint, add to list
    ffp_temp <- calc_fp_kljun(
      grid = grid,
      wd = site_ffp$wd[i],
      ustar = site_ffp$ustar[i],
      mo_length = site_ffp$mo_length[i],
      v_sigma = site_ffp$v_sigma[i],
      blh = site_ffp$blh[i],
      z = md$tower_height,
      zd = site_ffp$zd[i],
      zo = site_ffp$zo[i]
    )
    
    # Calculate AOI coverage
    ffp_phi_hh[i] <- sum(ffp_temp * aoi_grid)
    
    # Set output name as the timestamp
    ffp_temp_nm <- format(site_ffp$timestamp[i], "%Y-%m-%d-%H%M%S")
    names(ffp_phi_hh)[i] <- ffp_temp_nm
    
    # Write to file
    # - values are truncated to enable efficient storage as integers
    # - this can be suppressed by setting trunc = NA in write_matrix() 
    write_matrix(ffp_temp, file.path(ffp_out, "footprints", ffp_temp_nm))
    
    p$tick()
  }
  
  ffp_dir_size <- file.path(ffp_out, "footprints") %>% 
    list.files(full.names = TRUE) %>%
    magrittr::extract(1:50) %>%
    purrr::map_dbl(~ file.info(.)$size) %>% 
    mean() %>%
    magrittr::multiply_by(nrow(site_ffp)) %>%
    magrittr::divide_by(1e9) %>%
    round(2)
  
  cat(
    "\nFinished calculating footprints.\n", 
    "Matrix files can be found in:", file.path(ffp_out, "footprints"), "\n",
    "Total size of directory:", ffp_dir_size, "GB\n"
  )
  cat("Writing summary output...")
  
  # Write file with just footprint data so this doesn't have to be run again
  ffp_basics <- ffp_phi_hh %>% 
    tibble::enframe(name = "timestamp", value = "phi_ffp") %>% 
    dplyr::mutate(
      timestamp = lubridate::ymd_hms(timestamp, tz = md$tz_name)
    ) %>%
    tidyr::unnest(phi_ffp) %>%
    dplyr::right_join(dplyr::select(site, timestamp, ustar, mo_length)) %>%
    # Classify phi according to Gockede et al. 2008
    dplyr::mutate(
      phi_ffp_class = dplyr::case_when(
        # Homogenous measurements
        dplyr::between(phi_ffp, 0.95, 1.00) ~ 0L,
        # Representative measurements
        dplyr::between(phi_ffp, 0.80, 0.95) ~ 1L,
        # Acceptable measurements
        dplyr::between(phi_ffp, 0.50, 0.80) ~ 2L,
        # Disturbed measurements
        dplyr::between(phi_ffp, 0.00, 0.50) ~ 3L
      ),
      zeta = ((md$tower_height - md$displacement) / mo_length),
      ffp_valid = dplyr::if_else(
        magrittr::and(dplyr::between(ustar, 0.1, 2), zeta >= -15.5), 1L, 0L
      )
    ) %>% 
    dplyr::select(-ustar, -mo_length, -zeta)
  
  ffp_basics_out <- file.path(
    path_out, paste0("footprint_basics_ffp_", tag_out, ".csv")
  )
  write.csv(ffp_basics, ffp_basics_out, row.names = FALSE)
  
  # Create documentation for phi output
  ffp_basics_docu <- settings
  ffp_basics_docu <- append(
    settings, list(
      files = c(site_input, era_input, aoi_input),
      model_params = list(
        method = "K15", fetch = attr(grid, "fetch"), res = attr(grid, "res")
      )
    )
  )
  ffp_basics_docu_out <- stringr::str_replace(ffp_basics_out, ".csv", ".txt")
  # Save documentation
  sink(ffp_basics_docu_out)
  ffp_basics_docu 
  sink()
  
  cat("done.\n")
}


