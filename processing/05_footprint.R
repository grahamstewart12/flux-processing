### ============================================================================
# Footprint modeling ===========================================================
### ============================================================================

# Purpose: 

# References:

# Göckede, M., Foken, T., Aubinet, M., Aurela, M. A., Banza, J., Bernhofer, C., 
# et al. (2008). Quality control of CarboEurope flux data - Part 1: Coupling 
# footprint analyses with flux data quality assessment to evaluate sites in 
# forest ecosystems. Biogeosciences, 5(2), 433–450. 
# https://doi.org/10.5194/bg-5-433-2008

# Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple 
# two-dimensional parameterisation for Flux Footprint Prediction (FFP). 
# Geoscientific Model Development, 8(11), 3695–3713. 
# https://doi.org/10.5194/gmd-8-3695-2015


# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/R Projects/footprints")
library(progress)
library(lubridate, warn.conflicts = FALSE)
library(tidyverse)

# Load reference files
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")

# Load functions
source("~/Desktop/DATA/Flux/tools/engine/functions/dates_and_times.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/latest_version.R")
source("~/Desktop/DATA/Flux/tools/engine/functions/utilities.R")


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
path_in <- file.path(wd, "processing", "04_biomet_gapfill", "data")
  
# Input file - gap-filled biomet
data_input <- latest_version(path_in)
  
# Input file - processed ERA data for site location
era_input <- latest_version(
  file.path("~/Desktop", "DATA", "Flux", "JLL", "all", "era"), "era_proc"
)
# Input file - NOAA boundary layer height (/scripts/download_noaa_blh.R)
noaa_input <- file.path(
  "~/Desktop/DATA/Flux", "JLL", "all", "noaa", "noaa_blh_hh_interp.csv"
)
  
# AOI
aoi_input <- file.path(dirname(wd), paste0(settings$site, "_area_proc"))
  
# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)
  
# Set path for output files
path_out <- file.path(wd, "processing", "05_footprint")


### Load required input data ===================================================

cat("Importing data files.\n")

# Load the data
data <- readr::read_csv(
  data_input, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
)
# Force local time zone to align properly with external data
data <- dplyr::mutate(
  data, timestamp = lubridate::force_tz(timestamp, md$tz_name)
)

# Import site AOI polygon
aoi <- rgdal::readOGR(aoi_input, verbose = FALSE)

# Import ERA data, get variables
era <- readr::read_csv(
  era_input, col_types = readr::cols(.default = readr::col_guess()), 
  progress = FALSE
)
era <- era %>%
  dplyr::mutate(timestamp = lubridate::with_tz(timestamp, md$tz_name)) %>%
  dplyr::select(timestamp, blh_era = blh)

# Import NOAA data
noaa <- readr::read_csv(
  noaa_input, col_types = readr::cols(.default = readr::col_guess()), 
  progress = FALSE
)
noaa <- noaa %>%
  dplyr::mutate(timestamp = lubridate::with_tz(timestamp, md$tz_name)) %>%
  dplyr::select(timestamp, blh_noaa = blh)

# Add external data to main data frame
data <- data %>%
  # Roughness parameters
  dplyr::mutate(
    zd = md$displacement,
    zo = md$roughness_length
  ) %>%
  # ERA
  dplyr::left_join(era, by = "timestamp") %>%
  # Interpolate hourly -> half-hourly
  dplyr::mutate(blh_era = imputeTS::na_interpolation(blh_era)) %>%
  # NOAA (already interpolated)
  dplyr::left_join(noaa, by = "timestamp") %>%
  # Set primary BLH variable (ERA, since at finer spatial & temporal res)
  dplyr::mutate(blh = blh_era)


### Calculate 1-D footprints (fetch lengths) ===================================

if (control$fetch) {
  
  cat("Calculating half-hourly fetch lengths.\n")
  
  # Calculate half-hourly fetch lengths
  fetch <- fp_fetch(
    data, ws = "ws", ustar = "ustar", zeta = "zl", mo_length = "mo_length",
    z = md$tower_height, zd = md$displacement, zo = md$roughness_length, 
    percent = c(10, 30, 50, 70, 90, "peak"), method = control$fetch_model
  )
  
  # Rename with method ID tag
  fetch <- dplyr::rename_with(
    fetch, ~ stringr::str_c(.x, tolower(control$fetch_model))
  )
  
  fetch <- dplyr::bind_cols(dplyr::select(data, timestamp), fetch)
  
  fetch_out <- file.path(path_out, "fetch", paste0("fetch_", tag_out, ".csv"))
  readr::write_csv(fetch, fetch_out)
  
  # Create documentation for fetch output
  fetch_docu <- settings
  fetch_docu$files <- c(data_input)
  fetch_docu$method <- control$fetch_model
  fetch_docu_out <- stringr::str_replace(fetch_out, ".csv", ".txt")
  # Save documentation
  sink(fetch_docu_out)
  print(fetch_docu) 
  sink()
}


### Calculate 2-D footprints and coverage ======================================

if (control$fp) {
  
  # Select model function
  model_ref <- list(
    # H00 doesn't need any extra vars; timestamp satisfies dplyr::all_of()
    H00 = list(fun = calc_footprint_hsieh, vars = "timestamp"),
    KM01 = list(fun = calc_footprint_kormann, vars = "ws"),
    K15 = list(fun = calc_footprint_kljun, vars = "blh")
  )
  fp_fun <- purrr::pluck(model_ref, control$fp_model, "fun")
  
  # See references for description of models
  
  # Create new folder for output files
  fp_out <- file.path(
    path_out, "footprint", paste0("halfhourly_", control$fp_model, "_", tag_out)
  )
  dir.create(fp_out, showWarnings = FALSE)
  
  # Create sub-folders for grid and footprint matrices
  dir.create(file.path(fp_out, "grid"), showWarnings = FALSE)
  dir.create(file.path(fp_out, "footprints"), showWarnings = FALSE)
  
  # Select only necessary variables, remove missing data
  data_fp <- data %>% 
    # Filter out erroneous wind records
    dplyr::filter(qc_ws == 0, qc_wd == 0) %>%
    dplyr::select(
      timestamp, wd, ustar, mo_length, v_sigma, zd, zo,
      # Add model-specific vars
      dplyr::all_of(purrr::pluck(model_ref, control$fp_model, "vars"))
    ) %>%
    tidyr::drop_na()
  
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
  purrr::imap(
    purrr::list_modify(grid, aoi = aoi_grid), 
    ~ write_matrix(.x, file.path(fp_out, "grid", .y), trunc = NA)
  )
  
  # Calculate footprints, create raster images, write to file
  
  # Initialize loop
  phi <- vector("list", length = nrow(data_fp))
  p <- progress_info(nrow(data_fp))
  
  # Calculate footprint matrices
  for (i in 1:nrow(data_fp)) {
    
    # Calculate footprint, add to list
    fp_args <- list(
      grid = grid,
      wd = data_fp$wd[i],
      ustar = data_fp$ustar[i],
      mo_length = data_fp$mo_length[i],
      v_sigma = data_fp$v_sigma[i],
      blh = data_fp$blh[i],
      z = md$tower_height,
      zd = data_fp$zd[i],
      zo = data_fp$zo[i]
    )
    fp_temp <- rlang::exec(fp_fun, !!!fp_args)
    
    #fp_temp <- calc_fp_kljun(
    #  grid = grid,
    #  wd = site_ffp$wd[i],
    #  ustar = site_ffp$ustar[i],
    #  mo_length = site_ffp$mo_length[i],
    #  v_sigma = site_ffp$v_sigma[i],
    #  blh = site_ffp$blh[i],
    #  z = md$tower_height,
    #  zd = site_ffp$zd[i],
    #  zo = site_ffp$zo[i]
    #)
    
    # Calculate AOI coverage
    phi[i] <- sum(fp_temp * aoi_grid)
    
    # Set output name as the timestamp
    fp_temp_nm <- format(data_fp$timestamp[i], "%Y-%m-%d-%H%M%S")
    names(phi)[i] <- fp_temp_nm
    
    # Write to file
    # - values are truncated to enable efficient storage as integers
    # - this can be suppressed by setting trunc = NA in write_matrix() 
    write_matrix(fp_temp, file.path(fp_out, "footprints", fp_temp_nm))
    
    # if (!all(is.na(fp_temp))) fp_topo <- fp_topo + fp_temp
    
    p$tick()
  }
  
  dir_size <- file.path(fp_out, "footprints") %>% 
    list.files(full.names = TRUE) %>%
    magrittr::extract(1:50) %>%
    purrr::map_dbl(file.size) %>% 
    mean() %>%
    magrittr::multiply_by(nrow(data_fp)) %>%
    magrittr::divide_by(1e9) %>%
    round(2)
  
  cat(
    "\nFinished calculating footprints.\n", 
    "Matrix files can be found in:", file.path(fp_out, "footprints"), "\n",
    "Total size of directory:", dir_size, "GB\n"
  )
  cat("Writing summary output...")
  
  # Write file with just footprint data so this doesn't have to be run again
  fp_basics <- phi %>% 
    tibble::enframe(name = "timestamp", value = "phi") %>% 
    dplyr::mutate(
      timestamp = lubridate::ymd_hms(timestamp, tz = md$tz_name)
    ) %>%
    tidyr::unnest(phi) %>%
    dplyr::right_join(dplyr::select(data, timestamp, ustar, mo_length)) %>%
    # Classify phi according to Gockede et al. 2008
    dplyr::mutate(
      phi_class = dplyr::case_when(
        # Homogenous measurements
        dplyr::between(phi, 0.95, 1.00) ~ 0L,
        # Representative measurements
        dplyr::between(phi, 0.80, 0.95) ~ 1L,
        # Acceptable measurements
        dplyr::between(phi, 0.50, 0.80) ~ 2L,
        # Disturbed measurements
        dplyr::between(phi, 0.00, 0.50) ~ 3L
      ),
      # TODO allow for dynamic displacement here
      zeta = ((md$tower_height - md$displacement) / mo_length),
      # Indicate whether footprint is valid according to Kljun et al. 2015
      # - NOTE: for ustar Kljun only indicate > 0.1, but adding 2.0 upper limit
      # - Olson et al. 2004 (FLUXNET) suggest 6.0, but none >1 & <2
      # - this suggests spurious values above 2
      ffp_valid = dplyr::if_else(
        magrittr::and(dplyr::between(ustar, 0.1, 2), zeta >= -15.5), 1L, 0L
      )
    ) %>% 
    dplyr::select(-ustar, -mo_length, -zeta) %>%
    # Give phi vars suffix with model tag
    dplyr::rename_with(
      ~ stringr::str_c(.x, "_", tolower(control$fp_model)), c(phi, phi_class)
    )
  
  fp_basics_out <- file.path(
    path_out, "basics", 
    paste0("footprint_basics_", control$fp_model, "_", tag_out, ".csv")
  )
  readr::write_csv(fp_basics, fp_basics_out)
  
  # Create documentation for phi output
  fp_basics_docu <- settings
  fp_basics_docu <- append(
    settings, list(
      files = c(data_input, era_input, aoi_input),
      model_params = list(
        method = "K15", fetch = attr(grid, "fetch"), res = attr(grid, "res")
      )
    )
  )
  fp_basics_docu_out <- stringr::str_replace(fp_basics_out, ".csv", ".txt")
  # Save documentation
  sink(fp_basics_docu_out)
  print(fp_basics_docu) 
  sink()
  
  cat("done.\n")
}

# TODO footprint average module
# - to save memory, do this by cumulatively summing matrices
# - and then divide by n at the end
# - keep track of empty or NA matrices and don't add these

