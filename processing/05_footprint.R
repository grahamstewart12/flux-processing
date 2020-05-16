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

# Hsieh, C.-I., Katul, G., & Chi, T. (2000). An approximate analytical model for
# footprint estimation of scalar fluxes in thermally stratified atmospheric 
# flows. Advances in Water Resources, 23(7), 765–772. 
# https://doi.org/10.1016/S0309-1708(99)00042-1

# Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple 
# two-dimensional parameterisation for Flux Footprint Prediction (FFP). 
# Geoscientific Model Development, 8(11), 3695–3713. 
# https://doi.org/10.5194/gmd-8-3695-2015

# Kormann, R., & Meixner, F. X. (2001). An Analytical Footprint Model For 
# Non-Neutral Stratification. Boundary-Layer Meteorology, 99(2), 207–224. 
# https://doi.org/10.1023/A:1018991015119

# Leclerc, M. Y., & Foken, T. (2014). Footprints in Micrometeorology and 
# Ecology. Berlin Heidelberg: Springer-Verlag. 
# https://doi.org/10.1007/978-3-642-54545-0



# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/R Projects/footprints", quiet = TRUE)
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
    zm = md$tower_height,
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
  dplyr::mutate(blh = blh_era) %>%
  # Force time zone back to UTC
  dplyr::mutate(timestamp = lubridate::force_tz(timestamp, "UTC"))


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
  
  # Select model function and parameters
  model_ref <- list(
    H00 = list(
      fun = calc_footprint_hsieh, 
      vars = c("ustar", "mo_length", "v_sigma", "zo")
    ),
    KM01 = list(
      fun = calc_footprint_kormann, 
      vars = c("ws", "ustar", "mo_length", "v_sigma")
    ),
    K15 = list(
      fun = calc_footprint_kljun, 
      vars = c("ustar", "mo_length", "v_sigma", "blh", "zo")
    )
  )
  
  model_ref <- purrr::pluck(model_ref, control$fp_model)
  
  # Add WS to model vars if specified
  if (control$fp_model == "K15" & control$use_ws) {
    model_ref <- purrr::list_merge(model_ref, vars = "ws")
  }
  
  fp_fun <- purrr::pluck(model_ref, "fun")
  
  # (See references for description of models)
  
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
    # All models need zm, zd, and wd
    dplyr::select(timestamp, zm, zd, wd, dplyr::all_of(model_ref$vars)) %>%
    tidyr::drop_na()
  
  # Copy data into list (easier to access in footprint loop)
  data_fp_list <- data_fp %>%
    dplyr::select(-timestamp) %>%
    as.list()
  
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
    ~ write_matrix(
      .x, file.path(fp_out, "grid", .y), trunc = NA, compress = FALSE
    )
  )
  
  # Initialize loop for footprint calculation
  n_fp <- nrow(data_fp)
  phi <- vector("list", length = n_fp)
  fp_topo <- matrix(0, nrow = nrow(aoi_grid), ncol = ncol(aoi_grid))
  n_topo <- 0
  p <- progress_info(n_fp)
  
  # Calculate footprint matrices
  for (i in 1:n_fp) {
    
    # Get parameters from data list
    fp_args <- data_fp_list %>% 
      purrr::modify(i) %>% 
      # Splice in grid
      purrr::prepend(list(grid = grid))
    
    # Calculate footprint
    fp_temp <- rlang::exec(model_ref$fun, !!!fp_args)
    
    # Calculate AOI coverage
    phi[i] <- sum(fp_temp * aoi_grid)
    
    # Set output name as the timestamp
    fp_temp_nm <- format(data_fp$timestamp[i], "%Y-%m-%d-%H%M%S")
    names(phi)[i] <- fp_temp_nm
    
    # Write to file
    # - values are truncated to enable efficient storage as integers
    # - this can be suppressed by setting trunc = NA in write_matrix() 
    write_matrix(fp_temp, file.path(fp_out, "footprints", fp_temp_nm))
    
    # Add to footprint topology (unless it contains missing values)
    if (!any(is.na(fp_temp))) {
      fp_topo <- fp_topo + fp_temp
      n_topo <- n_topo + 1
    } 
    
    p$tick()
  }
  
  # Estimate total size of footprint output directory 
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
  
  # Calculate topology and write to file
  topo_out <- file.path(
    path_out, "topology", paste0("all_", control$fp_model, "_", tag_out)
  )
  fp_topo <- fp_topo / n_topo
  write_matrix(fp_topo, topo_out, trunc = NA, compress = FALSE)
  
  # Write file with footprint stats
  fp_stats <- phi %>% 
    tibble::enframe(name = "timestamp", value = "phi") %>% 
    dplyr::mutate(timestamp = lubridate::ymd_hms(timestamp)) %>%
    tidyr::unnest(phi) %>%
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
      )
    ) %>% 
    # Give phi vars suffix with model tag
    dplyr::rename_with(
      ~ stringr::str_c(.x, "_", tolower(control$fp_model)), c(phi, phi_class)
    )
  
  # Compute model-specific validity flags
  if (control$fp_model == "K15") {
    # Indicate whether footprint is valid according to Kljun et al. 2015
    fp_stats <- fp_stats %>%
      dplyr::right_join(
        dplyr::select(data, timestamp, ustar, mo_length, blh, zm, zd, zo), 
        by = "timestamp"
      ) %>%
      dplyr::arrange(timestamp) %>%
      dplyr::mutate(
        z = zm - zd,
        k15_valid = dplyr::case_when(
          ustar < 0.1 ~ 0L,
          z / mo_length < -15.5 ~ 0L,
          20 * zo >= z ~ 0L,
          z < 0.8 * blh ~ 0L,
          is.na(phi_k15) ~ NA_integer_,
          TRUE ~ 1L
        )
      ) %>%
      dplyr::select(-z, -zm, -zd, -zo, -ustar, -mo_length, -blh)
  }
  
  if (control$fp_model == "KM01") {
    # Indicate whether footprint is valid according to Kormann & Meixner 2001
    fp_stats <- fp_stats %>%
      dplyr::right_join(
        dplyr::select(data, timestamp, mo_length, zm, zd), by = "timestamp"
      ) %>%
      dplyr::arrange(timestamp) %>%
      dplyr::mutate(
        z = zm - zd,
        km01_valid = dplyr::if_else(abs(z / mo_length) > 3, 0L, 1L)
      ) %>%
      dplyr::select(-z, -zm, -zd, -mo_length)
  }
  
  fp_stats_out <- file.path(
    path_out, "stats", 
    paste0("footprint_stats_", control$fp_model, "_", tag_out, ".csv")
  )
  readr::write_csv(fp_stats, fp_stats_out)
  
  # Create documentation for phi output
  fp_stats_docu <- settings
  fp_stats_docu <- append(
    settings, list(
      files = c(data_input, era_input, aoi_input),
      model_params = list(
        model = control$fp_model, 
        roughness = if ("zo" %in% model_ref$vars) "zo" else "ws/ustar",
        fetch = attr(grid, "fetch"), 
        res = attr(grid, "res")
      )
    )
  )
  fp_stats_docu_out <- stringr::str_replace(fp_stats_out, ".csv", ".txt")
  # Save documentation
  sink(fp_stats_docu_out)
  print(fp_stats_docu) 
  sink()
  
  cat("done.\n")
}

