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
devtools::load_all("/Users/Graham/R Projects/footprints", quiet = TRUE)
library(progress)
library(lubridate, warn.conflicts = FALSE)
library(tidyverse)

# Load reference files
source("/Users/Graham/Desktop/DATA/Flux/tools/reference/site_metadata.R")

# Load functions
path_funs <- file.path(settings$dir, "Flux/tools/engine/functions")
source(file.path(path_funs, "dates_and_times.R"))
source(file.path(path_funs, "latest_version.R"))
source(file.path(path_funs, "utilities.R"))


### Helper functions ===========================================================

progress_info <- function(len) {
  
  progress_bar$new(
    total = len, clear = FALSE,
    format = paste0(
      "[:spin] Completed: :current | :percent  ", 
      "Elapsed: :elapsed  Remaining: :eta"
    )
  )
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path(settings$dir, "Flux", settings$site, settings$year)

paths <- list(
  # Input file - gap-filled biomet
  data = latest_version(
    file.path(wd, "processing", "04_biomet_gapfill", "data")
  ),
  # Input file - processed ERA data for site location
  era = latest_version(
    file.path(dirname(dirname(wd)), "JLL", "all", "era"), "era_proc"
  ),
  # Input file - NOAA boundary layer height (/scripts/download_noaa_blh.R)
  noaa = file.path(
    dirname(dirname(wd)), "JLL", "all", "noaa", "noaa_blh_hh_interp.csv"
  ),
  # Site delineation polygon
  delin = file.path(dirname(wd), "site_info", "delineation"),
  # Output files
  out = file.path(wd, "processing", "05_footprint")
)
  
# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)


### Load required input data ===================================================

cat("Importing data files.\n")

# Import site delineation polygon
delin <- sf::read_sf(paths$delin)

# Load the data
data <- readr::read_csv(
  paths$data, guess_max = 7000, 
  col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
)
# Force local time zone to align properly with external data
data <- dplyr::mutate(
  data, timestamp = lubridate::force_tz(timestamp, md$tz_name)
)

# Import ERA data, get variables
era <- readr::read_csv(
  paths$era, col_types = readr::cols(.default = readr::col_guess()), 
  progress = FALSE
)
era <- era %>%
  dplyr::mutate(timestamp = lubridate::with_tz(timestamp, md$tz_name)) %>%
  dplyr::select(timestamp, blh_era = blh)

# Import NOAA data
noaa <- readr::read_csv(
  paths$noaa, col_types = readr::cols(.default = readr::col_guess()), 
  progress = FALSE
)
noaa <- noaa %>%
  dplyr::mutate(timestamp = lubridate::with_tz(timestamp, md$tz_name)) %>%
  dplyr::select(timestamp, blh_noaa = blh)

# Add external data to main data frame
data <- data %>%
  # Roughness parameters
  dplyr::mutate(
    z = md$tower_height,
    zd = md$displacement,
    z0 = md$roughness_length
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
  
  # Select model function and parameters
  fetch_ref <- list(
    K15 = list(
      fun = calc_fetch_kljun, 
      vars = c("ustar", "mo_length", "blh", "z0")
    )
  )
  
  fetch_ref <- purrr::pluck(fetch_ref, control$fetch_model)
  
  # Add WS to model vars if specified
  if (control$fetch_model == "K15" & control$use_ws) {
    #control$fp_model <- "K15_ws"
    fetch_ref <- purrr::list_merge(fetch_ref, vars = "ws")
  }
  
  # Select only necessary variables, remove missing data
  data_fetch <- data %>% 
    # Filter out erroneous wind records
    dplyr::filter(qc_ws == 0) %>%
    dplyr::select(timestamp, z, zd, dplyr::all_of(fetch_ref$vars)) %>%
    tidyr::drop_na()
  
  # Copy data into list (easier to access in loop)
  data_fetch_list <- data_fetch %>%
    dplyr::select(-timestamp) %>%
    as.list()
  
  # Initialize loop for fetch calculation
  n_fetch <- nrow(data_fetch)
  fetch_list <- vctrs::vec_init(list(), n_fetch)
  p <- progress_info(n_fetch)
  
  # Apply one-dimensional footprint model to data
  for (i in 1:n_fetch) {
    
    # Get parameters from data list
    fetch_args <- data_fetch_list %>% 
      purrr::modify(i) %>% 
      # Splice in percents
      purrr::prepend(list(pct = control$fetch_pct, dx = control$fetch_dx))
    
    # Calculate fetch lengths
    fetch_list[[i]] <- rlang::exec(fetch_ref$fun, !!!fetch_args)
    
    p$tick()
  }
  
  fetch <- fetch_list %>%
    dplyr::bind_rows() %>%
    dplyr::bind_cols(dplyr::select(data_fetch, timestamp), .) %>%
    dplyr::right_join(dplyr::select(data, timestamp), by = "timestamp") %>%
    dplyr::arrange(timestamp)
  
  fetch_out <- file.path(
    paths$out, "fetch", 
    paste0("fetch_", control$fetch_model, "_", tag_out, ".csv")
  )
  readr::write_csv(fetch, fetch_out)
  
  # Create documentation for fetch output
  fetch_docu <- settings
  fetch_docu <- append(
    settings, list(
      files = c(paths$data, paths$era),
      model_params = list(
        model = control$fetch_model, 
        roughness = if ("z0" %in% fetch_ref$vars) "z0" else "ws/ustar",
        percents = control$fetch_pct
      )
    )
  )
  fetch_docu_out <- stringr::str_replace(fetch_out, ".csv", ".txt")
  # Save documentation
  sink(fetch_docu_out)
  print(fetch_docu) 
  sink()
  
  cat("\n")
}


### Calculate 2-D footprints and coverage ======================================

if (control$fp) {
  
  # Select model function and parameters
  model_ref <- list(
    H00 = list(
      fun = calc_footprint_hsieh, 
      vars = c("ustar", "mo_length", "v_sigma", "z0")
    ),
    KM01 = list(
      fun = calc_footprint_kormann, 
      vars = c("ws", "ustar", "mo_length", "v_sigma")
    ),
    K15 = list(
      fun = calc_footprint_kljun, 
      vars = c("ustar", "mo_length", "v_sigma", "blh", "z0")
    )
  )
  
  model_ref <- purrr::pluck(model_ref, control$fp_model)
  
  # Add WS to model vars if specified
  if (control$fp_model == "K15" & control$use_ws) {
    #control$fp_model <- "K15_ws"
    model_ref <- purrr::list_merge(model_ref, vars = "ws")
  }
  
  fp_fun <- purrr::pluck(model_ref, "fun")
  
  # (See references for description of models)
  
  # Create new folder for output files
  fp_out <- file.path(
    paths$out, "footprint", 
    paste0("halfhourly_", control$fp_model, "_", tag_out)
  )
  dir.create(fp_out, showWarnings = FALSE)
  
  # Create sub-folders for grid and footprint matrices
  dir.create(file.path(fp_out, "grid"), showWarnings = FALSE)
  dir.create(file.path(fp_out, "footprints"), showWarnings = FALSE)
  
  # Select only necessary variables, remove missing data
  data_fp <- data %>% 
    # Filter out erroneous wind records
    dplyr::filter(qc_ws == 0, qc_wd == 0) %>%
    # All models need z, zd, and wd
    dplyr::select(timestamp, z, zd, wd, dplyr::all_of(model_ref$vars)) %>%
    tidyr::drop_na()
  
  # Copy data into list (easier to access in footprint loop)
  data_fp_list <- data_fp %>%
    dplyr::select(-timestamp) %>%
    as.list()
  
  # Set up grid
  grid <- aoi_to_grid(delin, c(md$x_utm, md$y_utm), delta = 1)
  
  # Convert AOI to grid
  aoi_grid <- aoi_to_mask(delin, c(md$x_utm, md$y_utm), delta = 1)
  
  # Write grid to file
  purrr::imap(
    purrr::list_modify(grid, aoi = aoi_grid), 
    ~ write_matrix(
      .x, file.path(fp_out, "grid", .y), trunc = 0, compress = FALSE
    )
  )
  
  # Initialize loop for footprint calculation
  n_fp <- nrow(data_fp)
  phi <- vector("list", length = n_fp)
  fp_topo <- matrix(0, nrow = nrow(aoi_grid), ncol = ncol(aoi_grid))
  n_topo <- 0
  p <- progress_info(n_fp)
  
  cat("Calculating half-hourly footprints.\n")
  
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
    write_matrix(
      fp_temp, file.path(fp_out, "footprints", fp_temp_nm), trunc = 9
    )
    
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
    paths$out, "topology", paste0("all_", control$fp_model, "_", tag_out)
  )
  fp_topo <- fp_topo / n_topo
  write_matrix(fp_topo, topo_out, trunc = 0, compress = FALSE)
  
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
  if (stringr::str_detect(control$fp_model, "K15")) {
    
    # Indicate whether footprint is valid according to Kljun et al. 2015
    fp_stats <- fp_stats %>%
      dplyr::right_join(
        dplyr::select(data, timestamp, ustar, mo_length, blh, ws, z, zd, z0), 
        by = "timestamp"
      )
    
    if (control$use_ws) {
      fp_stats <- dplyr::mutate(
        fp_stats, z0 = calc_z0(ws, ustar, mo_length, z - zd)
      )
    }
    
    fp_stats <- fp_stats %>%
      dplyr::arrange(timestamp) %>%
      dplyr::mutate(
        zm = z - zd,
        # TODO fix this
        k15_valid = dplyr::case_when(
          ustar < 0.1 ~ 0L,
          zm / mo_length > -15.5 ~ 0L,
          # Using roughness sublayer limit as in latest FFP version, not paper
          12.5 * z0 >= zm ~ 0L,
          zm > 0.8 * blh ~ 0L,
          is.na(phi_k15) ~ NA_integer_,
          TRUE ~ 1L
        )
      ) %>%
      dplyr::select(-ustar, -mo_length, -ws)
  }
  
  if (control$fp_model == "KM01") {
    # Indicate whether footprint is valid according to Kormann & Meixner 2001
    fp_stats <- fp_stats %>%
      dplyr::right_join(
        dplyr::select(data, timestamp, mo_length, z, zd), by = "timestamp"
      ) %>%
      dplyr::arrange(timestamp) %>%
      dplyr::mutate(
        zm = z - zd,
        km01_valid = dplyr::if_else(abs(zm / mo_length) > 3, 0L, 1L)
      ) %>%
      dplyr::select(-z, -zm, -zd)
  }
  
  fp_stats_out <- file.path(
    paths$out, "stats", 
    paste0("footprint_stats_", control$fp_model, "_", tag_out, ".csv")
  )
  readr::write_csv(fp_stats, fp_stats_out)
  
  # Create documentation for phi output
  fp_stats_docu <- settings
  fp_stats_docu <- append(
    settings, list(
      files = c(paths$data, paths$era, paths$delin),
      model_params = list(
        model = control$fp_model, 
        roughness = if ("z0" %in% model_ref$vars) "z0" else "ws/ustar",
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

