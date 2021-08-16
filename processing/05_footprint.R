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
devtools::load_all("~/Projects/Flux/footprints", quiet = TRUE)
devtools::load_all("~/Projects/Flux/dscalr", quiet = TRUE)
library(progress)
library(lubridate, warn.conflicts = FALSE)
library(tidyverse)


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
md <- yaml::read_yaml(file.path("data", settings$site, "metadata.yml"))

# Set tag for file names
tag <- make_tag(md$site_code, settings$year)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path("data", settings$site, settings$year)
path_in <- file.path(wd, "04_biomet_gapfill", "data")
path_out <- file.path(wd, "05_footprint")

# Input file - gap-filled biomet
file_data <- file.path(wd, path_in, paste0("biomet_gf_", tag, ".csv"))

# Input file - site delineation polygon
file_delin <- file.path(dirname(wd), "site_info", "delineation")

# Input file - processed ERA data for site location
era_files <- list.files(
  "data/ERA", pattern = "era_proc.+.csv", recursive = TRUE, full.names = TRUE
)

# Load config for this script
config <- yaml::read_yaml(file.path(path_out, "config.yml"))


### Load required input data ===================================================

cat("Importing data files.\n")

# Import site delineation polygon
delin <- sf::read_sf(file_delin)

# Load the data
data <- read_csv(file_data, show_col_types = FALSE, progress = FALSE)
# Force local time zone to align properly with external data
data <- mutate(data, timestamp = force_tz(timestamp, md$tz$name))

# Load the ERA data library, get variables
era <- era_files %>%
  read_csv(show_col_types = FALSE, progress = FALSE) %>%
  mutate(timestamp = with_tz(timestamp, md$tz$name)) %>%
  select(timestamp, blh_era = blh)

# Add external data to main data frame
data <- data %>%
  # Roughness parameters
  mutate(
    z = md$tower$height,
    zd = md$displacement,
    z0 = md$roughness_length
  ) %>%
  # ERA
  left_join(era, by = "timestamp") %>%
  # Interpolate hourly -> half-hourly
  mutate(blh = imputeTS::na_interpolation(blh_era)) %>%
  # Force time zone back to UTC
  mutate(timestamp = force_tz(timestamp, "UTC"))


### Calculate 1-D footprints (fetch lengths) ===================================

if (config$fetch) {
  
  cat("Calculating half-hourly fetch lengths.\n")
  
  # Select model function and parameters
  fetch_ref <- list(
    K15 = list(
      fun = calc_fetch_kljun, 
      vars = c("ustar", "mo_length", "blh", "z0")
    ),
    KM01 = list(
      fun = calc_fetch_kormann, 
      vars = c("ws", "ustar", "mo_length", "z0")
    ),
    H00 = list(
      fun = calc_fetch_hsieh, 
      vars = c("mo_length", "z0")
    )
  )
  
  fetch_ref <- pluck(fetch_ref, config$fetch_model)
  
  # Add WS to model vars if specified
  if (config$fetch_model == "K15" & config$use_ws) {
    #config$fp_model <- "K15_ws"
    fetch_ref <- list_merge(fetch_ref, vars = "ws")
  }
  
  # Select only necessary variables, remove missing data
  data_fetch <- data %>% 
    # Filter out erroneous wind records
    filter(qc_ws == 0) %>%
    select(timestamp, z, zd, all_of(fetch_ref$vars)) %>%
    drop_na()
  
  # Copy data into list (easier to access in loop)
  data_fetch_list <- data_fetch %>%
    select(-timestamp) %>%
    as.list()
  
  # Initialize loop for fetch calculation
  n_fetch <- nrow(data_fetch)
  fetch_list <- vctrs::vec_init(list(), n_fetch)
  p <- progress_info(n_fetch)
  
  # Apply one-dimensional footprint model to data
  for (i in 1:n_fetch) {
    
    # Get parameters from data list
    fetch_args <- data_fetch_list %>% 
      modify(i) %>% 
      # Splice in percents
      prepend(list(pct = config$fetch_pct, dx = config$fetch_dx))
    
    # Calculate fetch lengths
    fetch_list[[i]] <- rlang::exec(fetch_ref$fun, !!!fetch_args)
    
    p$tick()
  }
  
  fetch <- fetch_list %>%
    bind_rows() %>%
    bind_cols(select(data_fetch, timestamp), .) %>%
    right_join(select(data, timestamp), by = "timestamp") %>%
    arrange(timestamp)
  
  fetch_out <- file.path(
    path_out, "fetch", paste0("fetch_", config$fetch_model, "_", tag, ".csv")
  )
  write_csv(fetch, fetch_out)
  
  # Create documentation for fetch output
  fetch_docu <- settings
  fetch_docu <- append(
    settings, list(
      files = c(file_data, era_files),
      model_params = list(
        model = config$fetch_model, 
        roughness = if ("z0" %in% fetch_ref$vars) "z0" else "ws/ustar",
        percents = config$fetch_pct
      )
    )
  )
  fetch_docu_out <- str_replace(fetch_out, ".csv", ".txt")
  # Save documentation
  sink(fetch_docu_out)
  print(fetch_docu) 
  sink()
  
  cat("\n")
}


### Calculate 2-D footprints and coverage ======================================

if (config$fp) {
  
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
  
  model_ref <- pluck(model_ref, config$fp_model)
  
  # Add WS to model vars if specified
  if (config$fp_model == "K15" & config$use_ws) {
    #config$fp_model <- "K15_ws"
    model_ref <- list_merge(model_ref, vars = "ws")
  }
  
  fp_fun <- pluck(model_ref, "fun")
  
  # (See references for description of models)
  
  # Create new folder for output files
  fp_out <- file.path(
    path_out, "footprint", paste0("halfhourly_", config$fp_model, "_", tag)
  )
  dir.create(fp_out, showWarnings = FALSE)
  
  # Create sub-folders for grid and footprint matrices
  dir.create(file.path(fp_out, "grid"), showWarnings = FALSE)
  dir.create(file.path(fp_out, "footprints"), showWarnings = FALSE)
  
  # Select only necessary variables, remove missing data
  data_fp <- data %>% 
    # Filter out erroneous wind records
    filter(qc_ws == 0, qc_wd == 0) %>%
    # All models need z, zd, and wd
    select(timestamp, z, zd, wd, all_of(model_ref$vars)) %>%
    drop_na()
  
  # Copy data into list (easier to access in footprint loop)
  data_fp_list <- data_fp %>%
    select(-timestamp) %>%
    as.list()
  
  # Set up grid
  grid <- aoi_to_grid(
    delin, c(md$tower$coords$east, md$tower$coords$north), delta = 1
  )
  
  # Convert AOI to grid
  aoi_grid <- aoi_to_mask(
    delin, c(md$tower$coords$east, md$tower$coords$north), delta = 1
  )
  
  # Write grid to file
  imap(
    list_modify(grid, aoi = aoi_grid), 
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
      modify(i) %>% 
      # Splice in grid
      prepend(list(grid = grid))
    
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
  dir_size <- fp_out %>% 
    file.path("footprints") %>% 
    fs::dir_info() %>% 
    summarize(size = sum(size)) %>% 
    pull(size)
  
  cat(
    "\nFinished calculating footprints.\n", 
    "Matrix files can be found in:", file.path(fp_out, "footprints"), "\n",
    "Total size of directory:", dir_size, "\n"
  )
  cat("Writing summary output...")
  
  # Calculate topology and write to file
  topo_out <- file.path(
    path_out, "topology", paste0("all_", config$fp_model, "_", tag)
  )
  fp_topo <- fp_topo / n_topo
  write_matrix(fp_topo, topo_out, trunc = 0, compress = FALSE)
  
  # Write file with footprint stats
  fp_stats <- phi %>% 
    enframe(name = "timestamp", value = "phi") %>% 
    mutate(timestamp = ymd_hms(timestamp)) %>%
    unnest(phi) %>%
    # Classify phi according to Gockede et al. 2008
    mutate(
      phi_class = case_when(
        # Homogenous measurements
        between(phi, 0.95, 1.00) ~ 0L,
        # Representative measurements
        between(phi, 0.80, 0.95) ~ 1L,
        # Acceptable measurements
        between(phi, 0.50, 0.80) ~ 2L,
        # Disturbed measurements
        between(phi, 0.00, 0.50) ~ 3L
      )
    ) %>% 
    # Give phi vars suffix with model tag
    rename_with(~ str_c(.x, "_", tolower(config$fp_model)), c(phi, phi_class))
  
  # Compute model-specific validity flags
  if (str_detect(config$fp_model, "K15")) {
    
    # Indicate whether footprint is valid according to Kljun et al. 2015
    fp_stats <- fp_stats %>%
      right_join(
        select(data, timestamp, ustar, mo_length, blh, ws, z, zd, z0), 
        by = "timestamp"
      )
    
    if (config$use_ws) {
      fp_stats <- mutate(fp_stats, z0 = calc_z0(ws, ustar, mo_length, z - zd))
    }
    
    fp_stats <- fp_stats %>%
      arrange(timestamp) %>%
      mutate(
        zm = z - zd,
        # TODO fix this
        k15_valid = case_when(
          ustar < 0.1 ~ 0L,
          zm / mo_length > -15.5 ~ 0L,
          # Using roughness sublayer limit as in latest FFP version, not paper
          12.5 * z0 >= zm ~ 0L,
          zm > 0.8 * blh ~ 0L,
          is.na(phi_k15) ~ NA_integer_,
          TRUE ~ 1L
        )
      ) %>%
      select(-ustar, -mo_length, -ws)
  }
  
  if (config$fp_model == "KM01") {
    # Indicate whether footprint is valid according to Kormann & Meixner 2001
    fp_stats <- fp_stats %>%
      right_join(
        select(data, timestamp, mo_length, z, zd), by = "timestamp"
      ) %>%
      arrange(timestamp) %>%
      mutate(
        zm = z - zd,
        km01_valid = if_else(abs(zm / mo_length) > 3, 0L, 1L)
      ) %>%
      select(-z, -zm, -zd)
  }
  
  fp_stats_out <- file.path(
    path_out, "stats", 
    paste0("footprint_stats_", config$fp_model, "_", tag, ".csv")
  )
  write_csv(fp_stats, fp_stats_out)
  
  # Create documentation for phi output
  fp_stats_docu <- settings
  fp_stats_docu <- append(
    settings, list(
      files = c(file_data, file_delin, era_files),
      model_params = list(
        model = config$fp_model, 
        roughness = if ("z0" %in% model_ref$vars) "z0" else "ws/ustar",
        fetch = attr(grid, "fetch"), 
        res = attr(grid, "res")
      )
    )
  )
  fp_stats_docu_out <- str_replace(fp_stats_out, ".csv", ".txt")
  # Save documentation
  sink(fp_stats_docu_out)
  print(fp_stats_docu) 
  sink()
  
  cat("done.\n")
}

# End