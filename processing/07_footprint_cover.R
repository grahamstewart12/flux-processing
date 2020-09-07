### ============================================================================
# Spatial analysis of flux area ================================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

rm(list = ls())
.rs.restartR()
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLN", # three letter site code
  year = 2019, # four digit year
  dir = "/Users/Graham/Desktop/DATA", # where all data can be found
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

control <- list(
  fp_model = "K15",
  features = c("ndvi", "ndwi", "z_max", "i_vp"),
  grid_width = 50,
  phi = 0.85,
  # Water level of flooded/dry threshold 
  flood_depth = -0.10
)

# Load the required packages
devtools::load_all("~/Projects/flux/footprints")
devtools::load_all("~/Projects/flux/dscalr")
library(progress)
library(lubridate)
library(tidyverse)

# Load reference files
source("/Users/Graham/Desktop/DATA/Flux/tools/reference/site_metadata.R")


### Helper functions ===========================================================

progress_info <- function(len) {
  
  progress_bar$new(
    total = len, 
    format = paste0(
      "[:spin] Completed: :current | :percent  ", 
      "Elapsed: :elapsed  Remaining: :eta"
    ),
    clear = FALSE
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
elev_corr <- md$rel_elev_well

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path(settings$dir, "Flux", settings$site, settings$year)
path_spatial <- file.path(dirname(wd), "analysis", "spatial")
path_in <- file.path(wd, "processing", "05_footprint")

paths <- list(
  # Half-hourly water level
  wtd = file.path(
    "/Users/Graham/Desktop", "DATA", "PT", "output", settings$site, 
    paste0(settings$site, "_", settings$year, ".csv")
  ),
  # Cover type raster
  class = file.path(path_spatial, "07_classification", "classes_rev.tif"),
  # Relative elevation raster
  elev = file.path(path_spatial, "08_relative_elevation", "rel_elev.tif"),
  # Other spatial features
  feat = file.path(path_spatial, "03_feature_extraction"),
  delin = file.path(dirname(wd), "site_info", "delineation"),
  # Footprint data
  fp = latest_version(
    file.path(path_in, "footprint"), paste0("halfhourly_", control$fp_model), ""
  ),
  phi = latest_version(
    file.path(path_in, "stats"), paste0("footprint_stats_", control$fp_model)
  ),
  # Output files
  out = file.path(wd, "processing", "07_footprint_cover", "data")
)

# Set tag for creating output file names
tag_out <- create_tag(settings$site, settings$year, settings$date)


### Load required input data ===================================================

# Import water level data
data_wtd <- readr::read_csv(
  paths$wtd, guess_max = 7000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Import footprint phi data
data_phi <- readr::read_csv(
  paths$phi, guess_max = 7000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Import footprint grid (including AOI grid)
grid <- paths$fp %>% 
  file.path("grid") %>%
  file.path(paste0(c("x", "y", "aoi"), ".txt")) %>% 
  purrr::map(read_matrix, trunc = 0) %>%
  rlang::set_names(c("x", "y", "aoi"))

# Import classified image
class <- raster::raster(paths$class)
class_key <- c("wtr" = 1, "veg" = 2, "frs" = 3)

# Import relative elevation
rel_elev <- raster::raster(paths$elev)
flood_key <- c("dry" = -1, "wet" = 1)

# Read all other spatial features
feat <- paths$feat %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE) %>%
  stringr::str_subset(
    stringr::str_c("/", control$feat, ".", collapse = "|")
  ) %>%
  purrr::map(raster::raster) %>%
  rlang::set_names(purrr::map(., names))


### Prepare input data =========================================================

# Map imagery onto same grid as footprint
class_grid <- class %>% 
  snap_to_grid(grid, c(md$x_utm, md$y_utm)) %>%
  with_matrix(~ as.integer(.x))
elev_grid <- snap_to_grid(rel_elev, grid, c(md$x_utm, md$y_utm))
feat_grid <- feat %>%
  purrr::map(~ snap_to_grid(.x, grid, c(md$x_utm, md$y_utm))) %>%
  purrr::prepend(list(rel_elev = elev_grid))

# Create key for combined cover types & flooding
classflood_key <- class_key %>% 
  purrr::map(~ .x * flood_key) %>% 
  purrr::imap(~ rlang::set_names(.x, stringr::str_c(.y, "_", names(.)))) %>%
  purrr::flatten() %>%
  purrr::simplify()

# Split site grid into objective quadrants
quad_grid <- with_matrix2(
  sign(grid$x), sign(grid$y), 
  ~ dplyr::case_when(
    .x ==  1 & .y ==  1 ~ 1,
    .x == -1 & .y ==  1 ~ 2,
    .x == -1 & .y == -1 ~ 3,
    .x ==  1 & .y == -1 ~ 4
  )
)
quad_key <- c(q1 = 1, q2 = 2, q3 = 3, q4 = 4)
quad_cover <- quad_grid %>% 
  with_matrix2(class_grid, ~ forcats::fct_cross(factor(.x), factor(.y))) %>% 
  with_matrix(~ as.integer(factor(.x))) 
quadclass_key <- quad_cover %>% 
  as.vector() %>% unique() %>% sort() %>%
  rlang::set_names(paste(
    rep(names(quad_key), each = length(class_key)), names(class_key), 
    sep = "_"
  ))

### Retrieve footprints, calculate cover =======================================

# Initialize loop

# Get list of footprint files
fp_files <- paths$fp %>% 
  file.path("footprints") %>%
  list.files(full.names = TRUE) %>%
  tibble::enframe(name = NULL, value = "file") %>%
  # Parse timestamps from file names
  dplyr::mutate(
    timestamp = file %>% 
      basename() %>% 
      stringr::str_sub(1, -5) %>% 
      lubridate::ymd_hms(),
    .before = 1
  ) %>%
  tibble::deframe()

fp_timestamps <- fp_files %>% names() %>% lubridate::ymd_hms()

n_fp <- length(fp_files)

# Empty list to hold cover data
cover <- vctrs::vec_init(list(), n_fp)

# Select water level data 
wtd <- data_wtd %>%
  dplyr::right_join(
    tibble::enframe(fp_timestamps, name = NULL, value = "timestamp"),
    by = "timestamp"
  ) %>%
  dplyr::pull(wtd_f)
p <- progress_info(n_fp)

for (i in 1:n_fp) {
  
  # Read and scale footprint matrix
  fp_temp <- read_matrix(fp_files[i], trunc = 9)

  # Mask integrated footprint to account for rapid expansion approaching 100%
  # - recommended: between 80% and 90% (Kljun et al. 2015)
  # - this means data should be filtered for phi >= p
  fp_mask <- mask_source_area(fp_temp, p = control$phi, mask_value = 0)
  fp_mask <- fp_mask * grid$aoi
  
  # Normalize integrated footprint to 1 (Tuovinen et al. 2019)
  # - the masked footprint is thus considered the 100% analytical footprint
  phi_sum <- sum(fp_mask)
  fp_norm <- fp_mask / phi_sum
  
  # Calculate water level over AOI
  # - as depth relative to land surface (negative = water below surface)
  # - give wtd a placeholder value if NA so that other covers can be summed
  aoi_wtd <- md$rel_elev_well + tidyr::replace_na(wtd[i], 0.1) - elev_grid
  # Classify flooded area (1 = wet, 0 = dry)
  aoi_flooded <- with_matrix(
    aoi_wtd, ~ dplyr::if_else(.x >= control$flood_depth, 1, 0)
  ) 
  
  # For cover classes, flooded is classified as (1 = wet, -1 = dry)
  aoi_cover <- class_grid * sign(aoi_wtd)
  
  # Calculate cover-type weights
  cover_wt <- summarize_cover(fp_norm, aoi_cover, levels = classflood_key)
  
  # Calculate quadrant weights
  quad_wt <- summarize_cover(fp_norm, quad_cover, levels = quadclass_key)
  
  # Calculate footprint-weighted WTD & flooded area
  wtd_wt <- summarize_cover(fp_norm, aoi_wtd, type = "numeric")
  flood_wt <- summarize_cover(fp_norm, aoi_flooded, type = "numeric")
  feat_wt <- feat_grid %>%
    purrr::map(~ summarize_cover(fp_norm, .x, type = "numeric")) %>%
    purrr::simplify()
  
  # Add all weights to list 
  cover[[i]] <- c(
    phi_mask = phi_sum, cover_wt, "wtd_wt" = wtd_wt, "flood_wt" = flood_wt, 
    feat_wt, quad_wt
  )

  p$tick()
}

# Gather everything into one data frame (this is very slow)
data_cover <- cover %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(timestamp = fp_timestamps) %>%
  dplyr::relocate(timestamp) %>%
  dplyr::right_join(dplyr::select(data_phi, timestamp), by = "timestamp") %>%
  dplyr::arrange(timestamp)

# Sum crosstabs
data_cover_all <- data_cover %>%
  dplyr::left_join(
    dplyr::select(data_wtd, timestamp, wtd = wtd_f), by = "timestamp"
  ) %>%
  dplyr::mutate(
    wtr = wtr_wet + wtr_dry,
    veg = veg_wet + veg_dry,
    frs = frs_wet + frs_dry,
    wet = wtr_wet + veg_wet + frs_wet,
    dry = wtr_dry + veg_dry + frs_dry,
    q1 = q1_wtr + q1_veg + q1_frs,
    q2 = q2_wtr + q2_veg + q2_frs,
    q3 = q3_wtr + q3_veg + q3_frs,
    q4 = q4_wtr + q4_veg + q4_frs,
    .after = 1
  ) %>%
  # Set wet/dry covers to NA if WTD was missing
  dplyr::mutate(dplyr::across(
    c(dplyr::ends_with("wet"), dplyr::ends_with("dry")), 
    ~ dplyr::if_else(is.na(wtd), NA_real_, .x)
  ))

# Peek distributions
data_cover_all %>%
  tidyr::pivot_longer(c(wet, dry)) %>%
  ggplot2::ggplot(ggplot2::aes(value, fill = name)) +
  ggplot2::geom_density(na.rm = TRUE, alpha = 0.6)

data_cover_all %>%
  dplyr::summarize(
    dplyr::across(c(-timestamp, -phi_mask), sum, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(dplyr::everything())

# Quadrant weights
ggplot2::ggplot() + 
  stars::geom_stars(
    data = quad_grid %>% 
      apply(2, rev) %>% t() %>%
      stars::st_as_stars() %>% 
      tibble::as_tibble(center = FALSE) %>% 
      dplyr::left_join(
        data_cover_all %>% 
          dplyr::filter(phi_mask > 0.8) %>%
          dplyr::select(q1, q2, q3, q4) %>% 
          dplyr::summarize(dplyr::across(.fns = ~ sum(.x, na.rm = TRUE))) %>% 
          tidyr::pivot_longer(dplyr::everything()) %>%
          dplyr::mutate(name = as.numeric(stringr::str_remove_all(name, "q"))),
        by = c("A1" = "name")
      ) %>% 
      dplyr::select(-A1) %>% 
      stars::st_as_stars()
  ) + 
  ggplot2::scale_fill_distiller(palette = "Spectral", trans = "log") + 
  ggplot2::coord_fixed()

# Write output to file
covers_out <- file.path(paths$out, paste0("footprint_cover_", tag_out, ".csv"))
readr::write_csv(data_cover_all, covers_out)


