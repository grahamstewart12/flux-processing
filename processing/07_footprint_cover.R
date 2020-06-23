### ============================================================================
# Spatial analysis of flux area ================================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLN", # three letter site code
  year = 2019, # four digit year
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

control <- list(
  fp_model = "K15",
  features = c("ndvi", "ndwi", "z_max", "i_vp", "rn_vp"),
  phi = 0.85
)

# Load the required packages
devtools::load_all("/Users/Graham/R Projects/footprints")
library(progress)
library(lubridate)
library(tidyverse)

# Load reference files
source("/Users/Graham/Desktop/DATA/Flux/tools/reference/site_metadata.R")

# Load functions
path_funs <- "/Users/Graham/Desktop/DATA/Flux/tools/engine/functions"
source(file.path(path_funs, "latest_version.R"))
source(file.path(path_funs, "utilities.R"))


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

cut_circle <- function(x, n_perim = 4, n_radius = 1) {
  
  center <- sf::st_centroid(x)
  
  circle_points <- x %>% 
    sf::st_cast("MULTIPOINT") %>% 
    sf::st_sfc() %>% 
    sf::st_cast("POINT") %>%
    sf::st_sf() 
  
  seg <- circle_points %>%
    dplyr::mutate(quad = rep(seq(1, n_perim), each = 360 / n_perim)) %>%
    dplyr::group_by(quad) %>%
    dplyr::slice((360 / n_perim / 2) + 1) %>% 
    dplyr::ungroup() %>%
    sf::st_cast("MULTIPOINT") %>% 
    dplyr::mutate(seg = purrr::map(geometry, ~ c(., center))) %>% 
    dplyr::pull(seg) %>% 
    sf::st_sfc() %>% 
    sf::st_cast("LINESTRING") %>%
    sf::st_cast("MULTILINESTRING") %>% 
    sf::st_line_merge() %>%
    sf::st_combine()
  
  perim_cut <- x %>%
    lwgeom::st_split(seg) %>% 
    sf::st_collection_extract()
  
  if (n_radius == 1) {
    return(sf::st_sf(perim_cut))
  }
  
  cookie_cutter <- function(x, y) {
    sf::st_collection_extract(lwgeom::st_split(x, y))
  }
  
  radius <- calc_radius(x)
  
  radii <- seq(0.5, radius, length.out = n_radius + 1)
  inside_radii <- radii[-c(1, length(radii))]
  
  radius_cut <- inside_radii %>%
    purrr::map(~ sf::st_buffer(center, .x, nQuadSegs = 360 / n_perim)) %>% 
    purrr::map(sf::st_cast, "LINESTRING") %>%
    purrr::reduce(cookie_cutter, .init = perim_cut)
  
  sector_names <- stringr::str_c(
    rep(c("E", "S", "W", "N"), each = n_radius),
    rep(1:n_radius, times = n_perim)
  )
  
  radius_cut %>%
    sf::st_sf() %>%
    dplyr::rename(geometry = 1) %>%
    dplyr::mutate(
      id = dplyr::row_number(), 
      sector = sector_names,
      .before = 1
    )
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- purrr::pluck(site_metadata, settings$site)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
wd <- file.path(
  "/Users/Graham/Desktop", "DATA", "Flux", settings$site, settings$year
)
path_spatial <- file.path(dirname(wd), "analysis", "spatial")
path_in <- file.path(wd, "processing", "05_footprint")

paths <- list(
  # Half-hourly water level
  wtd = file.path(
    "/Users/Graham/Desktop", "DATA", "PT", "output", settings$site, 
    paste0(settings$site, "_", settings$year, ".csv")
  ),
  # Cover type raster
  class = file.path(path_spatial, "07_classification", "classes_ens.tif"),
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
  paths$wtd, guess_max = 6000, 
  col_types = readr::cols(.default = readr::col_guess())
)

# Import footprint phi data
data_phi <- readr::read_csv(
  paths$phi, guess_max = 6000, 
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
feat_grid <- purrr::map(feat, ~ snap_to_grid(.x, grid, c(md$x_utm, md$y_utm)))

# Create key for combined cover types & flooding
key <- class_key %>% 
  purrr::map(~ .x * flood_key) %>% 
  purrr::imap(~ rlang::set_names(.x, stringr::str_c(.y, "_", names(.)))) %>%
  purrr::flatten() %>%
  purrr::simplify()

# Split site grid into objective quadrants
# - is it problematic to split radial quadrants, since not all same area??
# - I think yes
# quad <- sf::st_point(c(md$x_utm, md$y_utm)) %>%
#   # Important to set nQuadSegs to 90 so that each point = 1 degree
#   sf::st_buffer(500.5, nQuadSegs = 90) %>% 
#   cut_circle(n_perim = 4, n_radius = 10) %>% 
#   sf::st_set_crs(sf::st_crs(delin)) %>% 
#   stars::st_rasterize(dx = 1, dy = 1) %>%
#   # TODO make this not depend on 'class' being on the correct grid
#   stars::st_warp(stars::st_as_stars(class))
# which_quads <- quad %>% 
#   tibble::as_tibble() %>% 
#   dplyr::distinct(id) %>% 
#   dplyr::pull(id)
# # Create named key for quadrants
# quad_key <- c("e", "s", "w", "n") %>% 
#   rep(each = 10) %>% 
#   stringr::str_c(rep(rev(1:10), times = dplyr::n_distinct(.))) %>% 
#   rlang::set_names() %>%
#   purrr::map2(seq_along(.), ~ .y) %>%
#   purrr::keep(~ .x %in% which_quads) %>%
#   purrr::simplify()
# quad_grid <- as_matrix(quad)

# Split site grid into objective quadrants
quad_grid <- cut_grid(grid, width = 50)
# Create named key for quadrants
quad_key <- quad_grid %>%
  magrittr::multiply_by(grid$aoi) %>%
  as.vector() %>%
  vctrs::vec_unique() %>%
  purrr::discard(~ .x == 0) %>%
  rlang::set_names(stringr::str_c("q", .))


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
  # Join phi data
  dplyr::left_join(
    dplyr::select(data_phi, timestamp, phi = 2), by = "timestamp"
  ) %>%
  # Subset footprints with phi above threshold
  dplyr::filter(phi >= control$phi) %>%
  dplyr::select(-phi) %>%
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
  fp_temp <- read_matrix(fp_files[i], trunc = 9) # 6.61 ms

  # Mask integrated footprint to account for rapid expansion approaching 100%
  # - recommended: between 80% and 90% (Kljun et al. 2015)
  # - this means data should be first filtered for phi >= p
  fp_mask <- mask_source_area(fp_temp, p = control$phi, mask_value = 0) # 4.17 ms
  fp_mask <- fp_mask * grid$aoi # 0.192 ms
  
  # Normalize integrated footprint to 1 (Tuovinen et al. 2019)
  # - the masked footprint is thus considered the 100% analytical footprint
  fp_norm <- fp_mask / sum(fp_mask)
  
  # Calculate water level over AOI
  # - as depth relative to land surface (negative = water below surface)
  # - give wtd a placeholder value if NA so that other covers can be summed
  aoi_wtd <- md$rel_elev_well + tidyr::replace_na(wtd[i], 0.1) - elev_grid
  # Classify flooded area (1 = wet, 0 = dry)
  aoi_flooded <- with_matrix(aoi_wtd, ~ dplyr::if_else(.x >= 0, 1, 0)) 
  # aoi_flooded <- sign(aoi_wtd)
  
  # Can use cross_grids() if want to combine more than 2 matrices
  # - this is kind of awkward but flexible
  # - automatically detects if one matrix is all NA and drops it
  # - ensures that all cells are classified properly
  #aoi_cover <- cross_grids(class_grid, aoi_flooded) # 5.52 ms
  # For cover classes, flooded is classified as (1 = wet, -1 = dry)
  aoi_cover <- class_grid * sign(aoi_wtd)
  
  # Calculate cover-type weights
  cover_wt <- summarize_cover(fp_norm, aoi_cover, levels = key) # 3.79 ms
  
  # Calculate quadrant weights
  quad_wt <- summarize_cover(fp_norm, quad_grid, levels = quad_key) # 3.79 ms
  
  # Calculate footprint-weighted WTD & flooded area
  wtd_wt <- summarize_cover(fp_norm, aoi_wtd, type = "numeric") # 0.605 ms
  flood_wt <- summarize_cover(fp_norm, aoi_flooded, type = "numeric") # 0.605 ms
  feat_wt <- feat_grid %>%
    purrr::map(~ summarize_cover(fp_norm, .x, type = "numeric")) %>%
    purrr::simplify()
  
  # Add all weights to list 
  cover[[i]] <- c(
    cover_wt, "wtd_wt" = wtd_wt, "flood_wt" = flood_wt, feat_wt, quad_wt
  )

  p$tick()
}

# Gather everything into one data frame
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
    .after = 1
  ) %>%
  # Set wet/dry covers to NA if WTD was missing
  dplyr::mutate(dplyr::across(
    c(dplyr::ends_with("wet"), dplyr::ends_with("dry")), 
    ~ dplyr::if_else(is.na(wtd), NA_real_, .x)
  )) %>%
  dplyr::select(-wtd)

# Peek distributions
data_cover_all %>%
  ggplot2::ggplot(ggplot2::aes(frs)) +
  ggplot2::geom_density(na.rm = TRUE)

# Quadrant weights
ggplot2::ggplot() + 
  stars::geom_stars(
    data = quad_grid %>% 
      t() %>% 
      stars::st_as_stars() %>% 
      tibble::as_tibble(center = FALSE) %>% 
      dplyr::left_join(
        data_cover_all %>% 
          dplyr::select(dplyr::starts_with("q")) %>% 
          dplyr::summarize(dplyr::across(.fns = ~ sum(.x, na.rm = TRUE))) %>% 
          tidyr::pivot_longer(dplyr::everything()) %>% 
          dplyr::mutate(name = dplyr::recode(name, !!!quad_key)), 
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


