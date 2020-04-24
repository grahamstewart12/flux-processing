### ============================================================================
# Flux hotspots within footprint ===============================================
### ============================================================================

# For use with site JLR

# Purpose: 

# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/Desktop/RESEARCH/fluxtools")
library(openeddy)
library(lubridate)
library(tidyverse)


### ============================================================================
# Initialize script settings & documentation

# Set the desired working directory in RStudio interface
wd <- "~/Desktop/DATA/Flux/JLR"
wd_anal <- "~/Desktop/DATA/Flux/2019/JLR/analysis_data"
# This assumes that the subdirectory structure is already present

# Set the file names
# Half-hourly data
path_2018 <- latest_version(
  file.path(wd, "2018", "processing_data", "08_flux_qc", "output"), 
  "flux_qc_essentials"
)
path_2019 <- latest_version(
  file.path(wd, "2019", "processing_data", "08_flux_qc", "output"), 
  "flux_qc_essentials"
)

# Aerial image
rgbn_path <- file.path(
  wd, "analysis", "image_classification", "segm_feat", "JLR_NAIP_2018.tif"
)

# Classified image
img_path <- file.path(
  wd, "analysis", "image_classification", "output", "classified.tif"
)

# DEM
dem_path <- file.path(wd, "imagery", "JLR_dem.tif")

# Site area
aoi_path <- file.path(wd, "JLR_area_proc")


# Set the session information
# These settings will be used to determine which processing options are
# implemented in the script. It will also form part of the saved documentation,
# so options should be adjusted here rather than in the script directly.
settings <- list(
  # Session information
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  siteyear = "JLR", # three letter side code + two digit year
  session_time = lubridate::today(), # date script was run
  session_info = devtools::session_info(), # R session info: R, OS, packages
  # Processing settings
  files = c(), # names of eddypro/biomet files used
  note = "" # any additional notes
)

# Load metadata file
metadata <- read_metadata(file.path(wd, "site_info", "metadata_JLR.txt"))

# Set tag for creating output file names
tag_out <- paste0(settings$siteyear, "_", settings$session_time)


### Load required input data ===================================================

# Import flux data
site_2018 <- read_eddy(path_2018, stringsAsFactors = FALSE) %>%
  # Strip units and varnames - not needed here
  drop_all_attributes(c("units", "varnames")) %>%
  # Prevent integer & factor columns from messing with join
  mutate_if(is.integer, as.numeric) %>%
  # Parse timestamp
  mutate(timestamp = ymd_hms(timestamp, tz = "Etc/GMT+5"))
glimpse(site_2018) # check if loaded properly
site_2018 %>% 
  summarize(first(timestamp), last(timestamp)) # Starts/ends as expected?

site_2019 <- read_eddy(path_2019, stringsAsFactors = FALSE) %>% 
  # Strip units and varnames - not needed here
  drop_all_attributes(c("units", "varnames")) %>%
  # Prevent integer columns from messing with join
  mutate_if(is.integer, as.numeric) %>%
  # Parse timestamp
  mutate(timestamp = ymd_hms(timestamp, tz = "Etc/GMT+5")) 
glimpse(site_2019) # check if loaded properly
site_2019 %>% 
  summarize(first(timestamp), last(timestamp)) # Starts/ends as expected?

# Combine both years
site <- full_join(site_2018, site_2019)

# Subset flux data for footprint calculations
site_fp <- site %>%
  # Only good-quality daytime FCH4 data during growing season
  mutate(
    year = year(timestamp) %>% factor(),
    yday = yday(timestamp),
    fch4 = clean(fch4, qc_fch4_final),
    fc = clean(fc, qc_fc_final),
    le = clean(le, qc_le_final),
    h = clean(h, qc_h_final),
    tau = clean(tau, qc_tau_final),
    stab = atm_stability(
      mo_length, metadata$tower_height, metadata$displacement, 
      metadata$roughness_length
    ),
    grow_seas = between(month(timestamp), 4, 10)
  ) %>%
  filter(
    !is.na(fch4) | !is.na(fc) | !is.na(le) | !is.na(h) | !is.na(tau), 
    grow_seas, night == 0
  ) %>%
  select(-grow_seas)

# Check flux data
qplot(yday, fch4, color = year, data = site_fp)
qplot(yday, fc, color = year, data = site_fp)
qplot(yday, le, color = year, data = site_fp)
qplot(yday, h, color = year, data = site_fp)
qplot(yday, tau, color = year, data = site_fp)

# Import site imagery
dem <- raster::raster(dem_path)
img <- raster::raster(img_path)
img <- raster::setValues(img, factor(raster::values(img)))
rgbn <- raster::stack(rgbn_path)
rgbn <- raster::crop(rgbn, img)

# Import site area polygon
aoi <- rgdal::readOGR(aoi_path)
img %>%
  raster::rasterToPoints() %>% 
  data.frame() %>%
  rename(z = classified) %>%
  mutate(z = factor(z)) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = z)) +
  geom_path(aoi, mapping = aes(long, lat), alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  coord_fixed()
tower <- data.frame(x = metadata$x_utm, y = metadata$y_utm) %>% 
  sp::SpatialPoints(proj4string = raster::crs(img))
  #as("SpatialPointsDataFrame")


### Calculate footprints =======================================================

# All steps must be done in one loop due to memory constraints

# 1. grid

grid <- fp_grid(fetch = (metadata$tower_height - metadata$displacement) * 100)

# 2. footprints

# Initialize loop
len <- nrow(site_fp)
fp <- vector("list", length = len)
progress <- progress_estimated(len)

# Calculate footprint matrices
for (i in 1:len) {
  
  # Calculate footprint, add to list
  fp[[i]] <- fp_calculate(
    wd = site_fp$wd[i],
    ustar = site_fp$ustar[i],
    L = site_fp$mo_length[i],
    v_sd = site_fp$v_sigma[i],
    z = metadata$tower_height,
    zd = metadata$displacement,
    zo = metadata$roughness_length,
    grid = grid,
    model = "H00"
  )
  
  # Rasterize
  fp[[i]] <- fp_rasterize(
    fp[[i]], grid = grid, coords = c(metadata$x_utm, metadata$y_utm), 
    crs = raster::crs(img)
  )
  
  # Clip to 80%
  fp[[i]] <- fp_percent(fp[[i]], p = 0.80, fill = 0)
  
  # Crop to site extent
  fp[[i]] <- raster::crop(fp[[i]], aoi)
  
  # Set output name as the timestamp
  names(fp)[[i]] <- format(site_fp$timestamp[i], "%Y-%m-%d-%H%M%S")
  
  progress$tick()$print()
  
}

# "Weight" footprint probabilities by fluxes
fluxes <- c("fch4", "fc", "le", "h", "tau")
fp_flux <- vector("list", length = length(fluxes)) %>%
  lmap(list_modify, vector("list", length = 3)) %>%
  set_names(fluxes) %>%
  map(set_names, c("neutral", "stable", "unstable"))

for (i in seq_along(fp_flux)) {
  fp_flux[[i]]$neutral <- fp %>% 
    raster::stack() %>%
    raster::calc(function(x) x * site_fp[, fluxes[i]], forceapply = TRUE) %>%
    raster::subset(
      which(!is.na(site_fp[, fluxes[i]]) & site_fp$stab == "neutral")
    )
  fp_flux[[i]]$stable <- fp %>% 
    raster::stack() %>%
    raster::calc(function(x) x * site_fp[, fluxes[i]], forceapply = TRUE) %>%
    raster::subset(
      which(!is.na(site_fp[, fluxes[i]]) & site_fp$stab == "stable")
    )
  fp_flux[[i]]$unstable <- fp %>% 
    raster::stack() %>%
    raster::calc(function(x) x * site_fp[, fluxes[i]], forceapply = TRUE) %>%
    raster::subset(
      which(!is.na(site_fp[, fluxes[i]]) & site_fp$stab == "unstable")
    )
  progress_n(i)
}


# sum of the products of CH4 fluxes and footprint weights divided by the sum of footprint weights covering that grid cell
# sum(fch4 * fp) / sum(fp)

# Calculate summed footprint cells
#fp_sum <- fp %>% raster::stack() %>% raster::calc(fun = sum, na.rm = TRUE)
fp_flux_sum <- vector("list", length = length(fluxes)) %>%
  lmap(list_modify, vector("list", length = 3)) %>%
  set_names(fluxes) %>%
  map(set_names, c("neutral", "stable", "unstable"))

for (i in seq_along(fp_flux_sum)) {
  fp_flux_sum[[i]]$neutral <- fp_flux[[i]]$neutral %>%
    raster::calc(fun = sum, na.rm = TRUE)
  fp_flux_sum[[i]]$stable <- fp_flux[[i]]$stable %>%
    raster::calc(fun = sum, na.rm = TRUE)
  fp_flux_sum[[i]]$unstable <- fp_flux[[i]]$unstable %>%
    raster::calc(fun = sum, na.rm = TRUE)
  progress_n(i)
}


# Divide time-averaged flux by time-averaged footprint
fp_flux_top <- vector("list", length = length(fluxes)) %>%
  lmap(list_modify, vector("list", length = 3)) %>%
  set_names(fluxes) %>%
  map(set_names, c("neutral", "stable", "unstable"))

for (i in seq_along(fp_flux_top)) {
  fp_flux_top[[i]]$neutral <- fp %>% 
    raster::stack() %>%
    raster::subset(
      which(!is.na(site_fp[, fluxes[i]]) & site_fp$stab == "neutral")
    ) %>%
    raster::calc(fun = sum, na.rm = TRUE) %>%
    raster::overlay(fp_flux_sum[[i]]$neutral, fun = function(x, y) y / x)
  fp_flux_top[[i]]$stable <- fp %>% 
    raster::stack() %>%
    raster::subset(
      which(!is.na(site_fp[, fluxes[i]]) & site_fp$stab == "stable")
    ) %>%
    raster::calc(fun = sum, na.rm = TRUE) %>%
    raster::overlay(fp_flux_sum[[i]]$stable, fun = function(x, y) y / x)
  fp_flux_top[[i]]$unstable <- fp %>% 
    raster::stack() %>%
    raster::subset(
      which(!is.na(site_fp[, fluxes[i]]) & site_fp$stab == "unstable")
    ) %>%
    raster::calc(fun = sum, na.rm = TRUE) %>%
    raster::overlay(fp_flux_sum[[i]]$unstable, fun = function(x, y) y / x)
  progress_n(i)
}

fch4 <- fp %>% 
  raster::stack() %>% 
  raster::subset(which(!is.na(site_fp$fch4))) %>%
  raster::calc(fun = sum, na.rm = TRUE) %>% 
  raster::overlay(fp_fch4_sum, fun = function(x, y) y / x)
fc <- fp %>% 
  raster::stack() %>% 
  raster::subset(which(!is.na(site_fp$fc))) %>%
  raster::calc(fun = sum, na.rm = TRUE) %>% 
  raster::overlay(fp_fc_sum, fun = function(x, y) -y / x)
le <- fp %>% 
  raster::stack() %>% 
  raster::subset(which(!is.na(site_fp$le))) %>%
  raster::calc(fun = sum, na.rm = TRUE) %>% 
  raster::overlay(fp_le_sum, fun = function(x, y) y / x)
h <- fp %>% 
  raster::stack() %>% 
  raster::subset(which(!is.na(site_fp$h))) %>%
  raster::calc(fun = sum, na.rm = TRUE) %>% 
  raster::overlay(fp_h_sum, fun = function(x, y) y / x)
tau <- fp %>% 
  raster::stack() %>% 
  raster::subset(which(!is.na(site_fp$tau))) %>%
  raster::calc(fun = sum, na.rm = TRUE) %>% 
  raster::overlay(fp_tau_sum, fun = function(x, y) abs(y) / x)

# Plot the results
fp_flux_plots <- vector("list", length = length(fluxes)) %>%
  lmap(list_modify, vector("list", length = 3)) %>%
  set_names(fluxes) %>%
  map(set_names, c("neutral", "stable", "unstable"))

prepare_fp_df <- function(x) {
  x %>%
    raster::mask(aoi) %>%
    raster::mask(raster::buffer(tower, 3), inverse = TRUE) %>%
    raster::aggregate(fact = 2, fun = mean) %>%
    raster::clamp(
      lower = quantile(., 0.005, na.rm = TRUE), 
      upper = quantile(., 0.995, na.rm = TRUE)
    ) %>%
    raster::rasterToPoints() %>% 
    data.frame() %>%
    rename(z = layer)
}

plot_fp_flux <- function(x) {
  x %>%
    ggplot(aes(x, y)) +
    geom_raster(aes(fill = z)) +
    geom_path(data = aoi, mapping = aes(long, lat), alpha = 0.5) +
    geom_point(data = data.frame(tower), size = 3) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    coord_fixed()
}

for (i in seq_along(fp_flux_plots)) {
  fp_flux_plots[[i]]$neutral <- fp_flux_top[[i]]$neutral %>%
    prepare_fp_df() %>%
    plot_fp_flux() +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "neutral")
  fp_flux_plots[[i]]$stable <- fp_flux_top[[i]]$stable %>%
    prepare_fp_df() %>%
    plot_fp_flux() +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "stable")
  fp_flux_plots[[i]]$unstable <- fp_flux_top[[i]]$unstable %>%
    prepare_fp_df() %>%
    plot_fp_flux() +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "unstable")
}

fp_flux_plots$h$unstable


prepare_fp_rst <- function(x) {
  x %>%
    raster::mask(aoi) %>%
    raster::mask(raster::buffer(tower, 3), inverse = TRUE) %>%
    raster::aggregate(fact = 2, fun = mean) %>%
    raster::clamp(
      lower = quantile(., 0.005, na.rm = TRUE), 
      upper = quantile(., 0.995, na.rm = TRUE)
    )
}

leafsync::sync(
  mapview::viewRGB(raster::mask(rgbn, aoi), r = 1, g = 2, b = 3) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(prepare_fp_rst(fp_flux_top$le$unstable), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(raster::mask(dem, aoi), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(raster::mask(img, aoi), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  ncol = 2
)

leafsync::sync(
  mapview::viewRGB(raster::mask(rgbn, aoi), r = 1, g = 2, b = 3) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(prepare_fp_rst(fp_flux_top$fch4$unstable), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(prepare_fp_rst(fp_flux_top$fc$unstable), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(prepare_fp_rst(fp_flux_top$le$unstable), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(prepare_fp_rst(fp_flux_top$h$unstable), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  mapview::mapView(prepare_fp_rst(fp_flux_top$tau$unstable), legend = FALSE) + 
    mapview::mapView(aoi, fill = FALSE),
  ncol = 3
)

# Perhaps "long tails" of footprints are giving extra weight to grid cells that
# are further away from tower. 
# - solution: "clip" all footprints at e.g. 80%?
# - OR: only assign fch4 up to 80% (i.e. first 80% gets all the flux)?

# Write everything to file
filenames_out <- vector("list", length = length(fluxes)) %>%
  lmap(list_modify, vector("list", length = 3)) %>%
  set_names(fluxes) %>%
  map(set_names, c("neutral", "stable", "unstable"))

for (i in seq_along(filenames_out)) {
  filenames_out[[i]]$neutral <- paste0(
    "fp_", names(filenames_out)[i], "_neutral_", tag_out
  )
  filenames_out[[i]]$stable <- paste0(
    "fp_", names(filenames_out)[i], "_stable_", tag_out
  )
  filenames_out[[i]]$unstable <- paste0(
    "fp_", names(filenames_out)[i], "_unstable_", tag_out
  )
}

for (i in seq_along(filenames_out)) {
  # Rasters
  raster::writeRaster(fp_flux_top[[i]]$neutral, file.path(
    wd, "analysis", "flux_hotspots", paste0(filenames_out[[i]]$neutral, ".tif")
  ))
  raster::writeRaster(fp_flux_top[[i]]$stable, file.path(
    wd, "analysis", "flux_hotspots", paste0(filenames_out[[i]]$stable, ".tif")
  ))
  raster::writeRaster(fp_flux_top[[i]]$unstable, file.path(
    wd, "analysis", "flux_hotspots", paste0(filenames_out[[i]]$unstable, ".tif")
  ))
  # Plots
  ggsave(
    paste0(filenames_out[[i]]$neutral, ".pdf"), fp_flux_plots[[i]]$neutral, 
    path = file.path(wd, "analysis", "flux_hotspots"), width = 5, height = 5
  )
  ggsave(
    paste0(filenames_out[[i]]$stable, ".pdf"), fp_flux_plots[[i]]$stable, 
    path = file.path(wd, "analysis", "flux_hotspots"), width = 5, height = 5
  )
  ggsave(
    paste0(filenames_out[[i]]$unstable, ".pdf"), fp_flux_plots[[i]]$unstable, 
    path = file.path(wd, "analysis", "flux_hotspots"), width = 5, height = 5
  )
  progress_n(i)
}


