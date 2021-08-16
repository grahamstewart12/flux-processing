### ============================================================================
# Site information/metadata compilation & auxilliary data download =============
### ============================================================================

# Set the session information
# These settings will be used to determine which processing options are
# implemented in the script. It will also form part of the saved documentation,
# so options should be adjusted here rather than in the script directly.
settings <- list(
  # Session information
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  year = 2017, # four digit year
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# Purpose: 

# Load the required packages
devtools::load_all("~/Projects/flux/dscalr")
library(units)
library(lubridate)
library(tidyverse)

source("reference/site_metadata.R")


### Initialize script settings & documentation =================================

# Set directory where outputs are stored
# This assumes that the subdirectory structure is already present
wd <- file.path("data", "ERA", settings$year)

# Load metadata file
# md <- pluck(site_metadata, settings$site)

# Set tag for creating output file names
tag_out <- as.character(settings$year)


### ERA-5 hourly data ==========================================================

# Sources that helped in developing this script:
# - https://rpubs.com/boyerag/297592
# - https://www.r-bloggers.com/a-netcdf-4-in-r-cheatsheet/

# Load required packages for this section
library(ecmwfr)
library(ncdf4)

# Set token to local keychain
wf_set_key(
  user = "23795", 
  key = "b945c811-fc4f-44b8-a06a-056d3913df44", 
  service = "cds"
)

# Create data queries (MAX 11 vars per query)
request_land1 <- list(
  format = "netcdf",
  variable = c(
    "2m_temperature", "2m_dewpoint_temperature", "surface_pressure",
    "10m_u_component_of_wind", "10m_v_component_of_wind" 
  ),
  year = as.character(settings$year),
  month = str_pad(1:12, 2, "left", "0"),
  day = str_pad(1:31, 2, "left", "0"),
  time = str_c(0:23, "00", sep = ":") %>% str_pad(5, "left", "0"),
  dataset_short_name = "reanalysis-era5-land",
  target = paste0("era_land1_", tag_out, ".nc"),
  # These bounds include all sites
  area = "39.10/-75.80/39.05/-75.75" # North, West, South, East
)

request_land2 <- list_modify(
  request_land1,
  variable = c(
    "skin_temperature", "lake_mix_layer_temperature", 
    paste0("soil_temperature_level_", 1:3),
    paste0("volumetric_soil_water_layer_", 1:3)
  ),
  target = paste0("era_land2_", tag_out, ".nc")
)

# Single levels data
request_sl <- list_modify(
  request_land1,
  variable = c(
    "friction_velocity", "boundary_layer_height",
    "mean_total_precipitation_rate", 
    paste0(
      "mean_surface_", 
      c("downward_short", "downward_long", "net_short", "net_long"), 
      "_wave_radiation_flux"
    )
  ),
  product_type = "reanalysis",
  dataset_short_name = "reanalysis-era5-single-levels",
  target = paste0("era_sl_", tag_out, ".nc")
)

# Request the data - don't download yet, may take a while to process
file_land1 <- wf_request(request_land1, user = "23795", transfer = FALSE)
writeLines(file_land1$request_id, file.path(wd, "raw", "request_land1.txt"))
file_land2 <- wf_request(request_land2, user = "23795", transfer = FALSE)
writeLines(file_land2$request_id, file.path(wd, "raw", "request_land2.txt"))
file_sl <- wf_request(request_sl, user = "23795", transfer = FALSE)
writeLines(file_sl$request_id, file.path(wd, "raw", "request_sl.txt"))

# Download data when request is finished processing
# - check status: https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
wf_transfer(
  url = readLines(file.path(wd, "raw", "request_land1.txt")), 
  user = "23795",
  path = file.path(wd, "raw"), 
  filename = paste0("era_land1_", tag_out, ".nc"), 
  service = "cds"
)
wf_transfer(
  url = readLines(file.path(wd, "raw", "request_land2.txt")), 
  user = "23795",
  path = file.path(wd, "raw"), 
  filename = paste0("era_land2_", tag_out, ".nc"), 
  service = "cds"
)
wf_transfer(
  url = readLines(file.path(wd, "raw", "request_sl.txt")), 
  user = "23795",
  path = file.path(wd, "raw"), 
  filename = paste0("era_sl_", tag_out, ".nc"), 
  service = "cds"
)


# Import data, tidy, and combine

nc_files <- file.path(
  wd, "raw", paste0(c("era_land1_", "era_land2_", "era_sl_"), tag_out, ".nc")
)

# Create reference data frame with all variable attributes
era_vars <- nc_files %>%
  map(ncmeta::nc_atts) %>%
  set_names(c("land", "land", "sl")) %>%
  bind_rows(.id = "dataset") %>% 
  filter(!variable %in% c("longitude", "latitude", "time", "NC_GLOBAL")) %>% 
  mutate(value = flatten_chr(value)) %>% 
  pivot_wider(-id, names_from = name, values_from = value) %>% 
  type_convert(col_types = cols()) %>%
  select(variable, dataset, orig_units = units, everything(), -`_FillValue`)

era_dims <- nc_files %>%
  map(ncmeta::nc_dims)

# Create data frame of time series for all variables at study location
era_raw <- nc_files %>%
  map(stars::read_ncdf) %>%
  map(as_tibble) %>%
  map(select, -longitude, -latitude) %>%
  reduce(left_join, by = "time") %>%
  rename(timestamp = time) %>%
  # Recent data has two 'experiment versions' - need to remove one
  filter(across(any_of("expver"), ~ .x != 3)) %>%
  select(-any_of("expver"))

# Write unprocessed ERA data to .csv file
era_raw_out <- file.path(wd, "output", paste0("era_unproc_", tag_out, ".csv"))
write_csv(era_raw, era_raw_out)

# Create documentation for unprocessed ERA file
era_raw_docu <- append(
  settings, list(files = nc_files, var_info = as.list(era_vars))
)
# Save documentation
era_raw_docu_out <- str_replace(era_raw_out, ".csv", ".txt")
sink(era_raw_docu_out)
print(era_raw_docu)
sink()

# Harmonize units/formatting with EC data
# - based on Vuichard & Papale (2015) sec. 2.3.1
era_proc <- era_raw %>% 
  # Convert units
  mutate(
    # Temperatures from K to C
    across(c(t2m, d2m, stl1, stl2, stl3, lmlt, skt), set_units, degree_Celsius),
    # Precipitation from kg m-2 s-1 to mm
    mtpr = mtpr * set_units(3600, s) * set_units(1, mm/kg/m^-2),
    # Pressure from Pa to kPa
    sp = set_units(sp, kPa)
  ) %>%
  # Transform variables
  mutate(
    # Calculate outgoing radiation components
    msuwswrf = msdwswrf - msnswrf, 
    msuwlwrf = msdwlwrf - msnlwrf,
    # Calculate vapor pressure deficit
    vpd = bigleaf::e.to.VPD(
      bigleaf::Esat.slope(drop_units(d2m))$Esat, drop_units(t2m)
    ),
    vpd = set_units(vpd * 10, hPa)
  )

# Make sure time series looks correct
qplot(timestamp, t2m, data = drop_units(era_proc))
qplot(timestamp, msdwswrf, data = drop_units(era_proc))

# Write processed ERA data to .csv file
era_proc_out <- file.path(wd, "output", paste0("era_proc_", tag_out, ".csv"))
write_csv(era_proc, era_proc_out)

# Create documentation for processed ERA file
era_proc_docu <- append(
  settings, list(files = era_raw_out, var_info = as.list(era_vars))
)

# Save documentation
era_proc_docu_out <- str_replace(era_proc_out, ".csv", ".txt")
sink(era_proc_docu_out)
print(era_proc_docu)
sink()

# Finished
