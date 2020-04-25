### ============================================================================
# Run a script =================================================================
### ============================================================================

# Set the session information
# These settings will be used to determine which processing options are
# implemented in the script. It will also form part of the saved documentation,
# so options should be adjusted here rather than in the script directly.
rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLA", # three letter site code
  year = 2019, # four digit year
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# Set the script name (see below for options)
script <- "01_combine_eddypro"

# Run the script
source(file.path(
  "~/Desktop/DATA/Flux/tools/engine", "processing", paste0(script, ".R")
))


# Available scripts:
# 01_combine_eddypro
# 02_correct_eddypro
# 03_biomet_qc
# 04_biomet_gapfill
# 05_footprint
# 06_flux_qc
# 07_footprint_cover
# 08_daily_summary
