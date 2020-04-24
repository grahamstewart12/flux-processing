### ============================================================================
# Correlations between time series =============================================
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

# Import flux data - half-hourly
site_2018 <- read_eddy(path_2018, stringsAsFactors = FALSE) %>%
  # Strip units and varnames - not needed here
  drop_all_attributes(c("units", "varnames")) %>%
  # Prevent integer columns from messing with join
  mutate_if(is.integer, as.numeric) %>%
  # Parse timestamp
  mutate(timestamp = ymd_hms(timestamp, tz = "Etc/GMT+5"))
glimpse(site_2018) # check if data loaded correctly
site_2018 %>% 
  summarize(first(timestamp), last(timestamp)) # Starts/ends as expected?

site_2019 <- read_eddy(path_2019, stringsAsFactors = FALSE) %>% 
  # Strip units and varnames - not needed here
  drop_all_attributes(c("units", "varnames")) %>%
  # Prevent integer columns from messing with join
  mutate_if(is.integer, as.numeric) %>%
  # Parse timestamp
  mutate(timestamp = ymd_hms(timestamp, tz = "Etc/GMT+5")) 
glimpse(site_2019) # check if data loaded correctly
site_2019 %>% 
  summarize(first(timestamp), last(timestamp)) # Starts/ends as expected?

# Combine both years
site <- full_join(site_2018, site_2019)


### "Clean" data ===============================================================

site_c <- site %>%
  mutate(
    grow_seas = between(month(timestamp), 4, 10),
    fch4 = clean(fch4, qc_fch4_final),
    fc = clean(fch4, qc_fc_final),
    le = clean(le, qc_le_final),
    h = clean(h, qc_h_final),
    tau = clean(tau, qc_tau_final)
  )


### Half-hourly correlations ===================================================

site_c %>%
  ggplot(aes(ts_f, fch4)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE)

# All
site_c %>%
  do(
    ts_lm = lm(fch4 ~ ts_f, data = .),
    ts_nls = nls(
      fch4 ~ exp(a + b * ts_f), data = ., start = list(a = 0, b = 1)
    ),
    wtd_lm = lm(fch4 ~ ts_f, data = .)
  ) %>%
  mutate(nls_rsq = map(nls, ~modelr::rsquare(., site_c))) %>%
  pull(lm) %>% pluck(1) %>%
  summary()

# Growing season
site_c %>%
  filter(grow_seas) %>%
  do(
    data = .data,
    ts_lm = lm(fch4 ~ ts_f, data = .),
    ts_nls = nls(
      fch4 ~ exp(a + b * ts_f), data = ., start = list(a = 0, b = 1)
    ),
    wtd_lm = lm(fch4 ~ wtd_f, data = .),
    multi_lm = lm(fch4 ~ ts_f + wtd_f, data = .),
    multi_nls = nls(
      fch4 ~ exp(a + b * ts_f) + c * wtd_f, data = ., 
      start = list(a = 0, b = 1, c = 1)
    )
  ) %>%
  pivot_longer(-data, names_to = "model", values_to = "object") %>%
  transmute(
    model = model,
    glanced = map(object, broom::glance),
    tidied = map(object, broom::tidy),
    coef_1 = map_dbl(tidied, pluck, "estimate", 2),
    coef_2 = map_dbl(tidied, pluck, "estimate", 3, .default = NA_real_),
    p_1 = map_dbl(tidied, pluck, "p.value", 2),
    p_2 = map_dbl(tidied, pluck, "p.value", 3, .default = NA_real_),
    rsq = map2_dbl(object, data, modelr::rsquare)
  )

# All TS lm:  beta = 0.015171, r2 = 0.5455
# All TS nls: beta = 0.141950, r2 = 0.5964

# Growing season TS lm:  beta = 0.021247, r2 = 0.4325
# Growing season TS nls: beta = 0.148124, r2 = 0.4475

### Daily correlations =========================================================

daily <- site_c %>%
  mutate(date = date(timestamp)) %>%
  group_by(date, grow_seas) %>%
  summarize(
    fch4 = mean(fch4, na.rm = TRUE),
    ts = mean(ts_f, na.rm = TRUE),
    wtd = mean(wtd_f, na.rm = TRUE)
  ) %>%
  ungroup()

daily %>%
  ggplot(aes(ts, fch4)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

daily %>%
  do(
    lm = lm(fch4 ~ ts, data = .),
    nls = nls(fch4 ~ exp(a + b * ts), data = ., start = list(a = 0, b = 1))
  ) %>%
  mutate(nls_rsq = map(nls, ~modelr::rsquare(., daily))) %>%
  pull(lm) %>% pluck(1) %>%
  summary()

# Growing season
daily %>%
  filter(grow_seas) %>%
  do(
    lm = lm(fch4 ~ ts, data = .),
    nls = nls(fch4 ~ exp(a + b * ts), data = ., start = list(a = 0, b = 1))
  ) %>%
  mutate(nls_rsq = map(nls, ~modelr::rsquare(., filter(daily, grow_seas)))) %>%
  pull(lm) %>% pluck(1) %>%
  summary()

# TS lm:  beta = 0.021254, r2 = 0.6834
# TS nls: beta = 0.156733, r2 = 0.7535

# Growing season TS lm:  beta = 0.014306, r2 = 0.5971
# Growing season TS nls: beta = 0.149902, r2 = 0.6132

### Monthly correlations =======================================================

monthly <- site_c %>%
  mutate(month = month(timestamp), year = year(timestamp)) %>%
  group_by(month, year, grow_seas) %>%
  summarize(
    fch4 = mean(fch4, na.rm = TRUE),
    ts = mean(ts_f, na.rm = TRUE)
  ) %>%
  ungroup()

monthly %>%
  ggplot(aes(ts, fch4)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# All
monthly %>%
  do(
    lm = lm(fch4 ~ ts, data = .),
    nls = nls(fch4 ~ exp(a + b * ts), data = ., start = list(a = 0, b = 1))
  ) %>%
  mutate(nls_rsq = map(nls, ~modelr::rsquare(., monthly))) %>%
  pull(lm) %>% pluck(1) %>%
  summary()

# Growing season
monthly %>%
  filter(grow_seas) %>%
  do(
    lm = lm(fch4 ~ ts, data = .),
    nls = nls(fch4 ~ exp(a + b * ts), data = ., start = list(a = 0, b = 1))
  ) %>%
  mutate(nls_rsq = map(nls, ~modelr::rsquare(., filter(monthly, grow_seas)))) %>%
  pull(lm) %>% pluck(1) %>%
  summary()

# TS lm:  beta = 0.01384, r2 = 0.8005
# TS nls: beta = 0.14754, r2 = 0.8546

# Growing season TS lm:  beta = 0.02197, r2 = 0.7472
# Growing season TS nls: beta = 0.13891, r2 = 0.7277

### Seasonal correlations ======================================================

seasonal <- site_c %>%
  mutate(season = season(timestamp), year = year(timestamp)) %>%
  group_by(season, year) %>%
  summarize(
    fch4 = mean(fch4, na.rm = TRUE),
    ts = mean(ts_f, na.rm = TRUE)
  ) %>%
  ungroup()

seasonal %>%
  ggplot(aes(ts, fch4)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

seasonal %>%
  do(
    lm = lm(fch4 ~ ts, data = .),
    nls = nls(fch4 ~ exp(a + b * ts), data = ., start = list(a = 0, b = 1))
  ) %>%
  mutate(nls_rsq = map(nls, ~modelr::rsquare(., seasonal))) %>%
  pull(lm) %>% pluck(1) %>%
  summary()

# TS lm:  beta = 0.01434, r2 = 0.8436
# TS nls: beta = 0.12546, r2 = 0.8118
