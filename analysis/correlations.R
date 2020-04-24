### ============================================================================
# CH4 temperature-water level relationship =====================================
### ============================================================================

# For use with site JLR

# Set the session information
# These settings will be used to determine which processing options are
# implemented in the script. It will also form part of the saved documentation,
# so options should be adjusted here rather than in the script directly.
settings <- list(
  # Session information
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  session_time = lubridate::today(), # date script was run
  session_info = devtools::session_info() # R session info: R, OS, packages
)

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
wd <- file.path("~/Desktop", "DATA", "Flux", settings$site)
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
# Daily data
path_2018d <- latest_version(
  file.path(wd, "2018", "processing_data", "09_daily_averages", "output"), 
  "daily"
)
path_2019d <- latest_version(
  file.path(wd, "2019", "processing_data", "09_daily_averages", "output"), 
  "daily"
)

# Load metadata file
metadata <- read_metadata(
  file.path(wd, "site_info", paste0(settings$site, "_metadata", ".txt"))
)

# Set tag for creating output file names
tag_out <- paste0(settings$site, "_", settings$session_time)


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

# Import flux data - daily
daily_2018 <- read_eddy(path_2018d, stringsAsFactors = FALSE) %>%
  # Strip units and varnames - not needed here
  drop_all_attributes(c("units", "varnames")) %>%
  # Parse timestamp
  mutate(date = ymd(date), year = year(date)) %>%
  filter(year == 2018)

daily_2019 <- read_eddy(path_2019d, stringsAsFactors = FALSE) %>%
  # Strip units and varnames - not needed here
  drop_all_attributes(c("units", "varnames")) %>%
  # Parse timestamp
  mutate(date = ymd(date), year = year(date)) %>%
  filter(year == 2019)

# Combine both years
daily <- full_join(daily_2018, daily_2019) %>%
  mutate(
    year = factor(year),
    grow_seas = between(month(date), 4, 10)
  )


### "Clean" data ===============================================================

site <- site %>%
  mutate(
    fch4 = clean(fch4, qc_fch4_final),
    fc = clean(fch4, qc_fc_final)
    #le = clean(le, qc_le_final),
    #h = clean(h, qc_h_final),
    #tau = clean(tau, qc_tau_final)
  )

# Eliminate days with not enough half-hourly records
# - using threshold of 8 (i.e. need at least 4 hrs of data)
daily <- daily %>%
  mutate(
    fch4 = if_else(fch4_n < 8, NA_real_, fch4),
    fch4_avg = if_else(fch4_n < 8, NA_real_, fch4_avg),
    fch4_avg_day = if_else(fch4_n_day < 4, NA_real_, fch4_avg_day),
    fc = if_else(fc_n < 8 | fc_n_day < 2 | fc_n_night < 2, NA_real_, fc),
    fc_day = if_else(fc_n_day < 4, NA_real_, fc_day),
    fc_avg_day = if_else(fc_n_day < 4, NA_real_, fc_avg_day),
    fc_night = if_else(fc_n_night < 4, NA_real_, fc_night),
    fc_avg_night = if_else(fc_n_night < 4, NA_real_, fc_avg_night),
    # Other necessary variables
    wtd_diff = as.vector(diff(zoo::zoo(wtd_f_avg), na.pad = TRUE))
  )


### FCH4 distribution & transformations ========================================

# Distribution
daily %>%
  #filter(grow_seas == 1) %>%
  ggplot(aes(fch4_avg)) +
  geom_histogram(bins = 25)
# Does not appear that log transformation is needed


### Single variable relationships: Soil temperature ============================

# Soil temperature
daily %>%
  filter(grow_seas == 1) %>%
  ggplot(aes(ts_f_avg, fch4_avg, color = year)) +
  geom_point()

# Pearson correlation
cor(daily$ts_f_avg, daily$fch4_avg, use = "pairwise.complete.obs")

# Simple linear model
ts_lm <- lm(fch4_avg ~ ts_f_avg, data = daily)
summary(ts_lm)

# Nonlinear model 1: exp(a + b * x)
ts_nls1 <- nls(
  fch4_avg ~ exp(a + b * ts_f_avg), data = daily, start = list(a = 0, b = 1)
)
summary(ts_nls1)
modelr::rsquare(ts_nls1, daily)

# Nonlinear model 2: a * exp(b * x)
ts_nls2 <- nls(
  fch4_avg ~ a * exp(b * ts_f_avg), data = daily, start = list(a = 1, b = 1),
  control = nls.control(maxiter = 200)
)
summary(ts_nls2)
modelr::rsquare(ts_nls2, daily)


# Nonlinear model 3: a * b^((x - 10) / 10)
# Long2010
ts_nls2 <- nls(
  fch4_avg ~ a * b^((ts_f_avg - 10) / 10), data = daily, 
  start = list(a = 1, b = 1)
)
summary(ts_nls2)
modelr::rsquare(ts_nls2, daily)

# Seems like nonlinear models 1 & 2 are effectually the same
# - using model 1 from now on

# Points with different curves (650 x 450)
daily %>%
  #modelr::add_predictions(ts_lm, var = "lm") %>%
  modelr::add_predictions(ts_nls1, var = "nls1") %>%
  #modelr::add_predictions(ts_nls2, var = "nls2") %>%
  ggplot(aes(ts_f_avg, fch4_avg)) +
  geom_point(alpha = 0.3, size = 1.25) +
  #geom_line(aes(y = lm), size = 0.75, color = "blue", linetype = 1) +
  geom_line(aes(y = nls1), size = 1, color = "blue", linetype = 1) +
  labs(
    x = expression(Soil~~temp.~~(~degree~C)), 
    y = expression(CH[4]~~flux~~(mg~C~m^{-2}~s^{-1}))
  ) +
  theme_bw() +
  theme_agu

# Residuals over time (nls)***
daily %>%
  #mutate(day = yday(date)) %>%
  modelr::add_residuals(ts_nls1) %>%
  #filter(!is.na(fch4_avg))  %>%
  ggplot(aes(date, resid)) +
  geom_point(alpha = 0.4, size = 1.25) +
  facet_wrap(~ year, ncol = 1, scales = "free_x") +
  geom_hline(aes(yintercept = 0), alpha = 0.5) +
  geom_smooth(span = 0.1) +
  labs(x = NULL, y = "Soil temp. residuals") +
  scale_x_date(date_labels = "%b") +
  theme_bw() +
  theme_agu

# Residuals over time (nls, w/ other vars) (500 x 600)
daily %>%
  #mutate(day = yday(date)) %>%
  modelr::add_residuals(ts_nls1) %>%
  select(date, resid, wtd_f_avg, fc_avg_day) %>%
  filter(between(date, ymd("2018-04-01"), ymd("2019-12-01"))) %>%
  pivot_longer(-date, names_to = "var", values_to = "value") %>%
  mutate(
    var = recode_factor(
      var,
      resid = "Soil temp. residuals", wtd_f_avg = "Water table depth", 
      fc_avg_day = "Daytime NEE", .ordered = TRUE
    )
  ) %>%
  #filter(!is.na(fch4_avg))  %>%
  ggplot(aes(date, value)) +
  geom_point(alpha = 0.4, size = 1.25) +
  facet_wrap(~ var, ncol = 1, scales = "free") +
  geom_hline(aes(yintercept = 0), alpha = 0.5) +
  geom_smooth(span = 0.1) +
  labs(x = NULL, y = NULL) +
  #scale_x_date(date_labels = "%b") +
  theme_bw() +
  theme_agu

# Residuals vs TS (nls)
daily %>%
  modelr::add_residuals(ts_nls1) %>%
  ggplot(aes(ts_f_avg, resid)) +
  geom_point()


### Single variable relationships: Water table depth ===========================

# Water table depth
daily %>%
  filter(grow_seas == 1) %>%
  ggplot(aes(wtd_f_avg, fch4_avg, color = year)) +
  geom_point()


### Multi-variable relationships ===============================================

# TS residuals vs WTD
daily %>%
  modelr::add_residuals(ts_nls1) %>%
  filter(grow_seas, wtd_f_avg > -0.3) %>%
  ggplot(aes(wtd_f_avg, resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_smooth(span = 0.2) +
  theme_bw()

# TS residuals vs FC
daily %>%
  modelr::add_residuals(ts_nls1) %>%
  filter(grow_seas) %>%
  ggplot(aes(fc_avg_day, resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_smooth(span = 0.2) +
  theme_bw()

daily %>%
  filter(grow_seas == 1) %>%
  #mutate(ts_bin) %>%
  ggplot(aes(ts_f_avg, fch4_avg, color = year)) +
  geom_point()

daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  #mutate(wtd_class = cut_width(wtd_f_avg, width = 0.4, center = 0)) %>%
  mutate(
    wtd_class = cut_interval(wtd_f_avg, n = 2) %>% 
      fct_recode(
        "[-0.78,-0.16]  n=19" = "[-0.772,-0.158]", 
        "(-0.16,0.46]  n=258" = "(-0.158,0.456]"
      )
  ) %>%
  ggplot(aes(ts_f_avg, fch4_avg)) +
  #facet_wrap(~ year) +
  geom_point(size = 2) +
  geom_smooth(aes(color = wtd_class), method = "lm", se = FALSE) +
  #geom_smooth(
  #  aes(color = wtd_class), method = "nls", formula = y ~ a * exp(b * x), 
  #  se = FALSE, size = 1.25, method.args = list(start = list(a = 1, b = 1))
  #) +
  #scale_color_manual(
  #  values = RColorBrewer::brewer.pal(name = "YlGnBu", n = 3),
  #  name = expression(Water~~level~~(m))
  #) +
  labs(
    x = expression(Soil~~temp.~~(~degree~C)), 
    y = expression(CH[4]~~flux~~(mg~~C~~m^{-2}~~s^{-1})),
    color = "Water level (m)"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.25, 0.65),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  mutate(wtd_class = cut_width(wtd_f_avg, width = 0.5, center = 0.25)) %>%
  ggplot(aes(ts_f_avg, fch4_avg)) +
  geom_point(size = 2) +
  geom_smooth(aes(color = wtd_class), method = "lm", se = FALSE) +
  labs(
    x = expression(Soil~~temp.~~(~degree~C)), 
    y = expression(CH[4]~~flux~~(mg~~C~~m^{-2}~~s^{-1})),
    color = "Water level (m)"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.25, 0.65),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  mutate(wtd_class = cut_width(wtd_f_avg, width = 0.4, center = 0.2)) %>%
  ggplot(aes(ts_f_avg, fch4_avg)) +
  geom_point(size = 2) +
  geom_smooth(aes(color = wtd_class), method = "lm", se = FALSE) +
  labs(
    x = expression(Soil~~temp.~~(~degree~C)), 
    y = expression(CH[4]~~flux~~(mg~~C~~m^{-2}~~s^{-1})),
    color = "Water level (m)"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.25, 0.65),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

daily %>%
  filter(grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)) %>%
  #mutate(wtd_class = cut_width(wtd_f_avg, width = 0.4, center = 0)) %>%
  mutate(wtd_class = cut_interval(wtd_f_avg, n = 3)) %>%
  group_by(wtd_class) %>%
  summarize(n())

# FCH4 ~ TS
daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  do(lm = lm(fch4_avg ~ ts_f_avg + fc_avg_day + wtd_f_avg, data = .)) %>%
  pull(lm) %>% 
  pluck(1) %>%
  summary()
# rsq = 0.5582

# FCH4 ~ TS (NLS)
nls1 <- daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  do(nls = nls(
    fch4_avg ~ exp(a + b * ts_f_avg), data = ., start = list(a = 0, b = 1)
  )) %>%
  pull(nls) %>% 
  pluck(1)
summary(nls1)
modelr::rsquare(nls1, filter(daily,
  grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
))
# rsq = 0.6094068

# FCH4 ~ TS + WTD
daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  do(lm = lm(fch4_avg ~ ts_f_avg + wtd_f_avg, data = .)) %>%
  pull(lm) %>% 
  pluck(1) %>%
  summary()
# rsq = 0.6399

# FCH4 ~ TS + WTD (NLS)
nls2 <- daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
  ) %>%
  do(nls = nls(
    fch4_avg ~ exp(a + b * ts_f_avg) + c * wtd_f_avg + d, data = ., 
    start = list(a = 0, b = 1, c = 1, d = 0)
  )) %>%
  pull(nls) %>% 
  pluck(1)
summary(nls2)
modelr::rsquare(nls2, filter(daily, 
  grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg)
))
# rsq = 0.6620864

# FCH4 ~ TS (2018)
daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2018
  ) %>%
  do(lm = lm(fch4_avg ~ ts_f_avg, data = .)) %>%
  pull(lm) %>% 
  pluck(1) %>%
  summary()
# rsq = 0.5625

# FCH4 ~ TS (NLS, 2018)
nls3 <- daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2018
  ) %>%
  do(nls = nls(
    fch4_avg ~ exp(a + b * ts_f_avg), data = ., start = list(a = 0, b = 1)
  )) %>%
  pull(nls) %>% 
  pluck(1)
summary(nls3)
modelr::rsquare(nls3, filter(
  daily, grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
  year == 2018
))
# rsq = 0.5792944

# FCH4 ~ TS + WTD (2018)
daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2018
  ) %>%
  do(lm = lm(fch4_avg ~ ts_f_avg + wtd_f_avg, data = .)) %>%
  pull(lm) %>% 
  pluck(1) %>%
  summary()
# rsq = 0.5591

# FCH4 ~ TS + WTD (2018, NLS)
nls4 <- daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2018
  ) %>%
  do(nls = nls(
    fch4_avg ~ exp(a + b * ts_f_avg) + c * wtd_f_avg + d, data = ., 
    start = list(a = 0, b = 1, c = 1, d = 0)
  )) %>%
  pull(nls) %>% 
  pluck(1)
summary(nls4)
modelr::rsquare(nls4, filter(
  daily, grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
  year == 2018
))
# rsq = 0.586304

# FCH4 ~ TS (2019)
daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2019
  ) %>%
  do(lm = lm(fch4_avg ~ ts_f_avg, data = .)) %>%
  pull(lm) %>% 
  pluck(1) %>%
  summary()
# rsq = 0.597

# FCH4 ~ TS (NLS, 2019)
nls5 <- daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2019
  ) %>%
  do(nls = nls(
    fch4_avg ~ exp(a + b * ts_f_avg), data = ., start = list(a = 0, b = 1)
  )) %>%
  pull(nls) %>% 
  pluck(1)
summary(nls5)
modelr::rsquare(nls5, filter(
  daily, grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
  year == 2019
))
# rsq = 0.5952932

# FCH4 ~ TS + WTD (2019)
daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2019
  ) %>%
  do(lm = lm(fch4_avg ~ ts_f_avg + wtd_f_avg, data = .)) %>%
  pull(lm) %>% 
  pluck(1) %>%
  summary()
# rsq = 0.732

# FCH4 ~ TS + WTD (2019, NLS)
nls6 <- daily %>%
  filter(
    grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
    year == 2019
  ) %>%
  do(nls = nls(
    fch4_avg ~ exp(a + b * ts_f_avg) + c * wtd_f_avg + d, data = ., 
    start = list(a = 0, b = 1, c = 1, d = 0)
  )) %>%
  pull(nls) %>% 
  pluck(1)
summary(nls6)
modelr::rsquare(nls6, filter(
  daily, grow_seas == 1, !is.na(wtd_f_avg), !is.na(ts_f_avg), !is.na(fch4_avg),
  year == 2019
))
# rsq = 0.7391301


# Analysis from Olefeldt et al. 2017

# We used residuals from the stepwise linear model to assess potential time lags in relationships between water table position and ER and FCH4. Coefficients of determination were determined for linear correlations between model residuals and the net shift in water table position over a time period preceding a flux measurement. Time periods for lag effects ranging from 1 to 50 days were considered. Inter-annual timelag effects, i.e., the effect of the wetness of the preceding year on the current year fluxes, were assessed by linear correlations between the current year average model residual within each plot and its average water table position during the preceding year.
