### ============================================================================
# Flux quality control =========================================================
### ============================================================================


# Purpose: 

# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/Desktop/RESEARCH/fluxtools", quiet = TRUE)
library(openeddy)
library(lubridate)
library(tidyverse)

source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")
source("~/Desktop/DATA/Flux/tools/engine/flux_qc.R")

# Initialize script settings & documentation
script_init(settings, site_metadata)


### Load required input data ===================================================

# EddyPro
site <- read.csv(site_input, stringsAsFactors = FALSE)

# Add timestamp components
site <- site %>% 
  mutate(timestamp = ymd_hms(timestamp, tz = md$tz_name)) %>%
  add_time_comps()

# Footprint
#fp <- read.csv(fp_input, stringsAsFactors = FALSE)

# Add timestamp components
#fp <- mutate(fp, timestamp = ymd_hms(timestamp, tz = md$tz_name))

# Bind footprint data to main data frame
#site <- left_join(site, fp)

# FFP
ffp <- read.csv(ffp_input, stringsAsFactors = FALSE)

# Add timestamp components
ffp <- mutate(ffp, timestamp = ymd_hms(timestamp, tz = md$tz_name))

# Bind ffp data to main data frame
site <- left_join(site, ffp, by = "timestamp")

# Add 'daytime' variable based on ppfd_in_f
site <- mutate(site, daytime = sw_in_f > 20)


### Error levels ===============================================================

# Load RFlux functions
rflux <- list.files("~/Desktop/DATA/Flux/tools/RFlux", full.names = TRUE)
map(rflux, source)

# Create workset for RFlux input
ws_path <- file.path(
  dirname(path_out), paste0("rflux_workset_", tag_out, ".csv")
)
workset <- ecworkset(
  ep_path, qc_path, md_path, st_path, 
  dirname(tools::file_path_sans_ext(ws_path)), 
  basename(tools::file_path_sans_ext(ws_path))
)
ws_path <- latest_version(dirname(path_out), "rflux_workset")

# Pre-screen fluxes for periods when rain may have affected sensors
workset <- workset %>%
  bind_cols(select(site, p_rain_f)) %>%
  mutate_at(
    vars(H, LE, CO2flux, CH4flux, CO2str, CH4str),
    ~ if_else(!is.na(p_rain_f) & p_rain_f != 0, NA_real_, .)
  ) %>%
  select(-p_rain_f)

# Run the quality control procedure
clean <- cleanFlux(
  workset, md_path, path_output = NULL, FileName = NULL, 
  plotQC = FALSE, storage = TRUE
) 

clean <- clean %>%
  as_tibble() %>%
  mutate(timestamp = ymd_hm(TIMESTAMP_END, tz = md$tz_name)) %>%
  select(timestamp, everything()) %>%
  rename_all(tolower) %>%
  left_join(select(site, timestamp, daytime), by = "timestamp")

# Evidence of systematic error in raw data

# Percent fluxes removed due to severe error
clean %>% summarize(p = 1 - (length(na.omit(nee)) / length(na.omit(fc))))
clean %>% summarize(p = 1 - (length(na.omit(fch4)) / length(na.omit(fm))))

# Examine flags
clean %>%
  mutate(qc = factor(nee_data_flag)) %>%
  ggplot(aes(timestamp, fc, color = qc)) +
  facet_wrap(~ daytime) +
  geom_point() +
  ylim(-40, 40)
clean %>%
  mutate(qc = factor(fch4_data_flag)) %>%
  ggplot(aes(timestamp, fm, color = qc)) +
  facet_wrap(~ daytime) +
  geom_point() +
  ylim(-0.75, 1.25)

# Outlier detection
# - seasonal-trend decomposition based on LOESS

# Percent fluxes identified as outliers
clean %>% 
  filter(!is.na(nee_outlying_flag)) %>% 
  summarize(p = 1 - (n() - sum(nee_outlying_flag)) / n())
clean %>% 
  filter(!is.na(fch4_outlying_flag)) %>% 
  summarize(p = 1 - (n() - sum(fch4_outlying_flag)) / n())

# Examine flags
clean %>%
  mutate(
    ol = factor(nee_outlying_flag),
    fc = if_else(nee_data_flag == 2, NA_real_, fc), qc = factor(nee_data_flag)
  ) %>%
  ggplot(aes(timestamp, fc, color = ol, shape = qc)) +
  facet_wrap(~ daytime) +
  geom_point()
clean %>%
  mutate(
    ol = factor(fch4_outlying_flag),
    fm = if_else(fch4_data_flag == 2, NA_real_, fm), qc = factor(fch4_data_flag)
  ) %>%
  ggplot(aes(timestamp, fm, color = ol, shape = qc)) +
  facet_wrap(~ daytime) +
  geom_point()

# Examine data without outliers
clean %>%
  ggplot(aes(timestamp, nee)) +
  geom_point()
clean %>%
  ggplot(aes(timestamp, fch4)) +
  geom_point()

# Overall percent of fluxes remaining
clean %>% summarize(p = 1 - (length(na.omit(nee)) / length(na.omit(fc))))
clean %>% summarize(p = 1 - (length(na.omit(fch4)) / length(na.omit(fm))))

# Flag high uncertainty (i.e. low quality) in remaining fluxes
clean %>% 
  mutate(qc = factor(pull(site, h_randunc_hf) > 25)) %>%
  ggplot(aes(timestamp, h, color = qc)) +
  geom_point()
clean %>% 
  mutate(qc = factor(pull(site, le_randunc_hf) > 50)) %>%
  ggplot(aes(timestamp, le, color = qc)) +
  geom_point()
clean %>% 
  mutate(qc = factor(pull(site, fc_randunc_hf) > 5)) %>%
  ggplot(aes(timestamp, nee, color = qc)) +
  geom_point()
clean %>% 
  mutate(qc = factor(pull(site, fch4_randunc_hf) > 0.1)) %>%
  ggplot(aes(timestamp, fch4, color = qc)) +
  geom_point()
clean <- clean %>%
  mutate(
    h_rand_unc = pull(site, h_randunc_hf),
    le_rand_unc = pull(site, le_randunc_hf),
    nee_rand_unc = pull(site, fc_randunc_hf),
    fch4_rand_unc = pull(site, fch4_randunc_hf),
    
    h_unc_flag = apply_thr(h_rand_unc, c(25, 25)),
    le_unc_flag = apply_thr(le_rand_unc, c(50, 50)),
    nee_unc_flag = apply_thr(nee_rand_unc, c(5, 5)),
    fch4_unc_flag = apply_thr(fch4_rand_unc, c(0.1, 0.1))
  )

# Tokoro & Kuwae 2018 aquatic EC outlier detection
clean %>% 
  mutate(
    CO2_SIGMA = CO2_SIGMA / 400.1349, 
    H2O_SIGMA = H2O_SIGMA / 12315.13,
    SIGMA_MAX = map2_dbl(CO2_SIGMA, H2O_SIGMA, max)
  ) %>%
  ggplot(aes(SIGMA_MAX, NEE, color = NEE_RESID > 15)) +
  geom_point()
clean %>% 
  filter(!is.na(NEE)) %>% 
  arrange(desc(NEE_RESID)) %>% 
  select(NEE_RESID, CO2_SIGMA, H2O_SIGMA) %>% 
  mutate(
    co2_sigma = CO2_SIGMA / 400.1349, 
    H2O_SIGMA = H2O_SIGMA / 12315.13,
    SIGMA_MAX = map2_dbl(co2_sigma, H2O_SIGMA, max)
  ) %>%
  slice(1:15)

# Footprint coverage of AOI
clean %>%
  mutate(nee = clean(nee, nee_unc_flag), flag = factor(foot_flag)) %>%
  left_join(ffp) %>%
  ggplot(aes(timestamp, nee, color = flag)) +
  geom_point(size = 1)
clean %>%
  mutate(nee = clean(nee, nee_unc_flag)) %>%
  left_join(ffp) %>%
  ggplot(aes(timestamp, nee, color = phi_ffp > 0.80)) +
  geom_point(size = 1)
clean %>%
  mutate(fch4 = clean(fch4, fch4_unc_flag), flag = factor(foot_flag)) %>%
  left_join(ffp) %>%
  ggplot(aes(timestamp, fch4, color = flag)) +
  geom_point(size = 1)
clean %>%
  mutate(fch4 = clean(fch4, fch4_unc_flag)) %>%
  left_join(ffp) %>%
  ggplot(aes(timestamp, fch4, color = phi_ffp > 0.80)) +
  geom_point(size = 1)


clean <- clean %>%
  mutate(
    phi = pull(ffp, phi_ffp),
    phi_flag = apply_thr(phi, c(0.80, 0.80), flag = "lower")
  )

site %>% bind_cols(fch4_qc) %>% filter(qc_fch4_max != 2) %>% plot_dygraph(fch4)

clean_st %>%
  mutate(timestamp = ymd_hm(TIMESTAMP_END)) %>%
  bind_cols(select(site, phi, phi_ffp, night, fc_ss, fch4_randunc_hf)) %>%
  #filter(NEE_DATA_FLAG != 2, NEE_OUTLYING_FLAG != 1) %>%
  drop_na(NEE) %>% filter(phi_ffp > 0.80) %>% 
  group_by(night) %>% 
  summarize(
    low = quantile(NEE, c(0.005, 0.995), na.rm = TRUE)[1], 
    high = quantile(NEE, c(0.005, 0.995), na.rm = TRUE)[2]
  )


### Visual precheck ============================================================

# Display variables that can help identify problems with instruments
names <- c(
  "u_rot", "v_rot", "w_unrot", "w_rot", "sonic_temperature", "max_wind_speed",
  "tau", "ustar", "h", "le", "fc", "fch4",
  "u_var", "v_var", "w_var", "ts_var", "h2o_var", "co2_var", "ch4_var",
  "rand_err_tau", "rand_err_h", "rand_err_le", "rand_err_fc", "rand_err_fch4",
  "tau_scf", "h_scf", "le_scf", "co2_scf", "ch4_scf",
  "u_spikes", "v_spikes", "w_spikes", "ts_spikes", "co2_spikes", "ch4_spikes",
  "h2o_v_adv", "co2_v_adv", "ch4_v_adv",
  "co2_mixing_ratio", "h2o_mixing_ratio", "ch4_mixing_ratio", 
  "co2_time_lag", "h2o_time_lag", "ch4_time_lag",
  "x_peak_h00", "x_90perc_h00"
)


### Simplified QC strategy =====================================================

# Based on which tests/stats are actually effective at removing spurious values
# - not based on "standard procedure" because our sites are not "standard"
# - operating philosophy = glean as much information as possible from ecosystem

clean_1 <- site
site_qc <- select(site, timestamp)

# FC
# Choose theoretical flags 
site %>%
  # Storage correction, test possible flags
  mutate(
    fc = add_st(fc, st = sc_single), 
    flag = factor(stringr::str_sub(w_vm97_test, 5, 5))
  ) %>%
  # Clean data
  filter(
    !is.na(fc), phi_ffp > 0.85, ustar > 0.05, !p_rain > 0, fc_ss < 30, 
    fc_scf < 3, inst_li7500_agc_or_rssi >= 70, 
    stringr::str_sub(w_vm97_test, 5, 5) < 1
  ) %>%
  ggplot(aes(timestamp, fc, color = flag)) +
  geom_point(size = 1) +
  facet_wrap(~ night)
# MAYBE: 
# - **w_vm97_test[5] (skw/kur I think) - need to check JLR though before adding
# - w_vm97_test[6] (discontinuities)

# Set theoretical flags
site <- mutate(site, fc = add_st(fc, st = sc_single)) # add storage
site_qc <- site %>%
  transmute(
    qc_fc_ss = apply_thr(fc_ss, c(29, 29), flag = "higher"),
    qc_fc_agc = apply_thr(inst_li7500_agc_or_rssi, c(70, 70), flag = "lower"),
    qc_fc_phi = apply_thr(phi, c(0.7, 0.7), flag = "lower"),
    qc_fc_rain = apply_thr(p_rain, c(0, 0), flag = "higher"),
    qc_fc_skwkur = if_else(stringr::str_sub(w_vm97_test, 5, 5) == 1, 2L, 0L),
    qc_fc_ustar = apply_thr(ustar, c(0.05, 0.05), flag = "lower")
  ) %>%
  mutate_all(replace_na, 0) %>%
  mutate(qc_fc_comp = combn_QC(., names(.))) %>%
  bind_cols(site_qc, .)

# Empirical flags
# Spike detection
site_qc <- site %>%
  left_join(site_qc) %>%
  mutate(GR = sw_in_f) %>%
  filter(!is.na(GR)) %>%
  transmute(
    timestamp = timestamp,
    qc_fc_hhspikes = despikeLF(
      ., var = "fc", qc_flag = "qc_fc_comp",
      var_thr = NULL, light = "GR", night_thr = 20, z = 5.5
    ) %>% replace_na(0)
  ) %>%
  right_join(site_qc)

# Final flag
site_qc <- site_qc %>%
  mutate(
    qc_fc_final = combn_QC(., c("qc_fc_comp", "qc_fc_hhspikes"))
  )

# Check final cleaned time series
site %>%
  left_join(site_qc) %>%
  mutate(fc = clean(fc, qc_fc_final)) %>%
  ggplot(aes(timestamp, fc)) +
  geom_point()


# FCH4
# Convert nmol to umol
site <- mutate(
  site, 
  fch4 = fch4 / 1000, 
  sch4_single = sch4_single / 1000,
  ch4_mixing_ratio = ch4_mixing_ratio / 1000
)
site %>%
  # Look at storage correction, test possible flags
  mutate(
    fch4 = add_st(fch4, st = sch4_single), 
    flag = factor(stringr::str_sub(ch4_vm97_test, 5, 5))
  ) %>%
  # Clean data
  filter(
    !is.na(fch4), fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05,
    ch4_mixing_ratio < 5
  ) %>%
  ggplot(aes(timestamp, fch4, color = flag)) +
  geom_point()

# Set theoretical flags
site <- mutate(site, fch4 = add_st(fch4, st = sch4_single)) # add storage
site_qc <- site %>%
  transmute(
    qc_fch4_ss = apply_thr(fch4_ss, c(29, 29), flag = "higher"),
    qc_fch4_nsr = apply_thr(fch4_nsr, c(2, 2), flag = "higher"),
    qc_fch4_conc = apply_thr(ch4_mixing_ratio, c(5, 5), flag = "higher"),
    qc_fch4_phi = apply_thr(phi, c(0.7, 0.7), flag = "lower"),
    qc_fch4_ustar = apply_thr(ustar, c(0.05, 0.05), flag = "lower")
  ) %>%
  mutate_all(replace_na, 0) %>%
  mutate(qc_fch4_comp = combn_QC(., names(.))) %>%
  bind_cols(site_qc, .)

# Empirical flags
# Absolute limits from extreme statistical outliers
fch4_lims <- site %>%
  left_join(site_qc) %>%
  mutate(fch4 = clean(fch4, qc_fch4_comp)) %>%
  summarize(list(quantile(fch4, c(0.002, 0.998), na.rm = TRUE))) %>%
  pull()
fch4_lims
site_qc <- site %>%
  transmute(
    qc_fch4_lims = apply_thr(
      fch4, purrr::pluck(fch4_lims, 1), flag = "outside"
    ) %>% replace_na(0)
  ) %>%
  bind_cols(site_qc, .)

# Final flag
site_qc <- site_qc %>%
  mutate(
    qc_fch4_final = combn_QC(., c("qc_fch4_comp", "qc_fch4_lims"))
  )

# Check final cleaned time series
site %>%
  left_join(site_qc) %>%
  mutate(fch4 = clean(fch4, qc_fch4_final)) %>%
  ggplot(aes(timestamp, fch4)) +
  geom_point()

# Bind qc flags to main data frame
full <- left_join(site, site_qc)

# Save dataset
full_out <- paste0(
  wd, "/08_flux_qc/output/", "flux_qc_full_", tag_out, ".csv"
)
write_eddy(full, full_out)

# Gather & save essentials dataset
essentials <- site %>%
  select_at(vars(
    timestamp, sw_in_pot, night, h, le, fc, fh2o, fch4, u_unrot, v_unrot,
    w_unrot, ws, wd, wd_sigma, ustar, tke, mo_length, zl, bowen, tstar, t_sonic,
    ta_ep, pa_ep, rh_ep, air_density, specific_humidity, vpd_ep,
    ends_with("_mixing_ratio"), co2, h2o, ch4, ends_with("_sigma"), 
    starts_with("fetch_"), ends_with("_H00"), phi, ends_with("_f"),
    lw_in, lw_out, ppfd_in, rh, netrad, sw_out, g, swc, tw, tw_surf, tw_bot, 
    do_surf, do_bot, starts_with("qc_"), ends_with("_fqc"), ends_with("_fmeth"),
    h_ssitc_test, le_ssitc_test
  )) %>%
  left_join(select_at(site_qc, vars(timestamp, ends_with("final"))))

essentials_out <- paste0(
  wd, "/08_flux_qc/output/", "flux_qc_essentials_", tag_out, ".csv"
)
write_eddy(essentials, essentials_out)



### Quality checking - compare results =========================================

clean %>%
  bind_cols(site_qc) %>%
  mutate(
    fc1 = clean(fc, qc_fc_final),
    fc2 = clean(NEE)
  )
select(timestamp, fc, qc_fc_final)


### Quality checking - preliminary flags =======================================

# Set the reference list to determine which flags apply to which fluxes
combn_flags <- list(
  tau = "_tau|_sa|_all",
  h = "_h|_sa|_all",
  le = "_le|_h2o|_sa|_irga|_all",
  fc = "_fc|_co2|_sa|_irga|_all"
)

# Get a subsetted data frame with all recognizable QC information
site_qc <- site %>%
  select_at(vars(
    ends_with("_ssitc_test"), 
    ends_with("_ss_test"), 
    ends_with("_vm97_test"), # 800011199 
    # - spikes, ampres, dropout, abslim, skw/kur, discont, timelag, angle, steady
    ends_with("_itc"), ends_with("_ss"), # Foken stats used to calculate flags
    ends_with("_nsr"), # Mahrt 1998 Nonstationarity Ratios
    ends_with("_corrdiff"), # Correlation differences with/without repeat values
    ends_with("_zcd"), # Zero-Counts on Differenced variables
    ends_with("_kid"), # Kurtosis Index on Differenced variables
    ends_with("_hf")
  ))

# Extract the standardized tests/filters
# - these highly likely indicate problems with sensors, are used in all cases 
QC <- site %>%
  transmute(
    timestamp = timestamp,
    # Import EddyPro flags from main data frame
    #qc_tau = tau_ssitc_test, 
    #qc_le = le_ssitc_test, 
    #qc_h = h_ssitc_test, 
    #qc_fc = fc_ssitc_test,
    qc_tau_ss = apply_thr(tau_ss, c(30, 30)),
    qc_h_ss = apply_thr(h_ss, c(30, 30)),
    qc_le_ss = apply_thr(fh2o_ss, c(30, 30)),
    qc_fc_ss = apply_thr(fc_ss, c(30, 30)),
    # Excessive spectral correction factor
    qc_tau_scf = apply_thr(tau_scf, c(2, 3)),
    qc_le_scf = apply_thr(le_scf, c(2, 3)),
    qc_fc_scf = apply_thr(fc_scf, c(2, 3)),
    # Residual vertical wind component
    qc_all_wresid = apply_thr(w_unrot, c(0.35, 0.35)),
    # Add flags relating to irga components
    qc_irga_rssi = apply_thr(
      inst_li7500_agc_or_rssi, c(70, 70), flag = "lower"
    ),
    qc_co2_runs = flag_runs(co2_mixing_ratio),
    qc_co2_ppm = apply_thr(co2_mixing_ratio, c(300, 900), flag = "out"),
    qc_h2o_runs = flag_runs(h2o_mixing_ratio),
    # Add sonic anemometer equal length run flags
    qc_sa_t_runs = flag_runs(t_sonic),
    qc_sa_u_runs = flag_runs(u_unrot),
    qc_sa_v_runs = flag_runs(v_unrot),
    qc_sa_w_runs = flag_runs(w_unrot),
    # Add precipitation flags
    # - only accept gaps that were filled using nearby sites (NOT downscaled)
    # - replace NA values in p_rain with 0 so these do not get flagged
    p_rain = clean(p_rain_f, p_rain_fqc, 2) %>% replace_na(0),
    qc_all_precip = apply_thr(p_rain, c(0, 0))
  ) %>%
  select(-p_rain)
  # Attach automatically-extracted QC flags
  #bind_cols(extract_QC(site)) %>%
  # "sa_irga" just means it applies to both sensors, but this is implied
  #rename_at(vars(contains("sa_irga")), str_remove, "sa_")
str(QC) # Check the results

# Check results of QC flagging
lapply(select(QC, -timestamp), table, useNA = "always")

# Visual check on flags (should definitely write this plot into a function)
# - plot_qc_flags <- function(data, var, flags)
# - make sure none of the tests flagged every value, etc.
# fc
site %>%
  select(timestamp, fc) %>%
  mutate(fc = clean(fc, QC$qc_fc_ss), fc = clean(fc, QC$qc_irga_rssi)) %>%
  filter(!is.na(fc), site$phi > 0.7, site$ustar > 0.05) %>%
  left_join(select_at(QC, vars(timestamp, matches(combn_flags$fc)))) %>%
  #rename_at(vars(qc_sa_abslim, qc_sa_spikeshF), str_replace, "sa", "sa_sa") %>%
  coalesce_flags(-timestamp, -fc) %>%
  ggplot(aes(timestamp, fc, color = flag)) +
  geom_point() +
  ylim(-25, 25)
# le
site %>%
  select(timestamp, le) %>%
  filter(!is.na(le)) %>%
  left_join(select_at(QC, vars(timestamp, matches(combn_flags$le)))) %>%
  rename_at(vars(qc_sa_abslim, qc_sa_spikeshF), str_replace, "sa", "sa_sa") %>%
  coalesce_flags(-timestamp, -le) %>%
  ggplot(aes(timestamp, le, color = flag)) +
  geom_point()
# h
site %>%
  select(timestamp, h) %>%
  filter(!is.na(h)) %>%
  left_join(select_at(QC, vars(timestamp, matches(combn_flags$h)))) %>%
  coalesce_flags(-timestamp, -h) %>%
  ggplot(aes(timestamp, h, color = flag)) +
  geom_point() +
  ylim(-200, 300)

# Combine standard tests/filters into new flags

# Add key of combined flags to reference list
# - ThIS IS IMPORTANT: used to combine flags AND for reference later
combn_flags <- append(combn_flags, list(
  qc_tau_std = str_subset(names(QC), combn_flags$tau),
  qc_h_std = str_subset(names(QC), combn_flags$h),
  qc_le_std = str_subset(names(QC), combn_flags$le),
  qc_fc_std = str_subset(names(QC), combn_flags$fc)
))
combn_flags # check the list

# Create the combined standard flags
QC <- QC %>%
  mutate(
    qc_tau_std = combn_QC(., combn_flags$qc_tau_std, "qc_tau_std"),
    qc_h_std = combn_QC(., combn_flags$qc_h_std, "qc_h_std"),
    qc_le_std = combn_QC(., combn_flags$qc_le_std, "qc_le_std"),
    qc_fc_std = combn_QC(., combn_flags$qc_fc_std, "qc_fc_std")
  ) 

# Look at the overall effect of combined standard flags
lapply(select_at(QC, vars(ends_with("_std"))), table, useNA = "always")
site %>%
  left_join(QC) %>%
  filter(!is.na(fc)) %>%
  mutate(qc_fc_std = factor(qc_fc_std)) %>%
  ggplot(aes(timestamp, fc, color = qc_fc_std)) +
  geom_point() +
  ylim(-50, 50)
  
# Add flags to main data frame
site <- site %>% left_join(select_at(QC, vars(timestamp, ends_with("_std"))))


### Quality checking - other indicators of sensor problems/flux issues =========

# Test out any custom flags outside the standard set
# - e.g. precip lags, probable condensation, skewness/kurtosis, etc.

# Look at how sensors respond after precip events
# - vars to check: fc, co2_mixing_ratio, co2_var, inst_li7500_agc_or_rssi
# - tau, w_unrot, w_var, sonic_temperature (-1, 1), ts_var (0, 4)
site %>%
  mutate(var = clean(inst_li7500_agc_or_rssi, qc_fc_std)) %>%
  spread_lags("var", 24, 48, keep_vars = c("p_rain")) %>%
  filter(p_rain > 0) %>%
  mutate(event_id = dplyr::row_number()) %>%
  gather("n_after", "var", -p_rain, -event_id) %>%
  mutate(n_after = as.numeric(n_after)) %>%
  filter(!is.na(var)) %>%
  ggplot(aes(n_after, var)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1))

# Flux
site %>%
  mutate(fc = clean(fc, qc_fc_std)) %>%
  spread_lags("fc", 24, 48, keep_vars = c("p_rain")) %>%
  filter(p_rain > 0) %>%
  mutate(event_id = dplyr::row_number()) %>%
  gather("n_after", "fc", -p_rain, -event_id) %>%
  mutate(n_after = as.numeric(n_after)) %>%
  filter(!is.na(fc)) %>%
  group_by(n_after) %>%
  summarize(mean = mean(fc), sd = sd(fc))
# - erratic behavior from -1 to 9
# Gas concentration
site %>%
  mutate(co2 = clean(co2_mixing_ratio, qc_fc_std)) %>%
  spread_lags("co2", 24, 48, keep_vars = c("p_rain")) %>%
  filter(p_rain > 0) %>%
  mutate(event_id = dplyr::row_number()) %>%
  gather("n_after", "co2", -p_rain, -event_id) %>%
  mutate(n_after = as.numeric(n_after)) %>%
  filter(!is.na(co2)) %>%
  group_by(n_after) %>%
  summarize(mean = mean(co2), sd = sd(co2))
# - erratic behavior from -1 to 1
# Gas concentration variance
site %>%
  mutate(co2_var = clean(co2_meas_sigma, qc_fc_std)) %>%
  spread_lags("co2_var", 24, 48, keep_vars = c("p_rain")) %>%
  filter(p_rain > 0) %>%
  mutate(event_id = dplyr::row_number()) %>%
  gather("n_after", "co2_var", -p_rain, -event_id) %>%
  mutate(n_after = as.numeric(n_after)) %>%
  filter(!is.na(co2_var)) %>%
  group_by(n_after) %>%
  summarize(mean = mean(co2_var), sd = sd(co2_var))
# - erratic behavior from -2 to 3
# Sensor signal strength
site %>%
  mutate(rssi = clean(co2_signal_strength_7500_mean, qc_fc_std)) %>%
  spread_lags("rssi", 24, 48, keep_vars = c("p_rain")) %>%
  filter(p_rain > 0) %>%
  mutate(event_id = dplyr::row_number()) %>%
  gather("n_after", "rssi", -p_rain, -event_id) %>%
  mutate(n_after = as.numeric(n_after)) %>%
  filter(!is.na(rssi)) %>%
  group_by(n_after) %>%
  summarize(mean = mean(rssi), sd = sd(rssi))
# - erratic behavior from -2 to 6

# Test a possible precip lag filter
plotly::ggplotly(
  site %>%
    filter(month(timestamp) == 5) %>%
    mutate(
      fc = clean(fc, qc_fc_std),
      flag = apply_thr_around(p_rain_f, c(0, 0), 1, 1),
      flag = factor(flag)
    ) %>%
    filter(!is.na(fc)) %>%
    ggplot(aes(timestamp, fc, color = flag)) +
    geom_point()
)

site %>%
  mutate(
    fc = clean(fc, qc_fc_std),
    p_rain = clean(p_rain_f, p_rain_fqc, 2) %>% replace_na(0),
    qc_all_precip_range = apply_thr_around(p_rain_f, c(0, 0), 1, 1),
    # Skewness/kurtosis
    qc_co2_skw = apply_thr(co2_meas_skw, c(-2, 2), flag = "outside"),
    qc_co2_kur = apply_thr(co2_meas_kur, c(1, 8), flag = "outside"),
    qc_sa_w_skw = apply_thr(w_skw, c(-2, 2), flag = "outside"),
    qc_sa_w_kur = apply_thr(w_kur, c(1, 8), flag = "outside"),
    # Condensation/ high relative humidity
    qc_irga_cond = apply_thr(le, c(0, 0), flag = "lower"),
    qc_irga_hum = apply_thr(rh, c(99, 99))
  ) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(timestamp, fc, color = qc_irga_hum)) +
  facet_wrap(~ daytime) +
  geom_point() +
  ylim(-40, 40)

# Add flags to QC data frame, if applicable
QC <- site %>%
  transmute(
    timestamp = timestamp,
    qc_all_precip_range = apply_thr_around(p_rain_f, c(0, 0), 1, 1)
  ) %>%
  right_join(QC)

# Add key of combined flags to reference list
combn_flags <- append(combn_flags, list(
  qc_tau_prelim = c("qc_tau_std", "qc_all_precip_range"),
  qc_h_prelim = c("qc_h_std", "qc_all_precip_range"),
  qc_le_prelim = c("qc_le_std", "qc_all_precip_range"),
  qc_fc_prelim = c("qc_fc_std", "qc_all_precip_range")
))
combn_flags # check the list

# Create the combined preliminary flags
QC <- QC %>%
  mutate(
    qc_tau_prelim = combn_QC(., combn_flags$qc_tau_prelim, "qc_tau_prelim"),
    qc_h_prelim = combn_QC(., combn_flags$qc_h_prelim, "qc_h_prelim"),
    qc_le_prelim = combn_QC(., combn_flags$qc_le_prelim, "qc_le_prelim"),
    qc_fc_prelim = combn_QC(., combn_flags$qc_fc_prelim, "qc_fc_prelim")
  ) 

# Look at the overall effect of combined preliminary flags
lapply(select_at(QC, vars(ends_with("_prelim"))), table, useNA = "always")
site %>%
  left_join(QC) %>%
  filter(!is.na(fc)) %>%
  mutate(qc_fc_prelim = factor(qc_fc_prelim)) %>%
  ggplot(aes(timestamp, fc, color = qc_fc_prelim)) +
  geom_point() +
  ylim(-50, 50)

# Add flags to main data frame
site <- site %>% left_join(select_at(QC, vars(timestamp, ends_with("_prelim"))))

### Quality checking - flux interdependency and composite flags ================

# h_prelim is needed in interdep() only if irga = "open", otherwise not used
QC <- QC %>%
  mutate(
    timestamp = timestamp,
    qc_h_interdep = pull(
      interdep(qc_le_prelim, qc_h_prelim, "open"), qc_h_interdep
    ),
    qc_le_interdep = pull(
      interdep(qc_le_prelim, qc_h_prelim, "open"), qc_le_interdep
    ),
    qc_fc_interdep = pull(
      interdep(qc_le_prelim, qc_h_prelim, "open"), qc_fc_interdep
    )
  )

# Visual check on flags
# fc
site %>%
  left_join(select(QC, timestamp, qc_fc_interdep)) %>%
  mutate(
    fc = clean(fc, qc_fc_prelim),
    qc_fc_interdep = factor(qc_fc_interdep)
  ) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(timestamp, fc, color = qc_fc_interdep)) +
  facet_wrap(~ daytime) +
  geom_point() +
  ylim(-40, 40)
# le
site %>%
  left_join(select(QC, timestamp, qc_le_interdep)) %>%
  mutate(
    le = clean(le, qc_le_prelim),
    qc_le_interdep = factor(qc_le_interdep)
  ) %>%
  filter(!is.na(le)) %>%
  ggplot(aes(timestamp, le, color = qc_le_interdep)) +
  facet_wrap(~ daytime) +
  geom_point()
# h
site %>%
  left_join(select(QC, timestamp, qc_h_interdep)) %>%
  mutate(
    h = clean(h, qc_h_interdep),
    qc_h_interdep = factor(qc_h_interdep)
  ) %>%
  filter(!is.na(h)) %>%
  ggplot(aes(timestamp, h, color = qc_h_interdep)) +
  facet_wrap(~ daytime) +
  geom_point() +
  ylim(-200, 300)

# Combine interdep and preliminary flags

# Add key of combined flags to reference list
combn_flags <- append(combn_flags, list(
  qc_tau_comp = c("qc_tau_prelim"),
  qc_h_comp = c("qc_h_prelim", "qc_h_interdep"),
  qc_le_comp = c("qc_le_prelim", "qc_le_interdep"),
  qc_fc_comp = c("qc_fc_prelim", "qc_fc_interdep")
))
combn_flags # check the list

# Create the combined standard flags
QC <- QC %>%
  mutate(
    qc_tau_comp = combn_QC(., combn_flags$qc_tau_comp, "qc_tau_comp"),
    qc_h_comp = combn_QC(., combn_flags$qc_h_comp, "qc_h_comp"),
    qc_le_comp = combn_QC(., combn_flags$qc_le_comp, "qc_le_comp"),
    qc_fc_comp = combn_QC(., combn_flags$qc_fc_comp, "qc_fc_comp")
  )

# Look at the overall effect of combined preliminary flags
lapply(select_at(QC, vars(ends_with("_comp"))), table, useNA = "always")
site %>%
  left_join(QC) %>%
  filter(!is.na(fc)) %>%
  mutate(qc_fc_comp = factor(qc_fc_comp)) %>%
  ggplot(aes(timestamp, fc, color = qc_fc_comp)) +
  geom_point() +
  ylim(-40, 40)
# Look at cleaned data (at this point)
site %>%
  left_join(QC) %>%
  mutate(fc = clean(fc, qc_fc_comp)) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(timestamp, fc)) +
  geom_point() +
  ylim(-40, 40)

# Add flags to main data frame
site <- site %>% left_join(select_at(QC, vars(timestamp, ends_with("_comp"))))


### Storage correction of fluxes ===============================================

# Check storage flux diurnal cycles
# fc
site %>%
  mutate(
    date = date(timestamp),
    hour = decimal_hour(timestamp)
  ) %>%
  ggplot(aes(hour, sc_single)) +
  geom_line(aes(group = date), alpha = 0.1) +
  stat_summary(fun.y = mean, geom = "line", color = "cornflowerblue") +
  ylim(-2, 2)

# Add storage flux to the respective original flux, if applicable
# - correction  overwrites the original values of respective variables
# - original values are remapped to the names "[flux]_orig"
site <- site %>%
  mutate(
    h_orig = h, le_orig = le, fc_orig = fc,
    h = add_st(h, "h", sh_single),
    le = add_st(le, "le", sle_single),
    fc = add_st(fc, "fc", sc_single)
  )


### Low frequency flux despiking ===============================================

# Recode GR 
site <- mutate(site, GR = sw_in_f)

# Determine theoretically-acceptable ranges
# - taken as range of variable cleaned by composite flag allowing only QC == 0 
site %>% 
  mutate(fc = clean(fc, qc_fc_comp, 1)) %>% 
  ggplot(aes(timestamp, fc)) +
  geom_point()
plaus_lims <- site %>% 
  mutate(
    h = clean(h, qc_h_comp, 1),
    le = clean(le, qc_le_comp, 1),
    fc = clean(fc, qc_fc_comp, 1)
  ) %>% 
  summarize(
    h_ext = list(range(h, na.rm = TRUE)),
    le_ext = list(range(le, na.rm = TRUE)), 
    fc_ext = list(range(fc, na.rm = TRUE)),
    h_qnt = list(quantile(h, c(0.001, 0.999), na.rm = TRUE)),
    le_qnt = list(range(le, c(0.001, 0.999), na.rm = TRUE)), 
    fc_qnt = list(range(fc, c(0.001, 0.999), na.rm = TRUE))
  ) %>% 
  as.list() %>% purrr::flatten()
# Add general limits (from literature) & "eyeballed" limits to the list
plaus_lims <- append(plaus_lims, list(
  h_lit = c(-100, 500),
  le_lit = c(-100, 600), 
  fc_lit = c(-40, 40),
  h_man = c(-100, 300),
  le_man = c(-100, 500), 
  fc_man = c(-25, 30)
))

# Run spike detection algorithm
# (Not applied to tau)
QC <- site %>%
  # Remove missing part of the year
  filter(!is.na(GR)) %>%
  # Not using plausible limits for fc - too uncertain given number of spikes
  # Lowering z to 5.5 from 7 (see Papale2006)
  transmute(
    timestamp = timestamp,
    qc_h_spikesLF = despikeLF(
      ., var = "h", qc_flag = "qc_h_comp", name_out = "qc_h_spikesLF", 
      var_thr = plaus_lims$h_lit, light = "GR", night_thr = 20, z = 5.5
    ),
    qc_le_spikesLF = despikeLF(
      ., var = "le", qc_flag = "qc_le_comp", name_out = "qc_le_spikesLF", 
      var_thr = plaus_lims$le_lit, light = "GR", night_thr = 20, z = 5.5
    ),
    qc_fc_spikesLF = despikeLF(
      ., var = "fc", qc_flag = "qc_fc_comp", name_out = "qc_fc_spikesLF", 
      var_thr = NULL, light = "GR", night_thr = 20, z = 5.5
    )
  ) %>%
  right_join(QC)

# Check the results
lapply(select_at(QC, vars(ends_with("_spikesLF"))), table, useNA = "always")
site %>%
  left_join(QC) %>%
  filter(!is.na(fc)) %>%
  mutate(
    fc = clean(fc, qc_fc_comp), 
    qc_fc_spikesLF = factor(qc_fc_spikesLF)
  ) %>%
  ggplot(aes(timestamp, fc, color = qc_fc_spikesLF)) +
  geom_point()


### Fetch filter ===============================================================

# Two methods: (1) simple fetch filter, (2) footprint-weighted area

# Read in the boundary vector data
boundary <- read.csv(boundary_input)

# Test out both approaches, evaluate most appropriate
site %>% 
  mutate(
    fc = clean(fc, qc_fc_comp),
    flag = fetch_filter(., "x_70perc_h00", "wind_dir", boundary$dist),
    flag = apply_thr(phi, c(0.7, 0.7), flag = "lower"),
    flag = factor(flag)
  ) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(timestamp, fc, color = flag)) +
  geom_point()
# Phi appears to provide a more reliable filter

# NB: qc_all_fetch70 is actually not applied to tau
QC <- site %>% 
  transmute(
    timestamp = timestamp,
    qc_all_fetch70 = apply_thr(phi, c(0.7, 0.7), "qc_all_fetch70", "lower")
  ) %>%
  right_join(QC)

# Check the results
lapply(select_at(QC, vars(ends_with("_fetch70"))), table, useNA = "always")


### Insufficient turbulence ====================================================

# Look at overall turbulence relationship
site %>%
  left_join(QC) %>%
  mutate(
    fc = clean(fc, qc_fc_comp),
    fc = clean(fc, qc_all_fetch70),
    fc = clean(fc, qc_fc_spikesLF)
  ) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(ustar, fc)) +
  facet_wrap(~ daytime) +
  geom_point()

# Look at turbulence relationships across temperature classes
site %>%
  left_join(QC) %>%
  mutate(
    fc = clean(fc, qc_fc_comp),
    fc = clean(fc, qc_all_fetch70),
    fc = clean(fc, qc_fc_spikesLF)
  ) %>%
  filter(daytime == 0, !is.na(fc)) %>%
  mutate(TA_class = cut_number(TA_f, n = 6)) %>%
  ggplot(aes(ustar, fc)) +
  facet_wrap(~ TA_class, scales = "free") +
  geom_point() +
  stat_summary_bin(fun.y = mean, geom = "point", bins = 20, color = "orangered")

# Look at turbulence relationships across temperature classes & seasons
site %>%
  left_join(QC) %>%
  mutate(
    season = season(timestamp),
    fc = clean(fc, qc_fc_comp),
    fc = clean(fc, qc_all_fetch70),
    fc = clean(fc, qc_fc_spikesLF)
  ) %>%
  filter(daytime == 0, !is.na(fc)) %>%
  group_by(season) %>%
  mutate(TA_class = cut_number(TA_f, n = 6)) %>%
  ggplot(aes(ustar, fc)) +
  facet_wrap(season ~ TA_class, scales = "free") +
  geom_point() +
  stat_summary_bin(fun.y = mean, geom = "point", bins = 20, color = "orangered")

# Visualize potential turbulence flag
site %>%
  left_join(QC) %>%
  mutate(
    fc = clean(fc, qc_fc_comp),
    fc = clean(fc, qc_all_fetch70),
    fc = clean(fc, qc_fc_spikesLF),
    #zeta_flag = apply_thr(zeta, c(-0.1, 0.1), "outside"),
    flag = apply_thr(ustar, c(0.05, 0.75), flag = "outside"),
    flag = factor(flag)
  ) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(timestamp, fc, color = flag)) +
  geom_point()

QC <- site %>%
  transmute(
    timestamp = timestamp,
    qc_fc_ustar = apply_thr(ustar, c(0.05, 0.75), "qc_fc_ustar", "outside")
  ) %>%
  right_join(QC)

### QC Summary =================================================================

# NOTE: this is messy, but it does what I need it to do
# - should clean it up at some point though 
# - need to consider how best to bundle different pieces & retain flexibility

# Create final combined flags
# Add key of combined flags to reference list
combn_flags <- append(combn_flags, list(
  qc_tau_final = c("qc_tau_prelim"),
  qc_h_final = c("qc_h_comp", "qc_h_spikesLF", "qc_all_fetch70"),
  qc_le_final = c("qc_le_comp", "qc_le_spikesLF", "qc_all_fetch70"),
  qc_fc_final = c(
    "qc_fc_comp", "qc_fc_spikesLF", "qc_all_fetch70", "qc_fc_ustar"
  )
))
combn_flags # check the list

# Create the combined final flags
QC <- QC %>%
  mutate(
    qc_tau_final = combn_QC(., combn_flags$qc_tau_final, "qc_tau_final"),
    qc_h_final = combn_QC(., combn_flags$qc_h_final, "qc_h_final"),
    qc_le_final = combn_QC(., combn_flags$qc_le_final, "qc_le_final"),
    qc_fc_final = combn_QC(., combn_flags$qc_fc_final, "qc_fc_final")
  )

# Preview final flags
# h
site %>%
  left_join(QC) %>%
  mutate(
    clean = clean(h, qc_h_final),
    flag = factor(qc_h_final)
  ) %>%
  select(timestamp, orig = h, clean, flag) %>%
  gather("type", "h", -timestamp, -flag) %>%
  filter(!is.na(h)) %>%
  ggplot(aes(timestamp, h, color = flag)) +
  facet_wrap(~type, scales = "free") +
  geom_point(alpha = 0.5)
# le
site %>%
  left_join(QC) %>%
  mutate(
    clean = clean(le, qc_le_final),
    flag = factor(qc_le_final)
  ) %>%
  select(timestamp, orig = le, clean, flag) %>%
  gather("type", "le", -timestamp, -flag) %>%
  filter(!is.na(le)) %>%
  ggplot(aes(timestamp, le, color = flag)) +
  facet_wrap(~type, scales = "free") +
  geom_point(alpha = 0.5)
# fc
site %>%
  left_join(QC) %>%
  mutate(
    clean = clean(fc, qc_fc_final),
    flag = factor(qc_fc_final)
  ) %>%
  select(timestamp, orig = fc, clean, flag) %>%
  gather("type", "fc", -timestamp, -flag) %>%
  filter(!is.na(fc)) %>%
  ggplot(aes(timestamp, fc, color = flag)) +
  facet_wrap(~type, scales = "free") +
  geom_point(alpha = 0.5)

# Summarize results of QC flagging
QC %>%
  left_join(select(site, timestamp, tau, h, le, fc, sw_in_f)) %>%
  # Recode interdep flags to account for additivity
  mutate(
    qc_h_interdep = if_else(qc_h_comp == 2 & qc_h_prelim != 2, 2L, 0L),
    qc_le_interdep = if_else(qc_le_comp == 2 & qc_le_prelim != 2, 2L, 0L),
    qc_fc_interdep = if_else(qc_fc_comp == 2 & qc_fc_prelim != 2, 2L, 0L)
  ) %>%
  # Get vector of all flagged observations
  gather("flag", "value", -timestamp, -tau, -h, -le, -fc, -sw_in_f) %>%
  # Remove unflagged observations
  mutate(value = if_else(value > 1, 1, 0)) %>%
  #filter(value == 1) %>%
  group_by(timestamp) %>%
  # Add the number of flags applied for each flux at each timestamp
  mutate(
    tau_flags = sum(
      flag %in% expand_list(combn_flags, "qc_tau_final") & value == 1, 
      na.rm = TRUE
    ),
    h_flags = sum(
      flag %in% expand_list(combn_flags, "qc_h_final") & value == 1, 
      na.rm = TRUE
    ),
    le_flags = sum(
      flag %in% expand_list(combn_flags, "qc_le_final") & value == 1, 
      na.rm = TRUE
    ),
    fc_flags = sum(
      flag %in% expand_list(combn_flags, "qc_fc_final") & value == 1,
      na.rm = TRUE
    )
  ) %>%
  group_by(flag) %>%
  summarize(
    # Number of observations flagged by each test
    n = sum(value, na.rm = TRUE),
    n_tau = sum(value[!is.na(tau)], na.rm = TRUE),
    n_h = sum(value[!is.na(h)], na.rm = TRUE),
    n_le = sum(value[!is.na(le)], na.rm = TRUE),
    n_fc = sum(value[!is.na(fc)], na.rm = TRUE),
    n_fc_day = sum(value[!is.na(fc) & sw_in_f > 20], na.rm = TRUE),
    n_fc_night = sum(value[!is.na(fc) & sw_in_f < 20], na.rm = TRUE),
    # Percent of observations flagged by each test
    p_tau = n_tau / sum(!is.na(tau)),
    p_h = n_h / sum(!is.na(h)),
    p_le = n_le / sum(!is.na(le)),
    p_fc = n_fc / sum(!is.na(fc)),
    p_fc_day = n_fc_day / sum(!is.na(fc[sw_in_f > 20])),
    p_fc_night = n_fc_night / sum(!is.na(fc[sw_in_f < 20])),
    # Number of observations *uniquely* flagged by each test
    n_distinct_tau = sum(tau_flags == 1 & !is.na(tau)),
    n_distinct_h = sum(h_flags == 1 & !is.na(h)),
    n_distinct_le = sum(le_flags == 1 & !is.na(le)),
    n_distinct_fc = sum(fc_flags == 1 & !is.na(fc))
  ) %>%
  # Round proportions to two decimal places
  mutate_if(is.numeric, round, 2) %>%
  # Remove values for flags that do not apply to flux
  mutate_at(
    vars(contains("tau")), ~ifelse(str_detect(flag, combn_flags$tau), ., NA)
  ) %>%
  mutate_at(
    vars(contains("h")), ~ifelse(str_detect(flag, combn_flags$h), ., NA)
  ) %>%
  mutate_at(
    vars(contains("le")), ~ifelse(str_detect(flag, combn_flags$le), ., NA)
  ) %>%
  mutate_at(
    vars(contains("fc")), ~ifelse(str_detect(flag, combn_flags$fc), ., NA)
  )


# Next get number of observations uniquely flagged by each test
# - how to iterate this across other flux vars?
flag_stats <- site %>%
  # Remove missing flux values
  filter(!is.na(fc)) %>%
  select(timestamp) %>%
  # Attach QC flags for each observation
  left_join(select_at(QC, vars(timestamp, combn_flags$qc_fc_std))) %>%
  # Get vector of all flagged observations
  gather("flag", "value", -timestamp) %>%
  filter(value > 1) %>%
  select(-value) %>%
  # Remove observations that are flagged by multiple tests
  group_by(timestamp) %>%
  filter(n() == 1) %>%
  # Group by the flag name, summarize number of uniquely flagged observations
  group_by(flag) %>%
  summarize(unique = n_distinct(timestamp))


### QC for Ch4 flux data =======================================================

# Get a subsetted data frame with all recognizable QC information
site_qc <- site %>%
  select_at(vars(
    ends_with("_ssitc_test"), 
    ends_with("_ss_test"), 
    ends_with("_vm97_test"), # 800011199 
    # - spikes, ampres, dropout, abslim, skw/kur, discont, timelag, angle, steady
    ends_with("_itc"), ends_with("_ss"), # Foken stats used to calculate flags
    ends_with("_nsr"), # Mahrt 1998 Nonstationarity Ratios (NR = 2 crit. val.)
    ends_with("_corrdiff"), # Correlation differences with/without repeat values
    ends_with("_zcd"), # Zero-Counts on Differenced variables
    ends_with("_kid"), # Kurtosis Index on Differenced variables
    ends_with("_hf")
  ))

# ch4_tlag_actual, ch4_tlag_used, ch4_tlag_nominal, ch4_tlag_min, ch4_tlag_max
# ch4_meas_skw, ch4_meas_kur, w_skw, w_kur

qplot(timestamp, fch4, data = site)
qplot(timestamp, fch4 / 1000, color = factor(fch4_ssitc_test), data = site)
qplot(timestamp, fch4 / 1000, color = fch4_ss < 30, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = u_itc < 30, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = w_itc < 30, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = fch4_nsr < 2, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = fch4_corrdiff < 0.05, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = ch4_zcd, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = w_zcd, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = ch4_kid, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = w_kid, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = fch4_scf > 3, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = ch4_tlag_used, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = w_unrot > 0.35, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = ustar > 0.05, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = zl > 0.02, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = abs(ch4_meas_skw) > 2, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = ch4_meas_kur < 8, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = abs(w_skw) < 2, data = site, ylim = c(-0.5, 1))
qplot(timestamp, fch4 / 1000, color = w_kur < 8, data = site, ylim = c(-0.5, 1))

# QC scheme
site %>%
  filter(!is.na(fch4)) %>%
  separate(
    ch4_vm97_test, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9"), 
    1:8, remove = FALSE
  ) %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000, flag = factor(v9)) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  ggplot(aes(timestamp, fch4, color = flag)) +
  geom_point() +
  ylim(-0.1, 0.7)

# Get absolute limits
site %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  filter(!is.na(fch4)) %>%
  summarize(list(quantile(fch4, c(0.001, 0.999))))

# Check results
site %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000) %>%
  mutate(flag = between(fch4, -0.0212781, 0.5574030)) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  filter(!is.na(fch4)) %>%
  ggplot(aes(timestamp, fch4, color = flag)) +
  geom_point()
site %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  filter(!is.na(fch4), between(fch4, -0.0212781, 0.5574030)) %>%
  ggplot(aes(timestamp, fch4)) +
  geom_point()
site %>%
  mutate(date = date(timestamp)) %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  filter(!is.na(fch4), between(fch4, -0.0212781, 0.5574030)) %>%
  group_by(date) %>%
  summarize(fch4 = mean(fch4)) %>%
  ggplot(aes(date, fch4)) +
  #stat_summary(fun.data = mean_se, geom = "pointrange") +
  geom_point()
site %>%
  mutate(hour = decimal_hour(timestamp)) %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  filter(!is.na(fch4), between(fch4, -0.0212781, 0.5574030)) %>%
  ggplot(aes(hour, fch4)) +
  stat_summary(fun.data = mean_se, geom = "pointrange")
site %>%
  mutate(hour = decimal_hour(timestamp), month = month(timestamp)) %>%
  mutate(fch4 = (fch4 + sch4_single) / 1000) %>%
  filter(fch4_ss < 30, phi >= 0.7, fch4_nsr < 2, ustar > 0.05) %>%
  filter(!is.na(fch4), between(fch4, -0.0212781, 0.5574030)) %>%
  ggplot(aes(hour, fch4)) +
  facet_wrap(~month, scales = "free") +
  stat_summary(fun.data = mean_se, geom = "pointrange")

# Extract & combine the standardized tests/filters
QC_fch4 <- site %>%
  transmute(
    timestamp = timestamp,
    qc_fch4 = qc_fch4,
    # Low signal strength
    #qc_ch4_rssi = apply_thr(rssi_77_mean, c(15, 15), "qc_ch4_rssi", "lower"),
    # Mixing ratio outside reasonable limits
    qc_ch4_ppm = apply_thr(ch4_mixing_ratio, c(1.74, 1.74), flag = "lower"),
    # Excessive spectral correction factor
    qc_fch4_scf = apply_thr(fch4_scf, c(2, 3))
  ) %>%
  left_join(
    select_at(QC, vars(timestamp, qc_all_precip, qc_all_wresid, matches("_sa")))
  ) %>%
  mutate(
    qc_fch4_std = combn_QC(
      ., str_subset(names(.), "timestamp", negate = TRUE), "qc_fch4_std"
    )
  )

# Visual check on flags
# - make sure none of the tests flagged every value, etc.
site %>%
  select(timestamp, fch4) %>%
  filter(!is.na(fch4)) %>%
  left_join(QC_fch4) %>%
  rename_at(vars(qc_sa_abslim, qc_sa_spikeshF), str_replace, "sa", "sa_sa") %>%
  coalesce_flags(-timestamp, -fch4) %>%
  ggplot(aes(timestamp, fch4, color = flag)) +
  geom_point() +
  ylim(-1, 1.5)

# Fluxes with RSSI < 10 are already removed by EddyPro
# - check if there should be a higher minimum signal strength
# Flux
site %>% 
  left_join(QC_fch4) %>%
  mutate(
    fch4 = clean(fch4, qc_fch4_std)
    #rssi_77_mean = diff(zoo::zoo(rssi_77_mean), na.pad = TRUE)
  ) %>%
  filter(!is.na(fch4)) %>% 
  ggplot(aes(rssi_77_mean, fch4)) + 
  geom_point() + 
  #xlim(-10, 20) +
  ylim(-0.25, 0.75)
# Mixing ratio
site %>% 
  left_join(QC_fch4) %>%
  mutate(ch4 = clean(ch4_mixing_ratio, qc_fch4_std)) %>%
  filter(!is.na(ch4)) %>% 
  ggplot(aes(rssi_77_mean, ch4)) + 
  geom_point() + 
  ylim(2, 5)
# Not much evidence for a relationship with signal strength


# Other flags?
QC_fch4 <- site %>%
  transmute(
    # Probable condensation on ch4 sensor
    qc_ch4_cond = apply_thr(fh2o, c(0, 0), "qc_ch4_cond", "lower"),
    # Very low turbulence/stable conditions
    qc_fch4_ustar = apply_thr(ustar, c(0.05, 0.05), "qc_fch4_ustar", "lower"),
    # Footprint too far outside site area
    qc_fch4_fetch = apply_thr(phi, c(0.6, 0.6), "qc_fch4_fetch", "lower"),
    qc_ch4_var = apply_thr(sqrt(ch4_var), c(0.01, 0.01), "qc_ch4_var")
  )

# Look at effect of QC flags, adjust if needed
site %>%
  bind_cols(select(QC_fch4, -qc_fch4)) %>%
  mutate(
    fch4 = clean(fch4, qc_fch4),
    fch4 = clean(fch4, qc_ch4_rssi),
    fch4 = clean(fch4, qc_all_precip),
    fch4 = clean(fch4, qc_fch4_scf),
    fch4 = clean(fch4, qc_fch4_ustar),
    fch4 = clean(fch4, qc_ch4_ppm),
    fch4 = clean(fch4, qc_ch4_cond),
    fch4 = clean(fch4, qc_all_wresid),
    fch4 = clean(fch4, qc_sa_runs),
    #fch4 = clean(fch4, qc_fch4_fetch),
    fch4 = clean(fch4, qc_ch4_var),
    flag = factor(qc_fch4_fetch)
    #flag = factor(is.na(p_rain))
  ) %>%
  ggplot(aes(timestamp, fch4, color = flag)) +
  geom_point() +
  ylim(-0.5, 1)

# Preliminary flag
QC$qc_fch4_prelim <- combn_QC(QC_fch4, names(QC_fch4), "qc_fch4_prelim")

# Calculate reasonable flux limits based on quantiles
fch4_prelim <- clean(site$fch4, QC$qc_fch4_prelim)
quantile(fch4_prelim, c(0.0025, 0.9975), na.rm = TRUE)
QC$qc_fch4_spikesLF <- apply_thr(
  fch4_prelim, c(-0.1589996, 0.5131389), "qc_fch4_spikesLF", "outside"
)

# Final flag
QC$qc_fch4_comp <- combn_QC(
  QC, c("qc_fch4_prelim", "qc_fch4_spikesLF"), "qc_fch4_comp"
)
site$qc_fch4_comp <- QC$qc_fch4_comp

# Save dataset
site_out <- paste0(
  wd, "/08_flux_qc/output/", "flux_qc_", tag_out, ".csv"
)
write_eddy(site, site_out)

# remove random uncertainty > 0.2? (see Runkle2019)

# Runs in ch4 concentration
#QC$qc_ch4_runs <- flag_runs(site$ch4_mixing_ratio, "qc_ch4_runs")
# Very low turbulence/stable conditions
#QC$qc_fch4_stab <- apply_thr(abs(site$zeta), c(0.1, 0.1), "qc_fch4_stab")
#QC$qc_fch4_stab <- dplyr::if_else(site$PAR >= 20, 0L, QC$qc_fch4_stab)
#QC$qc_fch4_turb <- apply_thr(site$ustar, c(0.1, 0.1), "qc_fch4_turb", "lower")
# SD of mixing ratio too high
#QC$qc_ch4_var <- apply_thr(sqrt(site$ch4_var), c(0.005, 0.005), "qc_ch4_var")
# Another flag based on ts_var??
# -> QC$qc_ch4_ts <- apply_thr(site$ts_var, c(0.5, 0.5), "qc_ch4_ts")




