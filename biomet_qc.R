### ============================================================================
# Biomet data quality control ==================================================
### ============================================================================

# Purpose: 

# References:

# Estévez, J., Gavilán, P., Giráldez, J.V., 2011. Guidelines on validation 
# procedures for meteorological data from automatic weather stations. Journal of 
# Hydrology 402, 144–154. https://doi.org/10.1016/j.jhydrol.2011.02.031

# Olson, R.J., Holladay, S.K., Cook, R.B., Falge, E., Baldocchi, D., Gu, L., 
# 2004. FLUXNET: Database of fluxes, site characteristics, and flux-community 
# information. Oak Ridge National Laboratory, Oak Ridge, TN. 
# https://doi.org/10.2172/1184413

# ONEFlux Processing Pipeline
# https://github.com/AmeriFlux/ONEFlux/blob/master/oneflux_steps/qc_auto/src

# Input(s):

# Output(s):

# Load the required packages
devtools::load_all("~/Desktop/RESEARCH/fluxtools", quiet = TRUE)
library(openeddy)
library(lubridate)
library(tidyverse)

source("~/Desktop/DATA/Flux/tools/engine/biomet_qc_v2.R")
source("~/Desktop/DATA/Flux/tools/reference/var_attributes.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")

# Initialize script settings & documentation
script_init(settings, site_metadata)


### Biomet data initialization =================================================

# Load the Biomet file
biomet <- read.csv(biomet_input, stringsAsFactors = FALSE)

# Add timestamp components
biomet <- biomet %>% 
  mutate(timestamp = ymd_hms(timestamp, tz = md$tz_name)) %>%
  add_time_comps()

# Load Biomet file from closest site
biomet_aux <- read.csv(aux_input, stringsAsFactors = FALSE)
biomet_aux <- mutate(
  biomet_aux, 
  # Parse timestamp
  timestamp = ymd_hms(timestamp, tz = md$tz_name),
  # Remove records during system failure
  p_rain = clean(p_rain, flag_biomet_system(biomet_aux))
)

# Initialize QC data frame
qc_biomet <- select(biomet, timestamp)

# Warn if both ta or rh vars are the same
#biomet %>% filter(ta == ta_ep) %>% summarize(n())
#biomet %>% filter(rh == rh_ep & rh != 100) %>% summarize(n())


### Automatic statistical flags ================================================

# Compute all auto flags
auto_flags <- vars %>% 
  map(~ biomet_auto_flags(biomet, !!., p_rain)) %>%
  set_names(map_chr(vars, rlang::as_string))

# Validate flags for analagous sensors
# - can do this for sw_in vs sw_out because same unit NOT same sensor
auto_flags <- auto_flags %>%
  list_modify(sw_in2 = pluck(auto_flags, "sw_in")) %>%
  validate_flags(ta, ta_ep) %>%
  validate_flags(ppfd_in, sw_in) %>%
  validate_flags(sw_out, sw_in2) %>%
  list_modify(sw_in2 = NULL)

# Give flags QC names
for (i in seq_along(vars)) {
  auto_flags[[i]] <- rename_all(
    auto_flags[[i]], ~ str_c("qc_", names(auto_flags)[i], "_", .)
  )
}

# Plot all flags for each var
auto_plots <- map2(vars, auto_flags, ~ plot_flags(biomet, !!.x, .y))
#map(auto_plots, ~ .)

# Add combined flag to QC dataset
for (i in seq_along(vars)) {
  var_qc_name <- str_c("qc_", rlang::as_string(vars[[i]]), "_auto")
  qc_biomet[, var_qc_name] <- combine_flags(auto_flags[[i]], clean_value = 1)
}


### Flag theoretically implausible values ======================================

# Add p_rain to vars list
all_vars <- append(vars, list(expr(p_rain)))

for (i in seq_along(all_vars)) {
  var_qc_name <- str_c("qc_", rlang::as_string(all_vars[[i]]), "_plaus")
  qc_biomet[, var_qc_name] <- drop_attributes(apply_thr(
    pull(biomet, !!all_vars[[i]]),
    pluck(var_attrs, rlang::as_string(all_vars[[i]]), "limits"),
    flag = "outside"
  ))
}


### Special potential incoming radiation flags =================================

qc_biomet <- mutate(
  qc_biomet,
  # Check SW_IN against potential radiation
  qc_sw_in_pot = flag_rad(biomet$sw_in, biomet$sw_in_pot),
  # Check converted PPFD_IN against potential radiation
  qc_ppfd_in_pot = flag_rad(biomet$sw_in_ppfd, biomet$sw_in_pot)
)


### Special precipitation flags ================================================

p_rain_rh_lim <- biomet %>% 
  filter(p_rain > 0, !is.na(rh)) %>% 
  summarize(mean(rh) - 2 * sd(rh)) %>%
  pluck(1, 1) %>% signif(1)
p_rain_kt_lim <- biomet %>% 
  filter(p_rain > 0, !is.na(kt)) %>% 
  summarize(mean(kt) + 2 * sd(kt)) %>%
  pluck(1, 1) %>% signif(1)

qc_biomet <- mutate(
  qc_biomet,
  # Unlikely precipitation - soft flags
  # - non-rain events are not flagged because no plausible explanation for that
  # - Estevez et al. 2011 thresholds: rh = 80, kt = 0.5
  # - allowing empirical thresholds here due to site differences
  qc_p_rain_rh = if_else(biomet$p_rain > 0 & biomet$rh < p_rain_rh_lim, 1L, 0L),
  qc_p_rain_kt = if_else(biomet$p_rain > 0 & biomet$kt > p_rain_kt_lim, 1L, 0L),
  qc_p_rain_lik = combine_flags(qc_p_rain_rh, qc_p_rain_kt),
  # Did nearby site experience similar conditions? (+/- 1 lag)  
  qc_p_rain_aux = if_else(magrittr::and(
    biomet$p_rain > 0, around(biomet_aux$p_rain, ~ .x == 0, .and = TRUE)
  ), 1L, 0L),
  # Hard flag if unlikely precipitation AND closest site had no rain
  qc_p_rain_err = if_else(qc_p_rain_lik + qc_p_rain_aux > 1, 2L, 0L) %>%
    replace_na(0L)
)


### Difference between analagous measurement flags =============================

# TA, RH, SW_IN
qc_biomet <- bind_cols(
  qc_biomet,
  flag_var_diffs(biomet, ta, ta_ep, qc_biomet),
  # Only consider daytime differences for incoming rad
  flag_var_diffs(
    mutate_at(biomet, vars(sw_in, ppfd_in), ~ if_else(night_pot, NA_real_, .)), 
    sw_in, ppfd_in, qc_biomet
  )
)


### Combine flags and check results ============================================

for (i in seq_along(all_vars)) {
  var_name <- rlang::as_string(all_vars[[i]])
  var_flags <- qc_biomet %>% 
    select_at(vars(contains(str_c("_", var_name, "_")))) %>%
    mutate_all(~ if_else(. == 1, 0L, .))
  if (!str_detect(var_name, "_ep")) {
    var_flags <- select_at(var_flags, vars(-contains("_ep_")))
  }
  
  # Combine flags, flag isolated points
  biomet[, str_c("qc_", var_name)] <- flag_around(combine_flags(var_flags), 1)
}

# Plot of all flags for each var
flag_plots <- all_vars %>%
  map(~ plot_flags(biomet, !!., qc_biomet, geom = "line")) %>%
  set_names(str_c(map_chr(all_vars, rlang::as_label), "_flags"))
#map(flag_plots, ~ .)

# How many flagged?
biomet %>%
  select_at(vars(all_of(str_c("qc_", map_chr(all_vars, rlang::as_string))))) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  filter(value == 2) %>%
  summarize(n_flag = n())

# Plot overall flag for each var
qc_plots <- all_vars %>%
  map(~ plot_flags(biomet, !!., geom = "line")) %>%
  set_names(str_c(map_chr(all_vars, rlang::as_label), "_qc"))
#map(qc_plots, ~ .)

# Combine all plot lists
plots <- append(flag_plots, qc_plots)


### Save output ================================================================

# Save processed Biomet dataset
biomet_out <- file.path(path_out, str_c("biomet_qc_", tag_out, ".csv"))
write.csv(biomet, biomet_out, row.names = FALSE)

# Create documentation for processed Biomet output
biomet_docu <- append(
  settings, 
  list(
    files = c(biomet_input, aux_input),
    aux_site = md$closest_site[1],
    plaus_lims = modify(var_attrs, "limits"),
    p_rain_lik = list(kt_lim = p_rain_kt_lim, rh_lim = p_rain_rh_lim)
  )
)
biomet_docu_out <- str_replace(biomet_out, ".csv", ".txt")
# Save documentation
sink(biomet_docu_out)
biomet_docu
sink()

# Save Biomet QC flags datasets
auto_biomet_out <- str_replace(biomet_out, "_qc_", "_qc_auto_flags_")
write.csv(bind_cols(auto_flags), auto_biomet_out, row.names = FALSE)
qc_biomet_out <- str_replace(biomet_out, "_qc_", "_qc_flags_")
write.csv(qc_biomet, qc_biomet_out, row.names = FALSE)

# Save one pdf document with all diagnostic/summary plots
plot_path <- file.path(path_out, paste0("biomet_qc_plots_", tag_out, ".pdf"))
pdf(plot_path)
map(plots, ~ .)
dev.off()

# Finished
