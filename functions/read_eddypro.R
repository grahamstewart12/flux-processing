
read_eddypro <- function(file, timestamp = paste(date, time), ...) {
  
  timestamp <- rlang::enexpr(timestamp)
  
  # Default: skip = 0, units = TRUE
  data <- openeddy::read_eddy(file, ..., as.is = TRUE)
  
  orders <- c("YmdHM", "YmdHMS", "mdyHM", "mdyHMS")
  
  data %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    # This drops attributes
    #dplyr::mutate(dplyr::across(is.factor, as.character)) %>%
    dplyr::mutate(
      timestamp = lubridate::parse_date_time(!!timestamp, orders = orders)
    ) %>%
    # Only reason timestamp doesn't parse is if data is corrupt - remove
    dplyr::filter(!is.na(timestamp))
}

read_eddypro_settings <- function(file) {
  
  settings <- readr::read_lines(file)
  
  settings %>%
    # Separate names & values
    purrr::map(stringr::str_split, "=") %>%
    purrr::flatten() %>%
    # Remove list headers
    purrr::discard(~ length(.x) < 2) %>%
    # Set name as first value in each list item
    rlang::set_names(purrr::modify(., 1)) %>%
    # Remove name from list element
    purrr::modify(2) %>%
    purrr::map(readr::parse_guess)
}

read_ghg <- function(file, ext = c("data", "metadata"), biomet = FALSE, ...) {
  
  ext <- rlang::arg_match(ext)
  ext <- stringr::str_c(".", ext)
  
  if (biomet) {
    ext <- stringr::str_c("-biomet", ext)
    nl <- 5
    col_spec <- readr::cols_only(
      DATE = readr::col_date(),
      TIME = readr::col_character(),
      `LOGGERPOWER_1_1_1(other)` = readr::col_double(),
      `LOGGERTEMP_1_1_1(C)` = readr::col_double(),
      `LOGGERVIN_1_1_1(V)` = readr::col_double(),
      `LWIN_1_1_1(W/m^2)` = readr::col_double(),
      `LWOUT_1_1_1(W/m^2)` = readr::col_double(),
      `PPFD_1_1_1(umol/m^2/s^1)` = readr::col_double(),
      `RH_1_1_1(%)` = readr::col_double(),
      `RN_1_1_1(W/m^2)` = readr::col_double(),
      `SHF1_1_1_1(W/m^2)` = readr::col_double(),
      `SHF2_2_1_1(W/m^2)` = readr::col_double(),
      `SHF3_3_1_1(W/m^2)` = readr::col_double(),
      `SHFSENS1_1_1_1(other)` = readr::col_double(),
      `SHFSENS2_2_1_1(other)` = readr::col_double(),
      `SHFSENS3_3_1_1(other)` = readr::col_double(),
      `SWC1_1_1_1(m^3/m^3)` = readr::col_double(),
      `SWC2_2_1_1(m^3/m^3)` = readr::col_double(),
      `SWC3_3_1_1(m^3/m^3)` = readr::col_double(),
      `SWIN_1_1_1(W/m^2)` = readr::col_double(),
      `SWOUT_1_1_1(W/m^2)` = readr::col_double(),
      `TA_1_1_1(C)` = readr::col_double(),
      `TS1_1_1_1(C)` = readr::col_double(),
      `TS2_2_1_1(C)` = readr::col_double(),
      `TS3_3_1_1(C)` = readr::col_double(),
      `P_RAIN_1_1_1(mm)` = readr::col_double(),
      `VERSION_1_1_1(other)` = readr::col_double()
    )
  } else {
    nl <- 7
    col_spec <- readr::cols(
      DATAH = readr::col_skip(),
      Date = readr::col_date(),
      Time = readr::col_character(),
      CHK = readr::col_character(),
      .default = readr::col_double()
    )
  }
  
  name <- basename(file)
  file_con <- unz(file, stringr::str_replace(name, ".ghg", ext))
  #on.exit(closeAllConnections())
  
  out <- readr::read_tsv(
    file_con,
    col_types = col_spec,
    na = c("", "NA", "-9999"),
    skip = nl,
    progress = FALSE,
    ...
  )
  
  readr::stop_for_problems(out)
  
  out
}