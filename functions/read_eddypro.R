
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