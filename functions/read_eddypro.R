
read_eddypro <- function(file, timestamp = paste(date, time), ...) {
  
  timestamp <- rlang::enexpr(timestamp)
  
  # Default: skip = 0, units = TRUE
  data <- openeddy::read_eddy(file, ...)
  
  data %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(timestamp = lubridate::ymd_hm(!!timestamp)) %>%
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