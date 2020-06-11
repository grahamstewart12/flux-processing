
decimal_hour <- function(x) {
  lubridate::hour(x) + lubridate::minute(x) / 60
}

add_time_comps <- function(data, timestamp = timestamp, offset = 900) {
  
  timestamp <- rlang::enquo(timestamp)
  
  id_along <- function(x) {
    rle <- rle(as.vector(x))
    rep(seq_along(rle$lengths), times = rle$lengths)
  }
  
  data %>%
    dplyr::mutate(
      year = lubridate::year(!!timestamp - offset),
      month = lubridate::month(!!timestamp - offset),
      week = lubridate::week(!!timestamp - offset),
      date = lubridate::date(!!timestamp - offset),
      day = id_along(date),
      hour = decimal_hour(!!timestamp)
    )
}

remove_time_comps <- function(data, ...) {
  
  dots <- rlang::exprs(...)
  
  time_comps <- dplyr::select(
    data, timestamp, year, month, week, date, day, hour
  )
  
  if (!rlang::is_empty(dots)) time_comps <- dplyr::select(time_comps, !!!dots)
  
  time_names <- names(time_comps)
  
  dplyr::select(data, -dplyr::any_of(time_names))
}

# Generate Year-long Half-hourly Time Step Vector
create_timesteps <- function(year, dts = 48, tz = "UTC", shift_by = 720 / dts) {
  
  year <- as.numeric(year)
  
  if (!dts %in% c(24, 48)) {
    stop("Only implemented for 24 or 48 daily time steps.", call. = FALSE)
  }
  format <- "%Y-%m-%d-%H-%M"
  
  start <- paste(year, 1, 1, 0, shift_by, sep = "-")
  end <- paste(year + 1, 1, 1, 0, 30 - shift_by, sep = "-")
  
  # Timestamp vector with half-hourly timestamps
  out <- seq(
    strptime(start, format, tz), strptime(end, format, tz), (24 / dts * 60 * 60)
  )
  out
}

complete_time_seq <- function(data, time, time_diff, .first, .last) {
  
  # Not using yet--risky if timestamp doesn't start/end on start/end of year
  
  # Determine time var if not given
  if (missing(time)) {
    
    time_name <- data %>%
      purrr::map(is.POSIXct) %>%
      purrr::keep(~ .) %>%
      names() %>%
      vctrs::vec_slice(1)
    
    time <- rlang::sym(time_name)
  } else {
    
    time_name <- rlang::ensym(time) %>% rlang::as_string()
    time <- rlang::enquo(time)
  }
  
  time_seq <- dplyr::pull(data, !!time)
  if (missing(.first)) .first <- dplyr::first(time_seq)
  if (missing(.last)) .last <- dplyr::last(time_seq)
  
  # Determine interval if not given
  # TODO this should be its own function 'get_time_diff'
  if (missing(time_diff)) {
    time_diff <- time_seq %>%
      lubridate::int_diff() %>%
      lubridate::as.duration() %>%
      vctrs::vec_count() %>%
      vctrs::vec_slice(1) %>%
      purrr::pluck("key", 1)
  } else {
    time_diff <- lubridate::duration(time_diff)
  }
  
  # Create sequence vector
  full_seq <- tibble::tibble(
    !!time_name := seq(.first, .last, by = time_diff)
  )
  
  dplyr::left_join(full_seq, data, by = time_name)
}