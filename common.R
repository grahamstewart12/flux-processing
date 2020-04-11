
# X create_timestamps
# X create_tag
# X add_attr
# X latest version
# X decimal_hour
# X clean
# X combine_flags
# X flag_mahalanobis
# X coalesce_flags
# X add_time_comps
# X drop_attributes
# around
# X flag_around
# replace_real
# X fill_linear
# X detect_fog

# Find the latest version of an output file
latest_version <- function(dir, name, ext = ".csv", n = 1) {
  
  # Add path separator to dir (if needed)
  if (stringr::str_sub(dir, -1) != "/") dir <- stringr::str_c(dir, "/")
  # Add string separator to name (if needed)
  if (stringr::str_sub(name, -1) != "_") name <- stringr::str_c(name, "_")
  
  # Get files matching name and extension
  files <- list.files(dir) %>%
    stringr::str_subset(paste0("\\b", name, "[:upper:]")) %>%
    stringr::str_subset(ext)
  
  # Find file with the latest date tag
  ext_len <- stringr::str_length(ext)
  dates <- files %>%
    stringr::str_sub(-10 - ext_len, -1 - ext_len) %>%
    lubridate::ymd()
  date <- dates[which(dates == dplyr::last(sort(dates)))]
  out <- stringr::str_subset(files, as.character(date))
  
  # Choose longest file name if multiple files still remain
  # TODO allow n > 1
  if (length(out) > 1) out <- dplyr::first(out, stringr::str_length(out))
  paste0(dir, out)
  
}

create_tag <- function(site, year, date) {
  
  paste0(site, stringr::str_sub(year, -2), "_", date)
}

add_attr <- function(x, attr, values) {

  values <- tidyr::replace_na(values, "-")
  names <- colnames(x)
  l <- purrr::map2(names, values, ~ (x[[.x]] <-`attr<-`(x[[.x]], attr, .y)))
  names(l) <- names
  as.data.frame(l)
}

drop_attributes <- function(x, attrs = NULL) {
  if (is.null(attrs)) {
    attributes(x) <- NULL
    x
  } else {
    for (i in 1:length(attrs)) {
      attr(x, attrs[i]) <- NULL
    }
    x
  }
}

# Generate Year-long Half-hourly Time Step Vector
create_timesteps <- function(year, dts = 48, tz = "GMT", shift_by = 720 / dts) {
  
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

decimal_hour <- function(x) {
  lubridate::hour(x) + lubridate::minute(x) / 60
}

add_time_comps <- function(data, timestamp = timestamp, offset = 900) {
  
  timestamp <- rlang::enquo(timestamp)
  
  data %>%
    dplyr::mutate(
      year = lubridate::year(!!timestamp - offset),
      month = lubridate::month(!!timestamp - offset),
      week = lubridate::week(!!timestamp - offset),
      date = lubridate::date(!!timestamp - offset),
      day = tidy_rle(date)$id,
      hour = decimal_hour(!!timestamp)
    )
}

clean <- function(x, flag, value = 2, replace_with = NA_real_, na.as = 0) {
  
  if (is.character(value)) {
    flag <- dplyr::if_else(flag == value, 2, 0)
    value <- 2
  }
  
  flag <- tidyr::replace_na(flag, na.as)
  dplyr::if_else(flag >= value, replace_with, x)
}

combine_flags <- function(..., clean_value = 0) {
  
  dots <- rlang::list2(...)
  
  # Method for single data frame (assumes that all columns are flags)
  if (length(dots) == 1 & inherits(dots[[1]], "data.frame")) {
    dots <- as.list(dots[[1]])
  }
  
  dots <- dots %>%
    purrr::map(as.integer) %>%
    purrr::map(~ dplyr::if_else(. <= clean_value, NA_integer_, .))
  
  
  combined <- dplyr::coalesce(!!! dots)
  out <- tidyr::replace_na(combined, 0L)
  
  out
}

coalesce_flags <- function(data, ..., prefix = "qc_[[:alnum:]]+\\_") {
  
  data %>%
    tibble::as_tibble() %>%
    tidyr::gather("flag", "value", ...) %>%
    dplyr::mutate(
      flag = stringr::str_remove(flag, prefix),
      value = dplyr::if_else(value > 1, flag, NA_character_)
    ) %>%
    tidyr::spread(flag, value) %>%
    dplyr::mutate(
      flag = dplyr::coalesce(!!! dplyr::select(., ...)),
      flag = tidyr::replace_na(flag, "none")
    )
}

flag_around <- function(x, width = 1, flag_value = 2L, na.as = 0,
                        rule = c("between", "adjacent")) {
  # Useful when spike occurs
  rule <- match.arg(rule)
  x <- as.integer(x)
  flag_value <- as.integer(flag_value)
  na.as <- as.integer(na.as)
  out <- tidyr::replace_na(x, na.as)
  
  if (rule == "between") {
    for (i in 1:width) {
      out <- dplyr::if_else(
        dplyr::lag(out, i) == flag_value & dplyr::lead(out, i) == flag_value,
        flag_value, x
      )
    }
  }
  
  if (rule == "adjacent") {
    for (i in 1:width) {
      out <- dplyr::if_else(
        dplyr::lag(out, i) == flag_value | dplyr::lead(out, i) == flag_value,
        flag_value, x
      )
    }
  }
  
  out
}

flag_mahalanobis <- function(x, y, alpha = 0.001) {
  
  tbl <- tibble::tibble(
    x = dplyr::if_else(is.na(y), NA_real_, x),
    y = dplyr::if_else(is.na(x), NA_real_, y)
  )
  
  if (isTRUE(all.equal(tbl$x, tbl$y))) return(rep(0L, nrow(tbl)))
  
  # Get cutoff value (z) from alpha and degrees of freedom
  df <- ncol(tbl)
  z <- qchisq(p = 1 - alpha, df = df)
  
  means <- tbl %>%
    dplyr::summarize_all(mean, na.rm = TRUE) %>%
    purrr::flatten_dbl()
  
  cov <- cov(tbl, use = "pairwise.complete.obs")
  
  m_dist <- mahalanobis(tbl, means, cov)
  
  flag <- dplyr::if_else(m_dist > z, 2L, 0L) %>% tidyr::replace_na(0L)
  flag
}

fill_linear <- function(x, n_max = NULL) {
  
  # Option for returning original vector
  if (!is.null(n_max)) {
    if (is.na(n_max) | n_max == 0) return(x)
  }
  
  y <- seq_along(x)
  out <- approx(y, x, y)$y
  
  if (!is.null(n_max)) {
    .x <- tidyr::replace_na(x, 9999)
    rl <- rle(c(.x))
    n <- rep(rl$lengths, times = rl$lengths)
    out <- dplyr::if_else(n > n_max & is.na(x), NA_real_, out)
  }
  
  out
}

detect_fog <- function(rh, lw_in, lw_out, sw_in, p_rain) {
  
  # High humidity
  rh_ind <- rh > 90
  
  # Similar incident & outgoing long-wave radiation
  lw_ind <- (lw_out - lw_in) < 20
  
  # Rapid increases in shortwave radiation after day break
  sw_ind <- (sw_in - dplyr::lag(sw_in)) > 100
  
  # No precipitation
  p_ind <- p_rain == 0
  
  #inds <- unique(c(rh_ind, lw_ind, sw_ind, p_ind))
  
  out <- rep(0L, length(rh))
  out[which(rh_ind & lw_ind & sw_ind & p_ind)] <- 1L
  out
}
