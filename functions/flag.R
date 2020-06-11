
# Depends on functions/utilities.R
source("~/Desktop/DATA/Flux/tools/engine/functions/utilities.R")

flag_around <- function(x, width = 1, flag_value = 2L, na.as = 0L,
                        rule = c("between", "adjacent")) {
  
  # Useful for eliminating undetected isolated spikes
  
  rule <- rlang::arg_match(rule)
  and <- rule == "between"
  
  # Cast flag & na.as values to prototype of x
  flag_value <- vctrs::vec_cast(flag_value, vctrs::vec_ptype(x))
  na.as <- vctrs::vec_cast(na.as, vctrs::vec_ptype(x))
  
  # Get NA indicies so they can be applied to output vector
  na_loc <- which(is.na(x))
  out <- tidyr::replace_na(x, na.as)
  
  out <- dplyr::if_else(
    around(out, ~ .x == flag_value, .n = width, .and = and), flag_value, x
  )
  
  vctrs::vec_assign(out, na_loc, NA)
}

flag_distance <- function(x, y, alpha = 0.001) {
  
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

flag_repeats <- function(x, length = 2, exclude = NA) {
  rle <- tidy_rle(x)
  
  length <- length[1]
  # Don't flag anything if length is NULL or NA
  if (is.null(length) | is.na(length)) return(rep(0L, length(x)))
  
  if (!is.na(exclude)) {
    rle <- dplyr::mutate(
      rle, lengths = dplyr::if_else(values == exclude, 1L, lengths)
    )
  }
  
  lens <- dplyr::pull(rle, lengths)
  dplyr::if_else(lens >= length, 2L, 0L)
}

flag_thr <- function(x, thr, rule, flag_value = 2L) {
  
  thr <- vctrs::vec_sort(thr)
  rule <- rlang::arg_match(rule, c("higher", "lower", "outside", "between"))
  cln <- vctrs::vec_cast(0, vctrs::vec_ptype(flag_value))
  
  if (rule == "higher") {
    flag <- dplyr::if_else(x > thr, flag_value, cln)
  } else if (rule == "lower") {
    flag <- dplyr::if_else(x < thr, flag_value, cln)
  } else if (rule == "outside") {
    flag <- dplyr::if_else(x < thr[1] | x > thr[2], flag_value, cln)
  } else if (rule == "between") {
    flag <- dplyr::if_else(x > thr[1] & x < thr[2], flag_value, cln)
  }
  
  flag
}

combine_flags <- function(..., clean_value = 0L) {
  
  dots <- rlang::list2(...)
  
  # Method for single data frame (assumes that all columns are flags)
  if (length(dots) == 1 & inherits(dots[[1]], "data.frame")) {
    dots <- as.list(dots[[1]])
  }
  
  out <- purrr::lift_dl(pmax)(dots, na.rm = TRUE)
  
  #dots <- dots %>%
  #  purrr::map(as.integer) %>%
  #  purrr::map(~ dplyr::if_else(. <= clean_value, NA_integer_, .))
  
  
  #combined <- dplyr::coalesce(!!! dots)
  #out <- tidyr::replace_na(combined, 0L)
  
  out
}