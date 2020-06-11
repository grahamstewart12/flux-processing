
# Depends on functions/flag.R
source("~/Desktop/DATA/Flux/tools/engine/functions/flag.R")

clean <- function(x, ..., value = 2, replace_with = NA_real_, na.as = 0) {
  
  dots <- rlang::list2(...)
  
  # If multiple flags were passed, combine them
  if (length(dots) > 1) {
    flag <- combine_flags(!!!dots)
  } else {
    flag <- purrr::pluck(dots, 1)
  }
  
  if (is.character(value)) {
    flag <- dplyr::if_else(flag == value, 2, 0)
    value <- 2
  }
  
  flag <- tidyr::replace_na(flag, na.as)
  dplyr::if_else(flag >= value, replace_with, x)
}
