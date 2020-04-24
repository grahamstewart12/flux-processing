
clean <- function(x, flag, value = 2, replace_with = NA_real_, na.as = 0) {
  
  if (is.character(value)) {
    flag <- dplyr::if_else(flag == value, 2, 0)
    value <- 2
  }
  
  flag <- tidyr::replace_na(flag, na.as)
  dplyr::if_else(flag >= value, replace_with, x)
}