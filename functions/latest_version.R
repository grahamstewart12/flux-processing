
# Find the latest version of an output file
latest_version <- function(dir, name, ext = ".csv", n = 1) {
  
  # Get files matching extension
  files <- dir %>% list.files() %>% stringr::str_subset(ext)
  
  # If name is missing, assume all files in dir have the same name
  if (!missing(name)) {
    
    name <- as.list(name)
    
    # Add string separator to name (if needed)
    #if (!stringr::str_ends(name[1], "_")) name <- stringr::str_c(name, "_")
    name <- purrr::map(
      name, ~ dplyr::if_else(
        !stringr::str_ends(.x, "_"), stringr::str_c(.x, "_"), .x
      )
    )
    
    # Get files matching name
    #files <- stringr::str_subset(files, paste0("\\b", name, "[:upper:]"))
    files <- purrr::map(
      name, ~ stringr::str_subset(files, paste0("\\b", .x, "[:upper:]"))
    )
  } else {
    files <- list(files)
  }
  
  # Find file(s) with the latest date tag
  ext_len <- stringr::str_length(ext)
  
  #dates <- files %>%
  #  stringr::str_sub(-10 - ext_len, -1 - ext_len) %>%
  #  lubridate::ymd() %>%
  #  vctrs::vec_sort(direction = "desc") %>%
  #  vctrs::vec_slice(1:n)
  dates <- files %>%
    purrr::map(stringr::str_sub, -10 - ext_len, -1 - ext_len) %>%
    purrr::map(lubridate::ymd) %>%
    purrr::map(vctrs::vec_sort, direction = "desc") %>%
    purrr::map(vctrs::vec_slice, 1:n)
  
  #out <- stringr::str_subset(files, stringr::str_c(dates, collapse = "|"))
  out <- files %>%
    purrr::map2(
      dates, ~ stringr::str_subset(.x, stringr::str_c(.y, collapse = "|"))
    ) %>%
    purrr::map(~ file.path(dir, .x))

  # Return path to latest file(s)
  #out <- purrr::map(out, ~ file.path(dir, .x))
  purrr::simplify(out)
}
