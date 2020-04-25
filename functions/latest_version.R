
# Find the latest version of an output file
latest_version <- function(dir, name, ext = ".csv", n = 1) {
  
  # Get files matching extension
  files <- dir %>% list.files() %>% stringr::str_subset(ext)
  
  # If name is missing, assume all files in dir have the same name
  if (!missing(name)) {
    
    # Add string separator to name (if needed)
    if (stringr::str_sub(name, -1) != "_") name <- stringr::str_c(name, "_")
    
    # Get files matching name
    files <- stringr::str_subset(files, paste0("\\b", name, "[:upper:]"))
  }
  
  # Find file(s) with the latest date tag
  ext_len <- stringr::str_length(ext)
  
  dates <- files %>%
    stringr::str_sub(-10 - ext_len, -1 - ext_len) %>%
    lubridate::ymd() %>%
    vctrs::vec_sort(direction = "desc") %>%
    vctrs::vec_chop(list(1:n)) %>%
    purrr::pluck(1)
  
  out <- stringr::str_subset(files, stringr::str_c(dates, collapse = "|"))

  # Return path to latest file(s)
  file.path(dir, out)
}
