
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