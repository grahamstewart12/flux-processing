
create_tag <- function(site, year, date) {
  paste0(site, stringr::str_sub(year, -2), "_", date)
}

between2 <- function(x, range) {
  # Like dplyr::between, but accepts a range vector instead of two values
  
  # Checks
  len <- length(range)
  if (len > 2) stop("Range must be a vector of length 2.", call. = FALSE)
  if (range[2] < range[1]) range <- rev(range)
  
  dplyr::between(x, range[1], range[2])
}

create_around_mapper <- function(.x, .p, .n = 1, .and = FALSE) {
  
  sep <- if (.and) "&" else "|"
  logic <- rlang::f_text(.p)
  
  add <- ""
  for (i in 1:.n) {
    lag <- stringr::str_replace(logic, ".x", glue::glue("lag(.x, {i})"))
    lead <- stringr::str_replace(logic, ".x", glue::glue("lead(.x, {i})"))
    add <- glue::glue("{add} {sep} {lag} {sep} {lead}")
    
  }
  
  glue::glue("~ {logic} {add}") %>% as.formula() %>% purrr::as_mapper()
}

around <- function(.x, .p, .n = 1, .and = FALSE) {
  
  # e.g. around(x, ~ is.na(.x), .n = 1) is equivalent to:
  #   is.na(x) & is.na(lag(x, 1)) & is.na(lead(x, 1))
  around_mapper <- create_around_mapper(.x, .p, .n, .and)
  around_mapper(.x)
}

create_adjacent_mapper <- function(.x, .p, .n = 1, .and = FALSE, .fun) {
  
  rhs <- rlang::f_text(.p)
  seq <- 1:.n
  
  sep <- if (.and) " & " else " | "
  
  logic <- seq %>%
    purrr::map(
      ~ stringr::str_replace(rhs, ".x", glue::glue("{.fun}(.x, {.})"))
    ) %>%
    glue::glue_collapse(sep = sep)
  
  glue::glue("~ {logic}") %>% as.formula() %>% purrr::as_mapper()
}

before <- function(.x, .p, .n = 1, .and = FALSE) {
  
  mapper <- create_adjacent_mapper(.x, .p, .n, .and, .fun = "lag")
  mapper(.x)
}

after <- function(.x, .p, .n = 1, .and = FALSE) {
  
  mapper <- create_adjacent_mapper(.x, .p, .n, .and, .fun = "lead")
  mapper(.x)
}