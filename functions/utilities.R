
create_tag <- function(site, year, date) {
  
  paste0(site, stringr::str_sub(year, -2), "_", date)
}

tidy_rle <- function(x, explicit_na = FALSE, pull = NULL) {
  
  x <- as.vector(x)
  if (explicit_na) x <- tidyr::replace_na(x, -9999)
  
  rle <- rle(x)
  out <- tibble::tibble(
    id = rep(seq_along(rle$lengths), times = rle$lengths),
    lengths = rep(rle$lengths, times = rle$lengths),
    values = rep(rle$values, times = rle$lengths)
  )
  if (explicit_na) out$values <- dplyr::na_if(out$values, -9999)
  if (!is.null(pull)) out <- dplyr::pull(out, pull)
  
  out
}

between2 <- function(x, range) {
  # Like dplyr::between, but accepts a range vector instead of two values
  
  # Checks
  if (length(range) > 2) {
    stop("Range must be a vector of length 2.", call. = FALSE)
  } 
  if (range[2] < range[1]) range <- rev(range)
  
  dplyr::between(x, range[1], range[2])
}

compose_around_fn <- function(.x, .p, .n = 1, .and = FALSE, .env) {
  
  sep <- if (.and) "&" else "|"
  rhs <- rlang::f_text(.p)
  
  add <- ""
  for (i in 1:.n) {
    lag <- stringr::str_replace(rhs, ".x", glue::glue("lag(.x, {i})"))
    lead <- stringr::str_replace(rhs, ".x", glue::glue("lead(.x, {i})"))
    add <- glue::glue("{add} {sep} {lag} {sep} {lead}")
    
  }
  
  # Add 'around' predicates to original predicate
  glue::glue("~ {rhs} {add}") %>% 
    # Convert to function with specified environment
    as.formula(env = .env) %>% 
    rlang::as_function(env = .env)
}

around <- function(.x, .p, .n = 1, .and = FALSE) {
  
  # e.g. around(x, ~ is.na(.x), .n = 1) is equivalent to:
  #   is.na(x) & is.na(lag(x, 1)) & is.na(lead(x, 1))
  
  # Need to pass caller environment so symbols can be evaluated
  around_fn <- compose_around_fn(.x, .p, .n, .and, .env = rlang::caller_env())
  around_fn(.x)
}

compose_adjacent_fn <- function(.x, .p, .n = 1, .and = FALSE, .fun) {
  
  rhs <- rlang::f_text(.p)
  seq <- 1:.n
  
  sep <- if (.and) " & " else " | "
  
  logic <- seq %>%
    purrr::map(
      ~ stringr::str_replace(rhs, ".x", glue::glue("{.fun}(.x, {.})"))
    ) %>%
    glue::glue_collapse(sep = sep)
  
  glue::glue("~ {logic}") %>% as.formula() %>% rlang::as_function()
}

before <- function(.x, .p, .n = 1, .and = FALSE) {
  
  before_fn <- compose_adjacent_fn(.x, .p, .n, .and, .fun = "lag")
  before_fn(.x)
}

after <- function(.x, .p, .n = 1, .and = FALSE) {
  
  after_fn <- compose_adjacent_fn(.x, .p, .n, .and, .fun = "lead")
  after_fn(.x)
}