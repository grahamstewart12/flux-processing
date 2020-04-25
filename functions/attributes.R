
add_attr <- function(x, attr, values) {
  
  values <- tidyr::replace_na(values, "-")
  
  x %>% 
    colnames() %>% 
    rlang::set_names() %>%
    purrr::map2_df(values, ~ (x[[.x]] <-`attr<-`(x[[.x]], attr, .y)))
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

drop_all_attributes <- function(x, attrs = NULL) {
  
  purrr::modify(x, drop_attributes, attrs = attrs)
}