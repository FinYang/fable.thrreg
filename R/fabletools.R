
# This file contains (modified) internal function of fabletools
names_no_null <- function(x){
  names(x) %||% rlang::rep_along(x, "")
}

merge_named_list <- function(...){
  flat <- rlang::flatten(list(...))
  nm <- names_no_null(flat)
  lapply(split(flat, nm), function(x) rlang::flatten(unname(x)))
}
