
# This file contains (modified) internal function of fabletools


names_no_null <- function(x){
  names(x) %||% rlang::rep_along(x, "")
}

merge_named_list <- function(...){
  flat <- rlang::flatten(list(...))
  nm <- names_no_null(flat)
  lapply(split(flat, nm), function(x) rlang::flatten(unname(x)))
}


#' Recursively traverse an object
#'
#' Internal function in fabletools
#'
#' @param x The object to traverse
#' @param .f A function for combining the recursed components
#' @param .g A function applied to the object before recursion
#' @param .h A function applied to the base case
#' @param base The base case for the recursion
#'
#' @keywords internal
traverse <- function(x, .f = list, .g = identity, .h = identity, base = function(.x) rlang::is_syntactic_literal(.x) || rlang::is_symbol(.x)){
  # base case
  if(base(x))
    return(.h(x))
  # recursive case
  .f(lapply(.g(x), traverse, .f=.f, .g=.g, .h=.h, base=base), .h(x))
}
