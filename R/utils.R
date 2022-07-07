is_call_name <- function(x, name = NULL){
  is.call(x) && !is.null(name) && call_name(x) %in% name
}

quietly_squash <- function(x) quietly(squash)(x)$result

# @param include include only these names as call
# @param exclude exclude as call
find_leaf <- function(arg, branch_name = NULL, exclude = NULL, include = NULL) {
  tree <- new.env(parent=emptyenv())
  tree$i <- 1
  tree$leaves <- list()
  tree$j <- 1
  tree$calls <- list()
  pick_leaf <- function(x){
    tree$leaves[[tree$i]] <- x
    tree$i <- tree$i+1
  }
  pick_call <- function(x){
    tree$calls[[tree$j]] <- x
    tree$j <- tree$j+1
  }

  if(!is.null(branch_name)){
    fabletools::traverse(
      arg, .h=function(x) if(!is_call_name(x, branch_name)) pick_leaf(x),
      base = function(x) !is_call_name(x, branch_name)
    )
    return(
      tree$leaves[!sapply(tree$leaves, deparse) %in% branch_name]
    )
  } else if(!is.null(include)){
    include <- setdiff(include, exclude)
    fabletools::traverse(
      arg, .h=function(x) if(!is_call(x) || (!call_name(x) %in% include)) pick_leaf(x),
      base = function(x) if(is.call(x) && (call_name(x) %in% include)) {pick_call(x[[1]]); return(FALSE)} else TRUE
    )
    return(
      tree$leaves[!sapply(tree$leaves, deparse) %in% sapply(tree$calls, deparse)] %>%
        unique()
    )
  } else {
    fabletools::traverse(
      arg, .h=function(x) if(!is_call(x) || is_call_name(x, exclude)) pick_leaf(x),
      base = function(x) if(is.call(x) && (!call_name(x) %in% exclude)) {pick_call(x[[1]]); return(FALSE)} else TRUE
    )
    return(
      tree$leaves[!sapply(tree$leaves, deparse) %in% sapply(tree$calls, deparse)] %>%
        unique()
    )


  }

}

depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)

#' @examples
#' gsub_multir("a",1:3, "aaabb")
gsub_multir <- function(pattern, replacements, x){
  if(stringr::str_count(x, pattern) != length(replacements))
    abort("The number of matched patterns is not the same as the number of provided replacements.")
  purrr::reduce(replacements, function(x, r) stringi::stri_replace(x, r, regex = pattern), .init = x)

}


str_extract_perl <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string, perl = TRUE))
}
