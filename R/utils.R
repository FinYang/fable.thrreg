is_call_name <- function(x, name = NULL){
  is.call(x) && !is.null(name) && call_name(x) %in% name
}

has_call_name <- function(x, name = NULL){
  find_leaf(x, exclude = name) %>%
    sapply(is_call_name, name) %>%
    any()
}

Arithmetic <- c("+", "-", "*", "/", "^", "%%", "%/%")

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
    traverse(
      arg, .h=function(x) if(!is_call_name(x, branch_name)) pick_leaf(x),
      base = function(x) !is_call_name(x, branch_name)
    )
    return(
      tree$leaves[!sapply(tree$leaves, deparse) %in% branch_name]
    )
  } else if(!is.null(include)){
    include <- setdiff(include, exclude)
    traverse(
      arg, .h=function(x) if(!is_call(x) || (!call_name(x) %in% include)) pick_leaf(x),
      base = function(x) if(is.call(x) && (call_name(x) %in% include)) {pick_call(x[[1]]); return(FALSE)} else TRUE
    )
    return(
      tree$leaves[!sapply(tree$leaves, deparse) %in% sapply(tree$calls, deparse)] %>%
        unique()
    )
  } else {
    traverse(
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
gsub_multi <- function(pattern = NULL, replacements, x, strings = NULL, ...){
  if((is.null(pattern) && is.null(strings)) || (!is.null(pattern) && !is.null(strings)))
    abort("Either pattern or strings should be given")
  if(!is.null(pattern)){
    if(length(pattern) == 1) {
      if(stringr::str_count(x, pattern) != length(replacements))
        abort("The number of matched patterns is not the same as the number of provided replacements.")
      out <- purrr::reduce(replacements, function(x, r) sub(pattern, r, x, ...), .init = x)
    } else {
      abort("Only one pattern is supported.")
    }
  } else {
    if(length(strings) != length(replacements))
      abort("The number of given strings is not the same as the number of provided replacements.")
    out <- purrr::reduce2(replacements, strings, function(x, r, s) sub(s, r, x, fixed = TRUE, ...),
                          .init = x)
  }
  out

}


str_extract_perl <- function(string, pattern) {
  regmatches(string, gregexpr(pattern, string, perl = TRUE))
}


#' @param x numeric vector
#' @param value centre
#' @param num number of closest to find
#' @param exclude logical vector. Exclude these observations when finding closest
#' @return logical vector
closest <- function(x, value, num, exclude = logical(length(x))) {
  output <- logical(length(x))
  stopifnot(length(value) == 1 || length(x) == length(value))
  dist <- x-value

  dist[is.na(dist)] <- Inf
  dist[exclude] <- Inf
  output[dist == 0] <- TRUE
  output[dist<0] <- rank(abs(dist[dist<0]), ties.method = "min") <=num
  output[dist>0] <- rank(abs(dist[dist>0]), ties.method = "min") <=num
  output

}
