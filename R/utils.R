is_call_name <- function(x, name){
  is.call(x) && call_name(x) %in% name
}

quietly_squash <- function(x) quietly(squash)(x)$result

find_leaf <- function(arg, branch_name = NULL) {
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
  } else {
    fabletools::traverse(
      arg, .h=function(x) if(!is_call(x)) pick_leaf(x),
      base = function(x) if(is_call(x)) {pick_call(x[[1]]); return(FALSE)} else TRUE
    )
    return(
      tree$leaves[!sapply(tree$leaves, deparse) %in% sapply(tree$calls, deparse)] %>%
        unique()
    )


  }

}





