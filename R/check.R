check_gamma_order <- function(gamma_ind_pair){

  error_in_group <- gamma_ind_pair %>%
    sapply(function(ind_group){
      pairs <- matrix(ind_group, nrow = 2) %>%
        asplit(2)
      if(any(sapply(pairs, function(z) z != sort(z)))){
        abort("Please specify intercept first in a polynomial threshold, and put smaller id for intercept than slope.")
      }
      if(!identical(unlist(pairs), sort(unlist(pairs)))){
        TRUE
      } else {
        FALSE
      }
    })
  if(!identical(sort(temp <- unique(unlist(gamma_ind_pair))), temp) || any(error_in_group)){
    abort("Gamma id in parametric threshold indicates order of optimisation,
and within the same polynomial, the smaller one is the intercept.
For example, if gamma 1 and 2 are intercept and slope of one polynomial and gamma 3 and 4 are that of another,
gamma 1 and 2 will be estimated first.
You cannot put gamma 1 and 3 in one polynomial and 2 and 4 in another.")
  }
  NULL
}
