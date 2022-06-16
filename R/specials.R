#' Epaker kernel function
#'
#' Kernel function for grid windows
#'
#' @param x value to apply kernel to
#' @export
epaker <- function(x) {
  (3/4)*(1-x^2)*(abs(x)<=1)
}

specials_thrreg <- new_specials(
  delta = function(delta = FALSE){
    if(missing(delta))
      delta
  },
  gamma = function(var = NULL, kernel = epaker,
                   # n = 2^8,
                   min_points = NCOL(self$data) + 5, bw = sd(var)/10, max_iter = 10){
    var_name <- rlang::enexpr(var)
    bw <- rlang::enexpr(bw)
    if(!is.null(var_name)){
      var <- select(as_tibble(self$data), !!var_name) %>%
        unlist() %>%
        unname()
    }
    as.list(environment())
  },
  xreg = function(...) {
    return(select(as_tibble(self$data), !!!enexprs(...)))

    regs <- enexprs(...)
    if(any(sapply(regs[sapply(regs, is_call)], call_name) == "gamma"))
      stop("gamma() in the formula should interact with other terms using *.")
    parsed_regs <- map(
      regs,
      function(x){

        output <- traverse(
          x,
          .f = function(.x, ...) {
            .x <- flatten2(.x)
            fabletools:::merge_named_list(.x[[1]], .x[[2]])},
          .g = function(.x){
            map(as.list(get_expr(.x))[-1], expr)
          },
          .h = function(x){ # Base types
            x <- get_expr(x)
            if(!is_call(x)){
              list(list(x))
            }
            else{# Current call is a special function
              set_names(list(list(x)), call_name(x))
            }
          },
          base = function(.x){
            .x <- get_expr(.x)
            !is_call(.x) || !call_name(.x) %in% c("+", "*", "(")
            # !is_call(.x) || !call_name(.x) %in% c("+")
          }
        )
        flatten2(output)
      }
    )
    bare_xreg <- sapply(parsed_regs, length)==1

    if(all(bare_xreg)) {
      xregs_exog <- transmute(as_tibble(self$data), !!!unlist(parsed_regs[bare_xreg]))
      return(xregs_exog)
    }

    call_xreg <- c(parsed_regs[!bare_xreg],
                   list(list(unlist(parsed_regs[bare_xreg])))
    )

    xreg <- map(call_xreg, function(xcall) {
      bare_xreg <- !sapply(xcall, function(x) all(sapply(x, function(y) is_call(y) && call_name(y)=="gamma")))
      xcall[["xreg"]] <- expr(xreg(!!!xcall[[which(bare_xreg)]]))
      xcall[[which(bare_xreg)]] <- NULL
      xcall <- unlist(xcall)
      # xcall
      map(xcall, eval_tidy, env = caller_env(3))
    }
    )

    names(xreg) <- sapply(xreg, function(x) {
      output <- paste0(c(names(x)[names(x)!="xreg"], colnames(x$xreg)), collapse = "_")
      output[!grepl("gamma",output)] <- "xreg"
      output
    }
    )
    return(xreg)
  },
  .required_specials = c("delta"),
  .xreg_specials = "gamma"
)

flatten2 <- function(x) {
  if(is.list(x)){
    if(length(quietly(squash)(x)$result) > 1){
      if(length(x) == 1) {
        x <- flatten2(rlang::flatten(x))
      }
    } else {
      x <- quietly(squash)(x)$result
    }
  }
  x
}
