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
  one = function(arg){
    arg <- enexpr(arg)
    data <- break_interaction(arg, "*") %>%
      quietly_squash()
    bare_xreg <- !sapply(data, function(x) is_call_name(x, "ind"))
    data[[which(bare_xreg)]] <- expr(xreg(!!data[[which(bare_xreg)]]))

    data <- data %>%
      `names<-`(sapply(., deparse)) %>%
      map(eval_tidy, data = self$data, env = self$specials)
    list(data = data,
         expression = arg)
  },
  ind = function(x){
    x <- enexpr(x)
    out <- find_leaf(x, exclude = "gamma")
    bare_xreg <- !sapply(out, is_call_name, names(specials_thrreg))
    out$xreg <- expr(xreg(!!!out[which(bare_xreg)]))
    out[which(bare_xreg)] <- NULL
    names(out)[names(out)==""] <- sapply(out[names(out)==""], call_name)
    out <- map(out, eval_tidy, data = self$data, env = self$specials)
    out$expression <- x
    out$ind_expression <- find_leaf(x, include = c("<", "<=", ">", ">=")) %>%
      {.[!sapply(., function(y) is_call_name(y, "gamma"))]}

    out



  },
  gamma = function(id = NULL,
                   var = NULL,
                   grid = NULL,
                   kernel = epaker,
                   min_points = NCOL(self$data) + 5,
                   bw = sd(var)/10,
                   max_iter = 100){

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
    regs <- enexprs(...)

    if(any(sapply(regs[sapply(regs, is_call)], call_name) == "gamma"))
      stop("gamma() in the formula should used inside ind()")

    parsed_regs <- map(regs, break_interaction)
    bare_xreg <- sapply(parsed_regs, length)==1

    if(all(bare_xreg)) {
      xregs_exog <- reduce(unlist(parsed_regs), function(.x, .y) call2("+", .x, .y))
      return(
        list(xregs = select(as_tibble(self$data), !!!find_leaf(xregs_exog)),
             expression = xregs_exog)
      )
    }

    call_xreg <- c(parsed_regs[!bare_xreg],
                   list(list(unlist(parsed_regs[bare_xreg])))
    )
    xreg<- map(call_xreg, function(xcall) {
      bare_xreg <- !sapply(xcall, function(x) all(sapply(x, function(y) is_call_name(y,"ind"))))
      xcall[["xreg"]] <- expr(xreg(!!!xcall[[which(bare_xreg)]]))
      xcall[[which(bare_xreg)]] <- NULL
      xcall <- unlist(xcall)
      map(xcall, eval_tidy, data = self$data, env = self$specials)
    })
    return(xreg)
  },
  .required_specials = c("delta"),
  .xreg_specials = c("gamma", "ind")
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

break_interaction <- function(x, no_base = c("+", "*", "(")){
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
      !is_call(.x) || !call_name(.x) %in% no_base
      # !is_call(.x) || !call_name(.x) %in% c("+")
    }
  )
  flatten2(output)
}

inside_ind <- function(ind){
  inner <- ind[[2]]
  recur_list <- function(x){
    # browser()
    if(is.list(x)){
      return(lapply(x, recur_list))
    }
    if(is_call_name(x, c( "&", "<", "<=", ">", ">="))) {
      return(lapply(list(x[[1]], as.list(x[-1])), recur_list))
    }
    x
  }
  inner_ls <- recur_list(inner)

  if(depth(inner_ls)>4 || (depth(iner_ls)>2 && !identical(inner_ls[[1]], sym("&"))))
    stop("`ind` only support at most one `&` operation.")
  inner_ls[[2]] %>%
    lapply(function(x) `attr<-`(x[[2]], "compare", x[[1]]))
}
