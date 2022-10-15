

#' @importFrom dplyr lag
#' @export
dplyr::lag


train_thrreg <- function(.data, specials, ...){
  lag <- dplyr::lag
  if (length(tsibble::measured_vars(.data)) > 1) {
    stop("Only univariate response is supported by thrreg.")
  }
  gamma_fun_env <- new.env(parent=emptyenv())
  gamma_env <- new.env(parent=gamma_fun_env)
  gamma_env$id <- list()
  gamma_env$gamma <- list()
  gamma_env$gamma_special <- list()
  assign("gamma", function(id = var, var = NULL, ...){
    var <- enexpr(var)
    if(is.null(id)) {
      if(length(gamma_env$id)>0) abort("When there are more than one gamma, a numeric id needs to be given.")
      id <- max(unlist(gamma_env$id[vapply(gamma_env$id, is.numeric, logical())]) %||% 0, na.rm = TRUE) + 1
    }
    g <- paste0(".gamma_", id)
    gamma_env$id <- unique(c(gamma_env$id, id))
    gamma_env$gamma[[id]] <- g
    sym(g)
  }, envir = gamma_fun_env)


  eval_gamma_in_string <- function(string){
    gamma_pattern <- "gamma(\\((?:[^()]|(?1))*\\))"
    match_str <-
      str_extract_perl( string, gamma_pattern) %>%
      unlist()
    match_str %>%
      map(str2lang) %>%
      map(eval_tidy, env = gamma_fun_env) %>%
      gsub_multi(strings = match_str, replacements =  ., x = string)
  }

  replace_gamma_lang <- function(ind_expr){
    ind_expr_str <- paste0(deparse(ind_expr), collapse = "")
    str2lang(eval_gamma_in_string(ind_expr_str))

  }
  find_tibble <- function(x, layer){
    x %>%
      purrr::map_depth(layer, function(x) if(is_tibble(x)) x else NULL) %>%
      quietly_squash() %>%
      {.[!sapply(., is.null)]}
  }
  squash_tibble <- function(x, layers){
    lapply(layers, function(layer)
      find_tibble(x, layer)
    ) %>%
      quietly_squash() %>%
      purrr::reduce(function(x, y){
        bind_cols(x, select(y, !any_of(colnames(x)))) })
  }

  # END of functions

  rhs <- specials$xreg[[1]]
  ind_term <- rhs[sapply(rhs, length)>1]

  # gamma_terms <- unique(unlist(gamma_env$gamma))

  rhs_terms <- map(rhs, function(x){
    if("ind" %in% names(x)) {
      ind_expr <- x$ind$expression
      ind_expr <- replace_gamma_lang(ind_expr)
      id <- sapply(x$ind[names(x$ind) == "gamma"], function(y)y$id)
      if(is.null(unlist(id))) id <- gamma_env$id
      gamma_env$gamma_special[as.character(id)] <- x$ind[names(x$ind) == "gamma"]
      if(is_call_name(x$xreg$expression, "+")){

        expr(!!x$xreg$expression : I(as.numeric(!!ind_expr)))
      } else {
        expr(I(!!x$xreg$expression * (!!ind_expr)))
      }
    } else {
      x$xreg$expression
    }
  }
  )


  num_unique_ind <- rhs %>%
    map(function(x) {
      if("ind" %in% names(x)){
        x$xreg$expression
      }
    }) %>%
    .[!sapply(.,is.null)] %>%
    unique() %>%
    length()
  if(num_unique_ind > 1) {
    stop("The variables whose coefficients vary in different regimes need to be identical.")
  }


  if(is.numeric(unlist(gamma_env$id))){
    # reorder (multi_regime)
    temp_idx <- order(sapply(ind_term, function(ind_t) {
      mean(unlist(lapply(ind_t$ind[names(ind_t$ind) == "gamma"], getElement, "id")) %||% NA)
      }))
    ind_term[temp_idx] <- ind_term[order(temp_idx)]
    rhs_terms[sapply(rhs, function(x) "ind" %in% names(x))][temp_idx] <- rhs_terms[sapply(rhs, function(x) "ind" %in% names(x))][order(temp_idx)]
    gamma_env$gamma[temp_idx] <- gamma_env$gamma[order(temp_idx)]
    gamma_env$id[temp_idx] <- gamma_env$id[order(temp_idx)]
    gamma_env$gamma_special[temp_idx] <- gamma_env$gamma_special[order(temp_idx)]
  }

  model_data <- squash_tibble(rhs, 3) %>%
    bind_cols(select(squash_tibble(rhs[sapply(rhs, length)!=1], 4),
                     !any_of(colnames(.))))

  model_df <- .data %>%
    bind_cols(select(model_data, !any_of(colnames(.))))

  # When there is one term
  one_term <- specials$one[[1]]
  if(!is.null(one_term)){

    # if(kernel_est && grepl("gamma", deparse(one_term$expression))) {
    #   abort()
    #
    # }
    one_expr <- one_term$expression %>%
      deparse() %>%
      gsub("ind", "", .) %>%
      eval_gamma_in_string() %>%
      str2lang()
    one_data <- one_term$data
    one_data <- squash_tibble(one_data, 2:3)
    model_df <- model_df %>%
      bind_cols(select(one_data, !any_of(colnames(.))))
  } else {
    one_expr <- NULL
  }

  kernel_est <- FALSE
  if(any(sapply(gamma_env$gamma_special, function(x) !is.null(x$var)))){
    kernel_est <- TRUE
    multi_regime <- FALSE
    if(length(gamma_env$id)>1)
      abort("Only one gamma term can appear in the formula when using kernel estimation.")
  }



  resp <- measured_vars(.data)
  n <- nrow(model_data)
  # k <- ncol(select(rhs[[which(sapply(rhs, length)==1)]]$xreg$xregs, !any_of(resp)))+2
  # k <- length(rhs)
  if(length(tem <- which(sapply(rhs, length)==1)) == 0) {
    k <- length(rhs) + 1
  } else {
    k <- ncol(select(rhs[[tem]]$xreg$xregs, !any_of(resp))) + length(rhs) + 1
  }



  parametric <- FALSE
  # Finding the gamma grids
  if(!kernel_est){

    multi_regime <- FALSE
    if(length(unlist(gamma_env$gamma)) == length(ind_term)){
      get_grid_single_gamma <- function(ind_term, gamma_name, data_eval = NULL, weights_lgl = NULL){
        stopifnot(length(ind_term)==1)
        stopifnot(length(gamma_name)==1)
        stopifnot(is.null(weights_lgl) || is.logical(weights_lgl))
        all_gamma_grids <- map(ind_term, function(x){
          data_eva <- data_eval %||% x$ind$xreg$xregs
          if(is.null(weights_lgl)) weights_lgl <- rep(TRUE, nrow(data_eva))
          eval_tidy(x$ind$ind_expression[[1]], data = data_eva)[weights_lgl]
        }
        ) %>%
          split(unlist(gamma_name)) %>%
          lapply(function(x) x %>%
                   unlist() %>%
                   unique() %>%
                   sort() %>%
                   .[seq(k, length(.)-k, by = 1)] )
        all_gamma_grids <- all_gamma_grids %>%
          expand.grid() %>%
          split(row(.)) %>%
          lapply(unlist)

        all_gamma_grids
      }

      if(length(ind_term)==1) {
        # Single regime
        all_gamma_grids <- get_grid_single_gamma(ind_term, gamma_name = gamma_env$gamma)
      } else if(length(ind_term)>1) {
        num_gamma_each_ind <- 1
        multi_regime <- TRUE
      }

    } else if(length(unlist(gamma_env$gamma)) > length(ind_term)){
      # parametric ploy 2
      parametric <- TRUE

      if(length(unlist(gamma_env$gamma)) > 2*length(ind_term)){
        abort("Indicator function only support second order polynomial at most.")
      }



      gamma_ind_pair <- lapply(ind_term, function(x) sapply(x$ind[names(x$ind) == "gamma"],getElement, "id"))
      num_gamma_each_ind <- min(sapply(gamma_ind_pair, length))

      # Ensures the gamma id in order
      check_gamma_order(gamma_ind_pair)

      for(i in 2:length(gamma_ind_pair)){
        gamma_ind_pair[[i]] <- setdiff(gamma_ind_pair[[i]], gamma_ind_pair[[i-1]])
      }
      pair_idx <- numeric(length(unlist(gamma_env$gamma)))
      for(i in seq_along(gamma_ind_pair)){
        pair_idx[gamma_ind_pair[[i]]] <- i
      }

      given_grids <- gamma_env$gamma_special %>%
        lapply(getElement, "grid")
      if(any(sapply(given_grids, is.null))) abort("You need to provide grid to search for gamma.")
      all_gamma_grids <- given_grids %>%
        `names<-`(paste0(".gamma_", names(.))) %>%
        split(pair_idx) %>%
        lapply(
          do.call, what = expand.grid
        ) %>%
        lapply(function(x) lapply(asplit(x, 1), c))

      # for test
      based_var <- ind_term %>%
        map(function(x) x$ind$expression) %>%
        map(function(x) {
          find_leaf(x, include = c(">=", ">", "<=", "<", "&")) %>%
            {.[sapply(., has_call_name, "gamma")]} %>%
            map(function(x){
              find_leaf(x, include = Arithmetic) %>%
                {.[!sapply(., is_call_name, "gamma")]}
            }
            )
        }) %>%
        unlist() %>%
        unique()
      if(length(based_var)>1){
        abort("Variable to base the threshold should be identical.")
      }
      based_var <- based_var[[1]]
      if(length(ind_term)==1) {
        all_gamma_grids <- all_gamma_grids[[1]]

      } else {
        multi_regime <- TRUE
      }
    }
  }

  fm <- if(!is.null(one_term)){
    expr(I(!!sym(resp) - !!one_expr) ~ !!reduce(rhs_terms, function(.x, .y) call2("+", .x, .y)))
  } else {
    expr(!!sym(resp) ~ !!reduce(rhs_terms, function(.x, .y) call2("+", .x, .y)))
  }
  if(multi_regime){
    # find term with only one gamma
    temp_rhs_idx <- lapply(rhs_terms, find_leaf) %>%
      sapply(function(x) sum(sapply(x, function(y) grepl("\\.gamma", deparse(y))))) %>%
      `<=`(num_gamma_each_ind)
    temp_rhs <- rhs_terms[temp_rhs_idx]
    fm_top <- if(!is.null(one_term)){
      expr(I(!!sym(resp) - !!one_expr) ~ !!reduce(temp_rhs, function(.x, .y) call2("+", .x, .y)))
    } else {
      expr(!!sym(resp) ~ !!reduce(temp_rhs, function(.x, .y) call2("+", .x, .y)))
    }
  }


  get_fm_fun <- function(fm) {
    fm_str <- deparse(fm) %>%
      stringr::str_trim(side = "left") %>%
      paste0(collapse = "")
    purrr::walk(unname(unlist(gamma_env$gamma)), function(.x ){
      fm_str <<- gsub(.x, paste0("!!", .x), fm_str)
    })

    list_from_name <- function(name){
      `names<-`(vector("list", length(gamma_env$gamma)), unname(unlist(gamma_env$gamma)))
    }
    rlang::new_function(
      list_from_name(unname(unlist(gamma_env$gamma))),
      expr(str2lang(paste0(deparse(expr(!!str2lang(fm_str))), collapse = "")))
    )
  }



  get_ssr <- function(gammas, fm, weights = NULL, data = model_df){
    fm <- do.call(get_fm_fun(fm), as.list(gammas))
    weights <- enexpr(weights)
    data <- enexpr(data)

    if(length(gammas) == 2 || kernel_est){
      # test for identification
      test <- function(){
        if(is.null(weights)) weights <- rep(1, nrow(eval(data)))
        temp_env <- new.env()
        purrr::walk2(names(gammas), gammas,assign, envir = temp_env)

        # browser()
        # assign(".gamma_1", -2, temp_env)
        # assign(".gamma_2", 8.2, temp_env)

        temp <- ind_term[[1]]$ind$expression %>%
          replace_gamma_lang() %>%
          eval_tidy( data = eval(data)[eval(weights)>0,], env = temp_env)
        if((sum(temp, na.rm = TRUE) <= k ) || (sum(!temp, na.rm = TRUE) <= k) ) {
          return(TRUE)
        }
        FALSE
      }

      if(test()) return(Inf)
    }

    if(length(gammas) == 2 ){
      # test for intersection
      test2 <- function(){
        if(is.null(weights)) weights <- rep(1, nrow(eval(data)))

        ra <- eval_tidy(based_var, data = eval(data)[eval(weights)>0,]) %>%
          range()
        between(-gammas[[1]]/gammas[[2]], ra[[1]], ra[[2]])
      }

      # if(test2()) return(Inf)
    }


    fit <- eval(expr(lm(!!fm, data = !!data, model = FALSE, weights = !!weights)))
    if(anyNA(coef(fit))){
      # if(kernel_est)
      #   warning("Not all coefficient can be estimated. Consider increase bandwidth bw or raise number of min_points.")
      return(Inf)
    }
    sum(resid(fit)^2)
  }

  if(kernel_est){
    kernel <- gamma_env$gamma_special$fee
    model_df <- bind_cols(model_df,
                          tibble(!!kernel$var_name := kernel$var) %>%
                            select(!any_of(colnames(model_df))))
    kernel_grid <- kernel$var
    # kernel_grid <- do.call(seq, c(as.list(range(kernel$var)), list(length.out = kernel$n)))

    kernel_weight <- mapply(
      function(grid, bw) kernel$kernel((kernel$var - grid)/bw),
      grid = kernel_grid,
      bw = with(self$data, with(kernel, eval(bw))),
      SIMPLIFY = FALSE)
    if(is.null(kernel$min_points)) kernel$min_points <- 2*k

    which_drop <- kernel_weight %>%
      sapply(function(x) sum(x>0)<kernel$min_points)
    if(all(which_drop)) abort("Not enough data under current min_points.")
    kernel_grid_f <- function() kernel_grid[!which_drop]
    kernel_weight_f <- function() kernel_weight[!which_drop]

    coef_path <- list()
    gamma_path <- list()
    ssr_path <- list()

    # 1
    term_compared <- rhs[[1]]$ind$ind_expression[[1]]
    temp_var <- model_df %>%
      as_tibble() %>%
      dplyr::transmute(var = !!term_compared) %>%
      unlist() %>%
      unname()
    all_gamma_grids <- lapply(
      kernel_weight_f(),
      function(w) {
        stats::na.omit(unique(temp_var[w>0]))
      }
    )


    which_min <- lapply(seq_along(all_gamma_grids), function(i){
      # which_min <- lapply(seq_along(gamma_grids), function(i){
      i <- enexpr(i)
      eval(expr(sapply(
        all_gamma_grids[[!!i]] %>%
          lapply(`names<-`, gamma_env$gamma[[1]]),
        get_ssr,
        fm = fm,
        weights = kernel_weight_f()[[!!i]],
        data = model_df)
      ))
    }) %>%
      sapply(which.min)

    min_gamma <- mapply(function(.x, .i).x[.i], all_gamma_grids, which_min)




    gamma_col <- numeric(nrow(model_df))
    gamma_col[!which_drop] <- min_gamma
    gamma_col[which_drop] <- NA

    model_df_k <- model_df %>%
      mutate(.gamma=gamma_col)


    for(iter in seq_len(kernel$max_iter)){
      # !is.na(model_df_k$.gamma)
      # 2
      model_k <- eval(expr(lm(!!get_fm_fun(fm)(expr(.gamma)), data = model_df_k)))
      ssr_path[[iter]] <- sum((residuals(model_k))^2)

      coef_k <- coef_path[[iter]] <- coef(model_k)

      # 3
      model_frame <- model.frame(model_k)

      all_x_name <- intersect(
        names(coef_k),
        colnames(model_frame)
      )
      other_x_name <- all_x_name[!grepl("\\.gamma", all_x_name)]
      other_x_value <- model_frame[,other_x_name] %>%
        mapply(function(vec, coef) vec*coef,
               vec = ., coef = coef_k[other_x_name])
      resp_value <- getElement(model_frame, resp)

      item_before <- resp_value -
        coef_k[["(Intercept)"]] -
        rowSums(other_x_value)

      resp_1_value <- getElement(model_frame, deparse(rhs[[1]]$xreg$expression))
      item_times <- resp_1_value *
        coef_k[all_x_name[grepl("\\.gamma", all_x_name)]]

      find_multiplier_c <- function(gamma){
        temp <- resp_1_value < gamma
        if(sum(temp, na.rm = TRUE) < k || sum(!temp, na.rm = TRUE)>k) return(Inf)
        (item_before - item_times *temp)^2
      }

      temp_grids <- model_df %>%
        as_tibble() %>%
        dplyr::transmute( !!rhs[[1]]$ind$ind_expression[[1]]) %>%
        unlist(use.names =FALSE)
      # res_c <- lapply(kernel_weight, function(weights){
      temp_which_drop <- which_drop
      gamma_res <- res_c <- purrr::map2_dbl(
        kernel_weight_f(),
        seq_along(kernel_weight)[!which_drop],
        function(weights, ix){
          idx <- sapply(unique(temp_grids[weights>0]),
                        function(gamma){
                          a <- find_multiplier_c(gamma)
                          if(any(is.infinite(a))) return(Inf)
                          b <- weights[!`[<-`(which_drop, 1, TRUE)] # drop the first with na
                          # stopifnot(length(a) == length(b))
                          # if(length(a) != length(b)) {browser();abort()}
                          out <- sum(a*b)
                          stopifnot(!all(is.na(out)))
                          out
                        })
          if(all(is.infinite(idx))) {
            temp_which_drop[ix] <<- TRUE
            return(NA)
          }

          return(unique(temp_grids[weights>0])[which.min(idx)])
        })

      which_drop <- temp_which_drop

      # gamma_res <- sapply(res_c, function(x) x$minimum)
      gamma_path[[iter]] <- tibble(!!sym(kernel$var_name):=kernel_grid_f(), .gamma = na.omit(gamma_res))
      model_df_k <- model_df_k %>%
        select(-.gamma) %>%
        dplyr::left_join(distinct(gamma_path[[iter]]), by = as_string(kernel$var_name))

      if(iter>1){
        if(norm(coef_path[[iter]] - coef_path[[iter-1]], "2")< 1){
          break
        }
      }

    }

    path <- list(ssr_path = ssr_path,
                 coef_path = coef_path,
                 gamma_path = gamma_path)
    gamma_grid <- gamma_res
    .gamma <-  gamma_path[iter]

    interpo_func <- approxfun(getElement(.gamma[[1]], as_string(kernel$var_name)),
                              y = .gamma[[1]]$.gamma, rule = 2)
    df <- model_df %>%
      as_tibble() %>%
      left_join(distinct(.gamma[[1]]), by = as_string(kernel$var_name)) %>%
      mutate(.gamma = case_when(
        is.na(.gamma) ~ interpo_func(!!kernel$var_name),
        TRUE ~ .gamma
      ))

    mdl <- eval(expr(lm(!!get_fm_fun(fm)(expr(.gamma)), data = df)))



    .resid <- c(NA, residuals(mdl))
  } else {
    if(multi_regime){

      term_com <- unique(sapply(ind_term, function(x) x$ind$ind_expression))

      if(length(term_com)!=1)
        abort("Term to compare in ind terms should be identical.")

      term_com <- term_com[[1]]

      if(!parametric) {

        temp_weights_lgl <- rep(TRUE, nrow(model_df))
        temp_env <- new.env()
        for(i in seq_along(gamma_env$gamma)){

          all_gamma_grids <- get_grid_single_gamma(ind_term[i], gamma_name = gamma_env$gamma[i],
                                                   data_eval = model_df,
                                                   weights_lgl = temp_weights_lgl)
          gamma_in_fm <- fm_top %>%
            find_leaf() %>%
            sapply(deparse) %>%
            intersect(unlist(gamma_env$gamma))
          all_gamma_grids <- lapply(all_gamma_grids, function(x) {
            out <- rep(x, length(gamma_in_fm))
            names(out) <- gamma_in_fm
            out
          })
          ssr <- all_gamma_grids %>% sapply(get_ssr, fm_top, data = model_df, weights = as.numeric(temp_weights_lgl))
          # path <- all_gamma_grids %>% bind_rows() %>% mutate(ssr)
          .gamma <- unique(all_gamma_grids[[which.min(ssr)]])

          assign(gamma_env$gamma[[i]], .gamma, envir = temp_env)

          temp_weights_lgl <- term_com %>%
            eval_tidy(data = model_df) %>%
            {!closest(., .gamma, k+length(gamma_env$gamma), exclude = !temp_weights_lgl)} %>%
            `&`(temp_weights_lgl)


        }
        .gamma <- unlist(as.list(temp_env, all.names = TRUE))

        .gamma <- sort(.gamma, decreasing = TRUE)
        names(.gamma) <- sort(unlist(gamma_env$gamma))

        path <- NULL
      } else {
        #parametric

        temp_weights_lgl <- rep(TRUE, nrow(model_df))
        temp_env <- new.env()
        for(i in seq_along(gamma_ind_pair)){
          gamma_grids <- all_gamma_grids[[i]]

          gamma_in_fm <- fm_top %>%
            find_leaf() %>%
            sapply(deparse) %>%
            intersect(unlist(gamma_env$gamma))
          gamma_grids <- lapply(gamma_grids, function(x) {
            names(x) <- gamma_in_fm
            x
          })
          ssr <- gamma_grids %>% sapply(get_ssr, fm_top, data = model_df, weights = as.numeric(temp_weights_lgl))
          # path <- all_gamma_grids %>% bind_rows() %>% mutate(ssr)
          .gamma <- unique(gamma_grids[[which.min(ssr)]])

          walk2(paste0(".gamma_", gamma_ind_pair[[i]]), .gamma, assign, envir = temp_env)




          temp_weights_lgl <- term_com %>%
            eval_tidy(data = model_df) %>%
            closest(.gamma[[1]] + .gamma[[2]] * eval_tidy(based_var, data = model_df), ## assume quadratic
                    ( k+length(gamma_env$gamma))*10, exclude = !temp_weights_lgl) %>%
            `!`() %>%
            `&`(temp_weights_lgl)


        }
        .gamma <- unlist(as.list(temp_env, all.names = TRUE))

        path <- NULL
      }
    } else {

      ssr <- all_gamma_grids %>% sapply(get_ssr, fm)
      path <- all_gamma_grids %>% bind_rows() %>% mutate(ssr)
      .gamma <- all_gamma_grids[[which.min(ssr)]]
    }
    mdl <- eval(expr(lm(!!do.call(get_fm_fun(fm), as.list(.gamma)), data = model_df,
                        x=TRUE, y=TRUE)))

    .resid <- residuals(mdl)
  }

  structure(
    list(model = mdl,
         est = list(
           .gamma = .gamma,
           .path = path,
           .resid = .resid
         ),
         constrained = !is.null(one_term),
         multi_regime = multi_regime,
         kernel = kernel_est,
         data = self$data,
         num_regime = length(gamma_env$gamma)+1,
         fms = list(one_expr = one_term$expression,
                    ind_exprs= ind_term %>%
                      map(function(x) x$ind$expression))

    ),
    class = "fbl_thrreg"
  )

}

#' Estimate a Threshold Regression
#'
#' Searches through the grids generated from the specified variable to identify
#' the threshold level with minimised sum of squared error.
#' Supports estimation of (constrained) single regime (a.k.a with only one threshold),
#' multiple regimes (searched sequentially), variable dependent regime(s) with second order polynomial
#' (experimental; user supplied grids required),
#' and nonparametric function (experimental).
#'
#' @param formula Model specification (see "Specials" section).
#' @param ... Further arguments. Currently not used.
#'
#' @section Parameterisation:
#'
#' One regime threshold regression follows the form of
#'
#' \deqn{y_t = \beta_0 + \delta_0 y_{t-1} + \delta_1  y_{t-1}I(|y_{t-1}| \ge \gamma ) + controls + \epsilon_t}
#'
#' Where \eqn{I} is the indicator function.
#' The interface of THRREG allows user to specify the model as how it is written.
#' \eqn{ \delta_0} (or any other coefficient) can be optionally constrained to be 1 with special \code{one}.
#' The indicator can be evaluated with special \code{ind}.
#' See "Examples" and "Specials" section.
#'
#' @section Specials:
#'
#' The _specials_ define the regimes that \code{THRREG} will try to find.
#'
#'
#'
#' \subsection{one}{
#' The \code{one} special is used to constrain the coefficient of input variable to be 1.
#' When used on simple terms without indication function and threshold,
#' it is the same as [stats::offset()] but takes the calculation literally,
#' i.e. \code{-} is supported.
#' \preformatted{
#' one(arg)
#' }
#'
#' \tabular{ll}{
#'   \code{arg} \tab Regressor in the regression whose coefficient is contrained to be 1 \cr
#' }
#' }
#'
#' \subsection{ind and gamma}{
#' \code{gamma} defines the threshold location when the regimes are constant, or
#' the parameters that define the threshold location when the regimes depend on
#' other variables.
#' \code{ind} is used in conjunction with \code{gamma} to specify the indicator function
#' that seperates different regimes.
#' \preformatted{
#' gamma(id = NULL, var = NULL, grid = NULL)
#' }
#' \tabular{ll}{
#'   \code{id} \tab The identifier of parameter \eqn{\gamma}, namely the subscript.
#'   In the case of multiple constant regimes, this needs to be a numeric value
#'   indicating the relative position of the threshold ranked from large to small.
#'   For exmaple, threshold \code{gamma(1)} is above \code{gamma(2)} and so on.
#'   In the case of second order polynomial threshold, this needs to be a numeric value
#'   indicating the order of which the grid search is conducted, and the intercept in the
#'   polynomial should has a smaller id than the paired slope.
#'   For example, the threshold can be defined as \code{gamma(1) + gamma(2)*z}
#'   where \code{z} is the variable on which the threshold depends.\cr
#'   \code{var} \tab (Experimental) Variable in nonparametric estimation of the threshold. \cr
#'   \code{grid} \tab (Experimental) The grid of which to search for each parameter in the case of
#'   second order polynomial. For example \code{gamma(1, grid = 0:50) + gamma(2, grid = 0:20)*z}
#'   defines a threshold as a second order polynomial of variable \code{z}.
#'   \code{gamma(1)} can take any value in \code{0:50} and
#'   \code{gamma(1)} can take \code{0:20}.
#'   In constant regimes, the grids are the sample value of the comparing variable.  \cr
#' }
#'
#' \preformatted{
#' ind(expression)
#' }
#'
#' \tabular{ll}{
#'   \code{expression} \tab Input of the indicator function.
#'   Normally a comparison of transformation of response variable with the threshold defined by \code{gamma}.
#'   For example, for the one regime case this can be \code{abs(lag(y)) >= gamma(1)} \cr
#' }
#'
#' }
#' \subsection{xreg}{
#' Exogenous regressors can be included without explicitly using the \code{xreg()} special.
#' Interactions and other functionality behaves similarly to [stats::lm()].
#' }
#'
#' @return A model specification.
#'
#' @examples
#' # Simulate a time series
#' threshold_process <- function(lag1){
#'   if(abs(lag1) > 1) {
#'     # Above threshold of 1
#'     # Mean reverting
#'     # AR1 process with coeffcient 0.8
#'     lag1 * 0.8 + rnorm(1)
#'   } else {
#'     # Below threshold of 1
#'     # a unit root process
#'     lag1 + rnorm(1)
#'   }
#' }
#' set.seed(2222)
#' time_span <- 2000
#' y <- numeric(time_span)
#' for(i in 2:time_span) y[[i]] <- threshold_process(y[[i-1]])
#'
#' # Convert to tsibble
#' library(tsibble)
#' df <- tsibble(y = y, idx = seq_along(y), index = idx)
#'
#' # Fit the threshold regression
#' library(fable.thrreg)
#' fit <- df %>%
#'   model(thrreg = THRREG(y ~ offset(lag(y)) + lag(y)*ind(abs(lag(y)) >= gamma(1) )))
#' # Getting estimates
#' est <- tidy(fit)
#' est
#'
#' # Estimated threshold value
#' .gamma_1 <- est$estimate[est$term == ".gamma_1"]
#' # Plot
#' if(requireNamespace("ggplot2")){
#' autoplot(df, y) +
#'   ggplot2::geom_hline(yintercept = c(-1, 1)) +
#'   ggplot2::geom_hline(yintercept = c(-.gamma_1, .gamma_1), colour = "red")
#' }
#' @export
THRREG <- function(formula, ...){
  thrreg_model <- new_model_class("thrreg", train_thrreg, specials_thrreg)
  new_model_definition(thrreg_model, !!enquo(formula), ...)
}


#' @export
fitted.fbl_thrreg <- function(object, ...){
  c(NA, fitted(object$model, ...))
}

#' @export
residuals.fbl_thrreg <- function(object, ...){
  object$est[[".resid"]]
}

#' @export
tidy.fbl_thrreg <- function(x, ...){
  broom::tidy(x$model, ...) %>%
    # broom::tidy(x$model) %>%
    # mutate(estimate = `class<-`(estimate, class(x$est$.gamma))) %>%
    bind_rows(tibble(term = names(x$est$.gamma), estimate = as.numeric(x$est$.gamma)))

}

#' @export
report.fbl_thrreg <- function(object, ...){
  print(summary(object$model, ...))
}


#' @export
model_sum.fbl_thrreg <- function(x) {
  "thrreg"
}

