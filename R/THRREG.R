

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
    if(is.null(id)) id <- max(unlist(gamma_env$id[sapply(gamma_env$id, is.numeric)]) %||% 0, na.rm = TRUE) + 1
    g <- paste0(".gamma_", id)
    gamma_env$id <- c(gamma_env$id, id)
    gamma_env$gamma[[id]] <- g
    sym(g)
  }, env = gamma_fun_env)


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
      gamma_env$gamma_special[as.character(unlist(gamma_env$id))] <- x$ind[names(x$ind) == "gamma"]

      expr(I(!!x$xreg$expression * (!!ind_expr)))
    } else {
      x$xreg$expression
    }
  }
  )

  model_data <- squash_tibble(rhs, 3:4)

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
  }


  kernel_est <- FALSE
  if(any(sapply(gamma_env$gamma_special, function(x) !is.null(x$var)))){
    kernel_est <- TRUE
    if(length(gamma_env$id)>1)
      abort("Only one gamma term can appear in the formula when using kernel estimation.")
  }






  resp <- measured_vars(.data)
  n <- nrow(model_data)
  k <- ncol(select(rhs[[which(sapply(rhs, length)==1)]]$xreg$xregs, !any_of(resp)))+2




  # Finding the gamma grids
  if(!kernel_est){

    if(length(unlist(gamma_env$gamma)) == length(ind_term)){
      # Single regime
      parametric <- FALSE
      all_gamma_grids <- map(ind_term, function(x){
        eval_tidy(x$ind$ind_expression[[1]], data = x$ind$xreg$xregs)
      }
      ) %>%
        split(unlist(gamma_env$gamma)) %>%
        lapply(function(x) x %>%
                 unlist() %>%
                 unique() %>%
                 sort() %>%
                 .[seq(k, n-k, by = 1)] )
      all_gamma_grids <- all_gamma_grids %>%
        expand.grid() %>%
        split(row(.)) %>%
        lapply(unlist)

    } else if(length(unlist(gamma_env$gamma)) > length(ind_term)){
      # parametrix ploy 2
      parametric <- TRUE
      if(length(ind_term)>1) {
        abort("Only one indicator function allowed when using parametric gamma.")
      }

      ind_single <- ind_term[[1]]
      gamma_idx <- ind_single$ind$ind_expression %>%
        map(function(x) find_leaf(x, exclude = "gamma") %>%
              sapply(is_call_name, "gamma") %>%
              any()) %>%
        unlist() %>%
        which()

      x2 <- ind_single$ind$ind_expression[[gamma_idx]] %>%
        find_leaf() %>%
        sapply(deparse)

      given_grids <- gamma_env$gamma_special %>%
        lapply(getElement, "grid")
      if(any(sapply(given_grids, is.null))) abort("You need to provide grid to search for gamma.")
      all_gamma_grids <- given_grids %>%
        do.call(expand.grid, .) %>%
        `colnames<-`(unique(unlist(gamma_env$gamma))) %>%
        as_tibble() %>%
        # dplyr::filter(.gamma_2 < 270 - 12 * .gamma_1) %>%
        split(seq_len(nrow(.))) %>%
        lapply(unlist)

    }
  }

  fm <- if(!is.null(one_term)){
    expr(I(!!sym(resp) - !!one_expr) ~ !!reduce(rhs_terms, function(.x, .y) call2("+", .x, .y)))
  } else {
    expr(!!sym(resp) ~ !!reduce(rhs_terms, function(.x, .y) call2("+", .x, .y)))
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

    if(length(gammas) == 2){
      # test for identification
      test <- function(){
        temp_env <- new.env()
        purrr::walk2(names(gammas), gammas,assign, envir = temp_env)
        temp <- ind_single$ind$expression %>%
          replace_gamma_lang() %>%
          eval_tidy( data = eval(data), env = temp_env)
        if((sum(temp, na.rm = TRUE) <= k ) || (sum(!temp, na.rm = TRUE) <= k) ) {
          return(TRUE)
        }
        FALSE
      }

      if(test()) return(Inf)
    }


    fit <- eval(expr(lm(!!fm, data = !!data, model = FALSE, weights = !!weights)))
    if(anyNA(coef(fit))){
      if(kernel_est)
        warning("Not all coefficient can be estimated. Consider increase bandwidth bw or raise number of min_points.")
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

    which_drop <- kernel_weight %>%
      sapply(function(x) sum(x>0)<kernel$min_points)
    kernel_grid <- kernel_grid[!which_drop]
    kernel_weight <- kernel_weight[!which_drop]

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
      kernel_weight,
      function(w) {
        stats::na.omit(unique(temp_var[w>0]))
      }
    )


    which_min <- pbapply::pblapply(seq_along(all_gamma_grids), function(i){
      # which_min <- lapply(seq_along(gamma_grids), function(i){
      i <- enexpr(i)
      eval(expr(sapply(
        all_gamma_grids[[!!i]] %>%
          lapply(`names<-`, gamma_env$gamma[[1]]),
        get_ssr,
        fm = fm,
        weights = kernel_weight[[!!i]],
        data = model_df)
      ))
    }) %>%
      sapply(which.min)

    min_gamma <- mapply(function(.x, .i).x[.i], all_gamma_grids, which_min)


    model_df_k <- model_df %>%
      # mutate(!!rhs[[1]]$xreg$expression) %>%
      left_join(
        tibble(!!kernel$var_name:= kernel_grid, .gamma = min_gamma) %>%
          distinct(),
        by = as_string(kernel$var_name)
      )


    pb <- lazybar::lazyProgressBar(kernel$max_iter)
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
        (item_before -
           item_times *
           (resp_1_value < gamma))^2
      }

      temp_grids <- model_df %>%
        as_tibble() %>%
        dplyr::transmute( !!rhs[[1]]$ind$ind_expression[[1]]) %>%
        unlist(use.names =FALSE)

      # res_c <- lapply(kernel_weight, function(weights){
      gamma_res <- res_c <- sapply(kernel_weight, function(weights){
        idx <- sapply(unique(temp_grids[weights>0]),
                      function(gamma){
                        a <- find_multiplier_c(gamma)
                        b <- weights[!`[<-`(which_drop, 1, TRUE)] # drop the first with na
                        sum(a*b)
                      }) %>%
          which.min()
        unique(temp_grids[weights>0])[idx]
      })

      # gamma_res <- sapply(res_c, function(x) x$minimum)
      gamma_path[[iter]] <- tibble(!!sym(kernel$var_name):=kernel_grid, .gamma = gamma_res)
      model_df_k <- model_df_k %>%
        select(-.gamma) %>%
        dplyr::left_join(distinct(gamma_path[[iter]]), by = as_string(kernel$var_name))

      pb$tick()$print()
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
    if(FALSE){
      # if(parametric){

      .gamma <- c(0, 1) %>%
        `names<-`(unique(unlist(gamma_env$gamma))) %>%
        optim(get_ssr,fm = fm) %>%
        .$par


    } else {

      ssr <- all_gamma_grids %>% pbapply::pbsapply(get_ssr, fm)
      path <- all_gamma_grids %>% bind_rows() %>% mutate(ssr)
      .gamma <- all_gamma_grids[[which.min(ssr)]]
    }
    mdl <- eval(expr(lm(!!do.call(get_fm_fun(fm), as.list(.gamma)), data = model_df)))

    .resid <- residuals(mdl)
  }

  structure(
    list(model = mdl,
         est = list(
           .gamma = .gamma,
           .path = path,
           .resid = .resid
         )),
    class = "fbl_thrreg"
  )

}


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
    mutate(estimate = `class<-`(estimate, class(x$est$.gamma))) %>%
    bind_rows(tibble(term = names(x$est$.gamma), estimate = x$est$.gamma))

}

#' @export
report.fbl_thrreg <- function(object, ...){
  print(summary(object$model, ...))
}


#' @export
model_sum.fbl_thrreg <- function(x) {
  "thrreg"
}

