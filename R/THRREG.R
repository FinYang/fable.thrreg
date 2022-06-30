

#' @importFrom dplyr lag
#' @export
dplyr::lag


train_thrreg <- function(.data, specials, ...){
  lag <- dplyr::lag
  if (length(tsibble::measured_vars(.data)) > 1) {
    stop("Only univariate response is supported by thrreg.")
  }

  gamma_env <- new.env(parent=emptyenv())
  gamma_env$i <- 1
  gamma_env$gamma <- list()
  gamma <- function(id){
    g <- paste0(".gamma_", id)
    gamma_env$gamma[[gamma_env$i]] <- g
    gamma_env$i <- gamma_env$i + 1
    sym(g)
  }

  replace_gamma_lang <- function(ind_expr){
    gamma_pos <- which(sapply(as.list(ind_expr), function(y) is_call_name(y, "gamma")))
    ind_expr[[gamma_pos]] <- eval(ind_expr[[gamma_pos]])
    ind_expr
  }
  squash_tibble <- function(x, layers){
    lapply(layers, function(layer)
      find_tibble(x, layer)
    ) %>%
      quietly_squash() %>%
      purrr::reduce(function(x, y){
        bind_cols(x, select(y, !any_of(colnames(x)))) })
  }

  rhs <- specials$xreg[[1]]
  ind_term <- rhs[sapply(rhs, length)>1]



  gamma_terms <- unique(unlist(gamma_env$gamma))


  rhs_terms <- map(rhs, function(x){
    if("ind" %in% names(x)) {
      ind_expr <- x$ind$expression
      ind_expr <- replace_gamma_lang(ind_expr)

      expr(I(!!x$xreg$expression * (!!ind_expr)))
    } else {
      x$xreg$expression
    }
  }
  )

  find_tibble <- function(x, layer){
    x %>%
      purrr::map_depth(layer, function(x) if(is_tibble(x)) x else NULL) %>%
      quietly_squash() %>%
      {.[!sapply(., is.null)]}
  }



  model_data <- squash_tibble(rhs, 3:4)


  resp <- measured_vars(.data)
  n <- nrow(model_data)
  k <- ncol(select(rhs[[which(sapply(rhs, length)==1)]]$xreg$xregs, !any_of(resp)))+2

  model_df <- .data %>%
    bind_cols(select(model_data, !any_of(colnames(.))))

  one_term <- specials$one[[1]]
  if(!is.null(one_term)){
    one_expr <- one_term$expression %>%
      deparse() %>%
      gsub("ind", "", .) %>%
      gsub("gamma\\(([[:digit:]])*\\)", ".gamma_\\1", .) %>%
      str2lang()
    one_data <- one_term$data
    one_data <- squash_tibble(one_data, 2:3)
    model_df <- model_df %>%
      bind_cols(select(one_data, !any_of(colnames(.))))
  }


  gamma_grids <- map(ind_term, function(x){
    eval_tidy(x$ind$ind_expression[[1]], data = x$ind$xreg$xregs)
  }
  ) %>%
    split(unlist(gamma_env$gamma)) %>%
    lapply(function(x) x %>%
             unlist() %>%
             unique() %>%
             sort() %>%
             .[seq(k, n-k, by = 1)] )




  fm <- if(!is.null(one_term)){
    expr(I(!!sym(resp) - !!one_expr) ~ !!reduce(rhs_terms, function(.x, .y) call2("+", .x, .y)))
  } else {
    expr(!!sym(resp) ~ !!reduce(rhs_terms, function(.x, .y) call2("+", .x, .y)))
  }

  get_fm_fun <- function(fm) {
    fm_str <- deparse(fm) %>%
      stringr::str_trim(side = "left") %>%
      paste0(collapse = "")
    purrr::walk(names(gamma_grids), function(.x ){
      fm_str <<- gsub(.x, paste0("!!", .x), fm_str)
    })

    list_from_name <- function(name){
      `names<-`(list(rep(NULL, length(gamma_grids))), name)
    }
    rlang::new_function(
      list_from_name(names(gamma_grids)),
      expr(str2lang(paste0(deparse(expr(!!str2lang(fm_str))), collapse = "")))
    )
  }

  get_ssr <- function(gammas, fm, weights = NULL, data = model_df){
    gammas <- round(gammas, 6)
    fm <- do.call(get_fm_fun(fm), list(gammas))
    weights <- enexpr(weights)
    data <- enexpr(data)
    fit <- eval(expr(lm(!!fm, data = !!data, model = FALSE, weights = !!weights)))
    if(anyNA(coef(fit))){
      if(anyNA(coef(fit)[-3]))
        warning("Not all coefficient can be estimated. Consider increase bandwidth bw or raise number of min_points.")
      return(Inf)
    }
    sum(resid(fit)^2)
  }

  get_ssr2 <- function(gammas, fm, weights = NULL, data = model_df){
    fm <- do.call(get_fm_fun(fm), list(gammas))
    weights <- enexpr(weights)
    data <- enexpr(data)
    fit <- eval(expr(lm(!!fm, data = !!data, model = FALSE, weights = !!weights)))
    if(anyNA(coef(fit))){
      if(anyNA(coef(fit)[-3]))
        warning("Not all coefficient can be estimated. Consider increase bandwidth bw or raise number of min_points.")
      return(Inf)
    }
    sum(resid(fit)^2)
  }
  kernel_est <- FALSE
  if(kernel_est){
    model_df <- bind_cols(model_df,
                          tibble(!!kernel$var_name := kernel$var))
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

    # 1
    gamma_grids <- lapply(kernel_weight, function(w) stats::na.omit(unique(abs(y_1[w>0]))))

    which_min <- pbapply::pblapply(seq_along(gamma_grids), function(i){
      # which_min <- lapply(seq_along(gamma_grids), function(i){
      i <- enexpr(i)
      eval(expr(sapply(
        gamma_grids[[!!i]], get_ssr,
        fm = get_full_fm(gamma),
        weights = kernel_weight[[!!i]][kernel_weight[[!!i]]>0],
        data = model_df[kernel_weight[[!!i]]>0,])
      ))
    }) %>%
      sapply(which.min)

    min_gamma <- mapply(function(.x, .i).x[.i], gamma_grids, which_min)


    model_df_k <- tibble(!!kernel$var_name:= kernel_grid, .gamma = min_gamma) %>%
      distinct() %>%
      dplyr::inner_join(model_df, by = rlang::as_string(kernel$var_name)) %>%
      tidyr::drop_na(!!sym(resp_1))



    for(iter in seq_len(kernel$max_iter)){
      # !is.na(model_df_k$.gamma)
      # 2

      model_k <- eval(expr(lm(!!get_full_fm(expr(.gamma)), data = model_df_k)))

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
      resp_1_value <- getElement(model_frame,resp_1)
      item_before <- resp_value -
        coef_k[["(Intercept)"]] -
        rowSums(other_x_value)
      item_times <- resp_1_value *
        coef_k[all_x_name[grepl("\\.gamma", all_x_name)]]
      find_multiplier_c <- function(gamma){
        (item_before -
           item_times *
           (resp_1_value < gamma))^2
      }

      # res_c <- lapply(kernel_weight, function(weights){
      gamma_res <- res_c <- sapply(kernel_weight, function(weights){
        idx <- sapply(abs(unique(kernel$var[weights>0])),
                      function(gamma){
                        a <- find_multiplier_c(gamma)
                        b <- weights[!`[<-`(which_drop, 1, TRUE)] # drop the first with na
                        # if(length(a)!=length(b)) {
                        #   browser()
                        # }
                        sum(a*b)
                        # sum(find_multiplier_c(gamma) * weights[!which_drop])
                      }) %>%
          which.min()
        unique(kernel$var[weights>0])[idx]
      })

      # gamma_res <- sapply(res_c, function(x) x$minimum)
      gamma_path[[iter]] <- tibble(!!sym(kernel$var_name):=kernel_grid, .gamma = gamma_res)
      model_df_k <- model_df_k %>%
        select(-.gamma) %>%
        dplyr::left_join(distinct(gamma_path[[iter]]), by = as_string(kernel$var_name))

      if(iter>1){
        if(norm(coef_path[[iter]] - coef_path[[iter-1]], "2")< 1){
          break
        }
      }

    }

    path <- list(coef_path = coef_path,
                 gamma_path = gamma_path)
    gamma_grid <- gamma_res
    .gamma <-  gamma_path[iter]

    interpo_func <- approxfun(.gamma[[1]]$fee, y = .gamma[[1]]$.gamma, rule = 2)
    df <- model_df %>%
      left_join(.gamma[[1]], by = as_string(kernel$var_name)) %>%
      mutate(.gamma = case_when(
        is.na(.gamma) ~ interpo_func(!!kernel$var_name),
        TRUE ~ .gamma
      ))

    mdl <- eval(expr(lm(!!get_full_fm(expr(.gamma)), data = df)))




    .resid <- c(NA, residuals(mdl))
  } else {

    all_gamma_grids <- gamma_grids %>%
      expand.grid() %>%
      split(row(.)) %>%
      lapply(unlist)
    ssr <- all_gamma_grids %>%sapply(get_ssr, fm)
    # path <- tibble(gamma_grid, ssr)
    .gamma <- all_gamma_grids[[which.min(ssr)]]
    mdl <- eval(expr(lm(!!get_fm_fun(fm)(.gamma), data = model_df)))


    .resid <- residuals(mdl)
  }

  structure(
    list(model = mdl,
         est = list(
           .gamma = .gamma,
           # .path = path,
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
    bind_rows(tibble(term = ".gamma", estimate = x$est$.gamma))

}

#' @export
report.fbl_thrreg <- function(object, ...){
  print(summary(object$model, ...))
}


#' @export
model_sum.fbl_thrreg <- function(x) {
  "thrreg"
}

