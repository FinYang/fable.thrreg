




train_thrreg <- function(.data, specials, ...){
  if (length(tsibble::measured_vars(.data)) > 1) {
    stop("Only univariate response is supported by thrreg.")
  }

  # .data <- bitcoin %>%
  #   select(t, y)
  # xreg <- bitcoin %>%
  #   as_tibble() %>%
  #   select(vlm,vlm_K,vlt,vlt_K,Illiq,Illiq_K,trn,trn_K)
  # delta <- TRUE
  # kernel <- with(
  #   list(self = list(data = bitcoin)),
  #   (function(var = NULL, kernel = epaker,
  #             # n = 2^8,
  #             min_points = 30, bw = sd(var)/25, max_iter = 10){
  #     var_name <- rlang::enexpr(var)
  #     bw <- rlang::enexpr(bw)
  #     if(!is.null(var_name)){
  #       var <- dplyr::select(as_tibble(self$data), !!var_name) %>%
  #         unlist() %>%
  #         unname()
  #     }
  #     as.list(environment())
  #   })(fee))

  delta <- specials$delta[[1]]
  xreg <- specials$xreg[[1]]
  kernel <- specials$gamma[[1]]
  browser()

  kernel_est <- !is.null(kernel$var)
  if(kernel_est && !delta){
    delta <- TRUE
    warning("kernel estimation is only supported for unconstraint case. Setting delta = TRUE.")
  }


  n <- nrow(xreg)
  k <- ncol(xreg) + 2
  resp <- measured_vars(.data)
  resp_1 <- paste0(resp, "_1")

  # Find gamma grid
  y_1 <- dplyr::lag(getElement(.data, resp))
  gamma_grid <- y_1 %>%
    abs() %>%
    unique() %>%
    sort()
  gamma_grid <- gamma_grid[seq(k, n-k, by = 1)]

  model_df <- bind_cols(select(as_tibble(.data), !!sym(resp)),
                        !!sym(resp_1) := y_1,
                        xreg)
  get_full_fm <- function(.gamma){
    if(delta){
      fm <- expr(
        !!sym(resp) ~
          !!sym(resp_1) +
          I(!!sym(resp_1) * (abs(!!sym(resp_1)) <= !!.gamma)) +
          !!reduce(syms(colnames(xreg)), function(.x, .y) call2("+", .x, .y)))
    } else {
      fm <- expr(
        I(!!sym(resp) - !!sym(resp_1) * (abs(!!sym(resp_1)) <= !!.gamma)) ~
          I(!!sym(resp_1) -
              !!sym(resp_1) * (abs(!!sym(resp_1)) <= !!.gamma)) +
          !!reduce(syms(colnames(xreg)), function(.x, .y) call2("+", .x, .y)))
    }
    fm
  }


  get_ssr <- function(gamma, fm, weights = NULL, data = model_df){
    fm <- eval(enexpr(fm))
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
    ssr <- sapply(gamma_grid, get_ssr, fm = get_full_fm(gamma))
    path <- tibble(gamma_grid, ssr)
    .gamma <- gamma_grid[which.min(ssr)]
    mdl <- eval(expr(lm(!!get_full_fm(gamma_grid[which.min(ssr)]), data = model_df)))
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

