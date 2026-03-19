#' Estimate a Mixture via Efron-Style Deconvolution
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param x_thresh Numeric threshold used to define successes.
#' @param tau Support grid for the estimated mixture.
#' @param pDegree Degrees of freedom for the spline basis in [deconv_fast()].
#' @param c0 Penalty constant passed to [deconv_fast()].
#' @param weight_col Optional weight column name. Not supported for this method.
#'
#' @return A standardized mixture data frame with columns `theta`, `g`, and
#'   `cumul`.
#' @export
estimate_mixture_efron <- function(df,
                                   id_col = "id",
                                   val_col = "x",
                                   x_thresh,
                                   tau = seq(from = 0.005, to = 0.995, by = 0.005),
                                   pDegree = 7,
                                   c0 = 1,
                                   weight_col = NULL) {
  if (!is.null(weight_col)) {
    stop("`weight_col` is not supported for method = 'spline'.", call. = FALSE)
  }

  tau <- sort(unique(as.numeric(tau)))

  if (anyNA(tau) || any(tau < 0 | tau > 1)) {
    stop("`tau` must be numeric and lie in [0, 1].", call. = FALSE)
  }

  df <- .standardize_input_df(df, id_col = id_col, val_col = val_col)

  df_bin <- .aggregate_binomial_data(
    df = df,
    x_thresh = x_thresh,
    weight_col = NULL
  )

  x_input <- dplyr::select(df_bin, n, s)
  result <- deconv_fast(
    tau = tau,
    X = x_input,
    family = "Binomial",
    pDegree = pDegree,
    c0 = c0
  )

  stats_df <- data.frame(result$stats)

  .standardize_mixture_df(theta = stats_df$theta, g = stats_df$g)
}

#' Estimate a Mixture via NPMLE
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param x_thresh Numeric threshold used to define successes.
#' @param tau Support grid for the estimated mixture.
#' @param weight_col Optional column name containing id-level weights.
#' @param mixsqp_control Optional control list passed to [mixsqp::mixsqp()].
#'
#' @return A standardized mixture data frame with columns `theta`, `g`, and
#'   `cumul`.
#' @export
estimate_mixture_npmle <- function(df,
                                   id_col = "id",
                                   val_col = "x",
                                   x_thresh,
                                   tau = seq(from = 0, to = 1, by = 0.005),
                                   weight_col = NULL,
                                   mixsqp_control = list()) {
  tau <- sort(unique(as.numeric(tau)))

  if (anyNA(tau) || any(tau < 0 | tau > 1)) {
    stop("`tau` must be numeric and lie in [0, 1].", call. = FALSE)
  }

  df <- .standardize_input_df(
    df,
    id_col = id_col,
    val_col = val_col,
    weight_col = weight_col
  )

  df_bin <- .aggregate_binomial_data(
    df = df,
    x_thresh = x_thresh,
    weight_col = weight_col
  )

  if (is.null(mixsqp_control$verbose)) {
    mixsqp_control$verbose <- FALSE
  }

  logL <- outer(
    X = seq_len(nrow(df_bin)),
    Y = seq_along(tau),
    FUN = Vectorize(function(i, j) {
      stats::dbinom(
        x = df_bin$s[i],
        size = df_bin$n[i],
        prob = tau[j],
        log = TRUE
      )
    })
  )

  fit <- if (is.null(weight_col)) {
    do.call(
      mixsqp::mixsqp,
      list(
        L = logL,
        log = TRUE,
        control = mixsqp_control
      )
    )
  } else {
    do.call(
      mixsqp::mixsqp,
      list(
        L = logL,
        w = df_bin$weight,
        log = TRUE,
        control = mixsqp_control
      )
    )
  }

  .standardize_mixture_df(theta = tau, g = fit$x)
}

#' Estimate an Empirical Raw Mixture
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param x_thresh Numeric threshold used to define successes.
#' @param weight_col Optional column name containing id-level weights.
#'
#' @return A standardized mixture data frame with columns `theta`, `g`, and
#'   `cumul`.
#' @export
estimate_mixture_raw <- function(df,
                                 id_col = "id",
                                 val_col = "x",
                                 x_thresh,
                                 weight_col = NULL) {
  df <- .standardize_input_df(
    df,
    id_col = id_col,
    val_col = val_col,
    weight_col = weight_col
  )

  df_bin <- .aggregate_binomial_data(
    df = df,
    x_thresh = x_thresh,
    weight_col = weight_col
  )
  df_bin <- dplyr::mutate(df_bin, q = .data$s / .data$n)

  if (is.null(weight_col)) {
    res_raw <- dplyr::group_by(df_bin, q)
    res_raw <- dplyr::summarise(res_raw, total = dplyr::n(), .groups = "drop")
    res_raw <- dplyr::transmute(
      res_raw,
      theta = .data$q,
      g = .data$total / sum(.data$total)
    )
  } else {
    res_raw <- dplyr::group_by(df_bin, q)
    res_raw <- dplyr::summarise(
      res_raw,
      total = sum(.data$weight),
      .groups = "drop"
    )
    res_raw <- dplyr::transmute(
      res_raw,
      theta = .data$q,
      g = .data$total / sum(.data$total)
    )
  }

  if (!0 %in% res_raw$theta) {
    res_raw <- dplyr::bind_rows(data.frame(theta = 0, g = 0), res_raw)
  }

  .standardize_mixture_df(theta = res_raw$theta, g = res_raw$g)
}

#' Estimate a Beta Mixture Approximation
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param x_thresh Numeric threshold used to define successes.
#' @param tau Support grid for the estimated mixture.
#' @param weight_col Optional column name containing id-level weights.
#'
#' @return A standardized mixture data frame with columns `theta`, `g`, and
#'   `cumul`. The fitted beta parameters are attached in the `beta_mle`
#'   attribute when available.
#' @export
estimate_mixture_beta <- function(df,
                                  id_col = "id",
                                  val_col = "x",
                                  x_thresh,
                                  tau = seq(0.005, 0.995, by = 0.005),
                                  weight_col = NULL) {
  tau <- sort(unique(as.numeric(tau)))

  if (anyNA(tau) || any(tau < 0 | tau > 1)) {
    stop("`tau` must be numeric and lie in [0, 1].", call. = FALSE)
  }

  df <- .standardize_input_df(
    df,
    id_col = id_col,
    val_col = val_col,
    weight_col = weight_col
  )

  df_bin <- .aggregate_binomial_data(
    df = df,
    x_thresh = x_thresh,
    weight_col = weight_col
  )

  all_success <- all(df_bin$s == df_bin$n)
  all_failure <- all(df_bin$s == 0)

  if (all_success || all_failure) {
    g <- rep(0, length(tau))

    if (all_success) {
      idx <- which.min(abs(tau - 1))
      beta_info <- list(alpha = NA_real_, beta = NA_real_)
    } else {
      idx <- which.min(abs(tau - 0))
      beta_info <- list(alpha = NA_real_, beta = NA_real_)
    }

    g[idx] <- 1

    out <- .standardize_mixture_df(theta = tau, g = g)
    attr(out, "beta_mle") <- beta_info
    return(out)
  }

  fit_info <- .fit_betabinomial_vgam(df_bin, weight_col = weight_col)

  out <- .discretize_beta_on_tau(
    alpha = fit_info$alpha,
    beta = fit_info$beta,
    tau = tau
  )

  attr(out, "beta_mle") <- list(
    alpha = fit_info$alpha,
    beta = fit_info$beta
  )

  out
}

#' Estimate Mixtures Across a Threshold Grid
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param method Mixture estimation method.
#' @param x_grid Grid of thresholds to evaluate.
#' @param weight_col Optional column name containing id-level weights.
#' @param ... Additional arguments passed through to the chosen estimator.
#'
#' @return A named list of mixture estimates, one per threshold in `x_grid`.
#' @export
estimate_all_mixtures <- function(df,
                                  id_col = "id",
                                  val_col = "x",
                                  method = c("spline", "npmle", "raw", "beta"),
                                  x_grid = 1:10,
                                  weight_col = NULL,
                                  ...) {
  method <- match.arg(method)

  df <- .standardize_input_df(
    df,
    id_col = id_col,
    val_col = val_col,
    weight_col = weight_col
  )

  if (!is.null(weight_col) && method %in% c("spline")) {
    stop(
      "`weight_col` is only supported for methods 'npmle', 'raw', and 'beta'.",
      call. = FALSE
    )
  }

  df$id <- factor(df$id)

  res_list <- list()
  prev_k <- NULL
  prev_x_thresh <- NULL

  for (x_thresh in x_grid) {
    k <- as.integer(tapply(df$x <= x_thresh, df$id, sum))

    if (!is.null(prev_k) && identical(k, prev_k)) {
      res_list[[as.character(x_thresh)]] <- res_list[[as.character(prev_x_thresh)]]
    } else {
      mixture_est <- switch(
        method,
        spline = estimate_mixture_efron(
          df = df,
          id_col = "id",
          val_col = "x",
          x_thresh = x_thresh,
          ...
        ),
        npmle = estimate_mixture_npmle(
          df = df,
          id_col = "id",
          val_col = "x",
          x_thresh = x_thresh,
          weight_col = weight_col,
          ...
        ),
        raw = estimate_mixture_raw(
          df = df,
          id_col = "id",
          val_col = "x",
          x_thresh = x_thresh,
          weight_col = weight_col,
          ...
        ),
        beta = estimate_mixture_beta(
          df = df,
          id_col = "id",
          val_col = "x",
          x_thresh = x_thresh,
          weight_col = weight_col,
          ...
        )
      )

      res_list[[as.character(x_thresh)]] <- mixture_est
      prev_k <- k
      prev_x_thresh <- x_thresh
    }
  }

  res_list
}
