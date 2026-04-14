#' Estimate a Mixture via Efron-Style Deconvolution
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param x_thresh Numeric threshold used to define successes.
#' @param tau Support grid for the estimated mixture.
#' @param pDegree Degrees of freedom for the spline basis in [deconv_fast()].
#' @param c0 Penalty constant passed to [deconv_fast()].
#' @param weight_col Optional column name containing id-level weights.
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

  .estimate_mixture_efron_from_binomial(
    df_bin = df_bin,
    tau = tau,
    pDegree = pDegree,
    c0 = c0
  )
}

#' Estimate a Mixture via NPMLE
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param x_thresh Numeric threshold used to define successes.
#' @param tau Support grid for the estimated mixture.
#' @param weight_col Optional column name containing id-level weights.
#' @param backend Estimation backend: `"mixsqp"` (default) or `"REBayes"`.
#' @param mixsqp_control Optional control list passed to [mixsqp::mixsqp()].
#' @param rebayes_control Optional named list of additional arguments passed to
#'   [REBayes::Bmix()] when `backend = "REBayes"`.
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
                                   backend = c("mixsqp", "REBayes"),
                                   mixsqp_control = list(),
                                   rebayes_control = list()) {
  backend <- match.arg(backend)
  
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

  .estimate_mixture_npmle_from_binomial(
    df_bin = df_bin,
    tau = tau,
    backend = backend,
    mixsqp_control = mixsqp_control,
    rebayes_control = rebayes_control
  )
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

  .estimate_mixture_raw_from_binomial(df_bin = df_bin)
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

  .estimate_mixture_beta_from_binomial(df_bin = df_bin, tau = tau)
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

  precomputed <- .precompute_binomial_grid(
    df = df,
    x_grid = x_grid,
    weight_col = weight_col
  )

  .estimate_all_mixtures_from_precomputed(
    precomputed = precomputed,
    method = method,
    ...
  )
}

.estimate_mixture_efron_from_binomial <- function(df_bin,
                                                  tau = seq(from = 0.005, to = 0.995, by = 0.005),
                                                  pDegree = 7,
                                                  c0 = 1,
                                                  aStart = 1.0,
                                                  deconv_cache = NULL,
                                                  compute_stats = TRUE) {
  x_input <- dplyr::select(df_bin, n, s)

  if (is.null(deconv_cache)) {
    result <- deconv_fast(
      tau = tau,
      X = x_input,
      family = "Binomial",
      pDegree = pDegree,
      c0 = c0,
      aStart = aStart,
      compute_stats = compute_stats,
      obs_weights = if ("weight" %in% names(df_bin)) df_bin$weight else NULL
    )
  } else {
    logP <- .build_binomial_logP_matrix(
      n = df_bin$n,
      s = df_bin$s,
      tau = deconv_cache$tau
    )

    result <- deconv_fast(
      tau = deconv_cache$tau,
      family = "Binomial",
      y = 1,
      Q = deconv_cache$Q,
      P = logP,
      c0 = c0,
      aStart = aStart,
      compute_stats = compute_stats,
      obs_weights = if ("weight" %in% names(df_bin)) df_bin$weight else NULL
    )
  }

  stats_df <- data.frame(result$stats)
  out <- .standardize_mixture_df(theta = stats_df$theta, g = stats_df$g)
  attr(out, "mixture_type") <- "discrete"
  attr(out, "deconv_mle") <- result$mle

  out
}

.estimate_mixture_npmle_from_binomial <- function(df_bin,
                                                  tau = seq(from = 0, to = 1, by = 0.005),
                                                  backend = c("mixsqp", "REBayes"),
                                                  start = NULL,
                                                  mixsqp_control = list(),
                                                  rebayes_control = list()) {
  backend <- match.arg(backend)

  fit_g <- switch(
    backend,

    mixsqp = {
      if (is.null(mixsqp_control$verbose)) {
        mixsqp_control$verbose <- FALSE
      }

      collapsed <- .collapse_binomial_rows_for_mixsqp(df_bin)
      logL <- .build_binomial_logP_matrix(
        n = collapsed$n,
        s = collapsed$s,
        tau = tau
      )
      valid_cols <- colSums(is.finite(logL)) > 0

      if (!any(valid_cols)) {
        stop(
          paste(
            "No support points in `tau` yield positive likelihood for the observed binomial data.",
            "Try using a different `tau` grid."
          ),
          call. = FALSE
        )
      }

      if (sum(valid_cols) == 1L) {
        fit_g <- numeric(length(tau))
        fit_g[valid_cols] <- 1
        return(fit_g)
      }

      logL <- logL[, valid_cols, drop = FALSE]

      fit_args <- list(
        L = logL,
        w = collapsed$weight,
        log = TRUE,
        control = mixsqp_control
      )

      if (!is.null(start)) {
        fit_args$x0 <- .subset_mixsqp_start(start, valid_cols)
      }

      fit <- do.call(mixsqp::mixsqp, fit_args)
      fit_g <- numeric(length(tau))
      fit_g[valid_cols] <- fit$x
      fit_g
    },

    REBayes = {
      if (!requireNamespace("REBayes", quietly = TRUE)) {
        stop(
          paste(
            "The `REBayes` backend requires the optional package `REBayes`.",
            "Install it first, or use `backend = \"mixsqp\"`."
          ),
          call. = FALSE
        )
      }

      if (!requireNamespace("Rmosek", quietly = TRUE)) {
        stop(
          paste(
            "The `REBayes` backend requires a working MOSEK installation via `Rmosek`.",
            "Install/configure MOSEK and `Rmosek`, or use `backend = \"mixsqp\"`."
          ),
          call. = FALSE
        )
      }

      rebayes_args <- c(
        list(
          x = df_bin$s,
          k = df_bin$n,
          v = tau,
          collapse = TRUE
        ),
        if ("weight" %in% names(df_bin)) {
          list(weights = df_bin$weight / sum(df_bin$weight))
        } else {
          list()
        },
        rebayes_control
      )

      fit <- tryCatch(
        do.call(REBayes::Bmix, rebayes_args),
        error = function(e) {
          stop(
            paste(
              "The `REBayes` backend failed.",
              "This usually means MOSEK / `Rmosek` is not installed or configured correctly.",
              "Original error:", conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )

      fit$y
    }
  )

  out <- .standardize_mixture_df(theta = tau, g = fit_g)
  attr(out, "mixture_type") <- "discrete"
  out
}

.estimate_mixture_raw_from_binomial <- function(df_bin) {
  df_bin <- dplyr::mutate(df_bin, q = .data$s / .data$n)

  if ("weight" %in% names(df_bin)) {
    res_raw <- dplyr::group_by(df_bin, q)
    res_raw <- dplyr::summarise(
      res_raw,
      total = sum(.data$weight),
      .groups = "drop"
    )
  } else {
    res_raw <- dplyr::group_by(df_bin, q)
    res_raw <- dplyr::summarise(res_raw, total = dplyr::n(), .groups = "drop")
  }

  res_raw <- dplyr::transmute(
    res_raw,
    theta = .data$q,
    g = .data$total / sum(.data$total)
  )

  if (!0 %in% res_raw$theta) {
    res_raw <- dplyr::bind_rows(data.frame(theta = 0, g = 0), res_raw)
  }

  out <- .standardize_mixture_df(theta = res_raw$theta, g = res_raw$g)
  attr(out, "mixture_type") <- "discrete"
  out
}

.estimate_mixture_beta_from_binomial <- function(df_bin,
                                                 tau = seq(0.005, 0.995, by = 0.005),
                                                 start = NULL,
                                                 fallback_starts = NULL,
                                                 maxit = 200) {
  all_success <- all(df_bin$s == df_bin$n)
  all_failure <- all(df_bin$s == 0)

  if (all_success || all_failure) {
    g <- rep(0, length(tau))

    if (all_success) {
      idx <- which.min(abs(tau - 1))
    } else {
      idx <- which.min(abs(tau - 0))
    }

    g[idx] <- 1

    out <- .standardize_mixture_df(theta = tau, g = g)
    attr(out, "mixture_type") <- "beta"
    attr(out, "beta_mle") <- list(alpha = NA_real_, beta = NA_real_)
    return(out)
  }

  fit_info <- .fit_betabinomial_optim(
    df_bin,
    weight_col = if ("weight" %in% names(df_bin)) "weight" else NULL,
    start = start,
    fallback_starts = fallback_starts,
    maxit = maxit
  )

  out <- .discretize_beta_on_tau(
    alpha = fit_info$alpha,
    beta = fit_info$beta,
    tau = tau
  )

  attr(out, "mixture_type") <- "beta"
  attr(out, "beta_mle") <- list(
    alpha = fit_info$alpha,
    beta = fit_info$beta
  )

  out
}

.estimate_all_mixtures_from_precomputed <- function(precomputed,
                                                    method = c("spline", "npmle", "raw", "beta"),
                                                    ...) {
  method <- match.arg(method)

  res_list <- list()
  prev_k <- NULL
  prev_name <- NULL
  prev_beta_start <- NULL
  prev_npmle_start <- NULL
  prev_spline_start <- NULL
  spline_cache <- NULL
  spline_disable_warm_start <- FALSE

  if (identical(method, "spline")) {
    args <- list(...)
    tau <- if ("tau" %in% names(args)) args$tau else seq(from = 0.005, to = 0.995, by = 0.005)
    pDegree <- if ("pDegree" %in% names(args)) args$pDegree else 7
    spline_c0 <- if ("c0" %in% names(args)) args$c0 else 1
    spline_disable_warm_start <- is.numeric(spline_c0) &&
      length(spline_c0) == 1L &&
      !is.na(spline_c0) &&
      spline_c0 == 0
    spline_cache <- .make_spline_deconv_cache(
      n = precomputed$n,
      tau = tau,
      pDegree = pDegree
    )
  }

  for (idx in seq_along(precomputed$x_grid)) {
    x_thresh <- precomputed$x_grid[idx]
    k <- precomputed$s[, idx]
    name <- as.character(x_thresh)

    if (!is.null(prev_k) && identical(k, prev_k)) {
      res_list[[name]] <- res_list[[prev_name]]
      next
    }

    df_bin <- .binomial_df_from_precomputed(precomputed, idx)

    mixture_est <- if (identical(method, "beta")) {
      .estimate_mixture_beta_from_binomial(
        df_bin = df_bin,
        start = prev_beta_start,
        ...
      )
    } else if (identical(method, "spline")) {
      .estimate_mixture_efron_from_binomial(
        df_bin = df_bin,
        aStart = if (spline_disable_warm_start || is.null(prev_spline_start)) {
          1.0
        } else {
          prev_spline_start
        },
        deconv_cache = spline_cache,
        ...
      )
    } else {
      switch(
        method,
        npmle = .estimate_mixture_npmle_from_binomial(
          df_bin = df_bin,
          start = prev_npmle_start,
          ...
        ),
        raw = .estimate_mixture_raw_from_binomial(df_bin = df_bin)
      )
    }

    res_list[[name]] <- mixture_est
    prev_k <- k
    prev_name <- name

    if (identical(method, "beta")) {
      beta_mle <- attr(mixture_est, "beta_mle")
      if (!is.null(beta_mle) &&
          is.finite(beta_mle$alpha) &&
          is.finite(beta_mle$beta) &&
          beta_mle$alpha > 0 &&
          beta_mle$beta > 0) {
        prev_beta_start <- c(beta_mle$alpha, beta_mle$beta)
      } else {
        prev_beta_start <- NULL
      }
    } else if (identical(method, "spline")) {
      spline_mle <- attr(mixture_est, "deconv_mle")
      if (!spline_disable_warm_start &&
          !is.null(spline_mle) &&
          all(is.finite(spline_mle))) {
        prev_spline_start <- spline_mle
      } else {
        prev_spline_start <- NULL
      }
    } else if (identical(method, "npmle")) {
      if (all(is.finite(mixture_est$g)) && length(mixture_est$g) > 0) {
        prev_npmle_start <- mixture_est$g
      } else {
        prev_npmle_start <- NULL
      }
    }
  }

  res_list
}

.estimate_all_mixtures_for_bootstrap <- function(precomputed,
                                                 method = c("spline", "npmle", "raw", "beta"),
                                                 ...) {
  method <- match.arg(method)

  if (!identical(method, "spline")) {
    return(.estimate_all_mixtures_from_precomputed(
      precomputed = precomputed,
      method = method,
      ...
    ))
  }

  .estimate_all_mixtures_from_precomputed(
    precomputed = precomputed,
    method = method,
    compute_stats = FALSE,
    ...
  )
}

.build_binomial_logP_matrix <- function(n, s, tau) {
  n <- as.numeric(n)
  s <- as.numeric(s)
  tau <- as.numeric(tau)

  logP <- vapply(
    tau,
    FUN = function(prob) {
      stats::dbinom(s, size = n, prob = prob, log = TRUE)
    },
    FUN.VALUE = numeric(length(n))
  )

  matrix(logP, nrow = length(n), ncol = length(tau))
}

.collapse_binomial_rows_for_mixsqp <- function(df_bin) {
  has_weight <- "weight" %in% names(df_bin)
  key <- paste(df_bin$n, df_bin$s, sep = "\r")
  key_levels <- unique(key)
  group_idx <- match(key, key_levels)
  row_idx <- match(key_levels, key)

  out <- data.frame(
    n = df_bin$n[row_idx],
    s = df_bin$s[row_idx],
    weight = if (has_weight) {
      as.numeric(rowsum(df_bin$weight, group = group_idx, reorder = FALSE)[, 1])
    } else {
      as.numeric(tabulate(group_idx, nbins = length(key_levels)))
    }
  )

  out
}

.subset_mixsqp_start <- function(start, valid_cols) {
  start <- as.numeric(start)

  if (length(start) != length(valid_cols) || anyNA(start) || any(start < 0)) {
    return(NULL)
  }

  start <- start[valid_cols]

  if (length(start) <= 1L || sum(start) <= 0 || !all(is.finite(start))) {
    return(NULL)
  }

  start
}

.make_spline_deconv_cache <- function(n,
                                      tau = seq(from = 0.005, to = 0.995, by = 0.005),
                                      pDegree = 7,
                                      scale = TRUE) {
  tau <- sort(unique(as.numeric(tau)))
  Q <- splines::ns(tau, pDegree)

  if (scale) {
    Q <- scale(Q, center = TRUE, scale = FALSE)
    Q <- sweep(Q, 2, sqrt(colSums(Q * Q)), "/")
  }

  list(
    n = as.numeric(n),
    tau = tau,
    Q = as.matrix(Q)
  )
}
