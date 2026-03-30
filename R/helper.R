.validate_input_df <- function(df,
                               required_cols = c("id", "x"),
                               weight_col = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(weight_col) && !weight_col %in% names(df)) {
    stop(
      sprintf("`weight_col` = '%s' not found in `df`.", weight_col),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.standardize_input_df <- function(df,
                                  id_col = "id",
                                  val_col = "x",
                                  weight_col = NULL) {
  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("`id_col` must be a single non-empty string.", call. = FALSE)
  }

  if (!is.character(val_col) || length(val_col) != 1L || !nzchar(val_col)) {
    stop("`val_col` must be a single non-empty string.", call. = FALSE)
  }

  .validate_input_df(
    df,
    required_cols = c(id_col, val_col),
    weight_col = weight_col
  )

  out <- data.frame(
    id = df[[id_col]],
    x = df[[val_col]]
  )

  if (!is.null(weight_col)) {
    out$weight <- df[[weight_col]]
  }

  out
}

.validate_id_weights <- function(df, weight_col) {
  if (is.null(weight_col)) {
    return(invisible(TRUE))
  }

  tmp <- dplyr::group_by(df, id)
  tmp <- dplyr::summarise(
    tmp,
    n_weights = dplyr::n_distinct(.data[[weight_col]]),
    .groups = "drop"
  )

  bad_ids <- tmp$id[tmp$n_weights > 1]

  if (length(bad_ids) > 0) {
    stop(
      paste0(
        "Weights must be constant within id. ",
        "Found multiple weights for ", length(bad_ids), " id(s)."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.aggregate_binomial_data <- function(df, x_thresh, weight_col = NULL) {
  .validate_input_df(df, required_cols = c("id", "x"), weight_col = weight_col)
  .validate_id_weights(df, weight_col)

  if (!is.numeric(x_thresh) || length(x_thresh) != 1 || is.na(x_thresh)) {
    stop("`x_thresh` must be a single numeric value.", call. = FALSE)
  }

  out <- dplyr::group_by(df, id)

  if (is.null(weight_col)) {
    out <- dplyr::summarise(
      out,
      n = dplyr::n(),
      s = sum(.data$x <= x_thresh),
      .groups = "drop"
    )
  } else {
    out <- dplyr::summarise(
      out,
      n = dplyr::n(),
      s = sum(.data$x <= x_thresh),
      weight = dplyr::first(.data[[weight_col]]),
      .groups = "drop"
    )
  }

  if (nrow(out) == 0) {
    stop("No rows remain after aggregation.", call. = FALSE)
  }

  if (any(out$n <= 0)) {
    stop("Aggregated trial counts must be positive.", call. = FALSE)
  }

  if (any(out$s < 0 | out$s > out$n)) {
    stop("Aggregated successes must satisfy 0 <= s <= n.", call. = FALSE)
  }

  out
}

.precompute_binomial_grid <- function(df, x_grid, weight_col = NULL) {
  .validate_input_df(df, required_cols = c("id", "x"), weight_col = weight_col)
  .validate_id_weights(df, weight_col)

  x_grid <- sort(unique(as.numeric(x_grid)))

  if (length(x_grid) == 0 || anyNA(x_grid)) {
    stop("`x_grid` must be a non-empty numeric vector with no NA.", call. = FALSE)
  }

  by_id_x <- split(df$x, df$id, drop = TRUE)
  ids <- names(by_id_x)
  n <- as.integer(lengths(by_id_x))

  s_mat <- vapply(
    by_id_x,
    FUN = function(x) {
      findInterval(x_grid, sort.int(x, method = "quick"))
    },
    FUN.VALUE = numeric(length(x_grid))
  )
  s_mat <- matrix(
    s_mat,
    nrow = length(x_grid),
    ncol = length(by_id_x)
  )
  s_mat <- t(s_mat)

  weights <- NULL
  if (!is.null(weight_col)) {
    weight_df <- dplyr::group_by(df, id)
    weight_df <- dplyr::summarise(
      weight_df,
      weight = dplyr::first(.data[[weight_col]]),
      .groups = "drop"
    )
    weights <- weight_df$weight[match(ids, weight_df$id)]
  }

  list(
    x_grid = x_grid,
    ids = ids,
    n = n,
    s = s_mat,
    weight = weights
  )
}

.binomial_df_from_precomputed <- function(precomputed, idx) {
  if (!is.list(precomputed) ||
      !all(c("ids", "n", "s", "x_grid") %in% names(precomputed))) {
    stop(
      "`precomputed` must be a binomial grid object from `.precompute_binomial_grid()`.",
      call. = FALSE
    )
  }

  if (!is.numeric(idx) || length(idx) != 1L || is.na(idx) ||
      idx < 1 || idx > length(precomputed$x_grid)) {
    stop("`idx` must index one column of the precomputed grid.", call. = FALSE)
  }

  out <- data.frame(
    id = precomputed$ids,
    n = precomputed$n,
    s = precomputed$s[, idx],
    row.names = NULL
  )

  if (!is.null(precomputed$weight)) {
    out$weight <- precomputed$weight
  }

  out
}

.resample_precomputed_binomial_grid <- function(precomputed, sampled_idx) {
  if (!is.list(precomputed) ||
      !all(c("x_grid", "n", "s") %in% names(precomputed))) {
    stop(
      "`precomputed` must be a binomial grid object from `.precompute_binomial_grid()`.",
      call. = FALSE
    )
  }

  sampled_idx <- as.integer(sampled_idx)

  if (length(sampled_idx) == 0 || anyNA(sampled_idx) ||
      any(sampled_idx < 1 | sampled_idx > length(precomputed$n))) {
    stop("`sampled_idx` contains invalid row indices.", call. = FALSE)
  }

  out <- list(
    x_grid = precomputed$x_grid,
    ids = as.character(seq_along(sampled_idx)),
    n = precomputed$n[sampled_idx],
    s = precomputed$s[sampled_idx, , drop = FALSE],
    weight = if (is.null(precomputed$weight)) NULL else precomputed$weight[sampled_idx]
  )

  out
}

.standardize_mixture_df <- function(theta,
                                    g,
                                    sort_theta = TRUE,
                                    normalize = TRUE,
                                    add_cumul = TRUE) {
  theta <- as.numeric(theta)
  g <- as.numeric(g)

  if (length(theta) != length(g)) {
    stop("`theta` and `g` must have same length.", call. = FALSE)
  }

  if (anyNA(theta) || anyNA(g)) {
    stop("`theta` and `g` cannot contain NA.", call. = FALSE)
  }

  out <- data.frame(theta = theta, g = g)

  if (sort_theta) {
    out <- out[order(out$theta), , drop = FALSE]
  }

  if (normalize) {
    s <- sum(out$g)

    if (!is.finite(s) || s <= 0) {
      stop("Mixture weights must sum to a positive finite value.", call. = FALSE)
    }

    out$g <- out$g / s
  }

  if (any(out$g < -0.01)) {
    stop("Mixture weights contain materially negative values.", call. = FALSE)
  }

  out$g[out$g < 0] <- 0

  if (normalize) {
    out$g <- out$g / sum(out$g)
  }

  if (add_cumul) {
    out$cumul <- cumsum(out$g)
  }

  rownames(out) <- NULL
  out
}

.discretize_beta_on_tau <- function(alpha, beta, tau) {
  tau <- sort(unique(as.numeric(tau)))

  if (length(tau) == 0 || anyNA(tau)) {
    stop("`tau` must be a non-empty numeric vector with no NA.", call. = FALSE)
  }

  if (any(tau < 0 | tau > 1)) {
    stop("`tau` must lie in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(beta) || length(beta) != 1 || is.na(beta) || beta <= 0) {
    stop("`beta` must be a single positive number.", call. = FALSE)
  }

  n_tau <- length(tau)
  edges <- numeric(n_tau + 1)

  if (n_tau >= 2) {
    edges[2:n_tau] <- (tau[-n_tau] + tau[-1]) / 2
  }

  edges[1] <- 0
  edges[n_tau + 1] <- 1
  edges <- pmin(pmax(edges, 0), 1)
  edges <- cummax(edges)

  g <- stats::pbeta(edges[-1], shape1 = alpha, shape2 = beta) -
    stats::pbeta(edges[-length(edges)], shape1 = alpha, shape2 = beta)

  g[g < 0 & g > -1e-14] <- 0
  g[g < 0] <- 0

  .standardize_mixture_df(theta = tau, g = g)
}


.beta_mom_start <- function(df_bin, weight_col = NULL) {
  if (!all(c("n", "s") %in% names(df_bin))) {
    stop("`df_bin` must contain columns `n` and `s`.", call. = FALSE)
  }

  w <- if (is.null(weight_col)) {
    rep(1, nrow(df_bin))
  } else {
    as.numeric(df_bin[[weight_col]])
  }

  p_hat <- (df_bin$s + 0.5) / (df_bin$n + 1)
  w_sum <- sum(w)
  mu <- sum(w * p_hat) / w_sum
  var_p <- sum(w * (p_hat - mu)^2) / w_sum
  max_var <- max(mu * (1 - mu) - 1e-8, 1e-8)
  var_p <- min(max(var_p, 1e-8), max_var)

  phi <- mu * (1 - mu) / var_p - 1
  phi <- max(phi, 1e-4)

  c(
    alpha = max(mu * phi, 1e-4),
    beta = max((1 - mu) * phi, 1e-4)
  )
}

.fit_betabinomial_optim <- function(
    df_bin,
    weight_col = NULL,
    start = NULL,
    fallback_starts = NULL,
    maxit = 200,
    trace = FALSE
) {
  if (!all(c("n", "s") %in% names(df_bin))) {
    stop("`df_bin` must contain columns `n` and `s`.", call. = FALSE)
  }
  
  if (anyNA(df_bin$n) || anyNA(df_bin$s)) {
    stop("`df_bin$n` and `df_bin$s` must not contain NA values.", call. = FALSE)
  }
  
  if (any(df_bin$n < 0) || any(df_bin$s < 0) || any(df_bin$s > df_bin$n)) {
    stop("Require 0 <= s <= n for all rows.", call. = FALSE)
  }
  
  if (!is.null(weight_col)) {
    if (!weight_col %in% names(df_bin)) {
      stop(
        sprintf("`weight_col` = '%s' not found in `df_bin`.", weight_col),
        call. = FALSE
      )
    }
    
    weights <- as.numeric(df_bin[[weight_col]])
    
    if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
      stop(
        "Weights supplied to `.fit_betabinomial_vgam()` must be finite and non-negative.",
        call. = FALSE
      )
    }
  } else {
    weights <- NULL
  }

  if (is.null(start)) {
    start <- .beta_mom_start(
      df_bin,
      weight_col = if (is.null(weight_col)) NULL else weight_col
    )
  }

  start <- as.numeric(start)
  if (length(start) != 2L || anyNA(start) || any(!is.finite(start)) ||
      any(start <= 0)) {
    stop("`start` must contain two finite positive values.", call. = FALSE)
  }

  if (is.null(fallback_starts)) {
    mom <- .beta_mom_start(
      df_bin,
      weight_col = if (is.null(weight_col)) NULL else weight_col
    )
    fallback_starts <- rbind(
      mom,
      c(1, 1),
      c(0.5, 0.5),
      c(2, 2),
      c(5, 5),
      c(1, 5),
      c(5, 1)
    )
  }

  fallback_starts <- rbind(start, fallback_starts)
  fallback_starts <- unique(round(fallback_starts, 12))

  n_vec <- as.numeric(df_bin$n)
  s_vec <- as.numeric(df_bin$s)
  w_vec <- if (is.null(weights)) rep(1, length(n_vec)) else weights

  negloglik <- function(eta) {
    alpha <- exp(eta[1])
    beta <- exp(eta[2])

    ll_i <- lchoose(n_vec, s_vec) +
      lbeta(s_vec + alpha, n_vec - s_vec + beta) -
      lbeta(alpha, beta)

    -sum(w_vec * ll_i)
  }

  fit <- NULL
  best_value <- Inf
  attempts <- vector("list", nrow(fallback_starts))

  for (i in seq_len(nrow(fallback_starts))) {
    start_i <- as.numeric(fallback_starts[i, ])
    eta0 <- log(start_i)

    if (isTRUE(trace)) {
      message(
        sprintf(
          "Trying optim() beta-binomial fit with alpha = %g, beta = %g",
          start_i[1], start_i[2]
        )
      )
    }

    res <- tryCatch(
      stats::optim(
        par = eta0,
        fn = negloglik,
        method = "BFGS",
        control = list(maxit = maxit, trace = if (isTRUE(trace)) 1 else 0)
      ),
      error = function(e) e
    )

    if (inherits(res, "error")) {
      attempts[[i]] <- list(
        alpha_start = start_i[1],
        beta_start = start_i[2],
        ok = FALSE,
        error = conditionMessage(res)
      )
      next
    }

    eta_hat <- res$par
    alpha_hat <- exp(eta_hat[1])
    beta_hat <- exp(eta_hat[2])
    ok <- is.finite(res$value) && is.finite(alpha_hat) && is.finite(beta_hat)

    attempts[[i]] <- list(
      alpha_start = start_i[1],
      beta_start = start_i[2],
      ok = ok,
      convergence = res$convergence,
      value = res$value,
      message = res$message
    )

    if (ok && res$value < best_value) {
      best_value <- res$value
      fit <- res
    }

    if (ok && identical(res$convergence, 0L)) {
      fit <- res
      break
    }
  }

  if (is.null(fit)) {
    stop(
      paste0(
        "All optim() beta-binomial fits failed. Tried starting values: ",
        paste(
          sprintf("(alpha=%g, beta=%g)", fallback_starts[, 1], fallback_starts[, 2]),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  alpha <- exp(fit$par[1])
  beta <- exp(fit$par[2])

  list(
    fit = fit,
    alpha = alpha,
    beta = beta,
    attempts = attempts
  )
}

.lookup_cumul_at_alpha <- function(mixture_df, alpha) {
  if (!is.data.frame(mixture_df)) {
    stop("Each mixture estimate must be a data.frame.", call. = FALSE)
  }

  if (!all(c("theta", "g") %in% names(mixture_df))) {
    stop(
      "Each mixture estimate must contain columns `theta` and `g`.",
      call. = FALSE
    )
  }

  mix <- mixture_df
  if (!"cumul" %in% names(mix)) {
    mix <- .standardize_mixture_df(theta = mix$theta, g = mix$g)
  } else {
    mix <- mix[order(mix$theta), , drop = FALSE]
  }

  idx <- which(mix$theta <= alpha)
  if (length(idx) == 0) {
    return(0)
  }

  mix$cumul[max(idx)]
}

#' Compute Difference in Quantile Estimates
#'
#' Takes two quantile result objects and returns the estimated difference along
#' with delta-method standard errors and confidence intervals, assuming
#' independence.
#'
#' @param df1 A data frame containing columns `quantile`, `estimate`, and `se`,
#'   or a list with a `res` data frame in that format.
#' @param df2 A data frame containing columns `quantile`, `estimate`, and `se`,
#'   or a list with a `res` data frame in that format.
#' @param conf_level Confidence level for the interval. Default is 0.95.
#'
#' @return A data frame with columns:
#'   `quantile`, `estimate`, `se`, `ci_lo`, and `ci_hi`.
#' @export
quantile_diff <- function(df1, df2, conf_level = 0.95) {
  required_cols <- c("quantile", "estimate", "se")
  
  as_quantile_df <- function(x, arg_label, required_cols) {
    if (is.list(x) && !is.data.frame(x) && "res" %in% names(x)) {
      x <- x$res
    }
    
    if (!is.data.frame(x)) {
      stop(sprintf("%s must be a data frame or a list with `$res`.", arg_label),
        call. = FALSE
      )
    }
    
    if (!all(required_cols %in% names(x))) {
      stop(
        sprintf(
          "%s must contain `%s`.",
          arg_label,
          paste(required_cols, collapse = "`, `")
        ),
        call. = FALSE
      )
    }
    
    x
  }
  
  df1 <- as_quantile_df(df1, "`df1`", required_cols)
  df2 <- as_quantile_df(df2, "`df2`", required_cols)
  
  if (nrow(df1) != nrow(df2)) {
    stop("`df1` and `df2` must have the same number of rows.", call. = FALSE)
  }
  
  if (!identical(df1$quantile, df2$quantile)) {
    stop("`df1$quantile` and `df2$quantile` must match exactly.", call. = FALSE)
  }
  
  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha / 2)
  
  estimate <- df1$estimate - df2$estimate
  se <- sqrt(df1$se^2 + df2$se^2)
  
  out <- data.frame(
    quantile = df1$quantile,
    estimate = estimate,
    se = se,
    ci_lo = estimate - z * se,
    ci_hi = estimate + z * se
  )
  
  out
}
