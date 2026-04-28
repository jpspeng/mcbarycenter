#' Estimate the Distribution of an Alpha Quantile
#'
#' @param mixture_res A named list of mixture estimates.
#' @param alpha A single quantile level in `[0, 1]`.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_first_last Logical controlling endpoint treatment.
#' @param use_isotonic_dist Logical; whether to enforce monotonicity in the
#'   cumulative mixture curves across threshold values.
#'
#' @return A data frame with columns `x`, `pmf`, and `cdf`.
#' @export
est_dist_alpha <- function(mixture_res,
                           alpha,
                           use_midpoint = TRUE,
                           estimate_first_last = TRUE,
                           use_isotonic_dist = TRUE) {
  if (!is.list(mixture_res) || length(mixture_res) < 2) {
    stop("`mixture_res` must be a list of at least two mixture estimates.",
         call. = FALSE
    )
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha)) {
    stop("`alpha` must be a single non-missing numeric value.", call. = FALSE)
  }
  
  if (alpha < 0 || alpha > 1) {
    stop("`alpha` must lie in [0, 1].", call. = FALSE)
  }
  
  x_vals <- as.numeric(names(mixture_res))
  
  if (anyNA(x_vals)) {
    stop("`mixture_res` must be a named list with numeric names.", call. = FALSE)
  }
  
  ord <- order(x_vals)
  mixture_res <- mixture_res[ord]
  x_vals <- x_vals[ord]
  
  cumul_vals <- vapply(
    mixture_res,
    .lookup_cumul_at_alpha,
    numeric(1),
    alpha = alpha
  )
  
  if (use_isotonic_dist) {
    # enforce cumul_vals to be nonincreasing in x_vals
    iso_fit <- stats::isoreg(x = x_vals, y = -cumul_vals)
    cumul_vals <- -iso_fit$yf
  }
  
  if (estimate_first_last) {
    n_intervals <- length(mixture_res) - 1
  } else {
    if (length(mixture_res) < 3) {
      stop(
        "Need at least 3 mixture estimates when `estimate_first_last = FALSE`.",
        call. = FALSE
      )
    }
    n_intervals <- length(mixture_res) - 2
  }
  
  rows <- vector("list", n_intervals)
  
  for (i in seq_len(n_intervals)) {
    x1 <- x_vals[i]
    x2 <- x_vals[i + 1]
    
    x_temp <- if (use_midpoint) {
      (x1 + x2) / 2
    } else {
      x2
    }
    
    cumul1 <- cumul_vals[i]
    cumul2 <- cumul_vals[i + 1]
    
    rows[[i]] <- data.frame(
      x_temp = x_temp,
      delta = cumul1 - cumul2
    )
  }
  
  delta_df <- do.call(rbind, rows)
  
  if (!estimate_first_last) {
    delta_init <- 1 - cumul_vals[1]
    delta_last <- cumul_vals[length(cumul_vals) - 1]
    
    delta_df <- dplyr::bind_rows(
      delta_df,
      data.frame(x_temp = x_vals[1], delta = delta_init),
      data.frame(x_temp = x_vals[length(x_vals)], delta = delta_last)
    )
  }
  
  # guard against tiny negative values from floating point
  delta_df$delta <- pmax(delta_df$delta, 0)
  
  delta_sum <- sum(delta_df$delta)
  
  if (!is.finite(delta_sum) || delta_sum <= 0) {
    print(delta_df)
    stop(
      "Computed quantile weights are not positive; cannot form mean alpha quantile.",
      call. = FALSE
    )
  }
  
  dist_df <- dplyr::transmute(
    delta_df,
    x = .data$x_temp,
    pmf = .data$delta / delta_sum
  )
  dist_df <- dplyr::arrange(dist_df, .data$x)
  dist_df$cdf <- cumsum(dist_df$pmf)
  
  dist_df
}

#' Estimate Mean Alpha Quantiles from Mixture Distributions
#'
#' @param mixture_res A named list of mixture estimates.
#' @param alpha A single quantile level in `[0, 1]`.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_first_last Logical controlling endpoint treatment.
#' @param use_isotonic_dist Logical; whether to enforce monotonicity in the
#'   cumulative mixture curves across threshold values.
#'
#' @return A numeric scalar giving the estimated mean alpha quantile.
#' @export
est_mean_alpha_quantile <- function(mixture_res,
                                    alpha,
                                    use_midpoint = TRUE,
                                    estimate_first_last = TRUE,
                                    use_isotonic_dist = TRUE) {
  dist_df <- est_dist_alpha(
    mixture_res = mixture_res,
    alpha = alpha,
    use_midpoint = use_midpoint,
    estimate_first_last = estimate_first_last,
    use_isotonic_dist = use_isotonic_dist
  )
  
  sum(dist_df$x * dist_df$pmf)
}


#' Estimate Mean Quantiles Across an Alpha Grid
#'
#' @param mixture_res A named list of mixture estimates.
#' @param alpha_grid Grid of quantile levels.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_first_last Logical controlling endpoint treatment.
#' @param use_isotonic_dist Logical; whether to enforce monotonicity in the
#'   cumulative mixture curves across threshold values when constructing each
#'   alpha-quantile distribution.
#'
#' @return A data frame with columns `quantile` and `estimate`.
#' @export
est_all_quantiles <- function(mixture_res,
                              alpha_grid = seq(0.01, 0.99, by = 0.01),
                              use_midpoint = TRUE,
                              estimate_first_last = TRUE,
                              use_isotonic_dist = TRUE) {
  if (!is.numeric(alpha_grid) || anyNA(alpha_grid)) {
    stop("`alpha_grid` must be a numeric vector with no missing values.",
      call. = FALSE
    )
  }

  if (any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("`alpha_grid` must lie in [0, 1].", call. = FALSE)
  }

  alpha_grid <- sort(unique(alpha_grid))

  quantile_means <- vapply(
    alpha_grid,
    FUN = function(a) {
      est_mean_alpha_quantile(
        mixture_res = mixture_res,
        alpha = a,
        use_midpoint = use_midpoint,
        estimate_first_last = estimate_first_last,
        use_isotonic_dist = use_isotonic_dist
      )
    },
    FUN.VALUE = numeric(1)
  )

  data.frame(
    quantile = alpha_grid,
    estimate = quantile_means
  )
}

#' Monte Carlo Bootstrap for Quantile Estimation
#'
#' @param df A data frame.
#' @param id_col Column name identifying groups. Defaults to `"id"`.
#' @param val_col Column name containing the observed values. Defaults to `"x"`.
#' @param method Mixture estimation method.
#' @param x_grid Threshold grid used when `cutpoints` is `NULL`. Must be
#'   supplied if `cutpoints` is `NULL`.
#' @param cutpoints Optional number of evenly spaced thresholds. Either
#'   `cutpoints` or `x_grid` must be supplied.
#' @param bootstrap_samples Number of bootstrap resamples. When set to 1, only
#'   the quantile estimate is reported and uncertainty columns are returned as
#'   `NA`.
#' @param alpha_grid Grid of target quantile levels.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_endpoints Logical controlling endpoint treatment.
#' @param use_isotonic_dist Logical; whether to enforce monotonicity in the
#'   cumulative mixture curves across threshold values when constructing each
#'   alpha-quantile distribution.
#' @param use_isotonic_output Logical; if `TRUE`, apply isotonic regression to
#'   `res$estimate`, `res$ci_lo`, and `res$ci_hi`.
#' @param weight_col Optional column name containing id-level weights.
#' @param progress Logical; whether to display a bootstrap progress bar.
#' @param ci_level Confidence level.
#' @param ... Additional arguments passed to [estimate_all_mixtures()].
#'
#' @return A list with components `res`, `cov`, `mixtures`, `method`, and
#'   `data`. `res` is a data frame with columns `quantile`, `estimate`,
#'   `estimate_bs`, `se`, `ci_lo`, `ci_hi`, `pct_ci_lo`, and `pct_ci_hi`.
#'   `data` contains the original `id_col` and `val_col` columns from the
#'   input.
#' @export
mcbary <- function(df,
                id_col = "id",
                val_col = "x",
                method = c("spline", "npmle", "raw", "beta"),
                x_grid = NULL,
                cutpoints = NULL,
                bootstrap_samples = 100,
                alpha_grid = seq(0.01, 0.99, by = 0.01),
                use_midpoint = TRUE,
                estimate_endpoints = TRUE,
                use_isotonic_dist = TRUE,
                use_isotonic_output = FALSE,
                weight_col = NULL,
                progress = TRUE,
                ci_level = 0.95,
                ...) {
  method <- match.arg(method)
  data_out <- data.frame(
    id = df[[id_col]],
    val = df[[val_col]]
  )
  
  df <- .standardize_input_df(
    df,
    id_col = id_col,
    val_col = val_col,
    weight_col = weight_col
  )
  
  if (!is.numeric(bootstrap_samples) || length(bootstrap_samples) != 1 ||
      is.na(bootstrap_samples) || bootstrap_samples < 1) {
    stop("`bootstrap_samples` must be a single integer >= 1.", call. = FALSE)
  }
  
  bootstrap_samples <- as.integer(bootstrap_samples)
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 || is.na(ci_level) ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single numeric value in (0, 1).",
         call. = FALSE)
  }
  
  if (!is.logical(use_isotonic_dist) || length(use_isotonic_dist) != 1L ||
      is.na(use_isotonic_dist)) {
    stop("`use_isotonic_dist` must be a single TRUE/FALSE value.",
         call. = FALSE)
  }

  if (!is.logical(use_isotonic_output) || length(use_isotonic_output) != 1L ||
      is.na(use_isotonic_output)) {
    stop("`use_isotonic_output` must be a single TRUE/FALSE value.",
         call. = FALSE)
  }
  
  if (!is.logical(estimate_endpoints) || length(estimate_endpoints) != 1L ||
      is.na(estimate_endpoints)) {
    stop(
      "`estimate_endpoints` must be a single TRUE/FALSE value.",
      call. = FALSE
    )
  }
  
  if (is.null(cutpoints) && is.null(x_grid)) {
    stop(
      "Either `cutpoints` or `x_grid` must be supplied.",
      call. = FALSE
    )
  }
  
  if (!is.null(cutpoints)) {
    if (!is.numeric(cutpoints) || length(cutpoints) != 1 ||
        is.na(cutpoints) || cutpoints < 2) {
      stop("`cutpoints` must be a single numeric value >= 2.", call. = FALSE)
    }
    
    x_grid_overall <- seq(
      min(df$x),
      max(df$x),
      length.out = cutpoints
    )
  } else {
    x_grid_overall <- x_grid
  }

  use_precomputed_grid <- is.null(cutpoints)

  if (use_precomputed_grid) {
    precomputed_grid <- .precompute_binomial_grid(
      df = df,
      x_grid = x_grid_overall,
      weight_col = weight_col
    )

    overall_mixture_res <- .estimate_all_mixtures_from_precomputed(
      precomputed = precomputed_grid,
      method = method,
      ...
    )
  } else {
    overall_mixture_res <- estimate_all_mixtures(
      df = df,
      id_col = "id",
      val_col = "x",
      method = method,
      x_grid = x_grid_overall,
      weight_col = weight_col,
      ...
    )
  }
  
  overall_quantile_res <- est_all_quantiles(
    mixture_res = overall_mixture_res,
    alpha_grid = alpha_grid,
    use_midpoint = use_midpoint,
    estimate_first_last = estimate_endpoints,
    use_isotonic_dist = use_isotonic_dist
  )
  
  if (use_precomputed_grid) {
    n_ids <- length(precomputed_grid$n)
  } else {
    ids <- unique(df$id)
    by_id <- split(df, df$id)
  }
  
  res_matrix <- matrix(
    NA_real_,
    nrow = bootstrap_samples,
    ncol = length(alpha_grid)
  )
  
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  
  if (progress) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = bootstrap_samples,
      style = 3,
      file = stderr()
    )
    on.exit(close(pb), add = TRUE)
  }
  
  for (i in seq_len(bootstrap_samples)) {
    if (use_precomputed_grid) {
      sampled_idx <- sample.int(n_ids, size = n_ids, replace = TRUE)
      boot_precomputed <- .resample_precomputed_binomial_grid(
        precomputed = precomputed_grid,
        sampled_idx = sampled_idx
      )

      boot_mixture_res <- .estimate_all_mixtures_for_bootstrap(
        precomputed = boot_precomputed,
        method = method,
        ...
      )
    } else {
      sampled_ids <- sample(ids, size = length(ids), replace = TRUE)
      picked <- by_id[as.character(sampled_ids)]

      picked <- setNames(picked, seq_along(picked))

      df_boot <- dplyr::bind_rows(picked, .id = "boot_id")
      df_boot <- dplyr::mutate(df_boot, id = .data$boot_id)
      df_boot <- dplyr::select(df_boot, -boot_id)

      if (!is.null(cutpoints)) {
        x_grid_boot <- seq(
          min(df_boot$x),
          max(df_boot$x),
          length.out = cutpoints
        )
      } else {
        x_grid_boot <- x_grid
      }

      boot_mixture_res <- estimate_all_mixtures(
        df = df_boot,
        id_col = "id",
        val_col = "x",
        method = method,
        x_grid = x_grid_boot,
        weight_col = weight_col,
        ...
      )
    }
    
    boot_quantile_res <- est_all_quantiles(
      mixture_res = boot_mixture_res,
      alpha_grid = alpha_grid,
      use_midpoint = use_midpoint,
      estimate_first_last = estimate_endpoints,
      use_isotonic_dist = use_isotonic_dist
    )
    
    res_matrix[i, ] <- boot_quantile_res$estimate
    
    if (progress) {
      utils::setTxtProgressBar(pb, i)
      flush.console()
    }
  }
  
  if (bootstrap_samples == 1L) {
    estimate_bs <- rep(NA_real_, length(alpha_grid))
    se <- rep(NA_real_, length(alpha_grid))
    cov_mat <- matrix(NA_real_, nrow = length(alpha_grid), ncol = length(alpha_grid))
    pct_ci_lo <- rep(NA_real_, length(alpha_grid))
    pct_ci_hi <- rep(NA_real_, length(alpha_grid))
    ci_lo <- rep(NA_real_, length(alpha_grid))
    ci_hi <- rep(NA_real_, length(alpha_grid))
  } else {
    # Bootstrap summaries
    estimate_bs <- colMeans(res_matrix, na.rm = TRUE)
    se <- apply(res_matrix, 2, stats::sd, na.rm = TRUE)
    cov_mat <- stats::cov(res_matrix, use = "pairwise.complete.obs")
    
    # CI setup
    alpha_ci <- 1 - ci_level
    lo_prob <- alpha_ci / 2
    hi_prob <- 1 - alpha_ci / 2
    
    # Percentile CIs
    pct_ci_lo <- apply(
      res_matrix,
      2,
      stats::quantile,
      probs = lo_prob,
      na.rm = TRUE,
      names = FALSE
    )
    
    pct_ci_hi <- apply(
      res_matrix,
      2,
      stats::quantile,
      probs = hi_prob,
      na.rm = TRUE,
      names = FALSE
    )
    
    # Wald CIs
    z <- stats::qnorm(1 - alpha_ci / 2)
    ci_lo <- overall_quantile_res$estimate - z * se
    ci_hi <- overall_quantile_res$estimate + z * se
  }
  
  res <- data.frame(
    quantile = alpha_grid,
    estimate = overall_quantile_res$estimate,
    estimate_bs = estimate_bs,
    se = se,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    pct_ci_lo = pct_ci_lo,
    pct_ci_hi = pct_ci_hi
  )
  
  if (use_isotonic_output) {
    res$estimate <- stats::isoreg(x = res$quantile, y = res$estimate)$yf
    res$ci_lo <- stats::isoreg(x = res$quantile, y = res$ci_lo)$yf
    res$ci_hi <- stats::isoreg(x = res$quantile, y = res$ci_hi)$yf
  }

  if (any(diff(res$estimate) < 0, na.rm = TRUE)) {
    warning(
      "The estimated barycenter's quantile function is not monotone increasing.",
      call. = FALSE
    )
  }
  
  list(
    res = res,
    cov = cov_mat,
    mixtures = overall_mixture_res,
    method = method,
    data = data_out
  )
}
