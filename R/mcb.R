#' Estimate Mean Alpha Quantiles from Mixture Distributions
#'
#' @param mixture_res A named list of mixture estimates.
#' @param alpha A single quantile level in `[0, 1]`.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_first_last Logical controlling endpoint treatment.
#'
#' @return A numeric scalar giving the estimated mean alpha quantile.
#' @export
est_mean_alpha_quantile <- function(mixture_res,
                                    alpha,
                                    use_midpoint = TRUE,
                                    estimate_first_last = TRUE,
                                    use_isotonic = TRUE) {
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
  
  if (use_isotonic) {
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
    stop(
      "Computed quantile weights are not positive; cannot form mean alpha quantile.",
      call. = FALSE
    )
  }
  
  delta_df <- dplyr::mutate(
    delta_df,
    weighted_delta = .data$delta / delta_sum
  )
  
  # print(delta_df)
  
  sum(delta_df$weighted_delta * delta_df$x_temp)
}


#' Estimate Mean Quantiles Across an Alpha Grid
#'
#' @param mixture_res A named list of mixture estimates.
#' @param alpha_grid Grid of quantile levels.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_first_last Logical controlling endpoint treatment.
#'
#' @return A data frame with columns `quantile` and `estimate`.
#' @export
est_all_quantiles <- function(mixture_res,
                              alpha_grid = seq(0.01, 0.99, by = 0.01),
                              use_midpoint = TRUE,
                              estimate_first_last = TRUE) {
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
        estimate_first_last = estimate_first_last
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
#' @param x_grid Threshold grid used when `cutpoints` is `NULL`.
#' @param cutpoints Optional number of evenly spaced thresholds.
#' @param bootstrap_samples Number of bootstrap resamples.
#' @param alpha_grid Grid of target quantile levels.
#' @param use_midpoint Logical; whether to use interval midpoints between
#'   adjacent thresholds.
#' @param estimate_first_last Logical controlling endpoint treatment.
#' @param weight_col Optional column name containing id-level weights.
#' @param progress Logical; whether to display a bootstrap progress bar.
#' @param ci_level Confidence level.
#' @param ... Additional arguments passed to [estimate_all_mixtures()].
#'
#' @return A data frame with columns `quantile`, `estimate`, `se`, `ci_lo`,
#'   and `ci_hi`.
#' @export
mcb <- function(df,
                id_col = "id",
                val_col = "x",
                method = c("spline", "npmle", "raw", "beta"),
                x_grid = 1:10,
                cutpoints = NULL,
                bootstrap_samples = 100,
                alpha_grid = seq(0.01, 0.99, by = 0.01),
                use_midpoint = TRUE,
                estimate_first_last = TRUE,
                weight_col = NULL,
                progress = TRUE,
                ci_level = 0.95,
                ...) {
  method <- match.arg(method)
  
  df <- .standardize_input_df(
    df,
    id_col = id_col,
    val_col = val_col,
    weight_col = weight_col
  )
  
  if (!is.numeric(bootstrap_samples) || length(bootstrap_samples) != 1 ||
      is.na(bootstrap_samples) || bootstrap_samples < 2) {
    stop("`bootstrap_samples` must be a single integer >= 2.", call. = FALSE)
  }
  
  bootstrap_samples <- as.integer(bootstrap_samples)
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 || is.na(ci_level) ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single numeric value in (0, 1).",
         call. = FALSE)
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
  
  overall_mixture_res <- estimate_all_mixtures(
    df = df,
    id_col = "id",
    val_col = "x",
    method = method,
    x_grid = x_grid_overall,
    weight_col = weight_col,
    ...
  )
  
  overall_quantile_res <- est_all_quantiles(
    mixture_res = overall_mixture_res,
    alpha_grid = alpha_grid,
    use_midpoint = use_midpoint,
    estimate_first_last = estimate_first_last
  )
  
  ids <- unique(df$id)
  by_id <- split(df, df$id)
  
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
    
    boot_quantile_res <- est_all_quantiles(
      mixture_res = boot_mixture_res,
      alpha_grid = alpha_grid,
      use_midpoint = use_midpoint,
      estimate_first_last = estimate_first_last
    )
    
    res_matrix[i, ] <- boot_quantile_res$estimate
    
    if (progress) {
      utils::setTxtProgressBar(pb, i)
      flush.console()
    }
  }
  
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
  
  list(
    res = res,
    cov = cov_mat
  )
}