.betabinomial_score_hessian <- function(n, s, shape1, shape2) {
  n <- as.numeric(n)
  s <- as.numeric(s)

  if (length(n) == 0L || length(n) != length(s) || anyNA(n) || anyNA(s) ||
      any(!is.finite(n)) || any(!is.finite(s)) || any(n < 0) ||
      any(s < 0 | s > n)) {
    stop("`n` and `s` must be finite, equally sized, and satisfy 0 <= s <= n.",
         call. = FALSE)
  }
  if (!is.numeric(shape1) || length(shape1) != 1L || is.na(shape1) ||
      !is.finite(shape1) || shape1 <= 0 ||
      !is.numeric(shape2) || length(shape2) != 1L || is.na(shape2) ||
      !is.finite(shape2) || shape2 <= 0) {
    stop("`shape1` and `shape2` must be finite positive scalars.", call. = FALSE)
  }

  score_a <- digamma(s + shape1) - digamma(shape1) +
    digamma(shape1 + shape2) - digamma(n + shape1 + shape2)
  score_b <- digamma(n - s + shape2) - digamma(shape2) +
    digamma(shape1 + shape2) - digamma(n + shape1 + shape2)
  score <- cbind(shape1 = score_a, shape2 = score_b)

  hessian_aa <- trigamma(s + shape1) - trigamma(shape1) +
    trigamma(shape1 + shape2) - trigamma(n + shape1 + shape2)
  hessian_bb <- trigamma(n - s + shape2) - trigamma(shape2) +
    trigamma(shape1 + shape2) - trigamma(n + shape1 + shape2)
  hessian_ab <- trigamma(shape1 + shape2) -
    trigamma(n + shape1 + shape2)

  hessian <- array(
    0,
    dim = c(length(n), 2L, 2L),
    dimnames = list(NULL, c("shape1", "shape2"), c("shape1", "shape2"))
  )
  hessian[, 1L, 1L] <- hessian_aa
  hessian[, 2L, 2L] <- hessian_bb
  hessian[, 1L, 2L] <- hessian_ab
  hessian[, 2L, 1L] <- hessian_ab

  average_hessian <- apply(hessian, c(2L, 3L), mean)
  information <- -average_hessian
  information <- (information + t(information)) / 2
  information_eigenvalues <- if (all(is.finite(information))) {
    eigen(information, symmetric = TRUE, only.values = TRUE)$values
  } else {
    rep(NA_real_, 2L)
  }
  minimum_eigenvalue <- min(information_eigenvalues)
  condition_number <- if (all(is.finite(information_eigenvalues)) &&
      minimum_eigenvalue > 0) {
    max(information_eigenvalues) / minimum_eigenvalue
  } else {
    Inf
  }

  list(
    score = score,
    hessian = hessian,
    average_hessian = average_hessian,
    information = information,
    minimum_eigenvalue = minimum_eigenvalue,
    condition_number = condition_number,
    score_norm = sqrt(sum(colMeans(score)^2))
  )
}

.betabinomial_expected_information <- function(n, shape1, shape2) {
  n <- as.numeric(n)
  if (length(n) == 0L || anyNA(n) || any(!is.finite(n)) || any(n < 1) ||
      any(n != floor(n))) {
    stop("`n` must contain positive finite integer cluster sizes.",
         call. = FALSE)
  }
  if (!is.numeric(shape1) || length(shape1) != 1L || is.na(shape1) ||
      !is.finite(shape1) || shape1 <= 0 ||
      !is.numeric(shape2) || length(shape2) != 1L || is.na(shape2) ||
      !is.finite(shape2) || shape2 <= 0) {
    stop("`shape1` and `shape2` must be finite positive scalars.", call. = FALSE)
  }

  size_counts <- table(n)
  cluster_sizes <- as.integer(names(size_counts))
  size_weights <- as.numeric(size_counts) / length(n)
  conditional <- vector("list", length(cluster_sizes))
  information <- matrix(
    0,
    nrow = 2L,
    ncol = 2L,
    dimnames = list(c("shape1", "shape2"), c("shape1", "shape2"))
  )
  expected_score <- c(shape1 = 0, shape2 = 0)

  for (j in seq_along(cluster_sizes)) {
    m <- cluster_sizes[j]
    y <- 0:m
    log_probability <- lchoose(m, y) +
      lbeta(y + shape1, m - y + shape2) -
      lbeta(shape1, shape2)
    max_log_probability <- max(log_probability)
    scaled_probability <- exp(log_probability - max_log_probability)
    log_normalizer <- max_log_probability + log(sum(scaled_probability))
    probability <- exp(log_probability - log_normalizer)

    score <- .betabinomial_score_hessian(
      n = rep(m, length(y)),
      s = y,
      shape1 = shape1,
      shape2 = shape2
    )$score
    conditional_score <- colSums(score * probability)
    conditional_information <- crossprod(score * sqrt(probability))
    conditional_information <-
      (conditional_information + t(conditional_information)) / 2
    conditional_eigenvalues <- eigen(
      conditional_information,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    conditional_minimum_eigenvalue <- min(conditional_eigenvalues)
    conditional_condition_number <- if (conditional_minimum_eigenvalue > 0) {
      max(conditional_eigenvalues) / conditional_minimum_eigenvalue
    } else {
      Inf
    }

    information <- information + size_weights[j] * conditional_information
    expected_score <- expected_score + size_weights[j] * conditional_score
    conditional[[j]] <- list(
      cluster_size = m,
      frequency = as.integer(size_counts[j]),
      weight = size_weights[j],
      information = conditional_information,
      minimum_eigenvalue = conditional_minimum_eigenvalue,
      condition_number = conditional_condition_number,
      expected_score = conditional_score,
      expected_score_norm = sqrt(sum(conditional_score^2)),
      log_probability_normalizer = log_normalizer
    )
  }

  information <- (information + t(information)) / 2
  information_eigenvalues <- if (all(is.finite(information))) {
    eigen(information, symmetric = TRUE, only.values = TRUE)$values
  } else {
    rep(NA_real_, 2L)
  }
  minimum_eigenvalue <- min(information_eigenvalues)
  condition_number <- if (all(is.finite(information_eigenvalues)) &&
      minimum_eigenvalue > 0) {
    max(information_eigenvalues) / minimum_eigenvalue
  } else {
    Inf
  }

  list(
    information = information,
    minimum_eigenvalue = minimum_eigenvalue,
    condition_number = condition_number,
    expected_score = expected_score,
    expected_score_norm = sqrt(sum(expected_score^2)),
    conditional = conditional
  )
}

.beta_cdf_shape_gradient <- function(
    alpha,
    shape1,
    shape2,
    step = .Machine$double.eps^(1 / 3)
) {
  alpha <- as.numeric(alpha)
  if (length(alpha) == 0L || anyNA(alpha) || any(!is.finite(alpha)) ||
      any(alpha < 0 | alpha > 1)) {
    stop("`alpha` must be finite and lie in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(shape1) || length(shape1) != 1L || is.na(shape1) ||
      !is.finite(shape1) || shape1 <= 0 ||
      !is.numeric(shape2) || length(shape2) != 1L || is.na(shape2) ||
      !is.finite(shape2) || shape2 <= 0) {
    stop("`shape1` and `shape2` must be finite positive scalars.", call. = FALSE)
  }
  if (!is.numeric(step) || length(step) != 1L || is.na(step) ||
      !is.finite(step) || step <= 0) {
    stop("`step` must be a finite positive scalar.", call. = FALSE)
  }

  grad_a <- (
    stats::pbeta(alpha, shape1 = shape1 * exp(step), shape2 = shape2) -
      stats::pbeta(alpha, shape1 = shape1 * exp(-step), shape2 = shape2)
  ) / (2 * step * shape1)
  grad_b <- (
    stats::pbeta(alpha, shape1 = shape1, shape2 = shape2 * exp(step)) -
      stats::pbeta(alpha, shape1 = shape1, shape2 = shape2 * exp(-step))
  ) / (2 * step * shape2)

  cbind(shape1 = grad_a, shape2 = grad_b)
}

# Slower validation implementation obtained by differentiating the incomplete
# beta integral. Production inference uses `.beta_cdf_shape_gradient()`.
.beta_cdf_shape_gradient_integral <- function(alpha, shape1, shape2) {
  alpha <- as.numeric(alpha)
  if (any(alpha <= 0 | alpha >= 1)) {
    stop("The integral validation gradient requires `alpha` in (0, 1).",
         call. = FALSE)
  }

  one_gradient <- function(x) {
    cdf <- stats::pbeta(x, shape1 = shape1, shape2 = shape2)
    density_log_moment <- function(log_term) {
      integrand <- function(u) {
        theta <- x * u
        value <- x * stats::dbeta(theta, shape1 = shape1, shape2 = shape2) *
          log_term(theta)
        value[!is.finite(value)] <- 0
        value
      }
      stats::integrate(
        integrand,
        lower = 0,
        upper = 1,
        rel.tol = 1e-9,
        subdivisions = 1000L
      )$value
    }

    da <- density_log_moment(log) -
      cdf * (digamma(shape1) - digamma(shape1 + shape2))
    db <- density_log_moment(function(theta) log1p(-theta)) -
      cdf * (digamma(shape2) - digamma(shape1 + shape2))
    c(shape1 = da, shape2 = db)
  }

  out <- t(vapply(alpha, one_gradient, numeric(2L)))
  colnames(out) <- c("shape1", "shape2")
  out
}

.mcb_grid_jacobian <- function(
    cumul_vals,
    x_vals,
    use_midpoint,
    estimate_endpoints
) {
  cumul_vals <- as.numeric(cumul_vals)
  x_vals <- as.numeric(x_vals)
  if (length(cumul_vals) != length(x_vals) || length(x_vals) < 2L ||
      anyNA(cumul_vals) || anyNA(x_vals) ||
      any(!is.finite(cumul_vals)) || any(!is.finite(x_vals))) {
    stop("`cumul_vals` and `x_vals` must be equally sized finite vectors with at least two values.",
         call. = FALSE)
  }
  if (is.unsorted(x_vals, strictly = TRUE)) {
    stop("`x_vals` must be strictly increasing.", call. = FALSE)
  }
  if (!is.logical(use_midpoint) || length(use_midpoint) != 1L ||
      is.na(use_midpoint) ||
      !is.logical(estimate_endpoints) || length(estimate_endpoints) != 1L ||
      is.na(estimate_endpoints)) {
    stop("`use_midpoint` and `estimate_endpoints` must be single logical values.",
         call. = FALSE)
  }

  k <- length(x_vals)
  z <- if (use_midpoint) {
    (x_vals[-k] + x_vals[-1L]) / 2
  } else {
    x_vals[-1L]
  }

  if (estimate_endpoints) {
    denominator <- cumul_vals[1L] - cumul_vals[k]
    denominator_tol <- sqrt(.Machine$double.eps) *
      max(1, max(abs(cumul_vals)))
    if (!is.finite(denominator) || denominator <= denominator_tol) {
      stop(
        "The endpoint-estimated grid has a nonpositive or numerically zero total mass (`cumul_vals[1] - cumul_vals[K]`).",
        call. = FALSE
      )
    }

    gradient_numerator <- c(z[1L], diff(z), -z[length(z)])
    gradient_denominator <- c(1, rep(0, k - 2L), -1)
    numerator <- sum(z * (cumul_vals[-k] - cumul_vals[-1L]))
    estimate <- numerator / denominator
    gradient <- (gradient_numerator - estimate * gradient_denominator) /
      denominator
  } else {
    if (k < 3L) {
      stop(
        "Need at least 3 threshold values when `estimate_endpoints = FALSE`.",
        call. = FALSE
      )
    }

    z_used <- z[seq_len(k - 2L)]
    gradient <- numeric(k)
    gradient[1L] <- z_used[1L] - x_vals[1L]
    if (k > 3L) {
      gradient[2L:(k - 2L)] <- diff(z_used)
    }
    gradient[k - 1L] <- x_vals[k] - z_used[length(z_used)]
  }

  names(gradient) <- as.character(x_vals)
  gradient
}

.beta_mcb_analytic_inference <- function(
    precomputed,
    mixtures,
    alpha_grid,
    estimate,
    use_midpoint,
    estimate_endpoints,
    use_isotonic_dist,
    ci_level
) {
  if (!isTRUE(use_isotonic_dist)) {
    stop("Beta analytic inference requires `use_isotonic_dist = TRUE`.",
         call. = FALSE)
  }
  if (!is.list(precomputed) ||
      !all(c("x_grid", "ids", "n", "s") %in% names(precomputed))) {
    stop("`precomputed` must be a binomial grid object.", call. = FALSE)
  }

  x_grid <- precomputed$x_grid
  n_units <- length(precomputed$n)
  n_alpha <- length(alpha_grid)
  n_thresholds <- length(x_grid)
  if (n_units < 2L) {
    stop("Beta analytic inference requires at least two observational units.",
         call. = FALSE)
  }
  if (!any(precomputed$n > 1)) {
    stop(
      "Beta analytic inference requires at least one unit with more than one observation; the two Beta shape parameters are not identifiable when all unit sizes equal one.",
      call. = FALSE
    )
  }
  if (length(mixtures) != n_thresholds) {
    stop("The number of mixture estimates does not match the threshold grid.",
         call. = FALSE)
  }

  cumul_raw <- t(vapply(
    alpha_grid,
    function(alpha) {
      vapply(mixtures, .lookup_cumul_at_alpha, numeric(1L), alpha = alpha)
    },
    numeric(n_thresholds)
  ))
  cumul_isotonic <- cumul_raw
  for (j in seq_len(n_alpha)) {
    cumul_isotonic[j, ] <- -stats::isoreg(
      x = x_grid,
      y = -cumul_raw[j, ]
    )$yf
  }

  grid_jacobian <- t(vapply(
    seq_len(n_alpha),
    function(j) {
      .mcb_grid_jacobian(
        cumul_vals = cumul_isotonic[j, ],
        x_vals = x_grid,
        use_midpoint = use_midpoint,
        estimate_endpoints = estimate_endpoints
      )
    },
    numeric(n_thresholds)
  ))

  influence <- matrix(0, nrow = n_units, ncol = n_alpha)
  diagnostics <- vector("list", n_thresholds)
  information_diagnostics <- vector("list", n_thresholds)

  for (k in seq_len(n_thresholds)) {
    s_k <- precomputed$s[, k]
    all_failure <- all(s_k == 0)
    all_success <- all(s_k == precomputed$n)
    degenerate <- all_failure || all_success
    beta_fit <- attr(mixtures[[k]], "beta_fit", exact = TRUE)
    beta_mle <- attr(mixtures[[k]], "beta_mle", exact = TRUE)

    if (degenerate) {
      if_g_k <- matrix(0, nrow = n_alpha, ncol = n_units)
      shape1 <- NA_real_
      shape2 <- NA_real_
      observed_minimum_eigenvalue <- NA_real_
      observed_condition_number <- NA_real_
      expected_minimum_eigenvalue <- NA_real_
      expected_condition_number <- NA_real_
      expected_score_norm <- 0
      score_norm <- 0
      convergence <- if (is.null(beta_fit)) NA_integer_ else beta_fit$convergence
      objective <- if (is.null(beta_fit)) NA_real_ else beta_fit$objective
      at_boundary <- FALSE
      degeneracy <- if (all_success) "all_success" else "all_failure"
      information_diagnostics[[k]] <- list(
        observed = matrix(NA_real_, 2L, 2L),
        expected = matrix(NA_real_, 2L, 2L),
        expected_by_cluster_size = NULL,
        used = "fixed_degenerate"
      )
    } else {
      shape1 <- if (!is.null(beta_fit)) beta_fit$shape1 else beta_mle$alpha
      shape2 <- if (!is.null(beta_fit)) beta_fit$shape2 else beta_mle$beta
      if (is.null(shape1) || is.null(shape2) || !is.finite(shape1) ||
          !is.finite(shape2) || shape1 <= 0 || shape2 <= 0) {
        stop(
          sprintf(
            "Beta fit at threshold x = %s is nondegenerate but has missing or invalid shape estimates.",
            format(x_grid[k])
          ),
          call. = FALSE
        )
      }

      score_info <- .betabinomial_score_hessian(
        n = precomputed$n,
        s = s_k,
        shape1 = shape1,
        shape2 = shape2
      )
      observed_information <- score_info$information
      expected_info <- .betabinomial_expected_information(
        n = precomputed$n,
        shape1 = shape1,
        shape2 = shape2
      )
      expected_information <- expected_info$information
      if (any(!is.finite(expected_information))) {
        stop(
          sprintf(
            "Expected Beta Fisher information at threshold x = %s is nonfinite.",
            format(x_grid[k])
          ),
          call. = FALSE
        )
      }
      eig_values <- eigen(
        expected_information,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      singular_tol <- max(abs(eig_values)) * .Machine$double.eps * 100
      if (any(!is.finite(eig_values)) ||
          min(eig_values) <= singular_tol) {
        stop(
          sprintf(
            "Expected Beta Fisher information at threshold x = %s is nonfinite or numerically singular (minimum eigenvalue = %s).",
            format(x_grid[k]),
            format(min(eig_values))
          ),
          call. = FALSE
        )
      }

      gradient <- .beta_cdf_shape_gradient(
        alpha = alpha_grid,
        shape1 = shape1,
        shape2 = shape2
      )
      if_g_k <- gradient %*%
        solve(expected_information, t(score_info$score))
      observed_minimum_eigenvalue <- score_info$minimum_eigenvalue
      observed_condition_number <- score_info$condition_number
      expected_minimum_eigenvalue <- expected_info$minimum_eigenvalue
      expected_condition_number <- expected_info$condition_number
      expected_score_norm <- expected_info$expected_score_norm
      score_norm <- score_info$score_norm
      convergence <- if (is.null(beta_fit)) NA_integer_ else beta_fit$convergence
      objective <- if (is.null(beta_fit)) NA_real_ else beta_fit$objective
      at_boundary <- isTRUE(beta_fit$at_boundary)
      degeneracy <- NA_character_
      information_diagnostics[[k]] <- list(
        observed = observed_information,
        expected = expected_information,
        expected_by_cluster_size = expected_info$conditional,
        used = "expected"
      )
    }

    contribution <- sweep(
      if_g_k,
      MARGIN = 1L,
      STATS = grid_jacobian[, k],
      FUN = "*"
    )
    influence <- influence + t(contribution)
    diagnostics[[k]] <- data.frame(
      x = x_grid[k],
      shape1 = shape1,
      shape2 = shape2,
      degenerate = degenerate,
      degeneracy = degeneracy,
      convergence = convergence,
      objective = objective,
      minimum_eigenvalue = expected_minimum_eigenvalue,
      condition_number = expected_condition_number,
      observed_information_minimum_eigenvalue = observed_minimum_eigenvalue,
      observed_information_condition_number = observed_condition_number,
      expected_information_minimum_eigenvalue = expected_minimum_eigenvalue,
      expected_information_condition_number = expected_condition_number,
      score_norm = score_norm,
      expected_score_norm = expected_score_norm,
      at_parameter_boundary = at_boundary,
      pointwise_influence_norm = sqrt(sum(if_g_k^2)),
      isotonic_change = max(abs(cumul_isotonic[, k] - cumul_raw[, k])),
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- do.call(rbind, diagnostics)
  rownames(influence) <- as.character(precomputed$ids)
  colnames(influence) <- as.character(alpha_grid)
  influence_centered <- sweep(
    influence,
    MARGIN = 2L,
    STATS = colMeans(influence),
    FUN = "-"
  )
  cov_mat <- crossprod(influence_centered) /
    (n_units * (n_units - 1))
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  dimnames(cov_mat) <- list(as.character(alpha_grid), as.character(alpha_grid))
  se <- sqrt(pmax(diag(cov_mat), 0))
  z <- stats::qnorm(1 - (1 - ci_level) / 2)

  if (any(diagnostics$condition_number > 1e8, na.rm = TRUE)) {
    warning(
      "Some fitted expected Beta Fisher information matrices are poorly conditioned; inspect `fit$inference$diagnostics`.",
      call. = FALSE
    )
  }
  if (any(diagnostics$at_parameter_boundary, na.rm = TRUE)) {
    warning(
      "Some fitted Beta shape estimates are close to an optimizer boundary; analytic inference may be unreliable.",
      call. = FALSE
    )
  }
  if (any(diagnostics$score_norm > 1e-4, na.rm = TRUE)) {
    warning(
      "Some fitted Beta models have a non-negligible average score; inspect `fit$inference$diagnostics`.",
      call. = FALSE
    )
  }
  if (max(diagnostics$isotonic_change, na.rm = TRUE) > 1e-3) {
    warning(
      "Isotonic regression substantially changed at least one fitted Beta CDF curve; analytic inference does not differentiate through that projection.",
      call. = FALSE
    )
  }
  interior <- if (n_thresholds > 2L) 2L:(n_thresholds - 1L) else integer()
  if (length(interior) > 0L &&
      mean(diagnostics$degenerate[interior]) >= 0.5) {
    warning(
      "Many interior thresholds have degenerate all-success or all-failure Beta fits.",
      call. = FALSE
    )
  }

  list(
    se = se,
    ci_lo = estimate - z * se,
    ci_hi = estimate + z * se,
    cov = cov_mat,
    influence = influence,
    diagnostics = diagnostics,
    information = information_diagnostics,
    grid_jacobian = grid_jacobian,
    cumul_raw = cumul_raw,
    cumul_isotonic = cumul_isotonic
  )
}
