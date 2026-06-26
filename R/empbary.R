.trapz_integral <- function(x, y) {
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

.cumtrapz_integral <- function(x, y) {
  c(0, cumsum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2))
}

.presmooth_bandwidth <- function(values) {
  n <- length(values)
  scale_candidates <- c(stats::sd(values), stats::IQR(values) / 1.34)
  scale_candidates <- scale_candidates[is.finite(scale_candidates) & scale_candidates > 0]

  if (length(scale_candidates) == 0L) {
    value_range <- diff(range(values))
    if (is.finite(value_range) && value_range > 0) {
      scale_estimate <- value_range / 4
    } else {
      scale_estimate <- max(1, abs(values[[1]])) * 0.1
    }
  } else {
    scale_estimate <- min(scale_candidates)
  }

  max(2 * 0.9 * scale_estimate * n^(-1 / 5), sqrt(.Machine$double.eps))
}

.presmooth_grid <- function(split_vals, alpha_grid) {
  bandwidths <- vapply(
    split_vals,
    FUN.VALUE = numeric(1),
    FUN = .presmooth_bandwidth
  )
  pooled_values <- unlist(split_vals, use.names = FALSE)
  support_range <- range(pooled_values)
  max_bandwidth <- max(bandwidths)
  padding <- 4 * max_bandwidth

  if (!is.finite(diff(support_range)) || diff(support_range) == 0) {
    center <- support_range[[1]]
    support_range <- center + c(-padding, padding)
  } else {
    support_range <- c(
      support_range[[1]] - padding,
      support_range[[2]] + padding
    )
  }

  list(
    from = support_range[[1]],
    to = support_range[[2]],
    n = max(512L, 4L * length(alpha_grid))
  )
}

.presmooth_quantiles <- function(values, alpha_grid, grid_spec) {
  bandwidth <- .presmooth_bandwidth(values)
  density_fit <- stats::density(
    values,
    bw = bandwidth,
    kernel = "gaussian",
    from = grid_spec$from,
    to = grid_spec$to,
    n = grid_spec$n
  )

  density_area <- .trapz_integral(density_fit$x, density_fit$y)
  density_values <- density_fit$y / density_area
  cdf_values <- .cumtrapz_integral(density_fit$x, density_values)
  cdf_values <- cdf_values / tail(cdf_values, 1L)

  cdf_interp <- c(0, cdf_values, 1)
  x_interp <- c(
    density_fit$x[[1]],
    density_fit$x,
    density_fit$x[[length(density_fit$x)]]
  )
  keep <- !duplicated(cdf_interp)

  stats::approx(
    x = cdf_interp[keep],
    y = x_interp[keep],
    xout = alpha_grid,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

#' Empirical Barycenter by Grouped Quantiles
#'
#' Computes group-specific empirical quantiles over a shared grid and then
#' averages those quantiles across groups.
#'
#' @param df A data frame.
#' @param id_col A string giving the column name that identifies groups.
#' @param val_col A string giving the column name containing numeric values.
#' @param weight_col A string giving the column name containing group weights.
#'   Defaults to `NULL`, in which case the empirical barycenter is unweighted.
#' @param alpha_grid A numeric vector of probabilities passed to [stats::quantile()].
#'   Defaults to `seq(0.01, 0.99, 0.01)`.
#' @param presmooth Logical; if `TRUE`, estimate a Gaussian KDE for each group
#'   on a shared grid, convert that KDE to a CDF, and invert the CDF at
#'   `alpha_grid` before averaging. Defaults to `FALSE`.
#' @param quantile_type The `type` argument passed to [stats::quantile()].
#'   Defaults to `1`. Ignored when `presmooth = TRUE`.
#'
#' @return An S3 object of class `"empbary_result"` with:
#' \describe{
#'   \item{res}{A data frame with columns `quantile`, `estimate`, `se`,
#'   `ci_lo`, and `ci_hi`.}
#'   \item{cov}{An estimated covariance matrix for the quantile barycenter
#'   estimator, of dimension `length(alpha_grid) x length(alpha_grid)`.}
#'   \item{data}{The original input reduced to standardized `id`, `val`, and
#'   optional `weight` columns.}
#' }
#' @export
empbary <- function(
    df,
    id_col,
    val_col,
    weight_col = NULL,
    alpha_grid = seq(0.01, 0.99, 0.01),
    presmooth = FALSE,
    quantile_type = 1
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  
  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("`id_col` must be a single non-empty string.", call. = FALSE)
  }
  
  if (!is.character(val_col) || length(val_col) != 1L || !nzchar(val_col)) {
    stop("`val_col` must be a single non-empty string.", call. = FALSE)
  }
  
  if (!id_col %in% names(df)) {
    stop("`id_col` must name a column in `df`.", call. = FALSE)
  }
  
  if (!val_col %in% names(df)) {
    stop("`val_col` must name a column in `df`.", call. = FALSE)
  }
  
  if (!is.null(weight_col)) {
    if (!is.character(weight_col) || length(weight_col) != 1L || !nzchar(weight_col)) {
      stop("`weight_col` must be `NULL` or a single non-empty string.", call. = FALSE)
    }
    
    if (!weight_col %in% names(df)) {
      stop("`weight_col` must name a column in `df`.", call. = FALSE)
    }
    
    if (!is.numeric(df[[weight_col]])) {
      stop("`weight_col` must refer to a numeric column.", call. = FALSE)
    }
  }
  
  if (!is.numeric(df[[val_col]])) {
    stop("`val_col` must refer to a numeric column.", call. = FALSE)
  }
  
  data_out <- data.frame(
    id = df[[id_col]],
    val = df[[val_col]]
  )
  
  if (!is.null(weight_col)) {
    data_out$weight <- df[[weight_col]]
  }
  
  if (!is.numeric(alpha_grid) || length(alpha_grid) == 0L || anyNA(alpha_grid)) {
    stop(
      "`alpha_grid` must be a non-empty numeric vector without missing values.",
      call. = FALSE
    )
  }
  
  if (any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("`alpha_grid` values must lie between 0 and 1.", call. = FALSE)
  }
  
  if (!is.logical(presmooth) || length(presmooth) != 1L || is.na(presmooth)) {
    stop("`presmooth` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  
  if (!presmooth &&
      (!is.numeric(quantile_type) || length(quantile_type) != 1L ||
       is.na(quantile_type))) {
    stop("`quantile_type` must be a single numeric value.", call. = FALSE)
  }
  
  split_vals <- split(df[[val_col]], df[[id_col]])
  split_vals <- Filter(function(x) any(!is.na(x)), split_vals)
  split_vals <- lapply(split_vals, function(x) x[!is.na(x)])
  
  if (length(split_vals) == 0L) {
    stop(
      "`df` must contain at least one group with non-missing values.",
      call. = FALSE
    )
  }
  
  if (is.null(weight_col)) {
    group_weights <- rep(1, length(split_vals))
    names(group_weights) <- names(split_vals)
  } else {
    split_weights <- split(df[[weight_col]], df[[id_col]])
    invalid_weights <- vapply(
      split_weights,
      FUN.VALUE = logical(1),
      FUN = function(weights) {
        anyNA(weights) || length(unique(weights)) != 1L
      }
    )
    
    if (any(invalid_weights)) {
      bad_ids <- names(split_weights)[invalid_weights]
      stop(
        paste0(
          "`weight_col` must be non-missing and homogeneous within each `id`. ",
          "Problematic ids: ",
          paste(bad_ids, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
    
    group_weights <- vapply(
      split_weights[names(split_vals)],
      FUN.VALUE = numeric(1),
      FUN = function(weights) weights[[1]]
    )
    
    if (any(group_weights < 0)) {
      stop("`weight_col` must contain non-negative weights.", call. = FALSE)
    }
    
    if (all(group_weights == 0)) {
      stop(
        "`weight_col` must contain at least one positive weight among retained groups.",
        call. = FALSE
      )
    }
  }
  
  grid_spec <- NULL
  if (presmooth) {
    grid_spec <- .presmooth_grid(split_vals, alpha_grid)
  }
  
  quantile_matrix <- vapply(
    split_vals,
    FUN.VALUE = numeric(length(alpha_grid)),
    FUN = function(values) {
      if (presmooth) {
        .presmooth_quantiles(
          values = values,
          alpha_grid = alpha_grid,
          grid_spec = grid_spec
        )
      } else {
        stats::quantile(
          values,
          probs = alpha_grid,
          type = quantile_type,
          names = FALSE
        )
      }
    }
  )
  
  if (is.null(dim(quantile_matrix))) {
    quantile_matrix <- matrix(
      quantile_matrix,
      nrow = length(alpha_grid),
      dimnames = list(NULL, names(split_vals))
    )
  }
  
  n_groups <- ncol(quantile_matrix)
  normalized_weights <- group_weights / sum(group_weights)
  estimate <- as.vector(quantile_matrix %*% normalized_weights)
  
  if (n_groups > 1L && sum(normalized_weights > 0) > 1L) {
    centered_matrix <- sweep(quantile_matrix, 1L, estimate, FUN = "-")
    weight_sum_sq <- sum(normalized_weights^2)
    cov_mat <- (
      centered_matrix %*%
        diag(normalized_weights, nrow = length(normalized_weights)) %*%
        t(centered_matrix)
    ) * (weight_sum_sq / (1 - weight_sum_sq))
    se <- sqrt(diag(cov_mat))
  } else {
    cov_mat <- matrix(
      NA_real_,
      nrow = length(alpha_grid),
      ncol = length(alpha_grid)
    )
    se <- rep(NA_real_, length(alpha_grid))
  }
  
  z_value <- stats::qnorm(0.975)
  
  res <- data.frame(
    quantile = alpha_grid,
    estimate = estimate,
    se = se,
    ci_lo = estimate - z_value * se,
    ci_hi = estimate + z_value * se,
    row.names = NULL
  )
  
  structure(
    list(
    res = res,
    cov = cov_mat,
    data = data_out
    ),
    class = "empbary_result"
  )
}
