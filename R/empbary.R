#' Empirical Barycenter by Grouped Quantiles
#'
#' Computes group-specific empirical quantiles over a shared grid and then
#' averages those quantiles across groups.
#'
#' @param df A data frame.
#' @param id_col A string giving the column name that identifies groups.
#' @param val_col A string giving the column name containing numeric values.
#' @param alpha_grid A numeric vector of probabilities passed to [stats::quantile()].
#'   Defaults to `seq(0.01, 0.99, 0.01)`.
#' @param quantile_type The `type` argument passed to [stats::quantile()].
#'   Defaults to `3`.
#'
#' @return A list with:
#' \describe{
#'   \item{res}{A data frame with columns `quantile`, `estimate`, `se`,
#'   `ci_lo`, and `ci_hi`.}
#'   \item{cov}{An estimated covariance matrix for the quantile barycenter
#'   estimator, of dimension `length(alpha_grid) x length(alpha_grid)`.}
#'   \item{data}{The original input reduced to standardized `id` and `val`
#'   columns.}
#' }
#' @export
empbary <- function(
    df,
    id_col,
    val_col,
    alpha_grid = seq(0.01, 0.99, 0.01),
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
  
  if (!is.numeric(df[[val_col]])) {
    stop("`val_col` must refer to a numeric column.", call. = FALSE)
  }
  
  data_out <- data.frame(
    id = df[[id_col]],
    val = df[[val_col]]
  )
  
  if (!is.numeric(alpha_grid) || length(alpha_grid) == 0L || anyNA(alpha_grid)) {
    stop(
      "`alpha_grid` must be a non-empty numeric vector without missing values.",
      call. = FALSE
    )
  }
  
  if (any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("`alpha_grid` values must lie between 0 and 1.", call. = FALSE)
  }
  
  if (!is.numeric(quantile_type) || length(quantile_type) != 1L ||
      is.na(quantile_type)) {
    stop("`quantile_type` must be a single numeric value.", call. = FALSE)
  }
  
  split_vals <- split(df[[val_col]], df[[id_col]])
  split_vals <- Filter(function(x) any(!is.na(x)), split_vals)
  
  if (length(split_vals) == 0L) {
    stop(
      "`df` must contain at least one group with non-missing values.",
      call. = FALSE
    )
  }
  
  quantile_matrix <- vapply(
    split_vals,
    FUN.VALUE = numeric(length(alpha_grid)),
    FUN = function(values) {
      stats::quantile(
        values[!is.na(values)],
        probs = alpha_grid,
        type = quantile_type,
        names = FALSE
      )
    }
  )
  
  if (is.null(dim(quantile_matrix))) {
    quantile_matrix <- matrix(quantile_matrix, ncol = 1L)
  }
  
  n_groups <- ncol(quantile_matrix)
  estimate <- rowMeans(quantile_matrix)
  
  if (n_groups > 1L) {
    cov_mat <- stats::cov(t(quantile_matrix)) / n_groups
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
  
  list(
    res = res,
    cov = cov_mat,
    data = data_out
  )
}
