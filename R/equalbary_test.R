#' Test Equality of Two Barycenter Estimates
#'
#' Uses a Gaussian approximation with simulated supremum and integrated squared
#' statistics to compare two barycenter estimates.
#'
#' @param bary_res1 First barycenter result, containing `res$estimate` and
#'   `cov`.
#' @param bary_res2 Second barycenter result, containing `res$estimate` and
#'   `cov`.
#' @param B Number of Gaussian simulations.
#' @param seed Optional random seed.
#' @param variance_explained Number in `(0, 1]` giving the target cumulative
#'   variance explained for PCA truncation in the Hotelling test. Defaults to
#'   `0.95`.
#'
#' @return A list with components `sup_norm`, `int_sq_norm`, and `hotelling`.
#'   `sup_norm` and `int_sq_norm` each contain `pval`, `T_obs`, and `T_sim`.
#'   `hotelling` contains `pval`, `T_obs`, `df`, `mean_diff`, `cov_diff`,
#'   `variance_explained`, and `n_components`.
#' @export
equalbary_test <- function(bary_res1,
                           bary_res2,
                           B = 5000,
                           seed = NULL,
                           variance_explained = 0.95) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (is.null(bary_res1$data) || is.null(bary_res2$data)) {
    stop(
      "Both barycenter results must contain a `data` component.",
      call. = FALSE
    )
  }
  
  if (ncol(bary_res1$data) < 1 || ncol(bary_res2$data) < 1) {
    stop(
      "Both `data` components must contain an id column.",
      call. = FALSE
    )
  }

  est1 <- bary_res1$res$estimate
  est2 <- bary_res2$res$estimate
  alpha1 <- bary_res1$res$quantile
  alpha2 <- bary_res2$res$quantile
  cov1 <- as.matrix(bary_res1$cov)
  cov2 <- as.matrix(bary_res2$cov)
  
  if (length(est1) != length(est2) || !isTRUE(all.equal(alpha1, alpha2))) {
    stop(
      "`bary_res1` and `bary_res2` must use the same quantile grid.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) {
    stop("`B` must be a single integer >= 1.", call. = FALSE)
  }
  
  if (!is.numeric(variance_explained) ||
      length(variance_explained) != 1L ||
      is.na(variance_explained) ||
      variance_explained <= 0 ||
      variance_explained > 1) {
    stop(
      "`variance_explained` must be a single numeric value in (0, 1].",
      call. = FALSE
    )
  }
  
  B <- as.integer(B)
  
  n1 <- dplyr::n_distinct(bary_res1$data[[1]])
  n2 <- dplyr::n_distinct(bary_res2$data[[1]])

  diff_est <- est1 - est2
  cov_diff <- cov1 + cov2

  n <- n1 + n2
  Sigma <- n * cov_diff

  eig <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  
  cov_diff_sym <- (cov_diff + t(cov_diff)) / 2
  cov_diff_eig <- eigen(cov_diff_sym, symmetric = TRUE)
  cov_diff_eig$values[cov_diff_eig$values < 0] <- 0
  max_eig <- max(cov_diff_eig$values)
  tol <- if (max_eig > 0) {
    max_eig * .Machine$double.eps * length(diff_est)
  } else {
    .Machine$double.eps
  }
  keep_positive <- cov_diff_eig$values > tol
  
  if (!any(keep_positive)) {
    stop(
      "The combined covariance matrix is numerically singular.",
      call. = FALSE
    )
  }
  
  positive_values <- cov_diff_eig$values[keep_positive]
  cumulative_share <- cumsum(positive_values) / sum(positive_values)
  n_keep <- which(cumulative_share >= variance_explained)[1]
  keep_idx <- which(keep_positive)[seq_len(n_keep)]
  
  eigvec_keep <- cov_diff_eig$vectors[, keep_idx, drop = FALSE]
  eigval_keep <- cov_diff_eig$values[keep_idx]
  cov_diff_ginv <- eigvec_keep %*%
    diag(1 / eigval_keep, nrow = length(eigval_keep)) %*%
    t(eigvec_keep)
  hotelling_df <- length(keep_idx)
  variance_explained_used <- sum(eigval_keep) / sum(cov_diff_eig$values[keep_positive])

  trapz <- function(x, y) {
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }
  
  scaled_diff <- sqrt(n) * diff_est
  T_obs_sup <- max(abs(scaled_diff))
  T_obs_int <- trapz(alpha1, scaled_diff^2)
  T_obs_hotelling <- as.numeric(t(diff_est) %*% cov_diff_ginv %*% diff_est)

  Z <- matrix(rnorm(length(diff_est) * B), nrow = length(diff_est))
  G <- root %*% Z
  T_sim_sup <- apply(abs(G), 2, max)
  T_sim_int <- apply(G, 2, function(g) trapz(alpha1, g^2))

  pval_sup <- mean(T_sim_sup >= T_obs_sup)
  pval_int <- mean(T_sim_int >= T_obs_int)
  pval_hotelling <- stats::pchisq(
    T_obs_hotelling,
    df = hotelling_df,
    lower.tail = FALSE
  )

  list(
    sup_norm = list(
      pval = pval_sup,
      T_obs = T_obs_sup,
      T_sim = T_sim_sup
    ),
    int_sq_norm = list(
      pval = pval_int,
      T_obs = T_obs_int,
      T_sim = T_sim_int
    ),
    hotelling = list(
      pval = pval_hotelling,
      T_obs = T_obs_hotelling,
      df = hotelling_df,
      mean_diff = diff_est,
      cov_diff = cov_diff,
      variance_explained = variance_explained_used,
      n_components = hotelling_df
    )
  )
}
