#' Test Equality of Two Barycenter Estimates
#'
#' Uses a Gaussian approximation with simulated supremum statistics to compare
#' two barycenter estimates.
#'
#' @param bary_res1 First barycenter result, containing `res$estimate` and
#'   `cov`.
#' @param bary_res2 Second barycenter result, containing `res$estimate` and
#'   `cov`.
#' @param n1 Sample size associated with `bary_res1`.
#' @param n2 Sample size associated with `bary_res2`.
#' @param B Number of Gaussian simulations.
#' @param seed Optional random seed.
#'
#' @return A list with components `pval`, `T_obs`, and `T_sim`.
#' @export
equalbary_test <- function(bary_res1, bary_res2, n1, n2,
                           B = 5000, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  est1 <- bary_res1$res$estimate
  est2 <- bary_res2$res$estimate
  cov1 <- as.matrix(bary_res1$cov)
  cov2 <- as.matrix(bary_res2$cov)

  diff_est <- est1 - est2
  cov_diff <- cov1 + cov2

  n <- n1 + n2
  Sigma <- n * cov_diff

  eig <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0
  root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

  T_obs <- sqrt(n) * max(abs(diff_est))

  Z <- matrix(rnorm(length(diff_est) * B), nrow = length(diff_est))
  G <- root %*% Z
  T_sim <- apply(abs(G), 2, max)

  pval <- mean(T_sim >= T_obs)

  list(
    pval = pval,
    T_obs = T_obs,
    T_sim = T_sim
  )
}
