test_that("deconv_fast returns the expected top-level components", {
  tau <- seq(0.5, 2, length.out = 6)
  x <- c(1, 1, 2, 2, 3, 4, 1, 2, 3, 2)

  fit <- deconv_fast(
    tau = tau,
    X = x,
    family = "Poisson",
    n = 5,
    pDegree = 3
  )

  expect_true(is.list(fit))
  expect_true(all(c("mle", "Q", "P", "S", "cov", "cov.g", "stats") %in% names(fit)))
  expect_equal(nrow(fit$stats), length(tau))
  expect_true(all(c("theta", "g", "SE.g", "G", "SE.G", "Bias.g") %in% colnames(fit$stats)))
})

test_that("deconv_fast supports observation weights", {
  tau <- seq(0.5, 2, length.out = 6)
  x <- c(1, 1, 2, 2, 3, 4, 1, 2, 3, 2)

  fit_unweighted <- deconv_fast(
    tau = tau,
    X = x,
    family = "Poisson",
    n = 5,
    pDegree = 3
  )

  fit_weighted <- deconv_fast(
    tau = tau,
    X = x,
    family = "Poisson",
    n = 5,
    pDegree = 3,
    obs_weights = c(1, 3, 1, 1, 1)
  )

  expect_false(isTRUE(all.equal(fit_unweighted$stats[, "g"], fit_weighted$stats[, "g"])))
})

test_that("deconv_fast uses a stable Binomial log-likelihood path", {
  tau <- seq(0.1, 0.9, length.out = 6)
  x <- cbind(
    n = c(10, 12, 8, 15),
    s = c(2, 5, 1, 9)
  )

  fit <- deconv_fast(
    tau = tau,
    X = x,
    family = "Binomial",
    pDegree = 3
  )

  expect_null(fit$P)
  expect_true(is.matrix(fit$logP))
  expect_equal(dim(fit$logP), c(nrow(x), length(tau)))
  expect_equal(nrow(fit$stats), length(tau))
})
