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
