test_that("npmle mixsqp drops impossible tau columns without error", {
  df_bin <- data.frame(
    n = c(2, 2, 2),
    s = c(1, 1, 1)
  )

  fit_g <- .estimate_mixture_npmle_from_binomial(
    df_bin = df_bin,
    tau = c(0, 0.5, 1),
    backend = "mixsqp"
  )

  expect_equal(length(fit_g), 3)
  expect_equal(fit_g, c(0, 1, 0))
})

test_that("npmle mixsqp short-circuits when only one tau is feasible", {
  df_bin <- data.frame(
    n = c(3, 4, 5),
    s = c(0, 0, 0)
  )

  fit_g <- .estimate_mixture_npmle_from_binomial(
    df_bin = df_bin,
    tau = c(0, 1),
    backend = "mixsqp"
  )

  expect_equal(fit_g, c(1, 0))
})
