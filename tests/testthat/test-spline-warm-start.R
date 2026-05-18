test_that("spline method with c0 = 0 still produces quantile estimates", {
  data("sparse_dist_data", package = "mcbarycenter")

  fit <- mcbary(
    df = sparse_dist_data,
    id_col = "id",
    val_col = "value",
    method = "spline",
    cutpoints = 10,
    bootstrap_samples = 2,
    alpha_grid = 0.5,
    progress = FALSE,
    c0 = 0
  )

  expect_true(is.list(fit))
  expect_named(
    fit$res,
    c(
      "quantile", "estimate", "estimate_bs", "se", "ci_lo",
      "ci_hi", "pct_ci_lo", "pct_ci_hi"
    )
  )
  expect_equal(nrow(fit$res), 1)
  expect_true(is.finite(fit$res$estimate))
})
