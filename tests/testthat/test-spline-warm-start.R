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

test_that("spline method with endpoint masses works through warm-start cache", {
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
    tau = seq(0.1, 0.9, by = 0.2),
    pDegree = 3,
    c0 = 0.1,
    mass_at_endpoints = TRUE
  )

  expect_true(is.list(fit))
  expect_true(all(vapply(fit$mixtures, function(mix) 0 %in% mix$theta && 1 %in% mix$theta, logical(1))))
  expect_true(is.finite(fit$res$estimate))
})
