test_that("estimate_mixture_efron supports id-level weights", {
  df <- data.frame(
    id = c("a", "a", "b", "b", "c", "c"),
    x = c(0, 1, 0, 1, 1, 1),
    weight = c(1, 1, 4, 4, 1, 1)
  )

  fit_unweighted <- estimate_mixture_efron(
    df = df,
    id_col = "id",
    val_col = "x",
    x_thresh = 0.5,
    tau = seq(0.1, 0.9, by = 0.2),
    pDegree = 3,
    c0 = 0.1
  )

  fit_weighted <- estimate_mixture_efron(
    df = df,
    id_col = "id",
    val_col = "x",
    x_thresh = 0.5,
    tau = seq(0.1, 0.9, by = 0.2),
    pDegree = 3,
    c0 = 0.1,
    weight_col = "weight"
  )

  expect_named(fit_weighted, c("theta", "g", "cumul"))
  expect_equal(sum(fit_weighted$g), 1)
  expect_false(isTRUE(all.equal(fit_unweighted$g, fit_weighted$g)))
})

test_that("estimate_all_mixtures passes weights to spline method", {
  df <- data.frame(
    id = c("a", "a", "b", "b", "c", "c"),
    x = c(0, 1, 0, 1, 1, 1),
    weight = c(1, 1, 4, 4, 1, 1)
  )

  fit <- estimate_all_mixtures(
    df = df,
    id_col = "id",
    val_col = "x",
    method = "spline",
    x_grid = 0.5,
    weight_col = "weight",
    tau = seq(0.1, 0.9, by = 0.2),
    pDegree = 3,
    c0 = 0.1
  )

  expect_true(is.list(fit))
  expect_true("0.5" %in% names(fit))
  expect_named(fit[["0.5"]], c("theta", "g", "cumul"))
})
