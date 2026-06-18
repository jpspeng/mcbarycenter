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

test_that("estimate_mixture_efron can add explicit endpoint masses", {
  df <- data.frame(
    id = c("a", "a", "b", "b", "c", "c"),
    x = c(0, 1, 0, 1, 1, 1)
  )

  fit <- estimate_mixture_efron(
    df = df,
    id_col = "id",
    val_col = "x",
    x_thresh = 0.5,
    tau = seq(0.2, 0.8, by = 0.2),
    pDegree = 3,
    c0 = 0.1,
    mass_at_endpoints = TRUE
  )

  expect_true(0 %in% fit$theta)
  expect_true(1 %in% fit$theta)
  expect_equal(sum(fit$g), 1)
})

test_that("mixture estimators accept non-default weight column names", {
  df <- data.frame(
    grp = c("a", "a", "b", "b", "c", "c"),
    val = c(0, 1, 0, 1, 1, 1),
    wts = c(1, 1, 4, 4, 1, 1)
  )

  fit_raw <- estimate_mixture_raw(
    df = df,
    id_col = "grp",
    val_col = "val",
    x_thresh = 0.5,
    weight_col = "wts"
  )

  fit_beta <- estimate_mixture_beta(
    df = df,
    id_col = "grp",
    val_col = "val",
    x_thresh = 0.5,
    tau = seq(0.1, 0.9, by = 0.2),
    weight_col = "wts"
  )

  fit_all <- estimate_all_mixtures(
    df = df,
    id_col = "grp",
    val_col = "val",
    method = "raw",
    x_grid = 0.5,
    weight_col = "wts"
  )

  expect_named(fit_raw, c("theta", "g", "cumul"))
  expect_named(fit_beta, c("theta", "g", "cumul"))
  expect_true(is.list(fit_all))
  expect_true("0.5" %in% names(fit_all))
})
