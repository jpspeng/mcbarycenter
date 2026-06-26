test_that("empbary averages grouped quantiles and computes Wald intervals", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, 5, 7)
  )

  result <- empbary(
    df = df,
    id_col = "id",
    val_col = "value",
    alpha_grid = c(0, 0.5, 1),
    quantile_type = 3
  )

  expected_estimate <- c(3, 3, 5)
  expected_se <- c(2, 2, 2)
  z_value <- qnorm(0.975)

  expect_s3_class(result, "empbary_result")
  expect_identical(
    names(result),
    c("res", "cov", "data")
  )
  expect_equal(result$res$quantile, c(0, 0.5, 1))
  expect_equal(result$res$estimate, expected_estimate)
  expect_equal(result$res$se, expected_se)
  expect_equal(result$res$ci_lo, expected_estimate - z_value * expected_se)
  expect_equal(result$res$ci_hi, expected_estimate + z_value * expected_se)
  expect_equal(result$data, data.frame(id = df$id, val = df$value))
})

test_that("empbary drops groups with only missing values", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, NA, NA)
  )

  result <- empbary(
    df = df,
    id_col = "id",
    val_col = "value",
    alpha_grid = c(0, 0.5, 1)
  )

  expect_s3_class(result, "empbary_result")
  expect_equal(result$res$estimate, c(1, 1, 3))
  expect_true(all(is.na(result$res$se)))
  expect_equal(result$data, data.frame(id = df$id, val = df$value))
})

test_that("empbary computes a weighted barycenter when weight_col is supplied", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, 5, 7),
    weight = c(1, 1, 3, 3)
  )

  result <- empbary(
    df = df,
    id_col = "id",
    val_col = "value",
    weight_col = "weight",
    alpha_grid = c(0, 0.5, 1),
    quantile_type = 3
  )

  expect_s3_class(result, "empbary_result")
  expected_estimate <- c(4, 4, 6)
  expected_cov <- matrix(
    c(5, 5, 5,
      5, 5, 5,
      5, 5, 5),
    nrow = 3
  )
  expected_se <- rep(sqrt(5), 3)
  z_value <- qnorm(0.975)

  expect_equal(result$res$estimate, expected_estimate)
  expect_equal(result$res$se, expected_se)
  expect_equal(result$res$ci_lo, expected_estimate - z_value * expected_se)
  expect_equal(result$res$ci_hi, expected_estimate + z_value * expected_se)
  expect_equal(result$cov, expected_cov)
  expect_equal(result$data, data.frame(id = df$id, val = df$value, weight = df$weight))
})

test_that("empbary can presmooth grouped samples before averaging", {
  df <- data.frame(
    id = c("a", "a", "a", "b", "b", "b"),
    value = c(1, 1, 1, 3, 3, 3)
  )

  result <- empbary(
    df = df,
    id_col = "id",
    val_col = "value",
    alpha_grid = c(0.25, 0.5, 0.75),
    presmooth = TRUE
  )

  expect_s3_class(result, "empbary_result")
  expect_equal(result$res$quantile, c(0.25, 0.5, 0.75))
  expect_true(all(is.finite(result$res$estimate)))
  expect_true(all(diff(result$res$estimate) > 0))
  expect_equal(result$res$estimate[[2]], 2, tolerance = 1e-4)
  expect_true(all(is.finite(result$res$se)))
})

test_that("empbary presmoothing supports weighted barycenters", {
  df <- data.frame(
    id = c("a", "a", "a", "b", "b", "b"),
    value = c(1, 1, 1, 3, 3, 3),
    weight = c(1, 1, 1, 3, 3, 3)
  )

  result <- empbary(
    df = df,
    id_col = "id",
    val_col = "value",
    weight_col = "weight",
    alpha_grid = c(0.5),
    presmooth = TRUE
  )

  expect_equal(result$res$estimate[[1]], 2.5, tolerance = 1e-4)
})

test_that("empbary handles a single quantile level", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, 5, 7),
    weight = c(1, 1, 3, 3)
  )

  result <- empbary(
    df = df,
    id_col = "id",
    val_col = "value",
    weight_col = "weight",
    alpha_grid = 0.5,
    quantile_type = 3
  )

  expect_equal(result$res$quantile, 0.5)
  expect_equal(result$res$estimate[[1]], 4)
  expect_equal(result$res$se[[1]], sqrt(5))
})

test_that("empbary_result has a custom print method", {
  result <- structure(
    list(
      res = data.frame(quantile = c(0.25, 0.5), estimate = c(1, 2)),
      cov = diag(2),
      data = data.frame(id = c("a", "a"), val = c(1, 2))
    ),
    class = "empbary_result"
  )

  expect_output(print(result), "<empbary_result>")
  expect_output(print(result), "Quantiles: 2")
})

test_that("empbary requires weights to be homogeneous within id", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, 5, 7),
    weight = c(1, 2, 3, 3)
  )

  expect_error(
    empbary(
      df = df,
      id_col = "id",
      val_col = "value",
      weight_col = "weight"
    ),
    "`weight_col` must be non-missing and homogeneous within each `id`."
  )
})

test_that("empbary validates the presmooth flag", {
  df <- data.frame(
    id = c("a", "a"),
    value = c(1, 3)
  )

  expect_error(
    empbary(
      df = df,
      id_col = "id",
      val_col = "value",
      presmooth = NA
    ),
    "`presmooth` must be a single TRUE/FALSE value."
  )
})
