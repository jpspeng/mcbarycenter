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

  expect_equal(result$res$estimate, c(1, 1, 3))
  expect_true(all(is.na(result$res$se)))
  expect_equal(result$data, data.frame(id = df$id, val = df$value))
})
