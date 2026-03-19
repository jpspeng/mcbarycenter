test_that("empb averages grouped quantiles and computes Wald intervals", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, 5, 7)
  )

  result <- empb(
    df = df,
    id_col = "id",
    val_col = "value",
    grid = c(0, 0.5, 1),
    quantile_type = 3
  )

  expected_estimate <- c(3, 3, 5)
  expected_se <- c(2, 2, 2)
  z_value <- qnorm(0.975)

  expect_identical(
    names(result),
    c("quantile", "estimate", "se", "ci_lo", "ci_hi")
  )
  expect_equal(result$quantile, c(0, 0.5, 1))
  expect_equal(result$estimate, expected_estimate)
  expect_equal(result$se, expected_se)
  expect_equal(result$ci_lo, expected_estimate - z_value * expected_se)
  expect_equal(result$ci_hi, expected_estimate + z_value * expected_se)
})

test_that("empb drops groups with only missing values", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    value = c(1, 3, NA, NA)
  )

  result <- empb(
    df = df,
    id_col = "id",
    val_col = "value",
    grid = c(0, 0.5, 1)
  )

  expect_equal(result$estimate, c(1, 1, 3))
  expect_true(all(is.na(result$se)))
})
