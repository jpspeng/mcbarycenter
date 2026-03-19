test_that("quantile_graph returns a ggplot object", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3),
    ci_lo = c(0.8, 1.7, 2.5),
    ci_hi = c(1.2, 2.3, 3.5)
  )

  plot <- quantile_graph(df)

  expect_s3_class(plot, "ggplot")
})

test_that("quantile_graph validates required columns", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3)
  )

  expect_error(
    quantile_graph(df),
    "`df` must contain columns `quantile`, `estimate`, `ci_lo`, and `ci_hi`."
  )
})
