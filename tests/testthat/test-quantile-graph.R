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

test_that("quantile_graph can combine result frames with different extra columns", {
  df1 <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3),
    ci_lo = c(0.8, 1.7, 2.5),
    ci_hi = c(1.2, 2.3, 3.5)
  )

  df2 <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1.5, 2.5, 3.5),
    ci_lo = c(1.1, 2.1, 3.1),
    ci_hi = c(1.9, 2.9, 3.9),
    estimate_bs = c(1.4, 2.4, 3.4)
  )

  plot <- quantile_graph(df1, df2)

  expect_s3_class(plot, "ggplot")
})

test_that("quantile_graph uses input object names in the legend", {
  df1 <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3),
    ci_lo = c(0.8, 1.7, 2.5),
    ci_hi = c(1.2, 2.3, 3.5)
  )
  
  df2 <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1.5, 2.5, 3.5),
    ci_lo = c(1.1, 2.1, 3.1),
    ci_hi = c(1.9, 2.9, 3.9)
  )
  
  plot <- quantile_graph(df1, df2)
  
  expect_equal(levels(plot$data$.group), c("df1", "df2"))
})

test_that("quantile_graph validates required columns", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3)
  )

  expect_error(
    quantile_graph(df),
    "Argument 1 must contain `quantile`, `estimate`, `ci_lo`, `ci_hi`."
  )
})

test_that("graph_mixtures returns a combined ggplot", {
  fit <- list(
    mixtures = list(
      "1" = data.frame(theta = c(0, 0.5, 1), cumul = c(0.2, 0.7, 1)),
      "2" = data.frame(theta = c(0, 0.5, 1), cumul = c(0.1, 0.8, 1))
    )
  )
  
  plot <- graph_mixtures(fit, layout = "combined")
  
  expect_s3_class(plot, "ggplot")
  expect_equal(levels(plot$data$.mixture), c("1", "2"))
})

test_that("graph_mixtures returns individual ggplots with titled names", {
  fit <- list(
    mixtures = list(
      "0.5" = data.frame(theta = c(0, 1), cumul = c(0.3, 1)),
      "1.0" = data.frame(theta = c(0, 1), cumul = c(0.4, 1))
    )
  )
  
  plots <- graph_mixtures(fit, layout = "individual")
  
  expect_type(plots, "list")
  expect_named(plots, c("0.5", "1.0"))
  expect_s3_class(plots[[1]], "ggplot")
  expect_equal(plots[[1]]$labels$title, "X_0: 0.5")
  expect_equal(plots[[2]]$labels$title, "X_0: 1.0")
})
