test_that("graph_quantiles returns a ggplot object", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3),
    ci_lo = c(0.8, 1.7, 2.5),
    ci_hi = c(1.2, 2.3, 3.5)
  )

  plot <- graph_quantiles(df)

  expect_s3_class(plot, "ggplot")
})

test_that("graph_quantiles can combine result frames with different extra columns", {
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

  plot <- graph_quantiles(df1, df2)

  expect_s3_class(plot, "ggplot")
})

test_that("graph_quantiles uses input object names in the legend", {
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
  
  plot <- graph_quantiles(df1, df2)
  
  expect_equal(levels(plot$data$.group), c("df1", "df2"))
})

test_that("graph_quantiles validates required columns", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3)
  )

  expect_error(
    graph_quantiles(df),
    "Argument 1 must contain `quantile`, `estimate`, `ci_lo`, `ci_hi`."
  )
})

test_that("graph_quantiles allows missing confidence interval columns when disabled", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3)
  )

  plot <- graph_quantiles(df, show_ci = FALSE)

  expect_s3_class(plot, "ggplot")
})

test_that("graph_quantiles still requires confidence interval columns when enabled", {
  df <- data.frame(
    quantile = c(0, 0.5, 1),
    estimate = c(1, 2, 3)
  )

  expect_error(
    graph_quantiles(df, show_ci = TRUE),
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
  expect_equal(plot$labels$x, "alpha")
  expect_equal(plot$labels$y, "cdf")
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
  expect_equal(plots[[1]]$labels$x, "alpha")
  expect_equal(plots[[1]]$labels$y, "cdf")
})

test_that("graph_mixtures subsamples to about 10 by default and can show all", {
  mixtures <- lapply(seq_len(15), function(i) {
    data.frame(theta = c(0, 1), cumul = c(i / 20, 1))
  })
  names(mixtures) <- as.character(seq_len(15))
  fit <- list(mixtures = mixtures)
  
  plot_subset <- graph_mixtures(fit, layout = "combined")
  plot_all <- graph_mixtures(fit, layout = "combined", show_all = TRUE)
  
  expect_length(levels(plot_subset$data$.mixture), 10)
  expect_length(levels(plot_all$data$.mixture), 15)
})
