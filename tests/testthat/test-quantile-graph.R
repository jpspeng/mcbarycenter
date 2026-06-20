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

test_that("graph_quantiles accepts S3 barycenter results", {
  fit <- structure(
    list(
      res = data.frame(
        quantile = c(0, 0.5, 1),
        estimate = c(1, 2, 3),
        ci_lo = c(0.8, 1.7, 2.5),
        ci_hi = c(1.2, 2.3, 3.5)
      ),
      cov = diag(3),
      data = data.frame(id = c("a", "a"), val = c(1, 3))
    ),
    class = "empbary_result"
  )

  plot <- graph_quantiles(fit)

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

test_that("graph_individual_quantiles returns a ggplot object", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    val = c(1, 3, 2, 4)
  )

  plot <- graph_individual_quantiles(
    df,
    id_col = "id",
    val_col = "val",
    alpha_grid = c(0, 0.5, 1),
    interactive = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(levels(plot$data$id), c("a", "b"))
})

test_that("graph_individual_quantiles can thin the quantile grid for plotting", {
  df <- data.frame(
    id = rep(c("a", "b"), each = 20),
    val = seq_len(40)
  )

  plot <- graph_individual_quantiles(
    df,
    id_col = "id",
    val_col = "val",
    alpha_grid = seq(0, 1, length.out = 11),
    max_points = 5,
    interactive = FALSE
  )

  plotted_grid <- sort(unique(plot$data$quantile))

  expect_length(plotted_grid, 5)
  expect_equal(plotted_grid[c(1, 5)], c(0, 1))
  expect_true(is.unsorted(plotted_grid, strictly = TRUE) == FALSE)
})

test_that("graph_individual_quantiles validates required columns", {
  df <- data.frame(
    id = c("a", "a"),
    value = c(1, 2)
  )

  expect_error(
    graph_individual_quantiles(df, id_col = "missing", val_col = "value"),
    "`id_col` must name a column in `df`."
  )

  expect_error(
    graph_individual_quantiles(df, id_col = "id", val_col = "missing"),
    "`val_col` must name a column in `df`."
  )

  expect_error(
    graph_individual_quantiles(
      transform(df, val = value),
      id_col = "id",
      val_col = "val",
      max_points = 0
    ),
    "`max_points` must be NULL or a single positive whole number."
  )
})

test_that("graph_individual_quantiles requires interactive dependencies", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    val = c(1, 3, 2, 4)
  )

  skip_if(requireNamespace("plotly", quietly = TRUE))

  expect_error(
    graph_individual_quantiles(df, id_col = "id", val_col = "val"),
    "`plotly` must be installed when `interactive = TRUE`."
  )
})

test_that("graph_individual_quantiles returns a plotly htmlwidget", {
  df <- data.frame(
    id = c("a", "a", "b", "b"),
    val = c(1, 3, 2, 4)
  )

  skip_if_not_installed("plotly")
  skip_if_not_installed("htmlwidgets")

  plot <- graph_individual_quantiles(
    df,
    id_col = "id",
    val_col = "val",
    alpha_grid = seq(0, 1, length.out = 9),
    max_points = 4,
    interactive = TRUE
  )

  expect_s3_class(plot, "plotly")
  expect_s3_class(plot, "htmlwidget")
  expect_true(all(vapply(plot$x$data, function(trace) length(trace$x), integer(1)) == 4L))
})

test_that("graph_mixtures returns a combined ggplot", {
  fit <- structure(
    list(
      mixtures = list(
        "1" = data.frame(theta = c(0, 0.5, 1), cumul = c(0.2, 0.7, 1)),
        "2" = data.frame(theta = c(0, 0.5, 1), cumul = c(0.1, 0.8, 1))
      )
    ),
    class = "mcbary_result"
  )
  
  plot <- graph_mixtures(fit, layout = "combined")
  
  expect_s3_class(plot, "ggplot")
  expect_equal(levels(plot$data$.mixture), c("1", "2"))
  expect_equal(plot$labels$x, "alpha")
  expect_equal(plot$labels$y, "cdf")
})

test_that("graph_mixtures returns individual ggplots with titled names", {
  fit <- structure(
    list(
      mixtures = list(
        "0.5" = data.frame(theta = c(0, 1), cumul = c(0.3, 1)),
        "1.0" = data.frame(theta = c(0, 1), cumul = c(0.4, 1))
      )
    ),
    class = "mcbary_result"
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
  fit <- structure(list(mixtures = mixtures), class = "mcbary_result")
  
  plot_subset <- graph_mixtures(fit, layout = "combined")
  plot_all <- graph_mixtures(fit, layout = "combined", show_all = TRUE)
  
  expect_length(levels(plot_subset$data$.mixture), 10)
  expect_length(levels(plot_all$data$.mixture), 15)
})

test_that("graph_mixtures requires an mcbary_result object", {
  fit <- list(
    mixtures = list(
      "1" = data.frame(theta = c(0, 1), cumul = c(0.2, 1))
    )
  )

  expect_error(
    graph_mixtures(fit),
    "`mcb_res` must be an object of class `mcbary_result`."
  )
})
