test_that("mcbary returns mixtures, method, and filtered original data", {
  set.seed(1)
  input <- sample_data[1:40, ]
  input$group_id <- input$id
  input$extra <- seq_len(nrow(input))
  
  fit <- mcbary(
    df = input,
    id_col = "group_id",
    val_col = "value",
    method = "raw",
    cutpoints = 5,
    bootstrap_samples = 5,
    alpha_grid = c(0.25, 0.5),
    progress = FALSE
  )
  
  expect_true(is.list(fit))
  expect_true("mixtures" %in% names(fit))
  expect_true("method" %in% names(fit))
  expect_true("data" %in% names(fit))
  
  expect_equal(fit$method, "raw")
  expect_equal(
    fit$data,
    data.frame(
      id = input$group_id,
      val = input$value
    )
  )
  
  expect_true(is.list(fit$mixtures))
  expect_length(fit$mixtures, 5)
})

test_that("mcbary warns when the estimated quantile function is not monotone", {
  testthat::local_mocked_bindings(
    estimate_all_mixtures = function(...) list("1" = data.frame(
      theta = c(0, 1),
      g = c(1, 0),
      cumul = c(1, 1)
    )),
    est_all_quantiles = function(mixture_res, alpha_grid, ...) {
      data.frame(
        quantile = alpha_grid,
        estimate = c(2, 1)
      )
    }
  )
  
  input <- data.frame(
    id = c(1, 1, 2, 2),
    value = c(1, 2, 3, 4)
  )
  
  expect_warning(
    mcbary(
      df = input,
      id_col = "id",
      val_col = "value",
      method = "raw",
      x_grid = 1,
      bootstrap_samples = 2,
      alpha_grid = c(0.25, 0.75),
      progress = FALSE
    ),
    "The estimated barycenter's quantile function is not monotone increasing."
  )
})

test_that("mcbary can isotonicize estimate and confidence bounds", {
  testthat::local_mocked_bindings(
    estimate_all_mixtures = function(...) list("1" = data.frame(
      theta = c(0, 1),
      g = c(1, 0),
      cumul = c(1, 1)
    )),
    est_all_quantiles = function(mixture_res, alpha_grid, ...) {
      data.frame(
        quantile = alpha_grid,
        estimate = c(3, 1, 2)
      )
    }
  )
  
  input <- data.frame(
    id = c(1, 1, 2, 2),
    value = c(1, 2, 3, 4)
  )
  
  fit <- mcbary(
    df = input,
    id_col = "id",
    val_col = "value",
    method = "raw",
    x_grid = 1,
    bootstrap_samples = 2,
    alpha_grid = c(0.25, 0.5, 0.75),
    use_isotonic = TRUE,
    progress = FALSE
  )
  
  expect_true(all(diff(fit$res$estimate) >= 0))
  expect_true(all(diff(fit$res$ci_lo) >= 0))
  expect_true(all(diff(fit$res$ci_hi) >= 0))
})

test_that("mcbary requires cutpoints or x_grid", {
  input <- data.frame(
    id = c(1, 1, 2, 2),
    value = c(1, 2, 3, 4)
  )

  expect_error(
    mcbary(
      df = input,
      id_col = "id",
      val_col = "value",
      method = "raw",
      x_grid = NULL,
      cutpoints = NULL,
      bootstrap_samples = 2,
      alpha_grid = c(0.25, 0.75),
      progress = FALSE
    ),
    "Either `cutpoints` or `x_grid` must be supplied."
  )
})

test_that("est_dist_alpha returns x, pmf, and cdf and matches mean helper", {
  mixture_res <- list(
    "1" = data.frame(theta = c(0, 0.5, 1), g = c(0.2, 0.3, 0.5), cumul = c(0.2, 0.5, 1)),
    "2" = data.frame(theta = c(0, 0.5, 1), g = c(0.1, 0.4, 0.5), cumul = c(0.1, 0.5, 1)),
    "3" = data.frame(theta = c(0, 0.5, 1), g = c(0.05, 0.25, 0.7), cumul = c(0.05, 0.3, 1))
  )
  
  dist <- est_dist_alpha(
    mixture_res = mixture_res,
    alpha = 0.5,
    use_midpoint = TRUE,
    estimate_first_last = TRUE,
    use_isotonic = FALSE
  )
  
  expect_named(dist, c("x", "pmf", "cdf"))
  expect_equal(sum(dist$pmf), 1)
  expect_equal(tail(dist$cdf, 1), 1)
  expect_equal(
    est_mean_alpha_quantile(
      mixture_res = mixture_res,
      alpha = 0.5,
      use_midpoint = TRUE,
      estimate_first_last = TRUE,
      use_isotonic = FALSE
    ),
    sum(dist$x * dist$pmf)
  )
})
