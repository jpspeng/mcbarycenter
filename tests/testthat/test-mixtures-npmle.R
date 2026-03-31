test_that("mixsqp binomial collapse preserves grouped weights", {
  df_bin <- data.frame(
    n = c(5, 5, 5, 10, 10),
    s = c(2, 2, 4, 3, 3),
    weight = c(1.5, 2.5, 4, 1, 3)
  )

  collapsed <- .collapse_binomial_rows_for_mixsqp(df_bin)

  expect_equal(collapsed$n, c(5, 5, 10))
  expect_equal(collapsed$s, c(2, 4, 3))
  expect_equal(collapsed$weight, c(4, 4, 4))
})

test_that("mixsqp binomial collapse matches duplicate counts when unweighted", {
  df_bin <- data.frame(
    n = c(5, 5, 5, 10, 10, 10),
    s = c(2, 2, 4, 3, 3, 3)
  )

  collapsed <- .collapse_binomial_rows_for_mixsqp(df_bin)

  expect_equal(collapsed$n, c(5, 5, 10))
  expect_equal(collapsed$s, c(2, 4, 3))
  expect_equal(collapsed$weight, c(2, 1, 3))
})

test_that("npmle warm starts are forwarded across threshold-grid fits", {
  precomputed <- list(
    x_grid = c(0.1, 0.2),
    ids = c("a", "b"),
    n = c(2, 2),
    s = matrix(c(0, 1, 1, 2), nrow = 2, byrow = FALSE),
    weight = NULL
  )

  tracker <- new.env(parent = emptyenv())
  tracker$calls <- 0L
  tracker$starts <- list()
  old_fun <- get(".estimate_mixture_npmle_from_binomial", envir = asNamespace("mcbarycenter"))

  mock_fun <- function(df_bin,
                       tau = seq(from = 0, to = 1, by = 0.005),
                       backend = c("mixsqp", "REBayes"),
                       start = NULL,
                       mixsqp_control = list(),
                       rebayes_control = list()) {
    tracker$calls <- tracker$calls + 1L
    tracker$starts[[tracker$calls]] <- start
    g <- if (tracker$calls == 1L) c(0.2, 0.8) else c(1 / 3, 2 / 3)
    data.frame(
      theta = tau,
      g = g,
      cumul = cumsum(g)
    )
  }

  unlockBinding(".estimate_mixture_npmle_from_binomial", asNamespace("mcbarycenter"))
  assign(".estimate_mixture_npmle_from_binomial", mock_fun, envir = asNamespace("mcbarycenter"))
  lockBinding(".estimate_mixture_npmle_from_binomial", asNamespace("mcbarycenter"))

  on.exit({
    unlockBinding(".estimate_mixture_npmle_from_binomial", asNamespace("mcbarycenter"))
    assign(".estimate_mixture_npmle_from_binomial", old_fun, envir = asNamespace("mcbarycenter"))
    lockBinding(".estimate_mixture_npmle_from_binomial", asNamespace("mcbarycenter"))
  }, add = TRUE)

  tau <- c(0.25, 0.75)
  .estimate_all_mixtures_from_precomputed(
    precomputed = precomputed,
    method = "npmle",
    tau = tau
  )

  expect_equal(tracker$calls, 2L)
  expect_null(tracker$starts[[1]])
  expect_equal(tracker$starts[[2]], c(0.2, 0.8))
})
