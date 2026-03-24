test_that("equalbary_test returns sup and integrated squared norm results", {
  bary_res1 <- list(
    res = data.frame(
      quantile = c(0, 0.5, 1),
      estimate = c(1, 2, 3)
    ),
    cov = diag(c(0.1, 0.2, 0.3)),
    data = data.frame(
      id = c("a", "a", "b", "b"),
      value = c(1, 2, 3, 4)
    )
  )
  
  bary_res2 <- list(
    res = data.frame(
      quantile = c(0, 0.5, 1),
      estimate = c(1.5, 2, 2.5)
    ),
    cov = diag(c(0.2, 0.1, 0.2)),
    data = data.frame(
      subject = c("u", "u", "v", "v", "w", "w"),
      value = c(1, 1, 2, 2, 3, 3)
    )
  )
  
  out <- equalbary_test(bary_res1, bary_res2, B = 25, seed = 1)
  
  expect_named(out, c("sup_norm", "int_sq_norm", "hotelling"))
  expect_named(out$sup_norm, c("pval", "T_obs", "T_sim"))
  expect_named(out$int_sq_norm, c("pval", "T_obs", "T_sim"))
  expect_named(
    out$hotelling,
    c(
      "pval",
      "T_obs",
      "df",
      "mean_diff",
      "cov_diff",
      "variance_explained",
      "n_components"
    )
  )
  
  expect_length(out$sup_norm$T_sim, 25)
  expect_length(out$int_sq_norm$T_sim, 25)
  expect_true(is.numeric(out$sup_norm$T_obs))
  expect_true(is.numeric(out$int_sq_norm$T_obs))
  expect_true(is.numeric(out$hotelling$T_obs))
  expect_equal(out$hotelling$df, 3)
  expect_equal(out$hotelling$mean_diff, c(-0.5, 0, 0.5))
  expect_equal(out$hotelling$cov_diff, diag(c(0.3, 0.3, 0.5)))
  expect_equal(out$hotelling$variance_explained, 1)
  expect_equal(out$hotelling$n_components, 3)
  expect_gte(out$sup_norm$pval, 0)
  expect_lte(out$sup_norm$pval, 1)
  expect_gte(out$int_sq_norm$pval, 0)
  expect_lte(out$int_sq_norm$pval, 1)
  expect_gte(out$hotelling$pval, 0)
  expect_lte(out$hotelling$pval, 1)
})

test_that("equalbary_test supports PCA truncation for the Hotelling test", {
  bary_res1 <- list(
    res = data.frame(
      quantile = c(0, 0.5, 1),
      estimate = c(1, 2, 3)
    ),
    cov = diag(c(9, 1, 0)),
    data = data.frame(
      id = c("a", "a", "b", "b"),
      val = c(1, 2, 3, 4)
    )
  )
  
  bary_res2 <- list(
    res = data.frame(
      quantile = c(0, 0.5, 1),
      estimate = c(2, 1, 3)
    ),
    cov = matrix(0, nrow = 3, ncol = 3),
    data = data.frame(
      id = c("u", "u", "v", "v"),
      val = c(1, 2, 3, 4)
    )
  )
  
  out <- equalbary_test(
    bary_res1,
    bary_res2,
    B = 10,
    seed = 1,
    variance_explained = 0.9
  )
  
  expect_equal(out$hotelling$df, 1)
  expect_equal(out$hotelling$n_components, 1)
  expect_equal(out$hotelling$variance_explained, 0.9)
})

test_that("equalbary_test requires a shared quantile grid", {
  bary_res1 <- list(
    res = data.frame(quantile = c(0, 0.5, 1), estimate = c(1, 2, 3)),
    cov = diag(3),
    data = data.frame(id = c(1, 1, 2, 2), value = c(1, 2, 3, 4))
  )
  
  bary_res2 <- list(
    res = data.frame(quantile = c(0.25, 0.5, 0.75), estimate = c(1, 2, 3)),
    cov = diag(3),
    data = data.frame(id = c(3, 3, 4, 4), value = c(1, 2, 3, 4))
  )
  
  expect_error(
    equalbary_test(bary_res1, bary_res2, B = 10, seed = 1),
    "must use the same quantile grid"
  )
})
