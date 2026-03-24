test_that("shrink_dist computes raw and shrunk cdf values from x_vec", {
  mcb_obj <- list(
    mixtures = list(
      "1" = data.frame(theta = c(0.2, 0.8), g = c(0.5, 0.5)),
      "2" = data.frame(theta = c(0.3, 0.9), g = c(0.25, 0.75))
    ),
    data = data.frame(
      id = c("a", "a", "b", "b"),
      val = c(1, 2, 3, 4)
    )
  )
  
  out <- shrink_dist(x_vec = c(1, 3), mcb_obj = mcb_obj)
  
  expect_named(out, c("x", "cdf_raw", "cdf_shrunk"))
  expect_equal(out$x, c(1, 2))
  expect_equal(out$cdf_raw, c(0.5, 0.5))
  expect_true(all(out$cdf_shrunk >= 0 & out$cdf_shrunk <= 1))
})

test_that("shrink_dist can pull x_vec from mcb_obj$data using select_id", {
  mcb_obj <- list(
    mixtures = list(
      "1" = data.frame(theta = c(0.2, 0.8), g = c(0.5, 0.5))
    ),
    data = data.frame(
      id = c("a", "a", "b", "b"),
      val = c(1, 2, 3, 4)
    )
  )
  
  out <- shrink_dist(select_id = "a", mcb_obj = mcb_obj)
  
  expect_equal(out$x, 1)
  expect_equal(out$cdf_raw, 0.5)
})
