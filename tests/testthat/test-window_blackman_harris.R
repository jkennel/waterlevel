context("window_blackman_harris")

test_that("window_blackman_harris works", {
  expect_equal(window_blackman_harris(16)[1:8],
               c(0.0001, 0.0036, 0.0267, 0.1030, 0.2680, 0.5206, 0.7938, 0.9749),
               tolerance = 0.0001)
})
