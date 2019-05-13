context("window_blackman_nuttall")

test_that("window_blackman_harris works", {
  expect_equal(max(abs(window_blackman_harris(64) -
                       window_blackman_nuttall(64))), 0.0099, tolerance = 0.0001)
})
