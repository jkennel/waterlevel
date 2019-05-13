context("window_hann")

test_that("window_hann works", {
  expect_equal(window_hann(5), c(0.0, 0.5, 1.0, 0.5, 0.0))
})
