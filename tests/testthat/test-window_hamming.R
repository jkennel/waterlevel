context("window_hamming")

test_that("window_hamming works", {
  expect_equal(window_hamming(5), c(0.08, 0.54, 1.00, 0.54, 0.08))
})
