context("window_blackman")

test_that("window_blackman works", {
  expect_equal(window_blackman(5),
               c(-1.387779e-17,  0.34,  1.0,  0.34, -1.387779e-17))
})
