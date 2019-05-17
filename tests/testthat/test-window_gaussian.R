test_that("window_gaussian works", {
  expect_equal(length(window_gaussian(100, 0.3)), 100)
})
