test_that("brf_level_shift works", {
  
  expect_equal(nlevels(level_shift(1:100, c(2, 40, 86))), 4L)
  
})
