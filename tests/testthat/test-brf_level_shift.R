test_that("brf_level_shift works", {
  
  expect_equal(nlevels(set_level_shift(1:100, c(2, 40, 86))), 4L)
  
})
