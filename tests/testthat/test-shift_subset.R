test_that("shift_subset works", {
  
  expect_equal(check_lag(7, 2, 3), -1)
  expect_equal(check_lag(7, 3, 3), 0)
  expect_equal(check_lag(7, 3, -3), 6)
  expect_error(check_lag(7, 3, -5))
  
  expect_equal(length(shift_subset(1:7, 3, 1, 0)), 7)
  expect_equal(length(shift_subset(1:7, 3, 2, 0)), 4)
  expect_equal(length(shift_subset(1:6, 3, 2, 0)), 3)
  
  expect_equal(na.omit(shift_subset(1:6, 3, 1, 0)), 1:3, check.attributes = FALSE)
  expect_equal(na.omit(shift_subset(1:6, -3, 1, 0)), 4:6, check.attributes = FALSE)
  
  expect_equal(na.omit(shift_subset(1:6, 3, 2, 0)), 2, check.attributes = FALSE)
  
  expect_equal(na.omit(shift_subset(1:6, -3, 2, 0)), c(4,6), check.attributes = FALSE)
  expect_equal(na.omit(shift_subset(1:7, -3, 2, 0)), c(4,6), check.attributes = FALSE)
  expect_equal(na.omit(shift_subset(1:7, -3, 3, 0)), c(4, 7), check.attributes = FALSE)
  
  expect_equal(na.omit(shift_subset(1:7, -3, 2, 1)), c(5, 7), check.attributes = FALSE)
  expect_equal(na.omit(shift_subset(1:7, 3, 2, 1)), c(1, 3, 5), check.attributes = FALSE)
  
})
