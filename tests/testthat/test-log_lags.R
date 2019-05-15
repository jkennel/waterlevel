test_that("log_lags works", {
  
  expect_equal(log_lags(3, 100), c(1, 10, 100) - 1)
  
  suppressWarnings(expect_equal(log_lags(101, 100), 1:100 - 1))
  
  expect_warning(log_lags(101, 100))
  
})
