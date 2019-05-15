test_that("synthetic works", {
  exp1 <- -exp(-seq(0.01, 8, length.out = 60))
  exp2 <- -exp(-seq(0.01, 8, length.out = 43201))
  exp3 <- exp(-seq(0.01, 8, length.out = 600))
  
  exp1 <- exp1 /  sum(exp1)
  exp2 <- exp2 /  sum(exp2)
  exp3 <- exp3 / -sum(exp3)
  
  exp1 <- c(exp1, rep(0.0, length(exp2) - length(exp1)))
  exp3 <- c(exp3, rep(0.0, length(exp2) - length(exp3)))
  kern <- rev(exp1 + exp2 + exp3)
  
  resp_fun <- cumsum(rev(kern))*0.5
  
  dat <- synthetic(n = 1 * 86400,
                   seed = 1,
                   baro_kernel = kern)
  
  # length of the data is 86400 minus the kernel length, 
  # minus the length of ba_fft_sub*2
  expect_equal(nrow(dat), 86400-43200-4000)

})
