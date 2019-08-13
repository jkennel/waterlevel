test_that("convolve_fft works", {
  
  y <- rep(0.25, 4)
  x <- rnorm(100)
  
  expect_equal(waterlevel::convolve_fft(x,y), 
               as.numeric(stats::filter(x=x, filter = y, sides = 1)))
  
  
})
