test_that("pad_kernel works", {
  expect_equal(length(pad_kernel(kernel("modified.daniell", 3 %/% 2), 20)), 20)
  expect_equal(sum(pad_kernel(kernel("modified.daniell", 3 %/% 2), 20) == 0), 17)
  expect_error(pad_kernel(kernel("modified.daniell", rep(3 %/% 2, 20)), 10), 
               regexp = "x is shorter than kernel k")
  expect_error(pad_kernel(1:10, 10), 
               regexp = "'k' is not a kernel")
  
})

test_that("spec_pgram works", {
  expect_error(spec_pgram(1:100, kernel = 2))
  
  x <- rnorm(100)
  taper <- 0
  expect_equal(as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = TRUE, 
                                        taper = 0,
                                        fast = TRUE)))[2:51],
               spec.pgram(x, 
                          demean = TRUE, 
                          detrend = TRUE,
                          plot = FALSE,
                          taper = 0)$spec)
  
  expect_equal(as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = FALSE, 
                                        taper = 0)))[2:51],
               spec.pgram(x, 
                          demean = TRUE, 
                          detrend = FALSE,
                          plot = FALSE,
                          taper = 0)$spec)
  
  expect_equal(as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = TRUE, 
                                        taper = 0, 
                                        pad = 100)))[2:101],
               spec.pgram(x, 
                          demean = TRUE, 
                          detrend = TRUE,
                          plot = FALSE,
                          taper = 0,
                          pad = 1)$spec)
  
  
  x <- rnorm(2^17)
  taper <- 0
  
  expect_equal(as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = TRUE, 
                                        taper = 0)))[2:(length(x)/2+1)],
               spec.pgram(x, 
                          demean = TRUE, 
                          detrend = TRUE,
                          plot = FALSE,
                          taper = 0)$spec)
  
  
  expect_equal(as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = TRUE, 
                                        taper = 0,
                                        spans = 3)))[2:(length(x)/2+1)],
               spec.pgram(x, 
                          demean = TRUE, 
                          detrend = TRUE,
                          plot = FALSE,
                          taper = 0,
                          spans = 3)$spec)

  
  
  expect_equal(as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = TRUE, 
                                        taper = 0,
                                        kernel = kernel("modified.daniell", 3 %/% 2)))),
               as.numeric(Re(spec_pgram(x, 
                                        demean = TRUE, 
                                        detrend = TRUE, 
                                        taper = 0,
                                        spans = 3))))
  

})

