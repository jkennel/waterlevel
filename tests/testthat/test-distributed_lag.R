test_that("distributed_lag works", {
  
  library(splines)
  
  max_lag <- 28
  knots <- log_lags(6, 28)
  max_knot <- max(knots)
  x     <- 1:(60)
  
  # generate basis functions
  max_knot <- max(knots)
  n_knots  <- length(knots)
  
  
  # generate basis lag
  basis_lag <- ns(min(knots):max_knot,
                  knots = knots[-c(1, n_knots)],
                  Boundary.knots = range(knots),
                  intercept = TRUE)
  for (i in 1:(max_lag - 1)) {
    for (j in 0:(i-1)){
      n_sub <- i
      n_shift <- j
      
      dist_lag_par <- distributed_lag_parallel(rev(x),
                                                 t(as.matrix(basis_lag)),
                                                 max_knot - min(knots),
                                                 n_subset = n_sub,
                                                 n_shift = n_shift)
      dist_lag_fft <- cross_basis_fft(as.matrix(x), basis_lag)[seq(1, length(x)-n_shift, n_sub)+n_shift,]
      expect_equal(dist_lag_par, dist_lag_fft)
      
    }
  }
  
  expect_error(distributed_lag_parallel(rev(x),
                                        t(as.matrix(basis_lag)),
                                        max_knot - min(knots), 
                                        n_subset = 1, 
                                        n_shift = -1), class = 'std::range_error')
  expect_error(distributed_lag_parallel(rev(x),
                                        t(as.matrix(basis_lag)),
                                        max_knot - min(knots), 
                                        n_subset = 501, 
                                        n_shift = 500), class = 'std::range_error')
  expect_error(distributed_lag_parallel(rev(x),
                                        t(as.matrix(basis_lag)),
                                        max_knot - min(knots), 
                                        n_subset = 10, 
                                        n_shift = 500), class = 'std::range_error')
  
  expect_error(distributed_lag(1:100, knots = c(1,10,101)),
    regexp = 'The maximum knot cannot be larger than the number of elements in x')
  
  expect_equal(colnames(distributed_lag(1:100, knots = c(1,10), lag_name = 'baro')), 
                     c('distributed_lag_baro_1', 'distributed_lag_baro_10'))
  
  expect_equal(sum(is.na(distributed_lag(c(NA, 1:50, NA, 1:50), knots = c(1,10)))),
               40L)
  
  expect_equal(is.matrix(distributed_lag(1:10000, knots = c(1,6000))), TRUE)
  
})
