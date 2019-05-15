test_that("distributed_lag works", {
  
  library(splines)
  
  knots <- seq(0, 10, 2)
  max_knot <- max(knots)
  x     <- 1:(1000)
  
  # generate basis functions
  max_knot <- max(knots)
  n_knots  <- length(knots)
  
  
  # generate basis lag
  basis_lag <- ns(min(knots):max_knot,
                          knots = knots[-c(1, n_knots)],
                          Boundary.knots = range(knots),
                          intercept = TRUE)
  
  dist_lag_par <- distributed_lag_parallel(rev(x), 
                                           t(as.matrix(basis_lag)), 
                                           max_knot-min(knots) )
  
  dist_lag_fft <- cross_basis_fft(as.matrix(x), basis_lag)
  
  expect_equal(dist_lag_par, dist_lag_fft)
  
  expect_error(distributed_lag(1:100, knots = c(1,10,101)),
    regexp = 'The maximum knot cannot be larger than the number of elements in x')
  
  expect_equal(colnames(distributed_lag(1:100, knots = c(1,10), lag_name = 'baro')), 
                     c('distributed_lag_1_baro', 'distributed_lag_10_baro'))
  
  expect_equal(sum(is.na(distributed_lag(c(NA, 1:50, NA, 1:50), knots = c(1,10)))),
               40L)
  
  expect_equal(is.matrix(distributed_lag(1:10000, knots = c(1,6000))), TRUE)
  
})
