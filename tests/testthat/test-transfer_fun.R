test_that("transfer_fun works", {
  
  library(data.table)
  n <- 86400
  datetime <- as.POSIXct(1:n, origin = '1970-01-01')
  y <- rnorm(n)
  z <- rnorm(n)
  x <- y * 0.2 + z * 0.7
  dat <- data.table(x, y, z, datetime)
  
  tfp <- transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
               method = 'spec_pgram', spans = 3)
  
  expect_equal(Re(tfp$gain_x_y[1]), 0.2)
  expect_equal(Re(tfp$gain_x_z[1]), 0.7)
  
  
  tfw <- transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
                     method = 'spec_welch', n_subsets = 10)
  
  expect_equal(Re(tfw$gain_x_y[1]), 0.2)
  expect_equal(Re(tfw$gain_x_z[1]), 0.7)
  
  
  
  expect_error(transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
                            method = 'spectrum', n_subsets = 10), 
               regexp = 'spectrum method not yet implemented') 
  

  
})
