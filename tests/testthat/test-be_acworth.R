context("test-be_acworth")

test_that("be_acworth works", {
  
  # generate dataset
  require(earthtide)
  
  dt <- 3600
  ndays <- 1.9322736 * 363*2+1
  be <- 0.3
  s <- as.POSIXct('2000-01-01', tz = 'UTC')
  e <- as.POSIXct('2000-01-01', tz = 'UTC') + 86400 * ndays - 1
  x <- seq.POSIXt(s, e, dt)
  
  baro_1 <- sin(seq(0, 2*pi*ndays, length.out = length(x))) * 0.2
  baro_2 <- sin(seq(0, 2*pi*ndays*2, length.out = length(x))) * 0.1
  baro_3 <- sin(seq(0, 2*pi*ndays*1.9322736, length.out = length(x))) * 0.01
  ba <- baro_1 + baro_2 + baro_3
  
  
  dat <- data.table(ba, datetime = x)
  
  wg <- as.data.table(earthtide::eterna_wavegroups)
  wg <- wg[time == 'all', list(start, end)]
  dat[, et := calc_earthtide(datetime, 
                             latitude = 49.00937,
                             longitude = 8.40444,
                             astro_update = 600,
                             wave_groups = wg)[,2] * 3e-5]
  dat[, wl := -be * ba + et ]
  
  be_test <- be_acworth(dat, wl = 'wl', ba = 'ba', et = 'et', 
                        method = 'spec_pgram',
                        spans = c(3, 3, 3, 3, 3), inverse = TRUE)
  
  
  expect_equal(be_test, be, tolerance = 1e-2)
  
  be_test <- be_acworth(dat, wl = 'wl', ba = 'ba', et = 'et', 
                        method = 'spec_welch', window = window_hann,
                        n_subsets = 20, inverse = TRUE)
  
  
  expect_equal(be_test, be, tolerance = 1e-2)
  
  
  # check inverse
  dat[, wl := be * ba + et ]
  
  be_test <- be_acworth(dat, wl = 'wl', ba = 'ba', et = 'et', 
                        method = 'spec_pgram',
                        spans = c(3, 3, 3, 3, 3, 3), inverse = FALSE)
  
  
  expect_equal(be_test, be, tolerance = 1e-2)
  
  be_test <- be_acworth(dat, wl = 'wl', ba = 'ba', et = 'et', 
                        method = 'spec_welch', window = window_hann,
                        n_subsets = 50, inverse = TRUE)
  
  expect_equal(be_test, be, tolerance = 1e-2)
  
  
  expect_error(be_acworth(dat, wl = 'wl', ba = 'ba', et = 'et', 
                          method = 'spectrum', window = window_hann,
                          n_subsets = 50, inverse = TRUE), 
               regexp = 'spectrum method not yet implemented') 
  
  # plot(dat$et, type='l', xlim = c(0, 500))
  # plot(dat$ba, type='l', xlim = c(0, 500))
  # plot(dat$wl, type='l', xlim = c(0, 500))
  
})
