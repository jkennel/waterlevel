context("test-spec_welch")

test_that("spec_welch recovers input values", {
  library(data.table)
  dt <- 7200
  
  s <- as.POSIXct('2017-01-01', tz = 'UTC')
  e <- as.POSIXct('2017-01-01', tz = 'UTC') + 86400 * 3 * 5000-1
  x <- seq.POSIXt(s, e, dt)
  
  y <- harmonic(x, freq = c(1, 1.93, 2, 3))

  dat <- data.table(cbind(x, y[, 1:4]))
  names(dat) <- c('datetime', 'harm1', 'harm2', 'harm3', 'harm4')
  dat[, harm := harm1 + harm2 * 2 + harm3 * 0.25 + harm4 * 0.1]
  
  sp <- as.numeric(Re(spec_welch(as.matrix(dat$harm), 
                              window = window_hann)))
  
  # plot(sqrt(sp), type='l')
  
  n <- (length(sp))
  frequency <- seq.int(from = 1/n, by = 1/n, length.out = length(sp)) * 86400/dt
 
  wh1 <- which.min(abs(1-frequency))+1
  wh19 <- which.min(abs(1.93-frequency))+1
  wh2 <- which.min(abs(2-frequency))+1
  wh3 <- which.min(abs(3-frequency))+1
  expect_equal(sqrt(sp[wh1]), sqrt(sp[wh19])/2, tolerance = 1e-5, scale = 1.0)
  expect_equal(sqrt(sp[wh1])/4, sqrt(sp[wh2]), tolerance = 1e-5, scale = 1.0)
  expect_equal(sqrt(sp[wh1])/10, sqrt(sp[wh3]), tolerance = 1e-5, scale = 1.0)
  


  # # check that different dt gives similar results
  # dt <- 3600
  # s <- as.POSIXct('2017-01-01', tz = 'UTC')
  # e <- as.POSIXct('2017-01-01', tz = 'UTC') + 86400 * 3 * 5000-1
  # x <- seq.POSIXt(s, e, dt)
  # 
  # y <- harmonic(x, freq = c(1, 1.93, 2, 3))
  # 
  # dat <- data.table(cbind(x, y[, 1:4]))
  # names(dat) <- c('datetime', 'harm1', 'harm2', 'harm3', 'harm4')
  # dat[, harm := harm1 + harm2 * 2 + harm3 * 0.25 + harm4 * 0.1]
  # 
  # sp2 <- as.numeric(spec_welch(Re(as.matrix(dat$harm)), 
  #                              window = window_rectangular))
  # n <- 2 * (length(sp2))
  # frequency <- seq.int(from = 1/n, by = 1/n, length.out = length(sp2)) * 86400/dt
  # 
  # wh2 <- which.min(abs(1-frequency))
  # expect_equal(sqrt(sp[wh1]), sqrt(sp2[wh2]), tolerance = 1e-5)
  
  
  sp3 <- as.numeric(spec_welch(Re(as.matrix(dat$harm)), n_subsets = 10,
                               window = window_rectangular))
  n <- 2 * (length(sp3))
  frequency <- seq.int(from = 1/n, by = 1/n, length.out = length(sp3)) * 86400/dt
  
  wh3 <- which.min(abs(1-frequency))
  
  
  sp4 <- as.numeric(spec_welch(Re(as.matrix(dat$harm)), n_subsets = 5, 
                               window = window_rectangular))
  n <- 2 * (length(sp4))
  frequency <- seq.int(from = 1/n, by = 1/n, length.out = length(sp4)) * 86400/dt
  wh4 <- which.min(abs(1-frequency))
  
  expect_equal(sqrt(sp4[wh4]), sqrt(sp3[wh3]), tolerance = 1e-5)
  
  
  sp1 <- as.numeric(spec_welch(Re(as.matrix(dat$harm)),
                              window = window_rectangular),
                   detrend = FALSE, demean = TRUE)

  sp2 <- as.numeric(spec_welch(Re(as.matrix(dat$harm)),
                              window = window_rectangular),
                   detrend = FALSE, demean = FALSE)
  
  expect_equal(sp1, sp2)
  
  })


