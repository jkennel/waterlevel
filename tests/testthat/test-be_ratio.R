context("test-be_ratio")

test_that("be_ratio works", {
  library(data.table)
  datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
                         as.POSIXct("2016-01-05 12:00:00"), by='hour' )
  
  be = 0.4
  
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- -be * baro 
  dat <- data.table(baro, wl, datetime)
  expect_equal(be_ratio(dat, quant = 0.5, stat = median), be)
  
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- (1-be) * baro 
  dat <- data.table(baro, wl, datetime)
  expect_equal(be_ratio(dat, 
                        quant = 0.9,
                        stat = mean,
                        inverse = FALSE), 1-be)
  
  
})
