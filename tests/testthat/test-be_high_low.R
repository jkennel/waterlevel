context("test-be_high_low")

test_that("be_high_low works", {
  library(data.table)
  datetime <- seq.POSIXt(as.POSIXct("2016-01-01 00:00:00", tz = 'UTC'),
                         as.POSIXct("2016-01-06 00:00:00", tz = 'UTC'), by='hour')
  
  be = 0.4
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- -be * baro
  dat <- data.table(baro, wl, datetime)
  expect_equal(be_high_low(dat)$ratio,
               be_high_low(dat, time_group = 86400)$ratio)
  
  expect_equal(nrow(be_high_low(dat)), 6)
  
  #  expect_equal(be_high_low(dat)$ratio, be)
  
  
})
