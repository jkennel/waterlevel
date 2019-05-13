context("test-be_least_squares")

test_that("be_least_squares works", {
  library(data.table)
  datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00", tz = 'UTC'),
                         as.POSIXct("2016-01-05 12:00:00", tz = 'UTC'), by='sec' )
  
  be = 0.4
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- -be * baro
  dat <- data.table(baro, wl, datetime)
  expect_equal(be_least_squares(dat, dep='wl', ind='baro', inverse=TRUE), be, check.attributes = FALSE)
  
  
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- (1 - be) * baro
  dat <- data.table(baro, wl, datetime)
  
  expect_equal(be_least_squares(dat, dep='wl', ind='baro', inverse=FALSE), 1-be, check.attributes = FALSE)
  
  expect_equal(coefficients(be_least_squares(dat, dep='wl', ind='baro', return_model = TRUE, inverse=FALSE))[2], 1-be, check.attributes = FALSE)

})
