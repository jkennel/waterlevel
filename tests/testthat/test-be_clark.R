context("test-be_clark")

test_that("be_clark works", {
  datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00", tz = 'UTC'),
                         as.POSIXct("2016-01-05 12:00:00", tz = 'UTC'), by='sec' )
  
  be = 0.4
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- -be * baro
  dat <- data.table::data.table(baro, wl, datetime)
  expect_equal(be_clark(dat, dep='wl', ind='baro', lag_space=1, inverse=TRUE), be, check.attributes = FALSE)
  
  
  baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
  wl <- (1 - be) * baro
  dat <- data.table::data.table(baro, wl, datetime)
  
  expect_equal(be_clark(dat, dep='wl', ind='baro', lag_space=1, inverse=FALSE), 1-be, check.attributes = FALSE)
  expect_equal(coefficients(be_clark(dat, dep='wl', ind='baro', lag_space=1, return_model=TRUE,inverse=FALSE))[2], 1-be, check.attributes = FALSE)
  
  
})
