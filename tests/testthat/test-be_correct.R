context("be_correct")


test_that("be_correct works", {
  
  baro <- rnorm(1e3)
  wl <- (baro * 0.4) + 18
  dat <- data.table(baro, wl)
  expect_equal(dat$wl[2] + be_correct(dat, be = 0.4, inverse = FALSE)[2], 18, tolerance = 0.1)
  
  wl <- -baro * 0.4 + 18
  dat <- data.table(baro, wl)
  expect_equal(dat$wl[2] + be_correct(dat, be = 0.4, inverse = TRUE)[2], 18, tolerance = 0.1)
  
  # known mean is 0
  expect_equal(dat$wl[2] + be_correct(dat, be = 0.4, inverse = TRUE, known_mean = 0)[2], 18)
  
})


