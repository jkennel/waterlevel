test_that("konno_ohmachi works", {

  library(stats)
  library(RcppParallel)
  library(Rcpp)
  library(RcppArmadillo)

  n <- 2000
  x <- rnorm(n)
  f <- 1:n

  expect_equal(is.numeric(konno_ohmachi_serial(x, f, 10)), TRUE)
  expect_equal(is.numeric(konno_ohmachi_parallel(x, f, 10)), TRUE)

  # currently a problem when run with covr
  expect_equal(konno_ohmachi_parallel(x, f, 10),
               konno_ohmachi_serial(x, f, 10))
  
})

