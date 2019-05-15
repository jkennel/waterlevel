test_that("konno_ohmachi works", {
  
  n <- 3000
  x_in <- rnorm(n)
  x <- x_in
  ref_z <- seq(0.5, 2.0, 0.0005)

  f <- 1:n
  a <- konno_ohmachi_parallel(x, f, 20)
  
  x <- x_in
  ref_z <- seq(0.5, 2.0, 0.0005)
  f <- 1:n
  
  b <- konno_ohmachi_serial(x, f, 20)
  expect_equal(a, b)

})
