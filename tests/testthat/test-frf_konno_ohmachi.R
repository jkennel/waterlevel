test_that("konno_ohmachi works", {
  
  n <- 3000
  x_in <- rnorm(n)
  x <- x_in

  f <- 1:n
  a <- konno_ohmachi_parallel(x, f, 20)
  
  x <- x_in
  f <- 1:n
  
  b <- konno_ohmachi_serial(x, f, 20)
  expect_equal(a, b)

})
