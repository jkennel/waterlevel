test_that("predict_terms works", {
  
  # test single
  x <- 1:100
  y <- rnorm(100, sd = 0.01) + x/2
  dat <- data.frame(x, y)
  fit <- lm(y~x, dat)
  pt <- predict_terms(fit)
  pt <- pt + attr(pt,'constant')
  expect_equal(pt, predict(fit), check.attributes = FALSE)
  
  
  # test multiple
  x <- 1:100
  y <- rnorm(100, sd = 0.01) + x/2
  z <- rnorm(100, sd = 0.01) + x/4
  dat <- data.frame(x, y)
  fit <- lm(cbind(y,z)~x, dat)

  pt <- predict_terms(fit)
  pt <- pt + matrix(attr(pt,'constant'), byrow = TRUE, ncol=2, nrow = nrow(pt))
  
  expect_equal(pt, predict(fit),
               check.attributes = FALSE)
  
})
