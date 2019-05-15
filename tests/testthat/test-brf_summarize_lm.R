test_that("brf_summarize_lm works", {
  
  x <- 1:100
  y <- rnorm(100, sd = 0.01) + x/2
  fit <- lm(y~x)
  tmp <- summarize_lm(fit)
  
  expect_equal(tmp$r_squared, summary(fit)$r.squared)
  expect_equal(tmp$sigma, summary(fit)$sigma)
  expect_equal(tmp$AIC, AIC(fit))
  
})
