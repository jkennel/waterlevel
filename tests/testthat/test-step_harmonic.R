test_that("step_harmonic works", {

  library(recipes)
  library(earthtide)
  
  df <- tibble(x = 1, t = as.POSIXct(1:10, origin = '1970-01-01'))
  freqs <- c(1,2,3)
  # lags numeric data
  baked <- recipe(~ ., data = df) %>%
    step_harmonic(t, freq = freqs) %>%
    prep() %>%
    portion()
  
  expect_equal(baked$harmonic, 
               harmonic(df$t, freq = freqs),
               check.attributes = FALSE)
  
})
