test_that("portion works", {
  library(recipes)

  df <- tibble(x = 1, t = as.POSIXct(1:10, origin = '1970-01-01'))
  freqs <- c(1,2)
  
  # lags numeric data
  expect_error({recipe(~ ., data = df) %>%
    step_harmonic(t, freq = freqs) %>%
    portion()})
  
  expect_error({recipe(~ ., data = df) %>%
      step_harmonic(t, freq = freqs) %>%
      prep(retain = FALSE) %>%
      portion()})
  
})

