test_that("step_earthtide works", {
  library(recipes)
  library(earthtide)

  df <- tibble(x = 1, t = as.POSIXct(1:10, origin = '1970-01-01'))
  
  # lags numeric data
  baked <- recipe(~ ., data = df) %>%
    step_earthtide(t) %>%
    prep() %>%
    portion()
  
  expect_equal(baked$earthtide, 
               calc_earthtide(df$t, do_predict = FALSE, return_matrix = TRUE))
  
  
})
