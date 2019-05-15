test_that("step_lag_earthtide works", {

  library(recipes)
  library(earthtide)
  
  df <- tibble(x = 1, t = as.POSIXct(1:10, origin = '1970-01-01'))
  
  # lags numeric data
  baked <- recipe(~ ., data = df) %>%
    step_lag_earthtide(t, lag = 1) %>%
    prep() %>%
    portion()
  
  expect_equal(baked$lag_earthtide, 
               lag(calc_earthtide(df$t, 
                                  do_predict = TRUE, 
                                  return_matrix = TRUE), 1),
                check.attributes = FALSE)
  
  
})
