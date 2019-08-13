test_that("step_distributed_lag works", {
  
  library(recipes)

  df <- tibble(x = 1, t = 1:100)
  
  # lags numeric data
  baked <- recipe(~ ., data = df) %>%
    step_distributed_lag(t, knots = c(1, 3, 5)) %>%
    prep() %>%
    portion()
  
  expect_equal(baked$distributed_lag, 
               distributed_lag(df$t, knots = c(1,3,5)),
               check.attributes = FALSE)
  
  
  
  # lags numeric data
  baked <- recipe(~ ., data = df) %>%
    step_distributed_lag(t, knots = c(1, 3, 5), n_subset = 2) %>%
    prep() %>%
    portion()
  
  expect_equal(baked$distributed_lag, 
               distributed_lag(df$t, knots = c(1,3,5), n_subset = 2),
               check.attributes = FALSE)
  

  

})
