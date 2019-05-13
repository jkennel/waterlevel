test_that("step_lag_matrix works", {
  library(recipes)
  library(testthat)
  library(dplyr)
  library(rlang)
  
  context("step_lag")
  
  n <- 10
  start <- 1
  end <- 10
  
  set.seed(1)
  
  test_that("single feature",  {
    
    set.seed(27)
    df <- tibble(x = rnorm(n), t = 1:10)
    
    # lags numeric data
    baked <- recipe(~ ., data = df) %>%
      step_lag_matrix(t, lag = 2) %>%
      prep(df) %>%
      bake(df)
    expected <- mutate(df, lag_2_t = as.numeric(dplyr::lag(t, 2)))
    expect_equal(as.matrix(baked), as.matrix(expected),
                 check.attributes = FALSE)
    
    # lags date data
    baked <- recipe(~ ., data = df) %>%
      step_lag(t, lag = 2) %>%
      prep(df) %>%
      bake(df)
    expected <- mutate(df, lag_2_t = dplyr::lag(t, 2))
    expect_equal(baked, expected)
    
    
    # errors out on non-integer lag
    expect_error({
      prepped_rec <- recipe(~ ., data = df) %>%
        step_lag_matrix(x, lag = 0.5) %>%
        prep(df)
    })
  })
  
  test_that("multiple lags",  {
    
    set.seed(27)
    df <- tibble(x = rnorm(n), t = 1:10)
    
    baked_lag <- recipe(~ ., data = df) %>%
      step_lag(t, lag = c(1, 2)) %>%
      prep(df) %>%
      bake(df)
    baked_lag_matrix <- recipe(~ ., data = df) %>%
      step_lag_matrix(t, lag = c(1, 2)) %>%
      prep(df) %>%
      bake(df)
    
    expect_equal(as.matrix(baked_lag), as.matrix(baked_lag_matrix), check.attributes = FALSE)
  })
  
  test_that("negative lag",  {
    
    set.seed(27)
    n <- 10
    df <- tibble(x = rnorm(n), t = 1:n)
    
      
    baked_lag_matrix <- recipe(~ ., data = df) %>%
      step_lag_matrix(t, lag = -1) %>%
      prep(df) %>%
      bake(df)
    
    df$t <- rev(df$t)

    baked_lag <- recipe(~ ., data = df) %>%
      step_lag(t, lag = 1) %>%
      prep(df) %>%
      bake(df)

    baked_lag$lag_1_t <- rev(baked_lag$lag_1_t)
    expect_equal(ncol(baked_lag_matrix), 3)
    expect_equal(baked_lag$lag_1_t, baked_lag_matrix$lag_matrix_n1_t)
    
    
  })

  
})
