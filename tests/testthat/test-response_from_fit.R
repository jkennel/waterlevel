test_that("response_from_fit works", {
  library(earthtide)
  library(data.table)
  library(recipes)
  
  exp1 <- -exp(-seq(0.01, 8, length.out = 60))
  exp2 <- -exp(-seq(0.01, 8, length.out = 43201))
  exp3 <- exp(-seq(0.01, 8, length.out = 600))
  
  exp1 <- exp1 /  sum(exp1)
  exp2 <- exp2 /  sum(exp2)
  exp3 <- exp3 / -sum(exp3)
  
  exp1 <- c(exp1, rep(0.0, length(exp2) - length(exp1)))
  exp3 <- c(exp3, rep(0.0, length(exp2) - length(exp3)))
  kern <- rev(exp1 + exp2 + exp3)
  
  resp_fun <- cumsum(rev(kern))*0.5
 
  dat <- synthetic(sd_noise = 0e-2,
                   sd_noise_trend = 0e-12,
                   n = 2 * 86400,
                   linear_trend = 0e-12,
                   seed = 1,
                   baro_kernel = kern)[, list(datetime, baro, wl)]
  
  
  dt <- 60
  dat_sub <- dat[seq(1, nrow(dat), dt)]
  
  ba_lags <- seq(0, (43200+dt) / dt, 1)
  
  rec <- recipe(wl~., dat_sub) %>%                          
    step_lag_matrix(baro, lag = ba_lags) %>%                   
    prep() %>%                                                
    portion()
  
  
  fit <- lm(outcome~lag_matrix, rec,
            x = FALSE, y = FALSE, 
            tol = 1e-50, 
            na.action = na.exclude)
  
  
  resp <- response_from_fit(fit)
  expect_equal(tail(resp$value, 1), 0.5, tolerance = 1e-3)
  
  ba_lags <- log_lags(10, 43200)
  
  rec <- recipe(wl~., dat) %>%                          
    step_distributed_lag(baro, knots = ba_lags) %>%                   
    prep() %>%                                                
    portion()
  
  
  fit <- lm(outcome~distributed_lag, rec,
            x = FALSE, y = FALSE, 
            tol = 1e-50, 
            na.action = na.exclude)
  
  
  resp <- response_from_fit(fit)
  expect_equal(tail(resp$value, 1), 0.5, tolerance = 1e-3)
  
  
  # test harmonic
  library(tibble)
  library(data.table)
  
  x1 <- seq(0, 8 * pi, length.out = 100)
  x2 <- seq(0, 4 * pi, length.out = 100) 
  
  y <- sin(x1 + pi/6)* 0.2 + sin(x2 + pi/8) * 0.7
  

  df <- as_tibble(list(y = as.matrix(y),
                       harmonic_1 = cbind(sin_1 = sin(x1), cos_1 = cos(x1)),
                       harmonic_2 = cbind(sin_2 = sin(x2), cos_2 = cos(x2))))
  
  fit_2 <- lm(y~harmonic_1 + harmonic_2-1, df)
  resp <- response_from_fit(fit_2)
  
  expect_equal(resp[type == 'amplitude']$value, c(0.2, 0.7))
  
  
  
  library(tibble)
  library(data.table)
  library(recipes)
  data(transducer)

  # Select wavegroups
  wave_groups <- data.table(start = 1.9322736-1e-6, end = 1.9322736+1e-6)
  ns_df <- 7
  delta_t <- 120                        
  max_lag_baro <- 43200/delta_t 
  knots <- log_lags(15, max_lag_baro)
  
  dat_e <- recipe(wl~., transducer) %>%                          
    step_distributed_lag(baro, knots = knots) %>%   
    step_earthtide(datetime,
                   longitude = -118.67,
                   latitude = 34.23,
                   elevation = 550,
                   wave_groups = wave_groups, 
                   astro_update = 60) %>%      
    step_mutate(datetime_num = as.numeric(datetime)) %>%           
    step_ns(datetime_num, deg_free = ns_df, role = 'splines') %>%                        
    prep() %>%                                                 
    portion()   
  
  
  knots <- log_lags(15, max_lag_baro)
  dat_h <- recipe(wl~., transducer) %>%                          
    step_distributed_lag(baro, knots = knots) %>%   
    step_harmonic(datetime,freq = 1.9322736) %>%      
    step_mutate(datetime_num = as.numeric(datetime)) %>%           
    step_ns(datetime_num, deg_free = ns_df, role = 'splines') %>%                        
    prep() %>%                                                 
    portion()   
  

  fit_e <- lm(outcome~splines + earthtide + distributed_lag,  dat_e)
  fit_h <- lm(outcome~splines + harmonic + distributed_lag,  dat_h)

  resp_e <- response_from_fit(fit_e)
  resp_h <- response_from_fit(fit_h)
  
  # check the M2 signal amplitude with two methods
  expect_equal(resp_e[type=='amplitude']$value, 
               resp_h[type=='amplitude']$value,
               tol = 1e-7)

  
})
