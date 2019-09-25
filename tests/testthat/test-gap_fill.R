test_that("gap_fill works", {
  
  #library(waterlevel)
  library(data.table)
  library(tibble)
  library(recipes)
  
  data(transducer)

  transducer[21000:21500, wl := NA_real_]
  transducer[21200:nrow(transducer), wl := wl + 0.01]
  
  # transducer[21600:21700, wl := NA_real_]
  # transducer[21700:nrow(transducer), wl := wl + 0.03]

  tmp <- find_level_shift(transducer, dep_var = 'wl', 
                          time_var = 'datetime', 
                          time_interval = 120L)
  
  transducer[, level_shift := set_level_shift(datetime, tmp$midpoint)]

  ba_lags <- log_lags(12, 86400*2/120)
  et_lags <- seq(-3600, 21600, 1200) / 120
  
  rec <- recipe(wl~., transducer) %>%
    step_distributed_lag(baro, knots = ba_lags) %>%
    step_lag_earthtide(datetime,
                       latitude = 34.23,
                       longitude = -118.67,
                       lag = et_lags,
                       cutoff = 1e-6,
                       astro_update = 300) %>%
    step_mutate(datetime_num = as.numeric(datetime), role = 'datetime') %>%
    step_dummy(level_shift, one_hot = TRUE, role = 'level_shift') %>%
    step_naomit(has_role(match = "lag_earthtide")) %>%
    step_naomit(has_role(match = "distributed_lag")) %>%
    #step_dummy(level_shift, role = 'level_shift') %>%
    step_zv(has_role(match = "level_shift"))
  
  g <- gap_fill(transducer, tmp, rec, time_interval = 120, 
                buffer_start = 86400 * 6, buffer_end = 86400 * 4)
  
  a <- tail(g, 1)$coefs[[1]]

  expect_equal(diff(a[grep('level_shift', name)][['Estimate']]),
               0.01, tolerance = 0.0001, scale = 1.0)
  
  #gap_fill2(tmp, g)
  # s <- get_intercept_stats(g)
  # 
  # expect_equal(tail(s$min, 1), 0.01, tolerance = 0.0001)
  

  
})


