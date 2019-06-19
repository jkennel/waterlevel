test_that("gap_fill works", {
  
  library(waterlevel)
  library(data.table)

  data(transducer)

  transducer[21000:21500, wl := NA_real_]
  transducer[21200:nrow(transducer), wl := wl + 0.01]
  transducer[, level_shift := factor(1)]
  
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
    step_mutate(datetime_num = as.numeric(datetime)) %>%
    step_mutate(datetime_nu2m = as.numeric(datetime), role = 'datetime') %>%
    step_dummy(level_shift, role = 'level_shift') %>%
    step_zv(has_role(match = "level_shift"))
  
  g <- gap_fill(transducer, rec, time_interval = 120, 
                buffer_start = 86400 * 6, buffer_end = 86400 * 4)
  
  expect_equal(tail(g$level_shift, 1), 0.01, tolerance = 0.0001)
  
})
