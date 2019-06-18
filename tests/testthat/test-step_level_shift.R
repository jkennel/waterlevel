test_that("step_level_shift works", {
  library(recipes)
  library(waterlevel)  

  
  data(transducer)

  transducer[1000:1200, wl := NA_real_]
  
  transducer[20000:21000, wl := NA_real_]

  rec <- recipe(wl ~ .,
                data = transducer[, list(datetime, wl)])

  with_levels <- rec %>%
    step_level_shift(wl, datetime, time_interval = 120L) %>%
    prep() %>%
    juice()
  
  expect_equal(length(unique(with_levels$level_shift)), 3L)
  
})
