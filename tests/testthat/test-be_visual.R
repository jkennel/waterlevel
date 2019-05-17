context("be_visual_plot")

library(data.table)

be <- 0.43
x <- seq(0, 28*pi, pi / (12*12))

baro <- sin(x) + rnorm(length(x), sd = 0.04)
wl <- -sin(x) * be + rnorm(length(x), sd = 0.04)
dat <- data.table(datetime = as.POSIXct(x * 86400 / (2 * pi),
                                        origin = '1970-01-01', tz = 'UTC'),
                  wl = wl, baro = baro)



dat_be <- be_visual_data(dat)


test_that("be_visual runs", {
  
  expect_equal(class(be_visual(dat)), 
               c('plotly', 'htmlwidget'))
})

test_that("be_visual_plot runs", {
  expect_equal(class(be_visual_plot(dat_be)), c('plotly', 'htmlwidget'))
  
  expect_equal(class(be_visual_plot(dat_be, time = 'datetime')), 
               c('plotly', 'htmlwidget'))
  
  dat_be[, time := datetime]
  expect_equal(class(be_visual_plot(dat_be, time = 'time')), 
               c('plotly', 'htmlwidget'))
  
  dat_be[, other := datetime]
  expect_equal(class(be_visual_plot(dat_be, time = 'other')), 
               c('plotly', 'htmlwidget'))
  
  dat_be <- be_visual_data(dat, be_tests = seq(0, 1, 0.05))
  
  expect_equal(class(be_visual_plot(dat_be)), c('plotly', 'htmlwidget'))
  
})


context("be_visual_data")
library(data.table)
be <- 0.4
x <- seq(0, 28*pi, pi / (12*12))

baro <- sin(x) + rnorm(length(x), sd = 0.04)
wl <- -sin(x) * be + rnorm(length(x), sd = 0.04)
dat <- data.table(datetime = as.POSIXct(x * 86400 / (2 * pi),
                                        origin = '1970-01-01', tz = 'UTC'),
                  wl = wl, baro = baro)


test_that("be_visual_data works", {
  be_dat <- be_visual_data(dat)
  be_cor <- be_dat[, list(sse = sum(corrected^2)), by = be]
  
  expect_equal(be_cor[be==0.4]$sse, min(be_cor$sse))
})
