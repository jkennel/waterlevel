context("test-baro_elev")

test_that("baro_elev works", {
  expect_equal(
    baro_elev(
      elevation    = 2000,
      pressure_sea = 101325,
      temperature  = 20.0), 
    80258.8, tolerance = 0.1)
})
