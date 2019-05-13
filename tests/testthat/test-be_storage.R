context("be_storage")

# Batu 1998, pg 76
test_that("be_storage works", {
  expect_equal(be_storage(0.5, 0.32, 45, gamma = 9799.74, beta = 4.786e-10), 0.0001350765)
})
