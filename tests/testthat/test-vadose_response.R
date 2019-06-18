test_that("vadose_response works", {
  
  time <- seq(0, 43200, 1)
  D <- 0.1
  L <- 68
  precision <- 1e-8
  w <- vadose_response(time, D, L, precision)
  
  expect_gt(tail(w, 1), 0.1)
  expect_lt(head(w, 1), 1.0 + 1e-8)
})
