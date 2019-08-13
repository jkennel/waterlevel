test_that("transfer_fun works", {
  
  library(data.table)
  n <- 86400
  datetime <- as.POSIXct(1:n, origin = '1970-01-01')
  y <- rnorm(n)
  z <- rnorm(n)
  x <- y * 0.2 + z * 0.7
  
  dat <- data.table(x, y, z, datetime)
  
  tfp <- transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
               method = 'spec_pgram', spans = 7)
  
  expect_equal(Re(tfp$gain_x_y[1]), 0.2)
  expect_equal(Re(tfp$gain_x_z[1]), 0.7)
  
  expect_equal(all(abs(Re(tfp$phase_x_y)) < 1e-12), TRUE)
  expect_equal(all(abs(Re(tfp$phase_x_z)) < 1e-12), TRUE)
  
  

  tfw <- transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
                     method = 'spec_welch', n_subsets = 10)
  
  expect_equal(Re(tfw$gain_x_y[1]), 0.2)
  expect_equal(Re(tfw$gain_x_z[1]), 0.7)

  expect_equal(all(abs(Re(tfw$phase_x_y)) < 1e-12), TRUE)
  expect_equal(all(abs(Re(tfw$phase_x_z)) < 1e-12), TRUE)
  
  
  expect_error(transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
                            method = 'spectrum', n_subsets = 10), 
               regexp = 'spectrum method not yet implemented') 
  
  
  
  kern_y <- seq(0, 1, length.out = 100)
  kern_z <- 10^seq(-1, -10, length.out = 100)
  y <- convolve_fft(x, kern_y)
  z <- convolve_fft(x, kern_z)
  dat <- na.omit(data.table(x, y, z, datetime))
  
  tfp <- transfer_fun(dat, vars = c('x', 'y', 'z'), time = 'datetime',
                      method = 'spec_pgram', spans = 3, taper = 0.1)
  

  plot(Re(gain_x_y)~Re(frequency), tfp, type='l')
  points(Re(ba_efficiency)~Re(frequency), tfs, type='l', col = 'red')
  plot(Re(gain_x_z)~Re(frequency), tfp, type='l')
  points(Re(et_efficiency)~Re(frequency+1), tfs, type='l', col = 'red')
  
  
  
  tmp <- spec_pgram(as.matrix(dat[, list(x, y)]), spans = 3)
  rev(tail(tmp[,1,1], 10))
  head(tmp[,1,1], 11)[-1]
  head(tfs$Pxx, 10)
  all.equal(head(tmp[,1,1], 43201)[-1],
            head(tfs$Pxx, 43200))

  
  
})

library(xts)
library()
tfs <- transfer_smooth(dat, 'x', 'y', 'z', 'datetime', spans = 3)

transfer_smooth <- function(dat,
                            wl = 'wl',
                            baro = 'baro',
                            et = 'et',
                            datetime = 'datetime',
                            ...) {
  
  # replace with new method -->
  wl_xts <- xts(as.matrix(dat[, c(wl, baro, et), with = FALSE]), dat[[datetime]])
  spec <- cbind(setDT(cross_spectrum(as.ts(wl_xts[, c(baro, wl)]), ...))[, list(frequency = freq*86400, spec_ba = spec1, spec_wl = spec2, coh_ba_wl = coh, phase_ba_wl = phase, P_ba = Pxx, P_wl = Pyy, P_ba_wl = Pxy)],
                setDT(cross_spectrum(as.ts(wl_xts[, c(et, wl)]), ...))[, list(spec_et = spec1, coh_et_wl = coh, phase_et_wl = phase, P_et = Pxx, P_et_wl = Pxy)],
                setDT(cross_spectrum(as.ts(wl_xts[, c(baro, et)]), ...))[, list(coh_ba_et = coh, phase_ba_et = phase, P_ba_et = Pxy)])
  # <-- replace with new method
  
  spec <- spec[, list(frequency,
                      spec_wl, spec_ba, spec_et,
                      P_wl, P_ba, P_et,
                      P_ba_wl, P_et_wl, P_ba_et,
                      coh_ba_wl, coh_et_wl, coh_ba_et,
                      phase_ba_wl, phase_et_wl, phase_ba_et)]
  
  spec[, denominator := P_ba * P_et - P_ba_et * Conj(P_ba_et)]
  spec[, ba_transfer := (P_et * P_ba_wl - P_ba_et * P_et_wl) / denominator]
  spec[, et_transfer := (P_ba * P_et_wl - Conj(P_ba_et) * P_ba_wl) / denominator]
  spec[, denominator := NULL]
  
  spec[, ba_efficiency := Mod(ba_transfer)]
  spec[, et_efficiency := Mod(et_transfer)]
  
  spec[, ba_phase := atan2(Im(ba_transfer), Re(ba_transfer))]
  spec[, et_phase := atan2(Im(et_transfer), Re(et_transfer))]
  
  setnames(spec,
           c('P_wl', 'P_ba', 'P_et',
             'P_ba_wl', 'P_et_wl', 'P_ba_et',
             'coh_ba_wl', 'coh_et_wl', 'coh_ba_et',
             'phase_ba_wl', 'phase_et_wl', 'phase_ba_et'),
           c('Pxx', 'Pyy', 'Pzz',
             'Pxy', 'Pxz', 'Pyz',
             'coh_xy', 'coh_xz', 'coh_yz',
             'phase_xy', 'phase_xz', 'phase_yz')
  )
  return(spec)
}

