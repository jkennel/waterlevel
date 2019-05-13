#' synthetic
#'
#' This function is used for testing purposes
#'
#' @param sd_noise standard deviation of random noise to add (numeric)
#' @param sd_noise_trend standard deviation of noise to add to generate a 
#' trend (numeric)
#' @param linear_trend magnitude of linear trend in time (numeric)
#' @param n length of time series in seconds (integer)
#' @param seed random number seed for reproducibility (numeric)
#' @param baro_kernel vector values to convolve with barometric pressure
#'
#' @return data.table of results
#' @export
#' 
#' @importFrom stats rnorm
#'
synthetic <- function(sd_noise = 0.0002,
                      sd_noise_trend = 0.00003,
                      linear_trend = 1e-7,
                      n = 14 * 86400,
                      seed = NULL,
                      baro_kernel = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  datetime <- as.POSIXct('2016-08-24 08:46:40', tz = 'UTC') + 0:(n-1)
  
  # generate synthetic barometric pressure from fft coefs
  
  ba <- c(ba_fft_sub, rep(0.0, n - length(ba_fft_sub)))
  ba <- Re(fftw::IFFT(ba, scale=TRUE))
  ba[1:4000] <- NA_real_
  noise_baro <- rnorm(n, sd = 0.00018)
  baro <- ba * 2 - mean(ba, na.rm = TRUE) + noise_baro
  noise_wl <- 0#rnorm(n, sd = sd_noise) + cumsum(rnorm(n, sd = sd_noise_trend))
  
  if (is.null(baro_kernel)) {
    data('baro_kernel', envir=environment())
  }
  
  baro_response <- baro
  baro_response[!is.na(baro_response)] <- 0.5 * convolve_fft(x = na.omit(baro),
                                                             y = baro_kernel)
  
  lin_trend <- linear_trend * as.numeric(datetime)
  intercept <- 0
  wl <- baro_response + lin_trend + noise_wl + intercept
  
  return(na.omit(data.table(datetime, wl, baro,
                            baro_response, noise_baro,
                            noise_wl,
                            lin_trend)))
  
}

