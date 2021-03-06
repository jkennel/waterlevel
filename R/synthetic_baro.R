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
#' @param scale multiplier for barometric pressure (numeric)
#' @param baro_kernel vector values to convolve with barometric pressure
#'
#' @return data.table of synthetic water levels and barometric pressure
#' 
#' @export
#' 
#' @importFrom stats rnorm
#'
synthetic <- function(sd_noise = 0.0002,
                      sd_noise_trend = 0.00003,
                      linear_trend = 1e-7,
                      n = 14 * 86400,
                      seed = NULL,
                      scale = 0.5,
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
  
  
  baro_response <- baro
  baro_response[!is.na(baro_response)] <- scale * convolve_fft(x = na.omit(baro),
                                                             y = baro_kernel)
  
  lin_trend <- linear_trend * as.numeric(datetime)
  intercept <- 0
  wl <- baro_response + lin_trend + noise_wl + intercept
  
  return(na.omit(data.table(datetime, wl, baro,
                            baro_response, noise_baro,
                            noise_wl,
                            lin_trend)))
  
}



#' synthetic_wl
#'
#' This function is used for testing purposes
#'
#' @param baro barometric pressure (numeric vector)
#' @param datetime POSIXct dates
#' @param sd_noise standard deviation of random noise to add (numeric)
#' @param linear_trend magnitude of linear trend in time (numeric)
#' @param intercept 
#' @param kernel 
#' @param seed random number seed for reproducibility (numeric)
#'
#' @return data.table of synthetic water levels and barometric pressure
#' 
#' @export
#' 
#' @importFrom stats rnorm
#'
synthetic_wl <- function(baro,
                         datetime,
                         sd_noise = 0,
                         linear_trend = 0,
                         intercept = 0,
                         seed = NULL,
                         kernel = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # generate synthetic barometric pressure from fft coefs
  noise_wl <- rnorm(length(baro), sd = sd_noise)
  
  wl <- convolve_fft(x = baro,
                     y = kernel)
  
  lin_trend <- linear_trend * as.numeric(datetime)
  
  wl <- wl + lin_trend + noise_wl + intercept
  
  return(wl)
  
}


