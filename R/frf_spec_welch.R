#' spec_welch
#'
#' calculate the spectra and cross-spectra.  The scaling of individual series 
#' requires a closer look.
#' 
#' 
#' @param x \code{numeric matrix} univariate or multivariate time series
#' @param n_subsets number of subsets (integer)
#' @param overlap amount of overlap 0.0 < overlap < 1.0 (numeric)
#' @param window window function to apply (function). 
#'   The options are window_hann, window_hamming, window_blackman, 
#'   window_nuttall, window_blackman_nuttall, window_blackman_harris
#' @param demean should the mean be removed from x prior to calculation (logical)
#' @param detrend should the x be detrended to calculation (logical)
#' @param ... not used
#'
#'
#' @return array of spectra
#'
#' @export
#'
#'
spec_welch <- function (x,
                        n_subsets = 10,
                        overlap = 0.5,
                        window = window_hann,
                        demean = TRUE,
                        detrend = TRUE,
                        ...)
{
  
  n <- nrow(x)
  nser <- ncol(x)
  
  # determine subsets
  len_sub   <- floor(n / n_subsets)
  n_overlap <- floor(len_sub * overlap)
  sub_s     <- seq(1, n-len_sub, len_sub - n_overlap)
  sub_e     <- sub_s + len_sub-1
  
  # pad and apply window
  padding <- nextn(len_sub) - len_sub
  n       <- len_sub + padding
  win     <- window(len_sub)
  n_ffts  <- length(sub_s)
  nenbw   <- n * sum(win^2) / sum(win)^2
  
  #fs <- 3000
  # eq 21 from https://holometer.fnal.gov/GH_FFT.pdf
  scale   <- mean(win)^2 * nenbw * n_ffts * n 
  #scale   <- sum(win^2) * n_ffts * (n/4)
  #scale   <- sum(win)^2 * n_ffts * sqrt(n) * 2.0

  
  # pre-allocate matrices
  #n_out <- (n/2+1)
  n_out <- n
  pgram <- array(0i, dim = c(n_out, nser, nser))
  xfft  <- matrix(0i, nrow = n_out, ncol = nser)
  
  subs <- matrix(0.0, nrow = n, ncol = nser)
  subs_ind <- 1:len_sub
  
  for (i in 1:n_ffts) {
    
    # detrend and  demean
    subs <- x[sub_s[i]:sub_e[i],, drop = FALSE]
    if (detrend) {
      t <- 1L:len_sub - (len_sub + 1)/2
      sumt2 <- len_sub * (len_sub^2 - 1) / 12
      for (i in 1L:ncol(x))
        subs[, i] <- subs[, i] - mean(subs[, i]) - sum(subs[, i] * t) * t/sumt2
    } else if (demean) {
      x <- sweep(x, 2, colMeans(subs), check.margin=FALSE)
    }
    subs <- subs * win
    
    #subs <- subs
    #subs <- rbind(subs, matrix(0, nrow = padding, ncol = ncol(x)))
    
    
    for (j in 1:nser) {
      # xfft[, j] <- fftwtools::fftw(c(subs[, j], rep_len(0.0, padding)), HermConj = FALSE)
      # xfft[, j] <- fft(c(subs[, j], rep_len(0.0, padding)), HermConj = FALSE)[1:n_out]
      xfft[, j] <- fftw::FFT(c(subs[, j], rep_len(0.0, padding)))[1:n_out]
      
    }
    
    for (j in 1L:ncol(x)) {
      for (k in 1L:ncol(x)) {
        
        if(k >= j) {
          pgram[, j, k] <- pgram[, j, k] + xfft[, j] * Conj(xfft[, k])
        } else {
          pgram[, j, k] <- Conj(pgram[, k, j])
        }
        
        ## value at zero is invalid as mean has been removed, so interpolate:
        pgram[1, j, k] <- 0.5 * (pgram[2, j, k] + pgram[n_out, j, k])
      }
    }
    
    
  }
  
  #pgram <- pgram / n_ffts
  #pgram <- pgram[-1,,] / scale
  pgram <- pgram / scale
  
  
  return(pgram)
}



