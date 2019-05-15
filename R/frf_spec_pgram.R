#' pad_kernel
#'
#' Pad a smoothing kernel
#'  
#' @param k \code{tskernel} smoothing "tskernel" object.
#' @param n \code{integer} the total length
#'
#' @return vector of the padded kernel
#'
#' @importFrom stats is.tskernel
#' @export
#'
#'
pad_kernel <- function(k, n){
  
  if (!is.tskernel(k)) stop ("'k' is not a kernel")
  m <- k$m
  if (n <= 2L * m)
    stop ("x is shorter than kernel k")
  if (m == 0L)
    return (k)
  else {
    return(c(k[0L:m], rep_len(0, n - 2L * m - 1L), k[-m:-1L]))
  }
  
}







#' spec_pgram
#'
#' calculate the spectra and cross-spectra
#' 
#' @inheritParams stats::spec.pgram
#' 
#' @importFrom stats spec.taper
#' 
#' @return array of spectra
#'
#' @export
#'
#' @importFrom stats na.fail
#'
spec_pgram <-
  function (x, spans = NULL, kernel = NULL, taper = 0.1,
            pad = 0, fast = TRUE,
            demean = FALSE, detrend = TRUE,
            plot = TRUE, na.action = na.fail, ...) {
    
    ## Estimate spectral density from (smoothed) periodogram.
    series <- deparse(substitute(x))
    x <- na.action(as.ts(x))
    xfreq <- frequency(x)
    x <- as.matrix(x)
    N <- N0 <- nrow(x)
    nser <- ncol(x)
    
    if(!is.null(spans)) # allow user to mistake order of args
      kernel <- {
        if(is.tskernel(spans)) spans else
          kernel("modified.daniell", spans %/% 2)
      }
    
    if(!is.null(kernel) && !is.tskernel(kernel))
      stop("must specify 'spans' or a valid kernel")
    
    if (detrend) {
      t <- 1L:N - (N + 1)/2
      sumt2 <- N * (N^2 - 1)/12
      
      for (i in 1L:ncol(x))
        x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
      
    } else if (demean) {
      x <- sweep(x, 2, colMeans(x), check.margin=FALSE)
    }
    
    ## apply taper:
    x <- spec.taper(x, taper)
    ## to correct for tapering: Bloomfield (1976, p. 194)
    ## Total taper is taper*2
    # u2 <- (1 - (5/8) * taper * 2)
    # u4 <- (1 - (93/128) * taper * 2)
    
    if (pad > 0) {
      x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
      N <- nrow(x)
    }
    
    NewN <- if(fast) nextn(N) else N
    x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
    N <- nrow(x)
    
    Nspec <- floor(N/2)
    
    #freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)
    
    # use fftw package to calculate ffts for long series
    if (N > 1e5) {
      
      xfft <- matrix(NA_complex_, nrow = N, ncol = nser)
      for (i in 1:nser) {
        xfft[, i] <- fftw::FFT(x[, i])
      }
      
    } else {
      xfft <- mvfft(x)
    }
    
    pgram <- array(NA, dim = c(N, ncol(x), ncol(x)))
    
    # only calculate upper triangle vectors - jrk
    for (i in 1L:ncol(x)) {
      for (j in i:ncol(x)) { # N0 = #{non-0-padded}
        pgram[, i, j] <- xfft[, i] * Conj(xfft[, j]) / (N0 * xfreq)
        ## value at zero is invalid as mean has been removed, so interpolate:
        pgram[1, i, j] <- 0.5*(pgram[2, i, j] + pgram[N, i, j])
      }
    }
    
    if(!is.null(kernel)) {
      # fill in lower triangle with the complex conjugate to save compuation - jrk
      k_fft <- fftw::FFT(pad_kernel(kernel, N))
      
      for (i in 1L:ncol(x)) for (j in 1:ncol(x))
        if(i <= j) {
          if (N > 1e5) {
            pgram[, i, j] <- fftw::IFFT(fftw::FFT(pgram[, i, j]) * k_fft, scale = FALSE) / N
          } else {
            pgram[, i, j] <- fft(fft(pgram[, i, j]) * k_fft, inverse = TRUE) / N
          }
        } else {
          pgram[, i, j] <- Conj(pgram[, j, i])
        }
    }

    pgram <- pgram[2:(Nspec+1),,, drop=FALSE]

    return(pgram)

  }





#' #' kern_apply
#' #'
#' #' apply smoothing kernel using fftw
#' #'
#' #' @param x \code{numeric vector} an input vector to be smoothed.
#' #' @param k \code{tskernel} smoothing "tskernel" object.
#' #' @param circular \code{logical} a logical indicating whether the input sequence to be smoothed is treated as circular, i.e., periodic.
#' #' @param ... currently not used.
#' #'
#' #' @return smoothed periodogram
#' #' 
#' #' @importFrom stats fft
#' #' @importFrom stats mvfft
#' #' 
#' #' @export
#' #'
#' #'
#' kern_apply <- function (x, k, circular = FALSE, ...)
#' {
#'   if (!is.vector(x)) stop ("'x' is not a vector")
#'   if (!is.tskernel(k)) stop ("'k' is not a kernel")
#'   m <- k$m
#'   if (length(x) <= 2L*m)
#'     stop ("'x' is shorter than kernel 'k'")
#'   if (m == 0L)
#'     return (x)
#'   else
#'   {
#'     n <- length(x)
#'     w <- c(k[0L:m], rep_len(0, n-2L*m-1L), k[-m:-1L])
#'     if (n > 1e5) {
#'       y <- fftw::IFFT(fftw::FFT(x) * fftw::FFT(w), scale = FALSE) / n
#'     } else {
#'       y <- fft(fft(x)*fft(w), inverse = TRUE) / n
#'     }
#'     if (is.numeric(x)) y <- Re(y)
#'     if (circular)
#'       return (y)
#'     else
#'       return (y[(1L+m):(n-m)])
#'   }
#' }
#' 
#' #' xfft
#' #'
#' #' calculate the spectra and cross-spectra
#' #' 
#' #' @inheritParams stats::spec.pgram
#' #'
#' #' @return matrix of spectra
#' #' 
#' #' @importFrom stats as.ts
#' #' @importFrom stats frequency
#' #' @importFrom stats nextn
#' #' 
#' #' @export
#' #'
#' #'
#' xfft <-
#'   function (x, spans = NULL, kernel = NULL, taper = 0.1,
#'             pad = 0, fast = TRUE,
#'             demean = FALSE, detrend = TRUE,
#'             plot = TRUE, na.action = na.fail, ...) {
#'     
#'     ## Estimate spectral density from (smoothed) periodogram.
#'     series <- deparse(substitute(x))
#'     x <- na.action(as.ts(x))
#'     xfreq <- frequency(x)
#'     x <- as.matrix(x)
#'     N <- N0 <- nrow(x)
#'     nser <- ncol(x)
#'     
#'     # if(!is.null(spans)) # allow user to mistake order of args
#'     #   kernel <- {
#'     #     if(is.tskernel(spans)) spans else
#'     #       kernel("modified.daniell", spans %/% 2)
#'     #   }
#'     # 
#'     # if(!is.null(kernel) && !is.tskernel(kernel))
#'     #   stop("must specify 'spans' or a valid kernel")
#'     
#'     if (detrend) {
#'       t <- 1L:N - (N + 1)/2
#'       sumt2 <- N * (N^2 - 1)/12
#'       
#'       for (i in 1L:ncol(x))
#'         x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
#'       
#'     } else if (demean) {
#'       x <- sweep(x, 2, colMeans(x), check.margin=FALSE)
#'     }
#'     
#'     ## apply taper:
#'     #x <- spec.taper(x, taper)
#'     
#'     ## to correct for tapering: Bloomfield (1976, p. 194)
#'     ## Total taper is taper*2
#'     # u2 <- (1 - (5/8) * taper * 2)
#'     # u4 <- (1 - (93/128) * taper * 2)
#'     
#'     if (pad > 0) {
#'       x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
#'       N <- nrow(x)
#'     }
#'     
#'     NewN <- if(fast) nextn(N) else N
#'     x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
#'     N <- nrow(x)
#'     
#'     Nspec <- floor(N/2)
#'     
#'     #freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)
#'     
#'     # use fftw package to calculate ffts for long series
#'     if (N > 1e5) {
#'       
#'       xfft <- matrix(NA_complex_, nrow = N, ncol = nser)
#'       for (i in 1:nser) {
#'         xfft[, i] <- fftw::FFT(x[, i])
#'       }
#'       
#'     } else {
#'       xfft <- mvfft(x)
#'     }
#'     return(xfft[2:(Nspec+1),, drop=FALSE])
#'   }
