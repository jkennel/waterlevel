#https://en.wikipedia.org/wiki/Window_function


#' window_rectangular
#'
#' Rectangular window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_rectangular(100)
#'
window_rectangular <- function(n){
  
  rep(1.0, n)
  
}


#' window_hann
#'
#' Hann window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_hann(100)
#'
window_hann <- function(n){
  
  0.5 * (1 - cos((2 * pi * 0:(n - 1)) / (n - 1)))
  
}



#' window_hamming
#'
#' Hamming window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_hamming(100)
#'
window_hamming <- function(n){
  
  0.54 - (0.46 * cos((2 * pi * 0:(n - 1)) / (n - 1)))
  
}

#' window_blackman
#'
#' Blackman window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_blackman(100)
#'
window_blackman <- function(n){
  
  a0 <- 0.42
  a1 <- 0.5
  a2 <- 0.08
  k <- 0:(n - 1) / (n - 1)
  
  a0 - a1 * cos(2 * pi * k) +
    a2 * cos(4 * pi * k)
  
}

#' window_first_deriv
#'
#' First derivative window for FFT
#'
#' @param n length of signal (integer)
#' @param a0 \code{numeric} coefficient
#' @param a1 \code{numeric} coefficient
#' @param a2 \code{numeric} coefficient
#' @param a3 \code{numeric} coefficient
#'
#' @return window
#'
#' @export
#'
#' @examples
#' # nuttall window
#' window_first_deriv(100, 0.355768, 0.487396, 0.144232, 0.012604)
#'
window_first_deriv <- function(n,
                               a0,
                               a1,
                               a2,
                               a3) {
  
  k <- 0:(n - 1) / (n - 1)
  a0 - a1 * cos(2 * pi * k) +
    a2 * cos(4 * pi * k) -
    a3 * cos(6 * pi * k)
  
}

#' window_nuttall
#'
#' Nuttall window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_nuttall(100)
#'
window_nuttall <- function(n){
  
  a0 <- 0.355768
  a1 <- 0.487396
  a2 <- 0.144232
  a3 <- 0.012604
  
  window_first_deriv(n, a0, a1, a2, a3)
  
}

#' window_blackman_nuttall
#'
#' Blackman-Nuttall window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_blackman_nuttall(100)
#'
window_blackman_nuttall <- function(n){
  
  a0 <- 0.3635819
  a1 <- 0.4891775
  a2 <- 0.1365995
  a3 <- 0.0106411
  
  window_first_deriv(n, a0, a1, a2, a3)
  
}

#' window_blackman_harris
#'
#' Blackman-Harris window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_blackman_harris(100)
#'
window_blackman_harris <- function(n){
  
  a0 <- 0.35875
  a1 <- 0.48829
  a2 <- 0.14128
  a3 <- 0.01168
  
  window_first_deriv(n, a0, a1, a2, a3)
  
}



#' window_gaussian
#'
#' Gaussian window for FFT
#' https://en.wikipedia.org/wiki/Window_function
#'
#' @param n length of signal (integer)
#' @param sigma the standard deviation in periods (numeric)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_gaussian(100, 0.3)
#
window_gaussian <- function(n,
                            sigma) {
  
  exp(-0.5 * ((0:(n-1) - (n-1) / 2) / (sigma * (n-1) / 2))^2)
  
}

