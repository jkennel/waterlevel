#' distributed_lag
#'
#' Depending on the vector to lag and the maximum knot, distributed_lag will 
#' either use an FFT (no NA and large maximum knot), or a parallel method 
#' (NA, or small maximum knot).
#' 
#' @param x numeric vector to lag
#' @param knots specific knots for the lagging process
#' @param spline_fun spline function to use i.e. splines::ns, splines::bs
#' @param lag_name name of variable to lag
#'
#' @return matrix with distributed lag terms
#' 
#' @importFrom splines ns
#' 
#' @export
#' 
distributed_lag <- function(x, 
                            knots, 
                            spline_fun = ns,
                            lag_name = '') {

  # generate basis functions
  max_knot <- max(knots)
  n_knots  <- length(knots)
  
  
  # generate basis lag
  basis_lag <- spline_fun(min(knots):max_knot,
                          knots = knots[-c(1, n_knots)],
                          Boundary.knots = range(knots),
                          intercept = TRUE)
  
  # convolution - fft for large number of lags, otherwise use parallel version
  if(any(is.na(x)) | max_knot < 5000) {
    dist_lag_mat <- distributed_lag_parallel(rev(x), 
                                             t(as.matrix(basis_lag)), 
                                             max_knot - 1L)
  } else {
    dist_lag_mat <- cross_basis_fft(as.matrix(x), basis_lag)
  }
  
  colnames(dist_lag_mat) <- paste0("distributed_lag_", knots, '_', lag_name)
    
  return(dist_lag_mat)
}



#' convolve_fft
#'
#' @param x numeric vector
#' @param y convolution kernel
#' 
#' @return numeric vector result of convolution
#' @export
#'
convolve_fft <- function(x, y)
{
  n_in <- length(x)
  ny <- length(y)
  n1 <- ny - 1
  x <- c(rep.int(0, n1), x)
  n <- length(y <- c(y, rep.int(0, n_in - 1)))
  x <- fftw::IFFT(fftw::FFT(x) * Conj(fftw::FFT(y)), scale = FALSE) / n
  #x <- fftw::IFFT(fftw::FFT(x) * Conj(fftw::FFT(y)))
  x[1:n1] <- NA_real_
  return(Re(x)[1:n_in])
}


cross_basis_fft <- function(basisvar, basislag) 
{
  
  mat <- matrix(NA_real_, 
                nrow = nrow(basisvar), 
                ncol = ncol(basisvar) * ncol(basislag))
  
  s <- nrow(basisvar) + 1
  e <- nrow(basisvar) + nrow(basislag) -1
  
  for(v in seq(length=ncol(basisvar))) {
    for(l in seq(length=ncol(basislag))) {
      
      mat[, ncol(basislag)*(v-1)+l] <- 
        convolve_fft(basisvar[,v], rev(basislag[,l]))
      
    }
  }
  
  return(mat)
  
}
