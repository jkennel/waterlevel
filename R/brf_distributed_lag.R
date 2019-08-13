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
#' @param n_subset for distributed_lag parallel
#' @param n_shift for distributed_lag parallel
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
                            lag_name = '', 
                            n_subset = 1,
                            n_shift = 0) {

  # generate basis functions
  max_knot <- max(knots)
  n_knots  <- length(knots)
  x_len <- length(as.numeric(x))
  
  if(max_knot > x_len) {
    stop('The maximum knot cannot be larger than the number of elements in x')
  }
  
  # generate basis lag
  basis_lag <- spline_fun(min(knots):max_knot,
                          knots = knots[-c(1, n_knots)],
                          Boundary.knots = range(knots),
                          intercept = TRUE)
  
  # convolution - fft for large number of lags, otherwise use parallel version
  if(any(is.na(x)) | max_knot < 5000 | n_subset != 1 | n_shift != 0) {
    dist_lag_mat <- distributed_lag_parallel(rev(x), 
                                             t(as.matrix(basis_lag)), 
                                             max_knot-min(knots),
                                             n_subset = n_subset,
                                             n_shift = n_shift
                                             )
  } else {
    dist_lag_mat <- cross_basis_fft(as.matrix(x), basis_lag)
  }
  
  colnames(dist_lag_mat) <- paste0("distributed_lag_", lag_name, '_', knots)
    
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
  ny   <- length(y)
  n1   <- ny - 1
  x    <- c(rep.int(0, n1), x)
  n    <- length(y <- c(y, rep.int(0, n_in - 1)))
  x    <- fftw::IFFT(fftw::FFT(x) * Conj(fftw::FFT(y)), scale = FALSE) / n
  #x <- fftw::IFFT(fftw::FFT(x) * Conj(fftw::FFT(y)))
  x[1:n1] <- NA_real_
  return(Re(x)[1:n_in])
}



#' cross_basis_fft
#'
#' @param basisvar numeric vector
#' @param basislag lagging matrix
#' 
#' @return numeric vector result of convolution
#' @export
#'
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


# ==============================================================================
#' @title
#' fftw_convolve
#'
#' @description
#' Do convolution for a data series and matrix of equal length filters.
#' For long data series this function will be much faster than non-fft methods.
#'
#' @param x \code{numeric} vector or single column matrix the data series
#' @param y \code{numeric} matrix of filter(s)
#' @param normalize \code{logical} do you want the values normalized
#' @param align \code{character} what alignment for the convolution
#' @param dist_lag \code{logical} is this for a distributed lag?
#'
#' @return convolution of data series and filter
#'
#' @export
#'
#' @examples
#'
#' library(splines)
#'
#' knots <- round(10^seq(0.5, 4.635484, length.out = 9))
#' n_lags <- 86400*2
#'
#' bla <- ns(0:(n_lags), knots = knots,
#'           Boundary.knots = c(0, n_lags))
#' bv <- rnorm(86400*2 + 3600*24*4)
#' system.time({
#' dl <- fftw_convolve(bv, bla, normalize = FALSE)
#' })
#'
fftw_convolve <- function(x,
                          y,
                          normalize = TRUE,
                          align = 'center') {
  
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  n_x_in <- nrow(x)
  n_y    <- nrow(y)
  
  # reverse order
  y <- y[n_y:1, , drop = FALSE]
  
  if (n_y > n_x_in) {
    stop('The length of y should be less than or equal to x')
  }
  
  
  # set NA indices
  start <- 1:(ceiling(n_y / 2) - 1)
  end   <- (n_x_in - floor(n_y / 2) + 1):(n_x_in)
  
  
  if(n_x_in %% 2 == 1) {
    x     <- c(x, 0.0)
    n_x   <- length(x)
    x_pad <- n_x / 2
    sub   <- c((n_x + x_pad + 1):(n_x + 2 * x_pad), 1:(x_pad-1))
  } else {
    n_x   <- n_x_in
    x_pad <- n_x / 2
    sub   <- c((n_x + x_pad + 1):(n_x + 2 * x_pad), 1:(x_pad))
  }
  
  n_sub <- length(sub)
  
  # normalize filter to sum to 1
  if (normalize) {
    y <- t(t(y)  / colSums(y))
  }
  
  y_pad <- x_pad + ceiling((n_x - n_y) / 2)
  
  x_pad <- rep(0.0, x_pad)
  y_pad <- rep(0.0, y_pad)
  
  # do the FFT on the x vector
  f <- fftw::FFT(c(x_pad, x, x_pad))
  u <- matrix(NA_real_, nrow = n_sub, ncol = ncol(y))
  
  if (n_y %% 2 == 0) {
    
    warning('Values are shifted 0.5 units forward. Use odd number filter for better centering')
    for (i in 1:ncol(y)){
      u[, i] <- Re(fftw::IFFT(fftw::FFT(c(y_pad, y[,i], y_pad)) * f))[sub]
    }
    
  } else {
    for (i in 1:ncol(y)){
      u[,i] <- Re(fftw::IFFT(fftw::FFT(c(y_pad, y[,i], y_pad[-1])) * f))[sub]
    }
    
  }
  
  u[start,] <- NA_real_
  u[end,]   <- NA_real_
  
  # adjust output for alignment
  if (align == 'right') {
    for(i in 1:ncol(y)){
      u[,i] <- data.table::shift(u[,i], n = length(end), type = 'lag')
    }
  } else if (align == 'left') {
    for(i in 1:ncol(y)) {
      u[,i] <- data.table::shift(u[,i], n = length(start), type = 'lead')
    }
  }
  
  return(u)
}

