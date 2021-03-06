#' transfer_fun
#'
#' @param dat data that has the independent and dependent variables (data.table)
#' @param vars variables to include for transfer function calculation
#' @param time the name of the column that contains the POSIXct date and time
#' @param method either spec_pgram, spec_welch, or spec_multitaper
#' @param ... arguments to pass to method
#'
#' @return transfer function, gain and phase, coherency
#' 
#' @export
#'
transfer_fun <- function(dat, vars, time = 'datetime', method = 'spec_pgram', ...) {
  
  t_interval <- diff(as.numeric(dat[[time]][1:2]))
  
  if (method == 'spec_pgram') {
    
    pgram    <- spec_pgram(as.matrix(dat[, vars, with = FALSE]), ...)
    n_padded <- nrow(pgram)
    df       <- 1 / n_padded
    # frequency <- seq.int(from = df, by = df, 
    #                      length.out = floor(n_padded/2)) * 86400/t_interval
    frequency <- seq.int(from = 0, by = df, 
                         length.out = n_padded) * 86400/t_interval
    
  } else if (method == 'spec_welch') {
    
    pgram <- spec_welch(as.matrix(dat[, vars, with = FALSE]), ...)
    n_padded <- nrow(pgram)
    df       <- 1 / n_padded
    # frequency <- seq.int(from = df, by = df, 
    #                      length.out = floor(n_padded/2)) * 86400/t_interval
    frequency <- seq.int(from = 0, by = df, 
                         length.out = n_padded) * 86400/t_interval
  } else if (method == 'spec_multitaper') {
    
    x <- list(...)
    if(!'tapers' %in% names(x)) {
      tapers <- 7L
    } else {
      tapers <- x[['tapers']]
    }
    
    fftz <- apply(as.matrix(dat[, vars, with = FALSE]), 2, FFT)
    pgram <- resample_fft_parallel(fftz, tapers = tapers)$psd
    n_padded <- nrow(fftz)
    df       <- 1 / n_padded
    frequency <- seq.int(from = df, by = df, 
                         length.out = n_padded) * 2*86400/t_interval
    frequency <- frequency[1:nrow(pgram)]
  } else {
    stop(paste(method, 'method not yet implemented'))
  }
  
  if (length(vars) == 1) {
    warning('only one variable provided to transfer_fun, returning the spectral density')
    
    return(cbind(data.table(frequency = Re(frequency)),
                 spec = as.numeric(spec_from_pgram(pgram, method = method, ...))))
    
  }
  
  # solve for transfer function
  tf    <- solve_tf_parallel(pgram)
  
  # solve for gain
  gain  <- Mod(tf)
  
  # solve for phase
  phase <- atan2(Im(tf), Re(tf))
  
  # coherency
  coherency <- coh_phase(pgram)
  
  # spectrum
  spectrum <- spec_from_pgram(pgram, method = method, ...)
  
  tf_name <- c()
  gain_name <- c()
  phase_name <- c()
  for (i in 1:(length(vars) - 1)) {
    tf_name[i]    <- paste0('transfer_', vars[1], '_', vars[i+1])  
    gain_name[i]  <- paste0('gain_',     vars[1], '_', vars[i+1])  
    phase_name[i] <- paste0('phase_',    vars[1], '_', vars[i+1])  
  }
  
  
  coh_name <- c()
  coh_phase_name <- c()
  k <- 1
  # for (i in 1:(length(vars)-1)) {
    for (j in (2):(length(vars))) {
      coh_name[j-1]  <- paste0('coherency_', vars[1], '_', vars[j]) 
      coh_phase_name[j-1]  <- paste0('coherency_phase_', vars[1], '_', vars[j]) 
    }
  # }
  
  
  colnames(tf) <- tf_name
  colnames(gain) <- gain_name
  colnames(phase) <- phase_name
  colnames(coherency) <-c(coh_name, coh_phase_name)
  colnames(spectrum) <- vars
  
  # print(str(frequency))
  # print(str(tf))
  # print(str(gain))
  # print(str(phase))
  # print(str(coherency))
  tf <- as_tibble(tf)
  
  return(cbind(frequency = Re(frequency), tf, Re(gain), Re(phase), Re(coherency), Re(spectrum)))
}










#' coh_phase
#'
#' calculate phase and coherency from the spectra and cross-spectra
#'
#' @param pgram \code{numeric array} must be multivariate
#'
#'
#' @return list of coherency and phase
#'
#' @export
#'
#'
coh_phase <- function(pgram) {
  
  dims = dim(pgram)
  n_ser <- dims[2]
  n_r <- dims[1]
  
  if (n_ser == 1) {
    stop('Cannot calculate coherency for a single periodogram')
  }
  
  coh <- phase <- matrix(NA_real_, nrow = n_r, ncol = (n_ser - 1))
  
  for (i in 2L:(n_ser)) {
    
    #for (j in (i + 1):n_ser) {
      
      #ind <- i + (j - 1) * (j - 2)/2
      
      coh[, i-1]   <- Re(Mod(pgram[, 1, i])^2 / (pgram[, 1, 1] * pgram[, i, i]))
      phase[, i-1] <- Arg(pgram[, 1, i])
      
    #}
  }
  
  return(cbind(coh, phase))
}


#' spec_welch_tf_error
#'
#' @param coh coherency
#' @param amp amplitude
#' @param overlap overlap fraction from spec_welch
#' @param n_subsets number of subsets in spec welch
#' @param return either 'amp', or 'phase'
#'
#' @return error
#' 
#' @export
#'
spec_welch_tf_error <- function(coh, amp, overlap = 0.5, n_subsets = 10, return = 'amp') {
  
  # Hussein B4, B5, B6
  df <- n_subsets - ((n_subsets-1) * overlap)               # B6 text
  
  coh_error   <- sqrt(0.5 * (1 / df) * (( 1 / (coh)^2)-1))  # B6
  amp_error   <- coh_error * amp                            # B4
  phase_error <- coh_error                                  # B5
  
  if(return == 'phase') {
    return(phase_error)
  } else {
    return(amp_error)
  }

}