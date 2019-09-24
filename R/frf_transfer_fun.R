#' transfer_fun
#'
#' @param dat data that has the independent and dependent variables (data.table)
#' @param vars variables to include for transfer function calculation
#' @param time the name of the column that contains the POSIXct date and time
#' @param method either spec_pgram or spec_welch
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
    frequency <- seq.int(from = df, by = df, 
                         length.out = n_padded) * 86400/t_interval
    
  } else if (method == 'spec_welch') {
    
    pgram <- spec_welch(as.matrix(dat[, vars, with = FALSE]), ...)
    n <- (nrow(pgram) - 1)
    dt <- (n / 86400)
    frequency <- seq.int(from = 1/n, by = 1/n, 
                         length.out = nrow(pgram)) * 86400/t_interval
    
  } else {
    stop(paste(method, 'method not yet implemented'))
  }
  
  if (length(vars) == 1) {
    warning('only one variable provided to transfer_fun, returning the spectral density')
    
    if (method == 'spec_pgram') {
      
      dots <- list(...)
      if ('taper' %in% names(dots)) {
        u2 <- (1 - (5/8) * dots[['taper']] * 2)
      } else {
        u2 <- (1 - (5/8) * 0.1 * 2)
      }
      
      return(data.table(frequency = Re(frequency), 
                        phase     = Arg(pgram),
                        amplitude = Re(pgram) / u2))
    } else {
      return(data.table(frequency = Re(frequency), 
                        phase     = Arg(pgram),
                        amplitude = Re(pgram)))
    }
    
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
  spectrum <- spec_from_pgram(pgram)
  
  tf_name <- c()
  gain_name <- c()
  phase_name <- c()
  for (i in 1:(length(vars) - 1)) {
    tf_name[i]    <- paste0('transfer_', vars[1], '_', vars[i+1])  
    gain_name[i]  <- paste0('gain_',     vars[1], '_', vars[i+1])  
    phase_name[i] <- paste0('phase_',    vars[1], '_', vars[i+1])  
  }
  
  
  coh_name <- c()
  k <- 1
  for (i in 1:(length(vars)-1)) {
    for (j in (1+i):(length(vars))) {
      coh_name[k]  <- paste0('coherency_', vars[i], '_', vars[j]) 
      k <- k + 1
      coh_name[k]  <- paste0('coherency_phase_', vars[i], '_', vars[j]) 
      k <- k + 1
    }
  }
  
  
  colnames(tf) <- tf_name
  colnames(gain) <- gain_name
  colnames(phase) <- phase_name
  colnames(coherency) <- coh_name
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
  
  coh <- phase <- matrix(NA_real_, nrow = n_r, ncol = n_ser * (n_ser - 1) / 2)
  
  for (i in 1L:(n_ser - 1)) {
    
    for (j in (i + 1):n_ser) {
      
      ind <- i + (j - 1) * (j - 2)/2
      
      coh[, ind]   <- Re(Mod(pgram[, i, j])^2 / (pgram[, i, i] * pgram[, j, j]))
      phase[, ind] <- Arg(pgram[, i, j])
      
    }
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
#' @return
#' @export
#'
#' @examples
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