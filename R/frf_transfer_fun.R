
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
    
    pgram     <- spec_pgram(as.matrix(dat[, vars, with = FALSE]), spans = 3)
    n_padded  <- nextn(nrow(dat))
    df        <- 1 / n_padded
    frequency <- seq.int(from = df, by = df, 
                         length.out = floor(n_padded/2)) * 86400/t_interval
    
  } else if (method == 'spec_welch') {
    
    pgram <- spec_welch(as.matrix(dat[, vars, with = FALSE]), ...)
    n <- 2 * (nrow(pgram) - 1)
    dt <- (n / 86400)
    frequency <- seq.int(from = 1/n, by = 1/n, 
                         length.out = nrow(pgram)) * 86400/t_interval
    
  } else {
    stop('method must be either spec_welch or spec_pgram')
  }
  
  # solve for transfer function
  tf    <- solve_tf_parallel(pgram)
  
  # solve for gain
  gain  <- Mod(tf)
  
  # solve for phase
  phase <- atan2(Im(tf), Re(tf))
  
  # coherency
  coherency <- coh_phase(pgram)
  
  
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
      coh_name[k+length(vars)]  <- paste0('coherency_phase_', vars[i], '_', vars[j]) 
      k <- k + 1
    }
  }
  
  colnames(tf) <- tf_name
  colnames(gain) <- gain_name
  colnames(phase) <- phase_name
  colnames(coherency) <- coh_name
  
  return(as_tibble(cbind(frequency, tf, gain, phase, coherency)))
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
