#' frf_predict_terms
#'
#' @param frf frequency response function
#' @param x input variable dataset
#' @param time_col datetime
#' @param vars variables used in the transfer function call
#' @param coherency_cut values must have at least this coherency
#' @param limit gain must be lower than this value.  Often useful for Earth tides. 
#'
#' @return data.table of results
#' @export
#'
frf_predict_terms <- function(frf,
                              x, 
                              time_col = 'datetime',
                              vars = c('wl', 'baro'),
                              coherency_cut = 0.3,
                              limit = 2.0
                              ) {
  
  
  n <- nrow(x) %/% 2
  n_days   <- diff(range(as.numeric(x[[time_col]]))) / 86400
  n_cycles <- diff(as.numeric(x[[time_col]][1:2])) / 86400
  
  
  
  # get variable names
  gain_name <- c()
  phase_name <- c()
  coh_name <- c()
  for (i in 1:(length(vars) - 1)) {
    gain_name[i]  <- paste0('gain_',  vars[1], '_', vars[i+1])  
    phase_name[i] <- paste0('phase_', vars[1], '_', vars[i+1])  
    coh_name[i] <- paste0('coherency_', vars[1], '_', vars[i+1])  
  }
  
  
  
  freqs  <- seq(1 / n_days, 1 / (2 * n_cycles), length.out = n)
  
  res <- list()
  
  for (i in seq_along(gain_name)) {
    
    phase <- approx(x = (frf[["frequency"]]), 
                    y = frf[[phase_name[i]]], 
                    xout = (freqs))[['y']]
    gain  <- approx(x = (frf[["frequency"]]),
                    y = frf[[gain_name[i]]], 
                    xout = (freqs))[['y']]
    coh  <- approx(x = (frf[["frequency"]]),
                   y = frf[[coh_name[i]]], 
                   xout = (freqs))[['y']]
    
    gain_cut <- gain[which.min(abs(1.9322736-freqs))] * 2
    
    # replace NA values
    phase[is.na(phase)] <- 0
    gain[is.na(gain)] <- 0
    coh[is.na(coh)] <- 0
    
    # zero out low coherency values
    low_coh <- coh < coherency_cut
    # replace low coherency values with zero
    phase[low_coh] <- 0
    gain[low_coh] <- 0
    
    
    # zero out high gain values
    low_gain <- gain > limit
    # replace low coherency values
    phase[low_gain] <- 0
    gain[low_gain] <- 0
    
    
    tf_interp <- complex(modulus = gain, argument = -phase)
    
    if (length(tf_interp) %% 2 == 1) {
      tf_interp <- c(tf_interp[1],
                     tf_interp, Conj(rev(tf_interp)))
    } else {
      tf_interp <- c(tf_interp[1],
                     tf_interp, Conj(rev(tf_interp)[-1]))
    }
    
    res[[i]] <- Re(IFFT(tf_interp * FFT(x[[vars[[i+1]]]])))
  }
  
  names(res) <- vars[-1]
  
  res
}