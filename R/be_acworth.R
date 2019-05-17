#' be_acworth_eq_4
#'
#' @description Implementation of Acworth et. al. 2016 equation 4.
#' 
#' @param s2_gw \code{numeric} s2 component in the groundwater levels
#' @param s2_et \code{numeric} s2 component in the earth tides
#' @param s2_at \code{numeric} s2 component for atmospheric pressure
#' @param m2_gw \code{numeric} m2 component in the groundwater levels
#' @param m2_et \code{numeric} m2 component in the earth tides
#' @param d_phase \code{numeric} phase difference between Earth tide and atmospheric drivers s2_et and s2_at
#' @param inverse \code{logical} whether the barometric relationship is inverse (TRUE means that when the barometric pressure goes up the measured water level goes down (vented transducer, depth to water), FALSE means that when the barometric pressure goes up so does the measured pressure (non-vented transducer))
#'
#' @references Acworth, R. I., Halloran, L. J. S., Rau, G. C., Cuthbert, M. O., & Bernardi, T. L. (2016). An objective frequency-domain method for quantifying confined aquifer compressible storage using Earth and atmospheric tides. Geophysical Research Letters, 43 (November). \doi{https://doi.org/10.1002/2016GL071328}
#' 
#' @return barometric efficiency
#' @export
#'
#' @examples
#'
#' be_acworth_eq_4(s2_at = 7.461,
#'            s2_et=224.640,
#'            s2_gw=4.086,
#'            m2_gw = 0.471,
#'            m2_et = 492.526,
#'            d_phase=-56.709)
#' be_acworth_eq_4(s2_at = 6.164,
#'            s2_et=270.463,
#'            s2_gw=0.329,
#'            m2_gw = 0.225,
#'            m2_et = 551.572,
#'            d_phase=-71.726)
#' be_acworth_eq_4(s2_at = 5.897,
#'            s2_et=234.478,
#'            s2_gw=5.536,
#'            m2_gw = 0.773,
#'            m2_et = 558.075,
#'            d_phase=-70.393)
#'
be_acworth_eq_4 <- function(s2_gw,
                            s2_et,
                            s2_at,
                            m2_gw,
                            m2_et,
                            d_phase,
                            inverse = TRUE) {
  
  numer_term_2 <- s2_et * (cos(d_phase) * m2_gw / m2_et)
  
  if (inverse) {
    (s2_gw + numer_term_2) / s2_at
  } else {
    (s2_gw - numer_term_2) / s2_at
  }
  
}



#' be_acworth_eq_4
#'
#' @description Implementation of Acworth et. al. 2016.  This is basically the 
#' transfer function value for the semi-diurnal band.
#'
#' @param dat data that has the independent and dependent variables (data.table)
#' @param wl name of the water level column (character).
#' @param ba name of the barometric pressure column (character). 
#' @param et name of the Earth tide column (character). 
#' @param inverse \code{logical} whether the barometric relationship is inverse
#' (TRUE means that when the barometric pressure goes up the measured water
#' level goes down (vented transducer, depth to water), FALSE means that when 
#' the barometric pressure goes up so does the measured pressure 
#' (non-vented transducer)).
#' @param ... arguments to pass to transfer_acworth, spec_welch, spec_pgram
#'
#' @references Acworth, R. I., Halloran, L. J. S., Rau, G. C., Cuthbert, M. O., & Bernardi, T. L. (2016). An objective frequency-domain method for quantifying confined aquifer compressible storage using Earth and atmospheric tides. Geophysical Research Letters, 43 (November). \doi{https://doi.org/10.1002/2016GL071328}
#' @references Acworth, R. I., G. C. Rau, L. J. S. Halloran, and W. A. Timms (2017), Vertical groundwater stor1077 age properties and changes in confinement determined using hydraulic head response to atmospheric tides, Water Resources Research, 53(4), 2983–2997, \doi{https://doi.org/10.1002/2016WR020311}.
#' 
#' @return barometric efficiency
#' @export
#'
be_acworth <- function(dat, wl = 'wl', ba = 'ba', et = 'et', inverse = TRUE, ...) {
  
  spec <- transfer_acworth(dat, vars = c(wl, ba, et), 
                           time = 'datetime', ...)
  
  
  m2_freq <- 1.9322736
  s2_freq <- 2.0000000
  
  m2 <- spec[which.min(abs(spec$frequency-m2_freq)),]
  s2 <- spec[which.min(abs(spec$frequency-s2_freq)),]
  
  # print(m2)
  # print(s2)
  be_acworth_eq_4(s2_gw = s2$wl_amp, s2$et_amp, s2$ba_amp,
                  m2_et = m2$et_amp, m2_gw = m2$wl_amp, 
                  d_phase = s2$d_phase, inverse = inverse)
  
}

transfer_acworth <- function(dat, vars, time = 'datetime', 
                             method = 'spec_pgram', ...) {
  
  t_interval <- diff(as.numeric(dat[[time]][1:2]))
  
  if (method == 'spec_pgram') {
    
    
    spec_arr <- spec_pgram(as.matrix(dat[, vars, with = FALSE]), ...)
    spec <- matrix(NA_complex_, ncol = 3, nrow = nrow(spec_arr))
    n_padded <- nextn(nrow(dat))
    df <- 1 / n_padded
    frequency <- seq.int(from = df, by = df, 
                         length.out = floor(n_padded/2)) * 86400/t_interval
    
    
  } else if (method == 'spec_welch') {
    
    spec_arr <- spec_welch(as.matrix(dat[, vars, with = FALSE]), ...)
    spec <- matrix(NA_complex_, ncol = 3, nrow = nrow(spec_arr))
    n <- 2 * (nrow(spec) - 1)
    dt <- (n / 86400)
    frequency <- seq.int(from = 1/n, by = 1/n,
                         length.out = nrow(spec)) * 86400/t_interval
    
  } else {
    stop(paste(method, 'method not yet implemented'))
  }
  
  
  for(i in 1:length(vars)) { spec[,i] <- spec_arr[, i, i] }
  colnames(spec) <- vars
  
  amp   <-  sqrt(abs(Re(spec)))
  phase <-  Arg(spec_arr[, 2, 3])
  
  df <- data.table(frequency, amp, phase)
  setnames(df, c('frequency', 'wl_amp', 'ba_amp', 'et_amp', 'd_phase'))
  
  return(df)
  
}






