#' spec_from_pgram
#'
#' @param pgram the periodogram
#' @param method either spec_pgram or spec_welch
#' @param ... additional arguments
#'
#' @return spectrum from the pgram
#' @export
#'
spec_from_pgram <- function(pgram, method = 'spec_pgram',  ...) {
  
  dots <- list(...)
  
  if ('taper' %in% names(dots)) {
    u2 <- (1 - (5/8) * dots[['taper']] * 2)
  } else if (method == 'spec_pgram') {
    u2 <- (1 - (5/8) * 0.1 * 2)
  } else {
    u2 <- 1.0
  }
  
  m <- matrix(NA_real_, nrow = nrow(pgram), ncol=ncol(pgram))
  
  for(i in 1:ncol(pgram)) {
    m[,i] <- Re(pgram[, i, i]) / u2
  }
  
  return(m)
  
}
