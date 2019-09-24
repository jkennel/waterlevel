#' spec_from_pgram
#'
#' @param pgram the periodogram
#'
#' @return spectrum from the pgram
#' @export
#'
spec_from_pgram <- function(pgram) {
  
  m <- matrix(NA_real_, nrow = nrow(pgram), ncol=ncol(pgram))
  
  for(i in 1:ncol(pgram)) {
    m[,i] <- Re(pgram[, i, i])
  }
  
  return(m)
  
}