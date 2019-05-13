#' get_main_frequency
#'
#' Get the frequency of the wave with the maximum amplitude in a range.  This is
#' a placeholder until the earthtide package is updated
#'
#' @param start the starting frequency in cycles per day (numeric)
#' @param end the ending frequency in cycles per day (numeric)
#'
#' @return the main frequency between start and end
#' @export
#'
get_main_frequency <- function(start, end) {
  
  if(length(start) != length(end)) {
    stop('start and end must be the same length')  
  }
  
  if(any(start > end)) {
    stop('all values of start must be greater than the corresponding values of
         end')
  }
  
  mn <- c()
  for (i in seq_along(start)) {
    et     <- earthtide:::ksm04
    inds   <- et$frequency_cpd >= start[i] & et$frequency_cpd < end[i]
    et_sub <- et[inds, ]
    mn[i]  <- et_sub[which.max(et_sub$amplitude), 'frequency_cpd']
    
  }
  
  return(mn)
  
}