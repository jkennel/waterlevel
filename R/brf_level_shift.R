#' level_shift
#'
#' Create a factor column with the data subseted according to break points.
#' 
#' @param x vector of numeric
#' @param level_breaks break points in the data
#'
#' @return vector of type factor 
#' @export
#'
#' @examples
#' level_shift(1:100, c(2, 40, 86))
level_shift <- function(x,
                        level_breaks = NULL) {
  
  cut(x, 
      breaks = c(min(x, na.rm = TRUE),
                 level_breaks,
                 max(x, na.rm = TRUE)),
      include.lowest = TRUE)
}



#' find_level_shift
#'
#' @param x the data set (data.table)
#' @param dep_var the dependent variable name (character)
#' @param time_var the time variable name (character)
#' @param time_interval the delta t (numeric)
#'
#' @return data.table of gaps
#' @export
#'
find_level_shift <- function(x, dep_var = 'val', time_var = 'datetime', time_interval = 1) {
  
  start <- NULL
  end <- NULL
  midpoint <- NULL
  
  # remove times with no values
  tms <- x[!is.na(get(dep_var))][[time_var]]
  
  wh  <- which(diff(as.numeric(tms)) != time_interval)
  n   <- length(tms)
  
  subs <- list()
  
  if (length(wh) == 0) {
    subs <- data.table(start = tms[1], end = tms[n])
  } else {
    for (i in 1:(length(wh))) {
      subs[[i]] <- data.table(start = tms[(wh[i])], 
                              end   = tms[(wh[i])+1])
    }
  }
  
  subs <- rbindlist(subs)
  subs[, midpoint := as.POSIXct((as.numeric(start) + as.numeric(end)) / 2, 
                                origin = '1970-01-01', tz = 'UTC')]
  
  return(subs)
  
}



#' add_level_shift
#'
#' @param x the data set (data.table)
#' @param dep_var the time variable name (character)
#' @param time_var the time variable name (character)
#' @param time_interval the delta t (numeric)
#'
#' @return vector for grouping
#' @export
#'
add_level_shift <- function(x, dep_var = 'val', time_var = 'datetime', time_interval = 1) {
  
  ssss <- NULL
  eeee <- NULL
  level_shift <- NULL
  start <- NULL
  end <- NULL
  
  x <- setDT(x)
  # get the NA subsets
  gaps <- find_level_shift(x, dep_var, time_var, time_interval)
  
  x[, ssss := get(time_var)]
  x[, eeee := get(time_var)]
  
  mn <- min(x$datetime, na.rm = TRUE)
  mx <- max(x$datetime, na.rm = TRUE)
  
  shift_df <- data.table(start = c(mn, gaps$midpoint),
                         end   = c(gaps$midpoint, mx))
  
  # keep timezone info
  attr(shift_df$start, 'tzone') <- attr(mn, 'tzone')
  attr(shift_df$end, 'tzone') <- attr(mn, 'tzone')
  attr(x$ssss, 'tzone') <- attr(mn, 'tzone')
  attr(x$eeee, 'tzone') <- attr(mn, 'tzone')
  
  setkey(shift_df, start, end)
  setkey(x, ssss, eeee)
  
  
  x[, level_shift := sprintf("%04d", foverlaps(x, shift_df,
                                               type = 'within',
                                               which = TRUE,
                                               mult = 'first',
                                               nomatch = NULL))]
  
  as.factor(x$level_shift)
  
}





