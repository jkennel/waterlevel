#' set_level_shift
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
#' set_level_shift(1:100, c(2, 40, 86))
set_level_shift <- function(x,
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
  # subs[1, start = min(x$datetime, na.rm = TRUE)]
  # subs[nrow(subs), end = max(x$datetime, na.rm = TRUE)]
  subs[, midpoint := as.POSIXct((as.numeric(start) + as.numeric(end)) / 2, 
                                origin = '1970-01-01', tz = 'UTC')]
  
  return(subs)
  
}






#' gap_fill
#'
#' @param x data.table of water levels
#' @param recipe recipe to apply
#' @param dep_var the time variable name (character)
#' @param time_var the time variable name (character)
#' @param time_interval the delta t (numeric)
#' @param buffer_start how much buffer on each side of gap
#' @param buffer_end how much buffer on each side of gap
#'
#' @return data.table of predictions
#' @export
#'
gap_fill <- function(x, recipe, dep_var = 'wl', time_var = 'datetime',
                     time_interval = 1L, 
                     buffer_start = 86400*4,
                     buffer_end = 86400*4) {

  gaps <- find_level_shift(x, dep_var = dep_var, time_var = time_var,
                           time_interval = time_interval)

  gaps[, start := start - (buffer_start)]
  gaps[, end := end + (buffer_end)]

  x <- copy(x)

  g <- gaps[, get_shift(x, recipe, start, end),
            by = list(midpoint = as.POSIXct(midpoint, tz = 'UTC'))]

}

#' get_shift
#'
#' @param x data.table of water levels
#' @param recipe recipe to apply
#' @param start start subset time
#' @param end end subset tiem
#'
#' @return data.table of predictions
#' @export
#'
get_shift <- function(x, recipe, start, end) {
  
  y <- x[datetime %between% c(start, end)]

  dat <- recipe %>%
    prep(training = y) %>%
    portion()
  
  fit <- lm(outcome~distributed_lag + lag_earthtide + datetime + level_shift, 
            dat,  
            x = FALSE, y = FALSE, tol = 1e-50)
  
  print(summary(fit))
  
  pt <- as.data.table(predict(fit, dat, type = 'terms', terms = 'level_shift'))
  pt[, predict := predict(fit, dat)]
  pt[, level_shift := level_shift - level_shift[1]]
  pt[, datetime := as.POSIXct(dat$datetime, origin = '1970-01-01', tz = 'UTC')]
  
  pt
}



