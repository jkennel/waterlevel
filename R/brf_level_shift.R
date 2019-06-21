
# THESE FUNCTIONS ARE BRITTLE ---------------------------------------------



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
find_level_shift <- function(x,
                             dep_var = 'val', 
                             time_var = 'datetime', 
                             time_interval = 1) {

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
      
      start <- tms[(wh[i])]
      end   <- tms[(wh[i]) + 1]
      
      subs[[i]] <- data.table(start = start, 
                              end   = end,
                              start_val = x[get(time_var) == start][[dep_var]],
                              end_val = x[get(time_var) == end][[dep_var]])
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
  
  gaps[, get_shift(x, recipe, start, end),
       by = list(midpoint = as.POSIXct(midpoint, tz = 'UTC'))]
  
  
}


#' gap_fill2
#'
#' @param x gap data.table 
#' @param y interpolated values
#'
#' @return data.table with gap values adjusted for range
#' @export
#'
gap_fill2 <- function(x, y) {
  
  x <- x[y, on = 'midpoint']
  
  x[, list(interp = stretch_interp(start_val[1], end_val[1], predict_adj[datetime %between% c(start[1], end[1])]),
           datetime = datetime[datetime %between% c(start[1], end[1])]), 
       by = list(midpoint)]

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
  
  fit <- lm(outcome~distributed_lag + lag_earthtide + datetime + level_shift -1, 
            dat,  
            x = FALSE, y = FALSE, tol = 1e-50)
  
  print(summary(fit))
  
  pt <- as.data.table(predict(fit, dat, type = 'terms', terms = 'level_shift'))
  pt[, predict := predict(fit, dat)]
  pt[, level_shift := level_shift - level_shift[1]]
  pt[, datetime := as.POSIXct(dat$datetime, origin = '1970-01-01', tz = 'UTC')]
  pt[, predict_adj := predict - level_shift]
  pt
}



#' get_intercept_stats
#'
#' @param x data.table of level shifts from regression
#'
#' @return data.table of shifts
#' @export
#'
get_intercept_stats <- function(x) {
  
  x <- x[, list(shifts = unique(level_shift),
                 min_datetime = min(datetime),
                 max_datetime = max(datetime)), 
          by = list(midpoint)]
  x[, shift_diff := c(0.0, diff(shifts)), by = midpoint]
  mids <- unique(x[, list(shift_datetime = midpoint, end_toss = midpoint)])
  rngs <- unique(x[, list(start = min_datetime, end = max_datetime, midpoint)])
  setkey(mids, shift_datetime, end_toss)
  setkey(rngs, start, end)
  grps <- foverlaps(rngs, mids)[, list(start, end, shift_datetime, midpoint, type = 'reg')]
  grps <- grps[, rbind(data.table(shift_datetime = start[1], midpoint = midpoint[1], type = 'non-reg'), .SD), by = list(start, end)]
  grps <- grps[, list(shift_datetime, midpoint, type)]
  
  x <- cbind(x, shift_datetime = grps$shift_datetime)
  
  #return(x[midpoint == shift_datetime])
  
  x <- x[, list(min = min(shift_diff),
            max = max(shift_diff),
            mean = mean(shift_diff),
            n = .N),
     by = shift_datetime]
  
  x <- x[!(min == 0.0 & max == 0)]
  return(x)
}



#' add_level_shifts
#'
#' @param x datetimes
#' @param y shift vector
#'
#' @return vector of shifts
#' @export
#'
add_level_shifts <- function(x, y) {
  
  out <- rep(0.0, length(x))
  
  for(i in 1:nrow(y)) {
    wh <- which(x > y$shift_datetime[i])
    mn <- as.numeric(y[i, list(min, max, mean)])
    mn <- mn[which.min(abs(mn))]
    out[wh] <- out[wh] + mn
  }
  
  
  return(out)
  
}


#' stretch_interp
#'
#' @param start_val 
#' @param end_val 
#' @param values 
#'
#' @return shifted vector of values
#' @export
#' 
stretch_interp <- function(start_val = NA, 
                           end_val = NA,
                           values) {
  
  shifta <- start_val - values[1]
  shiftb <- end_val - values[length(values)]
  
  if(is.na(start_val)){
    shifta <- shiftb
  } 
  if(is.na(end_val)){
    shiftb <- shifta
  } 
  
  values <- values + seq(shifta, shiftb, length.out = length(values))
  
  return(values)
}


