
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
  
  br <- c(min(x, na.rm = TRUE),
          level_breaks,
          max(x, na.rm = TRUE))
  attr(br, "tzone") <- "UTC"
  cut(x, 
      breaks = br,
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
  
  # find locations of gaps
  wh  <- which(diff(as.numeric(tms)) != time_interval)
  n   <- length(tms)
  
  subs <- list()
  
  if (length(wh) == 0) {
    subs <- data.table(start = tms[1], end = tms[n])
  } else {
    for (i in 1:(length(wh))) {
      
      start <- tms[(wh[i])]
      end   <- tms[(wh[i]) + 1]
      
      # make sure it is even 
      if((as.numeric(start)-as.numeric(end)) %% 2 != 0) {
        end <- end + 1
      }
      
      subs[[i]] <- data.table(start     = start, 
                              end       = end,
                              start_val = x[get(time_var) == start][[dep_var]],
                              end_val   = x[get(time_var) == end][[dep_var]])
    }
  }
  
  subs <- rbindlist(subs)

  subs[, midpoint := as.POSIXct(round((as.numeric(start) + as.numeric(end)) / 2), 
                                origin = '1970-01-01', tz = 'UTC')]
  
  return(subs)
  
}





#' gap_fill
#'
#' @param x data.table of water levels
#' @param gaps data.table of gaps
#' @param recipe recipe to apply
#' @param dep_var the time variable name (character)
#' @param time_var the time variable name (character)
#' @param time_interval the delta t (numeric)
#' @param buffer_start how much buffer on each side of gap
#' @param buffer_end how much buffer on each side of gap
#' @param max_interp largest gap to interpolate
#' @param include_level_shift include level shift in adjustment
#'
#' @return data.table of predictions
#' 
#' @export
#'
gap_fill <- function(x, 
                     gaps, 
                     recipe, 
                     dep_var = 'wl',
                     time_var = 'datetime',
                     time_interval = 1L, 
                     buffer_start = 86400 * 4,
                     buffer_end = 86400 * 4,
                     max_interp = 86400 * 7,
                     include_level_shift = TRUE) {
  
  # gaps <- find_level_shift(x, dep_var = dep_var, time_var = time_var,
  #                          time_interval = time_interval)
  
  if (nrow(gaps) == 0) {
    stop('no gaps to fill')
    return(NULL)
  }
  
  gaps <- gaps[as.numeric(end) - as.numeric(start) < max_interp]
  
  if (nrow(gaps) == 0) {
    warning('no gaps to fill')
    return(NULL)
  }
  
  gaps[, start_reg := start - (buffer_start)]
  gaps[, end_reg   := end + (buffer_end)]
  
  x <- copy(x)
  
  gaps <- gaps[, get_fit_summary(x, recipe, start_reg, end_reg, start, end),
       by = list(midpoint = as.POSIXct(midpoint, tz = 'UTC'),
                 start_gap = as.POSIXct(start, tz = 'UTC'), 
                 end_gap = as.POSIXct(end, tz = 'UTC'),
                 start_val = start_val, 
                 end_val = end_val)]
  
  
  if(include_level_shift) {
    
    gaps <- get_intercept_stats(gaps)
  
    gaps[, start_val := start_val - cumsum(c(0, sh[-length(sh)]))]
    gaps[, end_val := end_val - cumsum(sh)]
    
  }
  
  gaps[, to_exclude := nrow(predict_adj[[1]]) == 0, by = midpoint]
  
  gaps[to_exclude == FALSE, `:=` (predict_adj =
                 list(data.table(datetime = predict_adj[[1]][['datetime']],
                                 adj = stretch_interp(start_val[1],
                                                      end_val[1],
                                                      predict_adj[[1]][['predict_adj']])))),
       by = midpoint]


  
}





#' get_fit_summary
#'
#' @param x data.table of water levels
#' @param recipe recipe to apply
#' @param start_reg start subset time
#' @param end_reg end subset time
#' @param start_gap start subset time
#' @param end_gap end subset time
#'
#' @return data.table summary
#' @export
#'
#'
#'
get_fit_summary <- function(x, recipe, start_reg, end_reg, start_gap, end_gap) {
  
  y   <- x[between(datetime, start_reg, end_reg)]
  
  dat <- recipe %>%
    prep(training = y) %>%
    portion()
  
  form <- formula_from_recipe(recipe = recipe)
  fit <- lm(form, 
            dat[!is.na(as.vector(dat$outcome)),],  
            x = FALSE, y = FALSE, tol = 1e-50)
  
  # get predictions without level shift
  dat_sub <- dat[is.na(as.vector(dat$outcome)),]
  pt  <- data.table(level_shift = as.numeric(predict(fit, dat_sub, type = 'terms', terms = 'level_shift')),
                    datetime = as.POSIXct(as.numeric(dat_sub[['datetime']]), origin = '1970-01-01', tz = 'UTC'))
  pt[, predict := predict(fit, dat_sub)]
  pt[, level_shift := level_shift - level_shift[1]]
  pt[, predict_adj := predict - level_shift]
  pt <- pt[between(datetime, start_gap, end_gap)]
  pt[, level_shift := NULL]
  pt[, predict := NULL]
  
  out <- summarize_lm(fit)
  #out[, `:=` (rank = fit$rank)]
  out[, start_reg := start_reg]
  out[, end_reg := end_reg]
  out[, `:=` (recipe = list(recipe))]
  #out[, `:=` (coef = list(summarize_coef(fit)))]
  out[, `:=` (formula = list(form))]
  out[, `:=` (terms = list(terms))]
  out[, `:=` (pivot = list(qr(fit)$pivot))]
  out[, `:=` (predict_adj = list(pt))]
  
  out
  
}




pree <- function(x, fit_dt) {
  
  dat <- recipe %>%
    prep(training = x) %>%
    portion()
  
  terms <- labels(terms(fit_dt$form))
  
  term_names <- lapply(select(dat, terms), colnames)
  
  lapply(term_names, function(x) {
    select
  })
  
  terms <- names(dat)
  
}






#' formula_from_recipe
#' guess formula from recipe
#'
#' @param recipe recipe to apply
#'
#' @return formula
#' @export
#'
formula_from_recipe <- function(recipe) {
  
  terms <- na.omit(
    vapply(recipe$steps, FUN = function(x) {
      as.character(x$role)
    }, FUN.VALUE = 'character'))
  terms <- terms[terms != 'exclude']
  rhs <- paste(terms, collapse = '+')
  
  if('level_shift' %in% terms) {
    add <- '-1'
    rhs <- paste0(rhs, add)
  }
  
  as.formula(paste0('outcome', '~', rhs))
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
  
  x[, list(interp = stretch_interp(start_val[1], 
                                   end_val[1], 
                                   predict_adj[between(datetime, start[1], end[1])]),
           datetime = datetime[between(datetime, start[1], end[1])]), 
       by = list(midpoint)]

}


#' get_shift
#'
#' @param x data.table of water levels
#' @param recipe recipe to apply
#' @param start start subset time
#' @param end end subset time
#'
#' @return data.table of predictions
#' @export
#'
get_shift <- function(x, recipe, start, end) {
  
  y <- x[between(datetime, start, end)]

  dat <- recipe %>%
    prep(training = y) %>%
    portion()
  
  fit <- lm(outcome~distributed_lag + lag_earthtide + datetime + level_shift -1, 
            dat,  
            x = FALSE, y = FALSE, tol = 1e-50)
  
  #print(summary(fit))
  
  pt <- as.data.table(predict(fit, dat, type = 'terms', terms = 'level_shift'))
  pt[, predict := predict(fit, dat)]
  pt[, level_shift := level_shift - level_shift[1]]
  pt[, datetime := as.POSIXct(dat$datetime, origin = '1970-01-01', tz = 'UTC')]
  pt[, predict_adj := predict - level_shift]
  pt
  
}


# x is the result of gap_fill
#' get_level_shift_coef
#'
#' @param x data.table of level shifts from regression
#'
#' @return data.table of shifts
#' @export
#'
get_level_shift_coef <- function(x) {
  
  co <- x[, (coefs[[1]]), by = list(start_reg, end_reg, start_gap, end_gap, midpoint)]
  # print(co, 100)
  co <- co[grep('level_shift', name)]
  setnames(co, 'Estimate', 'level_shift')
  
  # there shouldn't be any NA values but remove them if there are
  co <- co[!is.na(level_shift)]
  
  co 
  
}


# x is the result of gap_fill
#' level_shift_to_datetime
#'
#' @param x character string that includes the date
#'
#' @return datetime
#' @export
#'
level_shift_to_datetime <- function(x) {
  
  # Remove non-numeric
  x <- gsub("[^0-9.-]", "", x)
  datetime <- paste(scan(text = x, sep = '.', quiet = TRUE), collapse = ' ')
  
  as.POSIXct(datetime, 
             format = '%Y %m %d %H %M %S', 
             origin = '1970-01-01',
             tz = 'UTC')
  
}


#' get_intercept_stats
#'
#' @param x data.table of level shifts from regression
#'
#' @return data.table of shifts
#' @export
#'
get_intercept_stats <- function(x) {
  
  # get level_shift from regression
  y <- get_level_shift_coef(x)
  
  # the name is for grouping coefficients
  y <- y[, list(shifts = unique(level_shift),
                name), 
          by = list(midpoint)]
  
  
  y[, shift_diff := c(0.0, diff(shifts)), by = midpoint]
  
  # 
  # mids <- unique(x[, list(shift_datetime = midpoint, 
  #                         end_toss = midpoint)])
  # 
  # rngs <- unique(x[, list(start = min_datetime,
  #                         end = max_datetime, 
  #                         midpoint)])
  # 
  # 
  # setkey(mids, shift_datetime, end_toss)
  # setkey(rngs, start, end)
  # grps <- foverlaps(rngs, mids)[, list(start, end, shift_datetime, midpoint, type = 'reg')]
  # 
  # grps <- grps[, rbind(data.table(shift_datetime = start[1], midpoint = midpoint[1], type = 'non-reg'), .SD),
  #              by = list(start, end)]
  # 
  # grps <- grps[, list(shift_datetime, midpoint, type)]
  # 
  # x <- cbind(x, shift_datetime = grps$shift_datetime)
  
  #return(x[midpoint == shift_datetime])
  
  # the first value should always be 0 and can be removed
  
  y <- y[, .SD[-1], midpoint]
  
  small_mag <- function(x) { x[which.min(abs(x))]}
  
  y <- y[, list(min  = min(shift_diff),
                max  = max(shift_diff),
                mean = mean(shift_diff),
                sd   = sd(shift_diff),
                sh   = small_mag(c(range(shift_diff), mean(shift_diff))),
                n = .N,
                midpoint = level_shift_to_datetime(name)),
     by = list(name)]
  
  y <- y[!(min == 0.0 & max == 0)]
  
  y <- x[y, on = 'midpoint']
  
  return(y)
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
    wh <- which(x > y$midpoint[i])
    mn <- as.numeric(y[i, list(min, max, mean)])
    mn <- mn[which.min(abs(mn))]
    out[wh] <- out[wh] + mn
  }
  
  
  return(out)
  
}


#' stretch_interp
#'
#' @param start_val starting value to match
#' @param end_val end value to match
#' @param values set of values
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


# library(data.table)
# dat <- data.table(x = rnorm(100), y = rnorm(100), z = rnorm(100))
# dat[30:40, x := NA_real_]
# dat[50:70, z := NA_real_]
# dat[10:12, y := NA_real_]
# 
# fill_lm <- function(dat,
#                     fill_cols, 
#                     partial = FALSE) {
# 
#   dat <- copy(dat)
# 
#   if(length(fill_cols) < 2) {
#     stop('fill_cols must have length 2 or greater')
#   }
# 
#   for(i in seq_along(fill_cols)) {
# 
#     fill_col <- fill_cols[i]
#     ind_cols <- setdiff(fill_cols, fill_col)
# 
#     for(j in seq_along(ind_cols)) {
# 
#       for(k in seq_along(gaps)) {
# 
#         to_fill <- dat[is.na(get(fill_col))]
#         ind_col <- ind_cols[j]
#         n_na    <- sum(is.na(to_fill[[ind_col]]))
# 
#         if(n_na == 0) {
# 
#           form    <- as.formula(paste0(fill_col, '~', ind_col))
#           fit     <- lm(form, dat)
#           dat[is.na(get(fill_col)), (fill_col) := predict(fit, to_fill)]
#           break
#         }
#       }
#     }
#   }
# 
#   return(dat)
# 
# }


# library(data.table)
# library(recipes)
# library(waterlevel)
# 
# dt <- data.table(x = rnorm(100000), y = rnorm(100000), z = rnorm(100000), datetime = 1:100000)
# 
# recipe <- recipe(y+z~x, dt[1]) %>%
#   step_distributed_lag(x, knots = c(0, 2, 5)) %>%
#   step_mutate(test = x^2, role = 'squared')
# 
# tmp <- fit_gaps(dt, recipe = recipe)
# dtnew <- data.table(x = rnorm(1000), y = rnorm(1000), z = rnorm(1000))
# aa <- predict_gaps(dtnew, fits = tmp)


#' find_gaps
#'
#' @param x the data set (data.table)
#' @param dep_var the dependent variable name (character)
#' @param time_var the time variable name (character)
#' @param time_interval the delta t (numeric)
#'
#' @return data.table of gaps
#' @export
#'
find_gaps <- function(x,
                      dep_var = 'val', 
                      time_var = 'datetime', 
                      time_interval = 1) {
  
  midpoint <- NULL
  
  # remove times with no values
  tms <- x[!is.na(get(dep_var))][[time_var]]
  
  # find locations of gaps
  wh  <- which(diff(as.numeric(tms)) != time_interval)
  n   <- length(tms)
  
  subs <- list()
  
  if (length(wh) == 0) {
    subs <- data.table(start = tms[1], end = tms[n])
  } else {
    for (i in 1:(length(wh))) {
      
      start <- tms[(wh[i])]
      end   <- tms[(wh[i]) + 1]
      
      subs[[i]] <- data.table(start     = start, 
                              end       = end,
                              start_val = x[get(time_var) == start][[dep_var]],
                              end_val   = x[get(time_var) == end][[dep_var]])
    }
  }
  
  subs <- rbindlist(subs)
  
  subs[, midpoint := as.POSIXct((as.numeric(start) + as.numeric(end)) / 2, 
                                origin = '1970-01-01', tz = 'UTC')]
  
  return(subs)
  
}


#' fit_gaps
#'
#' @param x the data set (data.table)
#' @param recipe the recipe to apply
#' @param time_var the time variable name (character)
#'
#' @return data.table of fit results
#' @export
#'
fit_gaps <- function(x, recipe, time_var = 'datetime') {
  
  
  form <- formula_from_recipe(recipe)
  
  fit_dat <- recipe %>% 
    prep(training = x) %>% 
    portion()
  
  fit <- lm(formula = form, data = fit_dat,  
            x = FALSE, y = FALSE, tol = 1e-50)
  
  out <- summarize_lm(fit)
  #print(out)
  out[, `:=` (recipe      = list(recipe))]
  out[, `:=` (start_train = min(x$datetime))]
  out[, `:=` (end_train   = max(x$datetime))]
  out
  
}


#' predict_gaps
#'
#' @param x the data set (data.table)
#' @param fits data.table of fits
#' @param term_labels names of the group (character)
#'
#' @return data.table of predictions
#' @export
#'
predict_gaps <- function(x, fits, term_labels = NULL) {
  
  recipe      <- fits[['recipe']][[1]]
  
  if(is.null(term_labels)) {
    term_labels <- fits[['term_labels']][[1]]
  }
  
  coefs       <- fits[['coefs']][[1]]
  setkey(coefs, name)
  fit_dat <- recipe %>% 
    prep(training = x) %>% 
    portion()
  
  out <- matrix(NA_real_, 
                nrow = nrow(fit_dat), 
                ncol = length(term_labels))
  
  for (i in seq_along(term_labels)) {
    
    term_label <- term_labels[i]
    nms <- paste0(term_label, colnames(fit_dat[[term_label]]))
    
    out[, i] <- fit_dat[[term_label]] %*% as.matrix(coefs[nms, list(Estimate)])
  
  }

  out <- as.data.table(out)
  setnames(out, term_labels)
  
  out
  
}



