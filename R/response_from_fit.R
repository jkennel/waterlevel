#' response_from_fit
#'
#' @param fit model fit object, commonly lm
#'
#' @return tibble of model fit
#' 
#' @importFrom utils tail
#' @importFrom data.table as.data.table data.table
#' @importFrom stats coef
#' 
#' @export
#' 
response_from_fit <- function(fit){
  
  # hack for 'global variables NOTE
  variable = NULL
  
  term_group <- attr(attr(fit$model, 'terms'), "term.labels")
  term_subs  <- attr(fit$qr$qr, "assign")
  
  co  <- as.matrix(coef(fit))
  nms <- rownames(co)
  
  resp <- list()
  for (i in seq_along(term_group)) {
    
    tg <- term_group[i]
    nm <- nms[grep(tg, nms)]
    co_sub <- co[nms %in% nm, , drop = FALSE]
    
    if (grepl('distributed_lag', tg)) {
      resp[[i]] <- get_distributed_lag(co_sub, nm)
    } else if (grepl('lag_matrix', tg) | grepl('lag_earthtide', tg)) {
      resp[[i]] <- get_lag_matrix(co_sub, nm)
    } else if (grepl('earthtide', tg)) {
      resp[[i]] <- get_earthtide(co_sub, nm)
    } else if (grepl('harmonic', tg)) {
      resp[[i]] <- get_harmonic(co_sub, nm)
    } else {
      next
    }
    
    if(all(resp[[i]]$variable == 'V1')){
      resp[[i]][,  variable :=  tail(strsplit(nm, '_')[[1]], 1)]
    }
    
  }
  
  return(rbindlist(resp))
  
}



get_distributed_lag <- function(co, nm) {
  
  # hack for 'global variables NOTE
  type     = NULL
  value    = NULL
  variable = NULL
  
  
  nums <- lapply(strsplit(nm, '_'), function(x) {
    as.numeric(na.omit(as.numeric(x)))[1]
  })
  
  lags <- as.data.table(do.call(rbind, nums))
  setnames(lags, 'x')
  bl <- splines::ns(0:max(lags$x), 
                    knots = lags$x[-c(1, nrow(lags))],
                    Boundary.knots = c(0, max(lags$x)),
                    intercept = TRUE)
  lags <- data.table(x = 0:max(lags$x), bl %*% co)
  lags[, type := 'cumulative_response']
  lags <- melt(lags, id.vars = c('x', 'type'))
  lags[, value := cumsum(value), by = variable]
  
  return(lags)
}

get_lag_matrix <- function(co, nm) {
  
  # hack for 'global variables NOTE
  type = NULL
  
  nums <- lapply(strsplit(nm, '_'), function(x) {
    as.numeric(na.omit(as.numeric(x)))[1]
  })
  
  lags <- as.data.table(do.call(rbind, nums))
  setnames(lags, 'x')
  lags[, type := 'cumulative_response']
  lags <- cbind(lags, apply(co, 2, cumsum))
  lags <- melt(lags, id.vars = 1:2)
  
  return(lags)
}

get_harmonic <- function(co, nm) {
  
  nums <- lapply(strsplit(nm, '_'), function(x) {
    as.numeric(na.omit(as.numeric(x)))[1]
  })
  
  ranges <- unique(as.data.table(do.call(rbind,nums)))
  setnames(ranges, 'x')
  
  cc <- co[seq(1, nrow(co), 2), , drop = FALSE]
  ss <- co[seq(2, nrow(co), 2), , drop = FALSE]
  
  phase <- cbind(ranges, type = 'phase', atan2(ss, cc))
  amp <- cbind(ranges, type = 'amplitude', sqrt(ss^2 + cc^2))
  
  tmp <- rbind(melt(phase, id.vars = 1:2), melt(amp, id.vars = 1:2))
  
  
}

get_earthtide <- function(co, nm) {
  
  nums <- lapply(strsplit(nm, '_'), function(x) {
    as.numeric(na.omit(as.numeric(x)))[1:2]
  })
  
  ranges <- unique(do.call(rbind,nums))
  ranges <- data.table(x = get_main_frequency(ranges[,1], ranges[,2]))
  
  cc <- co[seq(1, nrow(co), 2), , drop = FALSE]
  ss <- co[seq(2, nrow(co), 2), , drop = FALSE]
  
  phase <- cbind(ranges, type = 'phase', atan2(ss, cc))
  amp   <- cbind(ranges, type = 'amplitude', sqrt(ss^2 + cc^2))
  
  tmp <- rbind(melt(phase, id.vars = 1:2), melt(amp, id.vars = 1:2))
  
}

