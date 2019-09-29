#' summarize_lm
#'
#' @param fit model fit from lm object
#' @param ... additional arguments
#'
#' @return data.table with fit results
#' @export 
#'
summarize_lm <- function (fit, ...) {
  UseMethod("summarize_lm", fit)
}

#' @export 
summarize_lm.mlm <- function(fit) {
  
  fit_summary <- summary(fit)
  
  if('mlm' %in% class(fit)) {
    len <- length(fit_summary)
  } else {
    len <- 1
  }
  
  out <- data.table(
    variable      = dimnames(fit$model$outcome)[[2]],
    form          = format(formula(fit)),
    n_obs         = rep(stats::nobs(fit), len),
    n_coef        = vapply(fit_summary, function(x) x$df[3], FUN.VALUE = 1),
    df            = vapply(fit_summary, function(x) x$df[2], FUN.VALUE = 2),
    sigma         = vapply(fit_summary, function(x) x$sigma, FUN.VALUE = 1.0),
    r_squared     = vapply(fit_summary, function(x) x$r.squared, FUN.VALUE = 1.0),
    adj_r_squared = vapply(fit_summary, function(x) x$adj.r.squared, FUN.VALUE = 1.0),
    AIC           = AIC_mlm(fit)$AIC
  )
  
  tt <- terms(fit)
  out[, `:=` (pivot = .(qr(fit)$pivot))]
  out[, `:=` (has_intercept =  attr(tt, "intercept") > 0L)]
  out[, `:=` (term_labels =  list(attr(tt, "term.labels")))]
  out[, `:=` (coefs = lapply(fit_summary, summarize_coef)) ]
  
  out
  
}


#' @export 
summarize_lm.lm <- function(fit) {
  
  fit_summary <- summary(fit)
  
  out <- data.table(
    variable      = dimnames(fit$model$outcome)[[2]],
    form          = Reduce(paste, deparse(formula(fit), width.cutoff = 500)),
    n_obs         = stats::nobs(fit),
    n_coef        = fit_summary$df[3],
    df            = fit_summary$df[2],
    sigma         = fit_summary$sigma,
    r_squared     = fit_summary$r.squared,
    adj_r_squared = fit_summary$adj.r.squared,
    AIC           = stats::AIC(fit)
  )
  
  tt <- terms(fit)
  out[, `:=` (pivot = .(qr(fit)$pivot))]
  out[, `:=` (has_intercept =  attr(tt, "intercept") > 0L)]
  out[, `:=` (term_labels =  list(attr(tt, "term.labels")))]
  out[, `:=` (coefs = list(summarize_coef(fit_summary))) ]
  
  out
  
}



#' summarize_coef
#'
#' @param fit model fit from lm object
#'
#' @return data.table with coefficient results
#' @export 
#'
summarize_coef <- function(fit) {
  
  co <- fit$coefficients
  
  if(is.vector(co)) {
    co <- data.table(
      name = names(co),
      co
    )
  } else {
    co <- data.table(
      name = rownames(fit$coefficients),
      co
    )
  }
  
  co
}


#' Computation of AIC for mlm objects
#'
#' Extends the \code{extractAIC} method from the \pkg{stats} package to handle 
#' multi-predictand linear models (objects of class mlm). FROM paleocar
#'
#' @param fit An object of class mlm.
#' @param scale The estimate of the error variance.
#' \code{scale = 0} indicates that it is to be 
#' estimated by maximum likelihood.
#' @param k Numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @return A list of length 2 giving
#' \itemize{
#'   \item{\code{df}  The 'equivalent degrees of freedom' for the fitted model \code{fit}.}
#'   \item{\code{AIC}  A vector of the (generalized) Akaike Information Criterion for the fits.}
#' }
#' @importFrom stats complete.cases deviance
AIC_mlm <- function(fit, scale = 0, k = 2){
  n <- nrow(fit$residuals)
  edf <- n - fit$df.residual
  RSS <- stats::deviance(fit)
  dev <- if (scale > 0) 
    RSS/scale - n
  else n * log(RSS/n)
  list(df=edf, AIC=dev + k * edf)
}

