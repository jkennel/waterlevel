#' summarize_lm
#'
#' @param fit model fit from lm object
#'
#' @return data.table with fit results
#' @export 
#'
summarize_lm <- function(fit) {
  
  fit_summary <- summary(fit)
  
  data.table(
    n_obs         = stats::nobs(fit),
    n_coef        = fit_summary$df[3],
    df            = fit_summary$df[2],
    sigma         = fit_summary$sigma,
    r_squared     = fit_summary$r.squared,
    adj_r_squared = fit_summary$adj.r.squared,
    AIC           = stats::AIC(fit),
    BIC           = stats::BIC(fit),
    logLik        = as.numeric(stats::logLik(fit))
  )
  
}


#' summarize_coef
#'
#' @param fit model fit from lm object
#'
#' @return data.table with coefficient results
#' @export 
#'
summarize_coef <- function(fit) {
  
  co <- summary(fit)$coefficients
  
  data.table(
    name = names(fit$coefficients),
    co
  )
 
}

