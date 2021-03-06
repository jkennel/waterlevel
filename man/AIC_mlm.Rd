% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brf_summarize_lm.R
\name{AIC_mlm}
\alias{AIC_mlm}
\title{Computation of AIC for mlm objects}
\usage{
AIC_mlm(fit, scale = 0, k = 2)
}
\arguments{
\item{fit}{An object of class mlm.}

\item{scale}{The estimate of the error variance.
\code{scale = 0} indicates that it is to be 
estimated by maximum likelihood.}

\item{k}{Numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
}
\value{
A list of length 2 giving
\itemize{
  \item{\code{df}  The 'equivalent degrees of freedom' for the fitted model \code{fit}.}
  \item{\code{AIC}  A vector of the (generalized) Akaike Information Criterion for the fits.}
}
}
\description{
Extends the \code{extractAIC} method from the \pkg{stats} package to handle 
multi-predictand linear models (objects of class mlm). FROM paleocar
}
