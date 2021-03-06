% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{lag_matrix}
\alias{lag_matrix}
\title{lag_matrix}
\usage{
lag_matrix(x, lags, n_subset = 1L, n_shift = 0L, var_name = "lag")
}
\arguments{
\item{x}{\code{numeric vector} to lag}

\item{lags}{\code{integer vector} with the lags}

\item{n_subset}{\code{integer} subset every n_subset values}

\item{n_shift}{\code{integer} amount to shift results}

\item{var_name}{\code{character} name for the generated matrix columns}
}
\value{
matrix with lagged values
}
\description{
lag data and subset the results
}
