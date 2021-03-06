% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frf_spec_welch.R
\name{spec_welch}
\alias{spec_welch}
\title{spec_welch}
\usage{
spec_welch(x, n_subsets = 10, overlap = 0.5, window = window_hann,
  demean = TRUE, detrend = TRUE, ...)
}
\arguments{
\item{x}{\code{numeric matrix} univariate or multivariate time series}

\item{n_subsets}{number of subsets (integer)}

\item{overlap}{amount of overlap 0.0 < overlap < 1.0 (numeric)}

\item{window}{window function to apply (function). 
The options are window_hann, window_hamming, window_blackman, 
window_nuttall, window_blackman_nuttall, window_blackman_harris}

\item{demean}{should the mean be removed from x prior to calculation (logical)}

\item{detrend}{should the x be detrended to calculation (logical)}

\item{...}{not used}
}
\value{
array of spectra
}
\description{
calculate the spectra and cross-spectra.  The scaling of individual series 
requires a closer look.
}
