% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frf_spec_pgram.R
\name{spec_pgram}
\alias{spec_pgram}
\title{spec_pgram}
\usage{
spec_pgram(x, spans = NULL, kernel = NULL, taper = 0.1, pad = 0,
  fast = TRUE, demean = FALSE, detrend = TRUE, plot = TRUE,
  na.action = na.fail, ...)
}
\arguments{
\item{x}{univariate or multivariate time series.}

\item{spans}{vector of odd integers giving the widths of modified
    Daniell smoothers to be used to smooth the periodogram.}

\item{kernel}{alternatively, a kernel smoother of class
    \code{"tskernel"}.}

\item{taper}{specifies the proportion of data to taper.  A split
    cosine bell taper is applied to this proportion of the data at the
    beginning and end of the series.}

\item{pad}{proportion of data to pad. Zeros are added to the end of
    the series to increase its length by the proportion \code{pad}.}

\item{fast}{logical; if \code{TRUE}, pad the series to a highly composite
    length.}

\item{demean}{logical. If \code{TRUE}, subtract the mean of the
    series.}

\item{detrend}{logical. If \code{TRUE}, remove a linear trend from
    the series. This will also remove the mean.}

\item{plot}{plot the periodogram?}

\item{na.action}{\code{NA} action function.}

\item{...}{graphical arguments passed to \code{plot.spec}.}
}
\value{
array of spectra
}
\description{
calculate the spectra and cross-spectra
}
