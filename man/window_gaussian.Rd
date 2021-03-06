% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frf_windows.R
\name{window_gaussian}
\alias{window_gaussian}
\title{window_gaussian}
\usage{
window_gaussian(n, sigma)
}
\arguments{
\item{n}{length of signal (integer)}

\item{sigma}{the standard deviation in periods (numeric)}
}
\value{
window
}
\description{
Gaussian window for FFT
https://en.wikipedia.org/wiki/Window_function
}
\examples{
window_gaussian(100, 0.3)
}
