% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/be_visual.R
\name{be_visual_data}
\alias{be_visual_data}
\title{be_visual_data}
\usage{
be_visual_data(dat, dep = "wl", ind = "baro", be_tests = seq(0, 1,
  0.1), inverse = TRUE)
}
\arguments{
\item{dat}{data that has the independent and dependent variables (data.table)}

\item{dep}{name of the dependent variable column (character).  This is 
typically the name for the column holding your water level data.}

\item{ind}{name of the independent variable column (character).  This is 
typically the name for the column holding your barometric pressure data.}

\item{be_tests}{vector of barometric efficiencies to test (between 0 and 1) (numeric)}

\item{inverse}{whether the barometric relationship is inverse
(TRUE means that when the barometric pressure goes up the measured water
level goes down (vented transducer, depth to water), FALSE means that when 
the barometric pressure goes up so does the measured pressure 
(non-vented transducer)) (logical).}
}
\value{
data.table of barometric efficiency compensated datasets
}
\description{
Generate dataset for comparing barometric efficiency
}
\examples{
library(data.table)
be <- 0.43
x <- seq(0, 28*pi, pi / (12*12))

baro <- sin(x) + rnorm(length(x), sd = 0.04)
wl <- -sin(x) * be + rnorm(length(x), sd = 0.04)
dat <- data.table(datetime = as.POSIXct(x * 86400 / (2 * pi),
                                        origin = '1970-01-01', tz = 'UTC'),
                  wl = wl, baro = baro)
be_visual_data(dat)

}
