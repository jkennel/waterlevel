% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/be_least_squares.R
\name{be_least_squares}
\alias{be_least_squares}
\title{be_least_squares}
\usage{
be_least_squares(dat, dep = "wl", ind = "baro", inverse = TRUE,
  return_model = FALSE)
}
\arguments{
\item{dat}{data that has the independent and dependent variables (data.table)}

\item{dep}{name of the dependent variable column (character).  This is 
typically the name for the column holding your water level data.}

\item{ind}{name of the independent variable column (character).  This is 
typically the name for the column holding your barometric pressure data.}

\item{inverse}{whether the barometric relationship is inverse
(TRUE means that when the barometric pressure goes up the measured water
level goes down (vented transducer, depth to water), FALSE means that when 
the barometric pressure goes up so does the measured pressure 
(non-vented transducer)) (logical).}

\item{return_model}{whether to return the lm model or just the barometric/loading
efficiency (logical).#'}
}
\value{
barometric efficiency calculated by least squares
}
\description{
Calculate the barometric efficiency by using least squares.
}
\examples{
library(data.table)
datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
                       as.POSIXct("2016-01-05 12:00:00"), by='hour' )
baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
wl <- -0.4 * baro
dat <- data.table(baro, wl, datetime)

be_least_squares(dat)

}
