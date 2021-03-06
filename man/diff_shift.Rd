% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/be_lag_difference.R
\name{diff_shift}
\alias{diff_shift}
\title{diff_shift}
\usage{
diff_shift(x, lag_space = 1)
}
\arguments{
\item{x}{vector to difference (numeric)}

\item{lag_space}{spacing for lags, useful for higher frequency monitoring
(integer).}
}
\value{
lagged differences
}
\description{
Calculate lagged differences padded with NA values
}
\examples{
diff_shift(1:100, 2)

}
