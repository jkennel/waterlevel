% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brf_level_shift.R
\name{gap_fill}
\alias{gap_fill}
\title{gap_fill}
\usage{
gap_fill(x, gaps, recipe, dep_var = "wl", time_var = "datetime",
  time_interval = 1L, buffer_start = 86400 * 4, buffer_end = 86400 *
  4, max_interp = 86400 * 7, include_level_shift = TRUE)
}
\arguments{
\item{x}{data.table of water levels}

\item{gaps}{data.table of gaps}

\item{recipe}{recipe to apply}

\item{dep_var}{the time variable name (character)}

\item{time_var}{the time variable name (character)}

\item{time_interval}{the delta t (numeric)}

\item{buffer_start}{how much buffer on each side of gap}

\item{buffer_end}{how much buffer on each side of gap}

\item{max_interp}{largest gap to interpolate}

\item{include_level_shift}{include level shift in adjustment}
}
\value{
data.table of predictions
}
\description{
gap_fill
}
