\name{est_dist_alpha}
\alias{est_dist_alpha}
\title{Estimate the Distribution of an Alpha Quantile}
\usage{
est_dist_alpha(
  mixture_res,
  alpha,
  use_midpoint = TRUE,
  estimate_first_last = TRUE,
  use_isotonic = TRUE
)
}
\arguments{
\item{mixture_res}{A named list of mixture estimates.}

\item{alpha}{A single quantile level in \code{[0, 1]}.}

\item{use_midpoint}{Logical; whether to use interval midpoints between
adjacent thresholds.}

\item{estimate_first_last}{Logical controlling endpoint treatment.}

\item{use_isotonic}{Logical; whether to enforce monotonicity in the
cumulative mixture curves across threshold values.}
}
\value{
A data frame with columns \code{x}, \code{pmf}, and \code{cdf}.
}
\description{
Estimate the distribution of an alpha quantile from a sequence of mixture
distributions.
}
