\name{est_all_quantiles}
\alias{est_all_quantiles}
\title{Estimate Mean Quantiles Across an Alpha Grid}
\usage{
est_all_quantiles(
  mixture_res,
  alpha_grid = seq(0.01, 0.99, by = 0.01),
  use_midpoint = TRUE,
  estimate_first_last = TRUE,
  use_isotonic_dist = TRUE
)
}
\arguments{
\item{mixture_res}{A named list of mixture estimates.}

\item{alpha_grid}{Grid of quantile levels.}

\item{use_midpoint}{Logical; whether to use interval midpoints between
adjacent thresholds.}

\item{estimate_first_last}{Logical controlling endpoint treatment.}

\item{use_isotonic_dist}{Logical; whether to enforce monotonicity in the
cumulative mixture curves across threshold values when constructing each
alpha-quantile distribution.}
}
\value{
A data frame with columns \code{quantile} and \code{estimate}.
}
\description{
Estimate mean quantiles across a grid of alpha values.
}
