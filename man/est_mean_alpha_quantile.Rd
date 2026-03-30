\name{est_mean_alpha_quantile}
\alias{est_mean_alpha_quantile}
\title{Estimate Mean Alpha Quantiles from Mixture Distributions}
\usage{
est_mean_alpha_quantile(
  mixture_res,
  alpha,
  use_midpoint = TRUE,
  estimate_first_last = TRUE,
  use_isotonic_dist = TRUE
)
}
\arguments{
\item{mixture_res}{A named list of mixture estimates.}

\item{alpha}{A single quantile level in \code{[0, 1]}.}

\item{use_midpoint}{Logical; whether to use interval midpoints between
adjacent thresholds.}

\item{estimate_first_last}{Logical controlling endpoint treatment.}

\item{use_isotonic_dist}{Logical; whether to enforce monotonicity in the
cumulative mixture curves across threshold values.}
}
\value{
A numeric scalar giving the estimated mean alpha quantile.
}
\description{
Estimate a mean alpha quantile from a sequence of mixture distributions.
}
