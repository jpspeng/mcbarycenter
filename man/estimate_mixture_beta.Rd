\name{estimate_mixture_beta}
\alias{estimate_mixture_beta}
\title{Estimate a Beta Mixture Approximation}
\usage{
estimate_mixture_beta(
  df,
  id_col = "id",
  val_col = "x",
  x_thresh,
  tau = seq(0.005, 0.995, by = 0.005),
  weight_col = NULL
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{x_thresh}{Numeric threshold used to define successes.}

\item{tau}{Support grid for the estimated mixture.}

\item{weight_col}{Optional column name containing id-level weights.}
}
\value{
A standardized mixture data frame with columns \code{theta}, \code{g}, and
\code{cumul}. The fitted beta parameters are attached in the \code{beta_mle}
attribute when available.
}
\description{
Estimate a discretized beta-mixture approximation on the supplied support grid.
}
