\name{estimate_mixture_efron}
\alias{estimate_mixture_efron}
\title{Estimate a Mixture via Efron-Style Deconvolution}
\usage{
estimate_mixture_efron(
  df,
  id_col = "id",
  val_col = "x",
  x_thresh,
  tau = seq(from = 0.005, to = 0.995, by = 0.005),
  pDegree = 7,
  c0 = 1,
  weight_col = NULL
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{x_thresh}{Numeric threshold used to define successes.}

\item{tau}{Support grid for the estimated mixture.}

\item{pDegree}{Degrees of freedom for the spline basis in \code{\link{deconv_fast}}.}

\item{c0}{Penalty constant passed to \code{\link{deconv_fast}}.}

\item{weight_col}{Optional column name containing id-level weights.}
}
\value{
A standardized mixture data frame with columns \code{theta}, \code{g}, and
\code{cumul}.
}
\description{
Estimate a mixture using the fast deconvolution approach.
}
