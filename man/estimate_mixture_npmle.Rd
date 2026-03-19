\name{estimate_mixture_npmle}
\alias{estimate_mixture_npmle}
\title{Estimate a Mixture via NPMLE}
\usage{
estimate_mixture_npmle(
  df,
  id_col = "id",
  val_col = "x",
  x_thresh,
  tau = seq(from = 0, to = 1, by = 0.005),
  weight_col = NULL,
  mixsqp_control = list()
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{x_thresh}{Numeric threshold used to define successes.}

\item{tau}{Support grid for the estimated mixture.}

\item{weight_col}{Optional column name containing id-level weights.}

\item{mixsqp_control}{Optional control list passed to \code{\link[mixsqp:mixsqp]{mixsqp()}}.}
}
\value{
A standardized mixture data frame with columns \code{theta}, \code{g}, and
\code{cumul}.
}
\description{
Estimate a mixture by nonparametric maximum likelihood.
}
