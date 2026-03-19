\name{estimate_mixture_raw}
\alias{estimate_mixture_raw}
\title{Estimate an Empirical Raw Mixture}
\usage{
estimate_mixture_raw(
  df,
  id_col = "id",
  val_col = "x",
  x_thresh,
  weight_col = NULL
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{x_thresh}{Numeric threshold used to define successes.}

\item{weight_col}{Optional column name containing id-level weights.}
}
\value{
A standardized mixture data frame with columns \code{theta}, \code{g}, and
\code{cumul}.
}
\description{
Estimate a raw empirical mixture from thresholded binomial summaries.
}
