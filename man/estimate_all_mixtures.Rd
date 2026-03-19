\name{estimate_all_mixtures}
\alias{estimate_all_mixtures}
\title{Estimate Mixtures Across a Threshold Grid}
\usage{
estimate_all_mixtures(
  df,
  id_col = "id",
  val_col = "x",
  method = c("spline", "npmle", "raw", "beta"),
  x_grid = 1:10,
  weight_col = NULL,
  ...
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{method}{Mixture estimation method.}

\item{x_grid}{Grid of thresholds to evaluate.}

\item{weight_col}{Optional column name containing id-level weights.}

\item{...}{Additional arguments passed through to the chosen estimator.}
}
\value{
A named list of mixture estimates, one per threshold in \code{x_grid}.
}
\description{
Estimate mixtures repeatedly over a grid of thresholds while reusing repeated
results when the thresholded counts do not change.
}
