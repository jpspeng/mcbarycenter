\name{empbary}
\alias{empbary}
\title{Empirical Barycenter by Grouped Quantiles}
\usage{
empbary(
  df,
  id_col,
  val_col,
  alpha_grid = seq(0.01, 0.99, 0.01),
  quantile_type = 1
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{A string giving the column name that identifies groups.}

\item{val_col}{A string giving the column name containing numeric values.}

\item{alpha_grid}{A numeric vector of probabilities passed to \code{\link[stats:quantile]{quantile()}}.
Defaults to \code{seq(0.01, 0.99, 0.01)}.}

\item{quantile_type}{The \code{type} argument passed to \code{\link[stats:quantile]{quantile()}}.
Defaults to \code{1}.}
}
\value{
A list with components \code{res}, \code{cov}, and \code{data}. \code{res}
is a data frame with columns \code{quantile}, \code{estimate}, \code{se},
\code{ci_lo}, and \code{ci_hi}. \code{data} contains standardized
\code{id} and \code{val} columns from the input.
}
\description{
Computes group-specific empirical quantiles over a shared grid and then
averages those quantiles across groups.
}
