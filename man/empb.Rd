\name{empb}
\alias{empb}
\title{Empirical Barycenter by Grouped Quantiles}
\usage{
empb(df, id_col, val_col, grid = seq(0, 1, 0.01), quantile_type = 3)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{A string giving the column name that identifies groups.}

\item{val_col}{A string giving the column name containing numeric values.}

\item{grid}{A numeric vector of probabilities passed to \code{\link[stats:quantile]{quantile()}}.
Defaults to \code{seq(0, 1, 0.01)}.}

\item{quantile_type}{The \code{type} argument passed to \code{\link[stats:quantile]{quantile()}}.
Defaults to \code{3}.}
}
\value{
A data frame with columns \code{quantile}, \code{estimate}, \code{se}, \code{ci_lo},
and \code{ci_hi}.
}
\description{
Computes group-specific empirical quantiles over a shared grid and then
averages those quantiles across groups.
}
