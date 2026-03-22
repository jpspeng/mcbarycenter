\name{quantile_diff}
\alias{quantile_diff}
\title{Compute Difference in Quantile Estimates}
\usage{
quantile_diff(df1, df2, conf_level = 0.95)
}
\arguments{
\item{df1}{A data frame containing columns \code{quantile}, \code{estimate}, and
\code{se}, or a list with a \code{res} data frame in that format.}

\item{df2}{A data frame containing columns \code{quantile}, \code{estimate}, and
\code{se}, or a list with a \code{res} data frame in that format.}

\item{conf_level}{Confidence level for the interval. Default is 0.95.}
}
\value{
A data frame with columns \code{quantile}, \code{estimate}, \code{se},
\code{ci_lo}, and \code{ci_hi}.
}
\description{
Takes two quantile result objects and returns the estimated difference along
with delta-method standard errors and confidence intervals, assuming
independence.
}
