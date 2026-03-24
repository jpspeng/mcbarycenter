\name{graph_quantiles}
\alias{graph_quantiles}
\title{Plot Estimated Quantiles with Confidence Bands}
\usage{
graph_quantiles(df)
}
\arguments{
\item{df}{A data frame containing \code{quantile}, \code{estimate}, \code{ci_lo}, and
\code{ci_hi} columns.}
}
\value{
A ggplot object.
}
\description{
Creates a ggplot2 line plot of the estimated quantile function with a
confidence ribbon.
}
