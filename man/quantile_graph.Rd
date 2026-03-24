\name{graph_quantiles}
\alias{graph_quantiles}
\title{Plot Estimated Quantiles with Confidence Bands}
\usage{
graph_quantiles(df)
}
\arguments{
\item{df}{A data frame containing \code{quantile} and \code{estimate}, and if confidence
bands are shown, also \code{ci_lo} and \code{ci_hi}.}
}
\value{
A ggplot object.
}
\description{
Creates a ggplot2 line plot of the estimated quantile function with a
confidence ribbon when available.
}
