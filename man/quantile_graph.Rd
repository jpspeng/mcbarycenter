\name{graph_quantiles}
\alias{graph_quantiles}
\title{Plot Estimated Quantiles with Confidence Bands}
\usage{
graph_quantiles(..., show_ci = TRUE)
}
\arguments{
\item{...}{One or more result objects. Each may be either a data frame
containing \code{quantile} and \code{estimate}, and if \code{show_ci = TRUE},
also \code{ci_lo} and \code{ci_hi}; or an object of class
\code{"empbary_result"} or \code{"mcbary_result"}.}

\item{show_ci}{Logical; if TRUE (default), show confidence ribbons.}
}
\value{
A ggplot object.
}
\description{
Creates a ggplot2 line plot of one or more estimated quantile functions with a
confidence ribbon when available.
}
