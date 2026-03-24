\name{graph_mixtures}
\alias{graph_mixtures}
\title{Plot Mixing Distributions from an MCB Result}
\usage{
graph_mixtures(mcb_res, layout = c("combined", "individual"))
}
\arguments{
\item{mcb_res}{A result list containing a \code{mixtures} component, such as the
output of \code{\link{mcbary}}.}

\item{layout}{Either \code{"combined"} to plot all mixing distributions on one
graph, or \code{"individual"} to return one graph per mixing distribution.}
}
\value{
If \code{layout = "combined"}, a ggplot object. If
\code{layout = "individual"}, a named list of ggplot objects.
}
\description{
Plots the mixing distributions stored in an MCB result using \code{theta} on
the x-axis and \code{cumul} on the y-axis.
}
