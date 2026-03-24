\name{graph_mixtures}
\alias{graph_mixtures}
\title{Plot Mixing Distributions from an MCB Result}
\usage{
graph_mixtures(
  mcb_res,
  layout = c("combined", "individual"),
  show_all = FALSE
)
}
\arguments{
\item{mcb_res}{A result list containing a \code{mixtures} component, such as the
output of \code{\link{mcbary}}.}

\item{layout}{Either \code{"combined"} to plot all mixing distributions on one
graph, or \code{"individual"} to return one graph per mixing distribution.}

\item{show_all}{Logical; if \code{TRUE}, show all mixing distributions. If
\code{FALSE} (default) and there are more than 10, plot an approximately
evenly spaced subset of 10.}
}
\value{
If \code{layout = "combined"}, a ggplot object. If
\code{layout = "individual"}, a named list of ggplot objects.
}
\description{
Plots the mixing distributions stored in an MCB result using \code{alpha} on
the x-axis and \code{cdf} on the y-axis.
}
