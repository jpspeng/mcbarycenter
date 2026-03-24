\name{shrink_dist}
\alias{shrink_dist}
\title{Empirical-Bayes Shrinkage of an Empirical CDF}
\usage{
shrink_dist(select_id = NULL, x_vec = NULL, mcb_obj)
}
\arguments{
\item{select_id}{Optional id value selecting one subject from
\code{mcb_obj$data}.}

\item{x_vec}{Optional numeric vector of observed values.}

\item{mcb_obj}{An \code{mcbary()} result object containing \code{data} and
\code{mixtures}.}
}
\value{
A data frame with columns \code{x}, \code{cdf_raw}, and \code{cdf_shrunk}.
}
\description{
Builds an empirical CDF from one selected subject or an explicit vector of
values, then applies empirical-Bayes shrinkage at each threshold using the
corresponding mixture distribution as a prior.
}
