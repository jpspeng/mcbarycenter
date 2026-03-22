\name{equalbary_test}
\alias{equalbary_test}
\title{Test Equality of Two Barycenter Estimates}
\usage{
equalbary_test(bary_res1, bary_res2, n1, n2, B = 5000, seed = NULL)
}
\arguments{
\item{bary_res1}{First barycenter result, containing \code{res$estimate} and
\code{cov}.}

\item{bary_res2}{Second barycenter result, containing \code{res$estimate} and
\code{cov}.}

\item{n1}{Sample size associated with \code{bary_res1}.}

\item{n2}{Sample size associated with \code{bary_res2}.}

\item{B}{Number of Gaussian simulations.}

\item{seed}{Optional random seed.}
}
\value{
A list with components \code{pval}, \code{T_obs}, and \code{T_sim}.
}
\description{
Uses a Gaussian approximation with simulated supremum statistics to compare
two barycenter estimates.
}
