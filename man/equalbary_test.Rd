\name{equalbary_test}
\alias{equalbary_test}
\title{Test Equality of Two Barycenter Estimates}
\usage{
equalbary_test(
  bary_res1,
  bary_res2,
  B = 5000,
  seed = NULL,
  variance_explained = NULL
)
}
\arguments{
\item{bary_res1}{First barycenter result, containing \code{res$estimate} and
\code{cov}.}

\item{bary_res2}{Second barycenter result, containing \code{res$estimate} and
\code{cov}.}

\item{B}{Number of Gaussian simulations.}

\item{seed}{Optional random seed.}

\item{variance_explained}{Optional number in \code{(0, 1]} giving the target
cumulative variance explained for PCA truncation in the Hotelling test. If
\code{NULL} (default), no PCA truncation is applied.}
}
\value{
A list with components \code{sup_norm}, \code{int_sq_norm}, and
\code{hotelling}. \code{sup_norm} and \code{int_sq_norm} each contain
\code{pval}, \code{T_obs}, and \code{T_sim}. \code{hotelling} contains
\code{pval}, \code{T_obs}, \code{df}, \code{mean_diff}, \code{cov_diff},
\code{variance_explained}, and \code{n_components}.
}
\description{
Uses a Gaussian approximation with simulated supremum and integrated squared
statistics to compare two barycenter estimates.
}
