\name{deconv_fast}
\alias{deconv_fast}
\title{Fast Deconvolution Estimator}
\usage{
deconv_fast(
  tau,
  X,
  y,
  Q,
  P,
  n = 40,
  family = c("Poisson", "Normal", "Binomial"),
  ignoreZero = TRUE,
  deltaAt = NULL,
  c0 = 1,
  scale = TRUE,
  pDegree = 5,
  aStart = 1,
  ...
)
}
\arguments{
\item{tau}{Support grid for the latent parameter.}

\item{X}{Observed data. For \code{family = "Binomial"}, this should be a two
column matrix or data frame containing \code{(n_i, s_i)}.}

\item{y}{Optional count vector. If \code{P} and \code{Q} are not supplied,
\code{y} is constructed from \code{X} when needed.}

\item{Q}{Optional design matrix for the exponential family prior.}

\item{P}{Optional sampling matrix.}

\item{n}{Number of bins/support points used when constructing \code{P}.}

\item{family}{Sampling family. One of \code{"Poisson"}, \code{"Normal"}, or
\code{"Binomial"}.}

\item{ignoreZero}{Logical; for Poisson models, whether to condition on
positive counts.}

\item{deltaAt}{Optional point mass location for the Normal family.}

\item{c0}{Penalty constant.}

\item{scale}{Logical; whether to center and scale the spline basis.}

\item{pDegree}{Degrees of freedom for the natural spline basis.}

\item{aStart}{Initial value for the optimizer. A scalar is recycled to the
number of columns of \code{Q}.}

\item{...}{Additional arguments passed to \code{\link[stats:nlm]{nlm()}}.}
}
\value{
A list containing the optimizer output, working matrices, covariance
estimates, summary statistics, and internal objective helpers.
}
\description{
A streamlined deconvolution routine adapted from Efron's
\pkg{deconvolveR} implementation, with vectorized construction of the
likelihood components.
}
