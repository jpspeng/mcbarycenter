\name{mcbary}
\alias{mcbary}
\title{Monte Carlo Bootstrap for Quantile Estimation}
\usage{
mcbary(
  df,
  id_col = "id",
  val_col = "x",
  method = c("spline", "npmle", "raw", "beta"),
  x_grid = 1:10,
  cutpoints = NULL,
  bootstrap_samples = 100,
  alpha_grid = seq(0.01, 0.99, by = 0.01),
  use_midpoint = TRUE,
  fix_endpoints = TRUE,
  use_isotonic = FALSE,
  weight_col = NULL,
  progress = TRUE,
  ci_level = 0.95,
  ci_method = c("percentile", "wald"),
  ...
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{method}{Mixture estimation method.}

\item{x_grid}{Threshold grid used when \code{cutpoints} is \code{NULL}.}

\item{cutpoints}{Optional number of evenly spaced thresholds.}

\item{bootstrap_samples}{Number of bootstrap resamples.}

\item{alpha_grid}{Grid of target quantile levels.}

\item{use_midpoint}{Logical; whether to use interval midpoints between
adjacent thresholds.}

\item{fix_endpoints}{Logical controlling endpoint treatment.}

\item{use_isotonic}{Logical; if \code{TRUE}, apply isotonic regression to
\code{res$estimate}, \code{res$ci_lo}, and \code{res$ci_hi}.}

\item{weight_col}{Optional column name containing id-level weights.}

\item{progress}{Logical; whether to display a bootstrap progress bar.}

\item{ci_level}{Confidence level.}

\item{ci_method}{Confidence interval method.}

\item{...}{Additional arguments passed to \code{\link{estimate_all_mixtures}}.}
}
\value{
A list with components \code{res}, \code{cov}, \code{mixtures},
\code{method}, and \code{data}. \code{res} is a data frame with columns
\code{quantile}, \code{estimate}, \code{se}, \code{ci_lo}, and
\code{ci_hi}. \code{data} contains the original \code{id_col} and
\code{val_col} columns from the input.
}
\description{
Estimate quantile curves with Monte Carlo bootstrap uncertainty.
}
