\name{mcbary}
\alias{mcbary}
\title{Marginal-Constructed Barycenter Estimation and Inference}
\usage{
mcbary(
  df,
  id_col = "id",
  val_col = "x",
  method = c("spline", "npmle", "raw", "beta"),
  x_grid = NULL,
  cutpoints = NULL,
  bootstrap_samples = 100,
  alpha_grid = seq(0.01, 0.99, by = 0.01),
  use_midpoint = TRUE,
  estimate_endpoints = TRUE,
  use_isotonic_dist = TRUE,
  use_isotonic_output = FALSE,
  weight_col = NULL,
  progress = TRUE,
  ci_level = 0.95,
  inference = c("bootstrap", "analytic", "none"),
  ...
)
}
\arguments{
\item{df}{A data frame.}

\item{id_col}{Column name identifying groups. Defaults to \code{"id"}.}

\item{val_col}{Column name containing the observed values. Defaults to \code{"x"}.}

\item{method}{Mixture estimation method.}

\item{x_grid}{Threshold grid used when \code{cutpoints} is \code{NULL}. Must be
supplied if \code{cutpoints} is \code{NULL}.}

\item{cutpoints}{Optional number of evenly spaced thresholds. Either
\code{cutpoints} or \code{x_grid} must be supplied.}

\item{bootstrap_samples}{Number of bootstrap resamples. When set to
\code{1}, only the quantile estimate is reported and uncertainty columns
are returned as \code{NA}.}

\item{alpha_grid}{Grid of target quantile levels.}

\item{use_midpoint}{Logical; whether to use interval midpoints between
adjacent thresholds.}

\item{estimate_endpoints}{Logical controlling endpoint treatment.}

\item{use_isotonic_dist}{Logical; whether to enforce monotonicity in the
cumulative mixture curves across threshold values when constructing each
alpha-quantile distribution.}

\item{use_isotonic_output}{Logical; if \code{TRUE}, apply isotonic regression to
\code{res$estimate}, \code{res$ci_lo}, and \code{res$ci_hi}.}

\item{weight_col}{Optional column name containing id-level weights.}

\item{progress}{Logical; whether to display a bootstrap progress bar.}

\item{ci_level}{Confidence level.}

\item{inference}{Inference method. \code{"bootstrap"} (the default) performs
the unit-level cluster bootstrap, \code{"analytic"} uses Beta-specific
unit-level influence functions, and \code{"none"} skips uncertainty
estimation. Analytic inference requires \code{method = "beta"}, no weights,
interior \code{alpha_grid} values, \code{use_isotonic_dist = TRUE}, and
\code{use_isotonic_output = FALSE}.}

\item{...}{Additional arguments passed to \code{\link{estimate_all_mixtures}}.}
}
\value{
A list with components \code{res}, \code{cov}, \code{mixtures},
\code{method}, \code{data}, \code{influence}, and \code{inference}. \code{res}
is a data frame with columns
\code{quantile}, \code{estimate}, \code{estimate_bs}, \code{se},
\code{ci_lo}, \code{ci_hi}, \code{pct_ci_lo}, and \code{pct_ci_hi}.
\code{data} contains the original \code{id_col} and \code{val_col} columns from
the input. \code{cov} is on the barycenter estimator scale. For analytic
inference, \code{influence} contains integrated unit-level influence functions;
otherwise it is \code{NULL}.
}
\description{
Estimate marginal-constructed barycenter quantile curves with bootstrap,
Beta-specific analytic, or no uncertainty estimation.
}
