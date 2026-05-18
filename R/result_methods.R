#' Print an MCB barycenter result
#'
#' @param x An object of class `"mcbary_result"`.
#' @param ... Unused.
#'
#' @method print mcbary_result
#' @export
print.mcbary_result <- function(x, ...) {
  n_quantiles <- if (is.data.frame(x$res)) nrow(x$res) else NA_integer_
  n_mixtures <- if (is.list(x$mixtures)) length(x$mixtures) else NA_integer_

  cat("<mcbary_result>\n")
  cat("  Method:", x$method, "\n")
  cat("  Quantiles:", n_quantiles, "\n")
  cat("  Mixtures:", n_mixtures, "\n")

  invisible(x)
}

#' Print an empirical barycenter result
#'
#' @param x An object of class `"empbary_result"`.
#' @param ... Unused.
#'
#' @method print empbary_result
#' @export
print.empbary_result <- function(x, ...) {
  n_quantiles <- if (is.data.frame(x$res)) nrow(x$res) else NA_integer_

  cat("<empbary_result>\n")
  cat("  Quantiles:", n_quantiles, "\n")

  invisible(x)
}
