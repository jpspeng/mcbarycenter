#' Plot Estimated Quantiles with Confidence Bands
#'
#' Creates a ggplot2 line plot of the estimated quantile function with a
#' confidence ribbon.
#'
#' @param df A data frame containing `quantile`, `estimate`, `ci_lo`, and
#'   `ci_hi` columns.
#'
#' @return A ggplot object.
#' @export
quantile_graph <- function(df) {
  required_cols <- c("quantile", "estimate", "ci_lo", "ci_hi")

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }

  if (!all(required_cols %in% names(df))) {
    stop(
      "`df` must contain columns `quantile`, `estimate`, `ci_lo`, and `ci_hi`.",
      call. = FALSE
    )
  }

  ggplot2::ggplot(df, ggplot2::aes(x = quantile, y = estimate)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      alpha = 0.2
    ) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Quantile", y = "Quantile level")
}
