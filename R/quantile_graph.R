#' Plot Estimated Quantiles with Optional Confidence Bands
#'
#' @param ... One or more data frames, each containing
#'   `quantile`, `estimate`, `ci_lo`, and `ci_hi`.
#' @param show_ci Logical; if TRUE (default), show confidence ribbons.
#'
#' @return A ggplot object.
#' @export
quantile_graph <- function(..., show_ci = TRUE) {
  dfs <- list(...)
  required_cols <- c("quantile", "estimate", "ci_lo", "ci_hi")
  
  if (length(dfs) == 0) {
    stop("At least one data frame must be provided.", call. = FALSE)
  }
  
  # Validate inputs
  for (i in seq_along(dfs)) {
    if (!is.data.frame(dfs[[i]])) {
      stop(sprintf("Argument %d is not a data frame.", i), call. = FALSE)
    }
    
    if (!all(required_cols %in% names(dfs[[i]]))) {
      stop(
        sprintf(
          "Data frame %d must contain `quantile`, `estimate`, `ci_lo`, and `ci_hi`.",
          i
        ),
        call. = FALSE
      )
    }
  }
  
  # Handle labels (use names if provided)
  labels <- names(dfs)
  if (is.null(labels) || any(labels == "")) {
    labels <- as.character(seq_along(dfs))
  }
  
  # Combine data
  df <- do.call(rbind, lapply(seq_along(dfs), function(i) {
    d <- dfs[[i]]
    d$.group <- factor(labels[i], levels = labels)
    d
  }))
  
  # Base plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = quantile,
      y = estimate,
      color = .group,
      group = .group
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Quantile",
      y = "Quantile level",
      color = "Series",
      fill = "Series"
    )
  
  # Add CI ribbons if requested
  if (show_ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = ci_lo,
          ymax = ci_hi,
          fill = .group
        ),
        alpha = 0.2,
        color = NA
      )
  }
  
  return(p)
}

#' Compute Difference in Quantile Estimates
#'
#' Takes two data frames containing quantile estimates and standard errors,
#' and returns the estimated difference along with delta-method standard errors
#' and confidence intervals, assuming independence.
#'
#' @param df1 A data frame containing columns `quantile`, `estimate`, and `se`.
#' @param df2 A data frame containing columns `quantile`, `estimate`, and `se`.
#' @param conf_level Confidence level for the interval. Default is 0.95.
#'
#' @return A data frame with columns:
#'   `quantile`, `estimate`, `se`, `ci_lo`, and `ci_hi`.
#' @export
quantile_diff <- function(df1, df2, conf_level = 0.95) {
  required_cols <- c("quantile", "estimate", "se")
  
  if (!is.data.frame(df1)) {
    stop("`df1` must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(df2)) {
    stop("`df2` must be a data frame.", call. = FALSE)
  }
  
  if (!all(required_cols %in% names(df1))) {
    stop(
      "`df1` must contain columns `quantile`, `estimate`, and `se`.",
      call. = FALSE
    )
  }
  if (!all(required_cols %in% names(df2))) {
    stop(
      "`df2` must contain columns `quantile`, `estimate`, and `se`.",
      call. = FALSE
    )
  }
  
  if (nrow(df1) != nrow(df2)) {
    stop("`df1` and `df2` must have the same number of rows.", call. = FALSE)
  }
  
  if (!identical(df1$quantile, df2$quantile)) {
    stop("`df1$quantile` and `df2$quantile` must match exactly.", call. = FALSE)
  }
  
  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha / 2)
  
  estimate <- df1$estimate - df2$estimate
  se <- sqrt(df1$se^2 + df2$se^2)
  
  out <- data.frame(
    quantile = df1$quantile,
    estimate = estimate,
    se = se,
    ci_lo = estimate - z * se,
    ci_hi = estimate + z * se
  )
  
  out
}