#' Plot Estimated Quantiles with Optional Confidence Bands
#'
#' @param ... One or more result objects. Each may be either a data frame
#'   containing `quantile`, `estimate`, `ci_lo`, and `ci_hi`, or a list with
#'   a `res` data frame in that format, such as the output of [empb()] or
#'   [mcb()].
#' @param show_ci Logical; if TRUE (default), show confidence ribbons.
#'
#' @return A ggplot object.
#' @export
quantile_graph <- function(..., show_ci = TRUE) {
  dfs <- list(...)
  required_cols <- c("quantile", "estimate", "ci_lo", "ci_hi")

  as_quantile_df <- function(x, arg_label, required_cols) {
    if (is.list(x) && !is.data.frame(x) && "res" %in% names(x)) {
      x <- x$res
    }

    if (!is.data.frame(x)) {
      stop(sprintf("%s must be a data frame or a list with `$res`.", arg_label),
        call. = FALSE
      )
    }

    if (!all(required_cols %in% names(x))) {
      stop(
        sprintf(
          "%s must contain `%s`.",
          arg_label,
          paste(required_cols, collapse = "`, `")
        ),
        call. = FALSE
      )
    }

    x
  }
  
  if (length(dfs) == 0) {
    stop("At least one data frame must be provided.", call. = FALSE)
  }
  
  for (i in seq_along(dfs)) {
    dfs[[i]] <- as_quantile_df(dfs[[i]], sprintf("Argument %d", i), required_cols)
  }
  
  labels <- names(dfs)
  if (is.null(labels) || any(labels == "")) {
    labels <- as.character(seq_along(dfs))
  }
  
  df <- do.call(rbind, lapply(seq_along(dfs), function(i) {
    d <- dfs[[i]]
    d$.group <- factor(labels[i], levels = labels)
    d
  }))
  
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
#' Takes two quantile result objects and returns the estimated difference along
#' with delta-method standard errors and confidence intervals, assuming
#' independence.
#'
#' @param df1 A data frame containing columns `quantile`, `estimate`, and `se`,
#'   or a list with a `res` data frame in that format.
#' @param df2 A data frame containing columns `quantile`, `estimate`, and `se`,
#'   or a list with a `res` data frame in that format.
#' @param conf_level Confidence level for the interval. Default is 0.95.
#'
#' @return A data frame with columns:
#'   `quantile`, `estimate`, `se`, `ci_lo`, and `ci_hi`.
#' @export
quantile_diff <- function(df1, df2, conf_level = 0.95) {
  required_cols <- c("quantile", "estimate", "se")

  as_quantile_df <- function(x, arg_label, required_cols) {
    if (is.list(x) && !is.data.frame(x) && "res" %in% names(x)) {
      x <- x$res
    }

    if (!is.data.frame(x)) {
      stop(sprintf("%s must be a data frame or a list with `$res`.", arg_label),
        call. = FALSE
      )
    }

    if (!all(required_cols %in% names(x))) {
      stop(
        sprintf(
          "%s must contain `%s`.",
          arg_label,
          paste(required_cols, collapse = "`, `")
        ),
        call. = FALSE
      )
    }

    x
  }

  df1 <- as_quantile_df(df1, "`df1`", required_cols)
  df2 <- as_quantile_df(df2, "`df2`", required_cols)
  
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
