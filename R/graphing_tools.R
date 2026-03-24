#' Plot Estimated Quantiles with Optional Confidence Bands
#'
#' @param ... One or more result objects. Each may be either a data frame
#'   containing `quantile` and `estimate`, and if `show_ci = TRUE`, also
#'   `ci_lo` and `ci_hi`; or a list with a `res` data frame in that format,
#'   such as the output of [empbary()] or [mcbary()].
#' @param show_ci Logical; if TRUE (default), show confidence ribbons.
#'
#' @return A ggplot object.
#' @export
graph_quantiles <- function(..., show_ci = TRUE) {
  dfs <- list(...)
  arg_exprs <- as.list(substitute(list(...)))[-1]
  required_cols <- c("quantile", "estimate")

  if (show_ci) {
    required_cols <- c(required_cols, "ci_lo", "ci_hi")
  }

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

  fill_bind_dfs <- function(dfs) {
    all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))

    filled <- lapply(dfs, function(d) {
      missing_cols <- setdiff(all_cols, names(d))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {
          d[[col]] <- NA
        }
      }

      d[all_cols]
    })

    do.call(rbind, filled)
  }
  
  if (length(dfs) == 0) {
    stop("At least one data frame must be provided.", call. = FALSE)
  }
  
  for (i in seq_along(dfs)) {
    dfs[[i]] <- as_quantile_df(dfs[[i]], sprintf("Argument %d", i), required_cols)
  }
  
  labels <- names(dfs)
  if (is.null(labels)) {
    labels <- rep("", length(dfs))
  }
  
  missing_labels <- labels == ""
  if (any(missing_labels)) {
    labels[missing_labels] <- vapply(
      arg_exprs[missing_labels],
      FUN = function(expr) paste(deparse(expr), collapse = " "),
      FUN.VALUE = character(1)
    )
  }
  
  df <- fill_bind_dfs(lapply(seq_along(dfs), function(i) {
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
      color = NULL,
      fill = NULL
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

#' Plot Mixing Distributions from an MCB Result
#'
#' @param mcb_res A result list containing a `mixtures` component, such as the
#'   output of [mcbary()].
#' @param layout Either `"combined"` to plot all mixing distributions on one
#'   graph, or `"individual"` to return one graph per mixing distribution.
#' @param show_all Logical; if `TRUE`, show all mixing distributions. If
#'   `FALSE` (default) and there are more than 10, plot an approximately evenly
#'   spaced subset of 10.
#'
#' @return If `layout = "combined"`, a ggplot object. If
#'   `layout = "individual"`, a named list of ggplot objects.
#' @export
graph_mixtures <- function(mcb_res,
                           layout = c("combined", "individual"),
                           show_all = FALSE) {
  layout <- match.arg(layout)
  
  if (!is.list(mcb_res) || is.null(mcb_res$mixtures)) {
    stop(
      "`mcb_res` must be a list with a `mixtures` component.",
      call. = FALSE
    )
  }
  
  mixtures <- mcb_res$mixtures
  
  if (!is.list(mixtures) || length(mixtures) == 0) {
    stop("`mcb_res$mixtures` must be a non-empty list.", call. = FALSE)
  }
  
  if (!is.logical(show_all) || length(show_all) != 1L || is.na(show_all)) {
    stop("`show_all` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  
  mixture_names <- names(mixtures)
  if (is.null(mixture_names) || any(mixture_names == "")) {
    mixture_names <- as.character(seq_along(mixtures))
  }
  
  for (i in seq_along(mixtures)) {
    mix <- mixtures[[i]]
    if (!is.data.frame(mix) || !all(c("theta", "cumul") %in% names(mix))) {
      stop(
        sprintf(
          "Each mixture must be a data frame containing `theta` and `cumul` (problem at `%s`).",
          mixture_names[i]
        ),
        call. = FALSE
      )
    }
  }
  
  if (!show_all && length(mixtures) > 10) {
    keep_idx <- unique(round(seq(1, length(mixtures), length.out = 10)))
    mixtures <- mixtures[keep_idx]
    mixture_names <- mixture_names[keep_idx]
  }
  
  if (layout == "combined") {
    plot_df <- do.call(rbind, lapply(seq_along(mixtures), function(i) {
      mix <- mixtures[[i]]
      mix$.mixture <- factor(mixture_names[i], levels = mixture_names)
      mix
    }))
    
    return(
      ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = theta,
          y = cumul,
          color = .mixture,
          group = .mixture
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::labs(
          x = "alpha",
          y = "cdf",
          color = "Mixture"
        )
    )
  }
  
  plots <- lapply(seq_along(mixtures), function(i) {
    mix <- mixtures[[i]]
    ggplot2::ggplot(
      mix,
      ggplot2::aes(x = theta, y = cumul)
    ) +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = "alpha",
        y = "cdf",
        title = paste0("X_0: ", mixture_names[i])
      )
  })
  
  names(plots) <- mixture_names
  plots
}
