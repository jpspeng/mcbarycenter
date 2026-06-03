#' Plot Estimated Quantiles with Optional Confidence Bands
#'
#' @param ... One or more result objects. Each may be either a data frame
#'   containing `quantile` and `estimate`, and if `show_ci = TRUE`, also
#'   `ci_lo` and `ci_hi`; or an object of class `"empbary_result"` or
#'   `"mcbary_result"`.
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
    if (inherits(x, "empbary_result") || inherits(x, "mcbary_result")) {
      x <- x$res
    }

    if (!is.data.frame(x)) {
      stop(
        sprintf(
          "%s must be a data frame or an `empbary_result`/`mcbary_result` object.",
          arg_label
        ),
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
      x = "Quantile level",
      y = "Quantile value",
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

#' Plot Individual Empirical Quantile Curves
#'
#' Computes empirical quantile curves for each individual in a raw data frame
#' and plots them either statically with ggplot2 or interactively with plotly.
#'
#' @param df A data frame.
#' @param id_col A string giving the column name that identifies individuals.
#' @param val_col A string giving the column name containing numeric values.
#' @param alpha_grid A numeric vector of probabilities passed to [stats::quantile()].
#'   Defaults to `seq(0.01, 0.99, 0.01)`.
#' @param quantile_type The `type` argument passed to [stats::quantile()].
#'   Defaults to `1`.
#' @param interactive Logical; if `TRUE`, return an interactive plotly object.
#'   If `FALSE`, return a ggplot object.
#'
#' @return A ggplot object if `interactive = FALSE`; otherwise a plotly htmlwidget.
#' @export
graph_individual_quantiles <- function(
    df,
    id_col,
    val_col,
    alpha_grid = seq(0.01, 0.99, 0.01),
    quantile_type = 1,
    interactive = TRUE
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }

  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("`id_col` must be a single non-empty string.", call. = FALSE)
  }

  if (!is.character(val_col) || length(val_col) != 1L || !nzchar(val_col)) {
    stop("`val_col` must be a single non-empty string.", call. = FALSE)
  }

  if (!id_col %in% names(df)) {
    stop("`id_col` must name a column in `df`.", call. = FALSE)
  }

  if (!val_col %in% names(df)) {
    stop("`val_col` must name a column in `df`.", call. = FALSE)
  }

  if (!is.numeric(df[[val_col]])) {
    stop("`val_col` must refer to a numeric column.", call. = FALSE)
  }

  if (!is.numeric(alpha_grid) || length(alpha_grid) == 0L || anyNA(alpha_grid)) {
    stop(
      "`alpha_grid` must be a non-empty numeric vector without missing values.",
      call. = FALSE
    )
  }

  if (any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("`alpha_grid` values must lie between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(quantile_type) || length(quantile_type) != 1L || is.na(quantile_type)) {
    stop("`quantile_type` must be a single numeric value.", call. = FALSE)
  }

  if (!is.logical(interactive) || length(interactive) != 1L || is.na(interactive)) {
    stop("`interactive` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  split_vals <- split(df[[val_col]], df[[id_col]])
  split_vals <- Filter(function(x) any(!is.na(x)), split_vals)

  if (length(split_vals) == 0L) {
    stop(
      "`df` must contain at least one group with non-missing values.",
      call. = FALSE
    )
  }

  quantile_df <- do.call(
    rbind,
    lapply(names(split_vals), function(id) {
      values <- split_vals[[id]]
      data.frame(
        id = id,
        quantile = alpha_grid,
        estimate = stats::quantile(
          values[!is.na(values)],
          probs = alpha_grid,
          type = quantile_type,
          names = FALSE
        ),
        row.names = NULL
      )
    })
  )

  quantile_df$id <- factor(quantile_df$id, levels = unique(quantile_df$id))

  if (!interactive) {
    return(
      ggplot2::ggplot(
        quantile_df,
        ggplot2::aes(
          x = quantile,
          y = estimate,
          color = id,
          group = id
        )
      ) +
        ggplot2::geom_line(alpha = 0.5, linewidth = 0.4) +
        ggplot2::labs(
          x = "Quantile level",
          y = "Quantile value",
          color = "ID"
        )
    )
  }

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "`plotly` must be installed when `interactive = TRUE`.",
      call. = FALSE
    )
  }

  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop(
      "`htmlwidgets` must be installed when `interactive = TRUE`.",
      call. = FALSE
    )
  }

  ids <- levels(quantile_df$id)
  palette <- grDevices::hcl.colors(length(ids), palette = "Dynamic")
  names(palette) <- ids

  p <- plotly::plot_ly()

  for (id in ids) {
    d <- quantile_df[quantile_df$id == id, , drop = FALSE]

    p <- plotly::add_trace(
      p,
      data = d,
      x = ~quantile,
      y = ~estimate,
      type = "scatter",
      mode = "lines",
      name = id,
      customdata = rep(id, nrow(d)),
      hovertemplate = paste0(
        "ID: ", id,
        "<br>Quantile: %{x:.2f}",
        "<br>Value: %{y:.3f}",
        "<extra></extra>"
      ),
      line = list(
        color = palette[[id]],
        width = 1
      ),
      opacity = 0.5,
      showlegend = TRUE
    )
  }

  p <- plotly::layout(
    p,
    xaxis = list(title = "Quantile level"),
    yaxis = list(title = "Quantile value"),
    title = list(text = "Click a line to highlight its ID"),
    hovermode = "closest"
  )

  htmlwidgets::onRender(
    p,
    "
    function(el, x) {
      var fadedOpacity = 0.5;
      var fadedWidth = 1;
      var activeOpacity = 1.0;
      var activeWidth = 3;

      function resetTraces(gd) {
        for (var i = 0; i < gd.data.length; i++) {
          Plotly.restyle(gd, {
            opacity: fadedOpacity,
            'line.width': fadedWidth
          }, [i]);
        }
      }

      el.on('plotly_click', function(evt) {
        if (!evt || !evt.points || evt.points.length === 0) {
          return;
        }

        var traceIndex = evt.points[0].curveNumber;
        var id = evt.points[0].data.name;

        resetTraces(el);

        Plotly.restyle(el, {
          opacity: activeOpacity,
          'line.width': activeWidth
        }, [traceIndex]);

        Plotly.relayout(el, {
          title: {text: 'Selected ID: ' + id}
        });
      });

      el.on('plotly_doubleclick', function(evt) {
        resetTraces(el);
        Plotly.relayout(el, {
          title: {text: 'Click a line to highlight its ID'}
        });
      });
    }
    "
  )
}

#' Plot Mixing Distributions from an MCB Result
#'
#' @param mcb_res An object of class `"mcbary_result"`, such as the output of
#'   [mcbary()].
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
  
  if (!inherits(mcb_res, "mcbary_result")) {
    stop(
      "`mcb_res` must be an object of class `mcbary_result`.",
      call. = FALSE
    )
  }

  if (is.null(mcb_res$mixtures)) {
    stop("`mcb_res$mixtures` must be a non-empty list.", call. = FALSE)
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
