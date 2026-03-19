.validate_input_df <- function(df,
                               required_cols = c("id", "x"),
                               weight_col = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(weight_col) && !weight_col %in% names(df)) {
    stop(
      sprintf("`weight_col` = '%s' not found in `df`.", weight_col),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.standardize_input_df <- function(df,
                                  id_col = "id",
                                  val_col = "x",
                                  weight_col = NULL) {
  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("`id_col` must be a single non-empty string.", call. = FALSE)
  }

  if (!is.character(val_col) || length(val_col) != 1L || !nzchar(val_col)) {
    stop("`val_col` must be a single non-empty string.", call. = FALSE)
  }

  .validate_input_df(
    df,
    required_cols = c(id_col, val_col),
    weight_col = weight_col
  )

  out <- data.frame(
    id = df[[id_col]],
    x = df[[val_col]]
  )

  if (!is.null(weight_col)) {
    out$weight <- df[[weight_col]]
  }

  out
}

.validate_id_weights <- function(df, weight_col) {
  if (is.null(weight_col)) {
    return(invisible(TRUE))
  }

  tmp <- dplyr::group_by(df, id)
  tmp <- dplyr::summarise(
    tmp,
    n_weights = dplyr::n_distinct(.data[[weight_col]]),
    .groups = "drop"
  )

  bad_ids <- tmp$id[tmp$n_weights > 1]

  if (length(bad_ids) > 0) {
    stop(
      paste0(
        "Weights must be constant within id. ",
        "Found multiple weights for ", length(bad_ids), " id(s)."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.aggregate_binomial_data <- function(df, x_thresh, weight_col = NULL) {
  .validate_input_df(df, required_cols = c("id", "x"), weight_col = weight_col)
  .validate_id_weights(df, weight_col)

  if (!is.numeric(x_thresh) || length(x_thresh) != 1 || is.na(x_thresh)) {
    stop("`x_thresh` must be a single numeric value.", call. = FALSE)
  }

  out <- dplyr::group_by(df, id)

  if (is.null(weight_col)) {
    out <- dplyr::summarise(
      out,
      n = dplyr::n(),
      s = sum(.data$x <= x_thresh),
      .groups = "drop"
    )
  } else {
    out <- dplyr::summarise(
      out,
      n = dplyr::n(),
      s = sum(.data$x <= x_thresh),
      weight = dplyr::first(.data[[weight_col]]),
      .groups = "drop"
    )
  }

  if (nrow(out) == 0) {
    stop("No rows remain after aggregation.", call. = FALSE)
  }

  if (any(out$n <= 0)) {
    stop("Aggregated trial counts must be positive.", call. = FALSE)
  }

  if (any(out$s < 0 | out$s > out$n)) {
    stop("Aggregated successes must satisfy 0 <= s <= n.", call. = FALSE)
  }

  out
}

.standardize_mixture_df <- function(theta,
                                    g,
                                    sort_theta = TRUE,
                                    normalize = TRUE,
                                    add_cumul = TRUE) {
  theta <- as.numeric(theta)
  g <- as.numeric(g)

  if (length(theta) != length(g)) {
    stop("`theta` and `g` must have same length.", call. = FALSE)
  }

  if (anyNA(theta) || anyNA(g)) {
    stop("`theta` and `g` cannot contain NA.", call. = FALSE)
  }

  out <- data.frame(theta = theta, g = g)

  if (sort_theta) {
    out <- out[order(out$theta), , drop = FALSE]
  }

  if (normalize) {
    s <- sum(out$g)

    if (!is.finite(s) || s <= 0) {
      stop("Mixture weights must sum to a positive finite value.", call. = FALSE)
    }

    out$g <- out$g / s
  }

  if (any(out$g < -1e-12)) {
    stop("Mixture weights contain materially negative values.", call. = FALSE)
  }

  out$g[out$g < 0] <- 0

  if (normalize) {
    out$g <- out$g / sum(out$g)
  }

  if (add_cumul) {
    out$cumul <- cumsum(out$g)
  }

  rownames(out) <- NULL
  out
}

.discretize_beta_on_tau <- function(alpha, beta, tau) {
  tau <- sort(unique(as.numeric(tau)))

  if (length(tau) == 0 || anyNA(tau)) {
    stop("`tau` must be a non-empty numeric vector with no NA.", call. = FALSE)
  }

  if (any(tau < 0 | tau > 1)) {
    stop("`tau` must lie in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(beta) || length(beta) != 1 || is.na(beta) || beta <= 0) {
    stop("`beta` must be a single positive number.", call. = FALSE)
  }

  n_tau <- length(tau)
  edges <- numeric(n_tau + 1)

  if (n_tau >= 2) {
    edges[2:n_tau] <- (tau[-n_tau] + tau[-1]) / 2
  }

  edges[1] <- 0
  edges[n_tau + 1] <- 1
  edges <- pmin(pmax(edges, 0), 1)
  edges <- cummax(edges)

  g <- stats::pbeta(edges[-1], shape1 = alpha, shape2 = beta) -
    stats::pbeta(edges[-length(edges)], shape1 = alpha, shape2 = beta)

  g[g < 0 & g > -1e-14] <- 0
  g[g < 0] <- 0

  .standardize_mixture_df(theta = tau, g = g)
}


.fit_betabinomial_vgam <- function(
    df_bin,
    weight_col = NULL,
    start_grid = NULL,
    maxit = 100,
    trace = FALSE
) {
  if (!requireNamespace("VGAM", quietly = TRUE)) {
    stop("Package `VGAM` must be installed.", call. = FALSE)
  }
  
  if (!all(c("n", "s") %in% names(df_bin))) {
    stop("`df_bin` must contain columns `n` and `s`.", call. = FALSE)
  }
  
  if (anyNA(df_bin$n) || anyNA(df_bin$s)) {
    stop("`df_bin$n` and `df_bin$s` must not contain NA values.", call. = FALSE)
  }
  
  if (any(df_bin$n < 0) || any(df_bin$s < 0) || any(df_bin$s > df_bin$n)) {
    stop("Require 0 <= s <= n for all rows.", call. = FALSE)
  }
  
  if (!is.null(weight_col)) {
    if (!weight_col %in% names(df_bin)) {
      stop(
        sprintf("`weight_col` = '%s' not found in `df_bin`.", weight_col),
        call. = FALSE
      )
    }
    
    weights <- as.numeric(df_bin[[weight_col]])
    
    if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
      stop(
        "Weights supplied to `.fit_betabinomial_vgam()` must be finite and non-negative.",
        call. = FALSE
      )
    }
  } else {
    weights <- NULL
  }
  
  # Default starting values to try.
  # You can expand or replace this grid from the caller.
  if (is.null(start_grid)) {
    start_grid <- expand.grid(
      ishape1 = c(0.25, 0.5, 1, 2, 5, 10),
      ishape2 = c(0.25, 0.5, 1, 2, 5, 10),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  if (!all(c("ishape1", "ishape2") %in% names(start_grid))) {
    stop(
      "`start_grid` must be a data.frame with columns `ishape1` and `ishape2`.",
      call. = FALSE
    )
  }
  
  start_grid$ishape1 <- as.numeric(start_grid$ishape1)
  start_grid$ishape2 <- as.numeric(start_grid$ishape2)
  
  if (anyNA(start_grid$ishape1) || anyNA(start_grid$ishape2) ||
      any(!is.finite(start_grid$ishape1)) || any(!is.finite(start_grid$ishape2)) ||
      any(start_grid$ishape1 <= 0) || any(start_grid$ishape2 <= 0)) {
    stop(
      "All `ishape1` and `ishape2` values in `start_grid` must be finite and > 0.",
      call. = FALSE
    )
  }
  
  fit <- NULL
  last_error <- NULL
  attempts <- vector("list", nrow(start_grid))
  
  for (i in seq_len(nrow(start_grid))) {
    ishape1_i <- start_grid$ishape1[i]
    ishape2_i <- start_grid$ishape2[i]
    
    if (isTRUE(trace)) {
      message(
        sprintf(
          "Trying VGAM::vglm with ishape1 = %g, ishape2 = %g",
          ishape1_i, ishape2_i
        )
      )
    }
    
    res <- tryCatch(
      {
        fit_i <- VGAM::vglm(
          cbind(s, n - s) ~ 1,
          family = VGAM::betabinomialff(
            ishape1 = ishape1_i,
            ishape2 = ishape2_i
          ),
          data = df_bin,
          weights = weights,
          control = VGAM::vglm.control(maxit = maxit)
        )
        
        list(ok = TRUE, fit = fit_i, error = NULL)
      },
      error = function(e) {
        list(ok = FALSE, fit = NULL, error = conditionMessage(e))
      }
    )
    
    attempts[[i]] <- list(
      ishape1 = ishape1_i,
      ishape2 = ishape2_i,
      ok = res$ok,
      error = res$error
    )
    
    if (isTRUE(res$ok)) {
      fit <- res$fit
      break
    } else {
      last_error <- res$error
    }
  }
  
  if (is.null(fit)) {
    tried_txt <- paste(
      sprintf(
        "(ishape1=%g, ishape2=%g)",
        start_grid$ishape1, start_grid$ishape2
      ),
      collapse = ", "
    )
    
    stop(
      paste0(
        "All VGAM::vglm fits failed for betabinomialff. Tried starting values: ",
        tried_txt,
        if (!is.null(last_error)) paste0(". Last error: ", last_error) else ""
      ),
      call. = FALSE
    )
  }
  
  cf <- stats::coef(fit)
  
  if (length(cf) != 2) {
    stop(
      paste0(
        "Expected exactly 2 coefficients from intercept-only betabinomialff fit, got ",
        length(cf),
        ". Coefficient names were: ", paste(names(cf), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # betabinomialff(): eta1 = log(alpha), eta2 = log(beta)
  eta_alpha <- unname(cf[1])
  eta_beta  <- unname(cf[2])
  
  alpha <- exp(eta_alpha)
  beta  <- exp(eta_beta)
  
  # Replace Inf with max double
  if (!is.finite(alpha)) {
    warning("alpha overflowed; capping to .Machine$double.xmax")
    alpha <- .Machine$double.xmax
  }
  if (!is.finite(beta)) {
    warning("beta overflowed; capping to .Machine$double.xmax")
    beta <- .Machine$double.xmax
  }
  
  list(
    fit = fit,
    alpha = alpha,
    beta = beta,
    attempts = attempts
  )
}

.lookup_cumul_at_alpha <- function(mixture_df, alpha) {
  if (!is.data.frame(mixture_df)) {
    stop("Each mixture estimate must be a data.frame.", call. = FALSE)
  }

  if (!all(c("theta", "g") %in% names(mixture_df))) {
    stop(
      "Each mixture estimate must contain columns `theta` and `g`.",
      call. = FALSE
    )
  }

  mix <- mixture_df
  if (!"cumul" %in% names(mix)) {
    mix <- .standardize_mixture_df(theta = mix$theta, g = mix$g)
  } else {
    mix <- mix[order(mix$theta), , drop = FALSE]
  }

  idx <- which(mix$theta <= alpha)
  if (length(idx) == 0) {
    return(0)
  }

  mix$cumul[max(idx)]
}
