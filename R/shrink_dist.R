#' Empirical-Bayes Shrinkage of an Empirical CDF
#'
#' @param select_id Optional id value selecting one subject from `mcb_obj$data`.
#' @param x_vec Optional numeric vector of observed values.
#' @param mcb_obj An `mcbary()` result object containing `data` and `mixtures`.
#'
#' @return A data frame with columns `x`, `cdf_raw`, and `cdf_shrunk`.
#' @export
shrink_dist <- function(select_id = NULL, x_vec = NULL, mcb_obj) {
  if (missing(mcb_obj) || !is.list(mcb_obj)) {
    stop("`mcb_obj` must be a list.", call. = FALSE)
  }
  
  if (is.null(mcb_obj$mixtures) || !is.list(mcb_obj$mixtures) ||
      length(mcb_obj$mixtures) == 0) {
    stop("`mcb_obj` must contain a non-empty `mixtures` list.", call. = FALSE)
  }
  
  has_select_id <- !is.null(select_id)
  has_x_vec <- !is.null(x_vec)
  
  if (has_select_id == has_x_vec) {
    stop("Provide exactly one of `select_id` or `x_vec`.", call. = FALSE)
  }
  
  if (has_select_id) {
    if (is.null(mcb_obj$data) || !is.data.frame(mcb_obj$data) ||
        !all(c("id", "val") %in% names(mcb_obj$data))) {
      stop(
        "`mcb_obj$data` must be a data frame containing `id` and `val`.",
        call. = FALSE
      )
    }
    
    x_vec <- mcb_obj$data$val[mcb_obj$data$id == select_id]
    
    if (length(x_vec) == 0) {
      stop("`select_id` was not found in `mcb_obj$data$id`.", call. = FALSE)
    }
  }
  
  if (!is.numeric(x_vec) || anyNA(x_vec)) {
    stop("`x_vec` must be a numeric vector with no missing values.", call. = FALSE)
  }
  
  if (length(x_vec) == 0) {
    stop("`x_vec` must contain at least one value.", call. = FALSE)
  }
  
  mixtures <- mcb_obj$mixtures
  x_grid <- as.numeric(names(mixtures))
  
  if (anyNA(x_grid)) {
    stop("`mcb_obj$mixtures` must be a named list with numeric names.", call. = FALSE)
  }
  
  ord <- order(x_grid)
  mixtures <- mixtures[ord]
  x_grid <- x_grid[ord]
  
  n_obs <- length(x_vec)
  
  rows <- lapply(seq_along(mixtures), function(i) {
    mixture_df <- mixtures[[i]]
    
    if (!is.data.frame(mixture_df) || !all(c("theta", "g") %in% names(mixture_df))) {
      stop(
        sprintf(
          "Each element of `mcb_obj$mixtures` must contain `theta` and `g` (problem at `%s`).",
          names(mixtures)[i]
        ),
        call. = FALSE
      )
    }
    
    x_val <- x_grid[i]
    s <- sum(x_vec <= x_val)
    cdf_raw <- s / n_obs
    
    prior <- .standardize_mixture_df(
      theta = mixture_df$theta,
      g = mixture_df$g,
      add_cumul = FALSE
    )
    log_post_unnorm <- stats::dbinom(
      x = s,
      size = n_obs,
      prob = prior$theta,
      log = TRUE
    ) + log(prior$g)
    log_post_unnorm <- log_post_unnorm - max(log_post_unnorm)
    post_weights <- exp(log_post_unnorm)
    post_weights <- post_weights / sum(post_weights)
    cdf_shrunk <- sum(prior$theta * post_weights)
    
    data.frame(
      x = x_val,
      cdf_raw = cdf_raw,
      cdf_shrunk = cdf_shrunk
    )
  })
  
  do.call(rbind, rows)
}
