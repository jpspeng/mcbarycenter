#' Fast Deconvolution Estimator
#'
#' A streamlined deconvolution routine adapted from Efron's
#' \pkg{deconvolveR} implementation, with vectorized construction of the
#' likelihood components.
#'
#' @param tau Support grid for the latent parameter.
#' @param X Observed data. For `family = "Binomial"`, this should be a two
#'   column matrix or data frame containing `(n_i, s_i)`.
#' @param y Optional count vector. If `P` and `Q` are not supplied, `y` is
#'   constructed from `X` when needed.
#' @param Q Optional design matrix for the exponential family prior.
#' @param P Optional sampling matrix.
#' @param n Number of bins/support points used when constructing `P`.
#' @param family Sampling family. One of `"Poisson"`, `"Normal"`, or
#'   `"Binomial"`.
#' @param ignoreZero Logical; for Poisson models, whether to condition on
#'   positive counts.
#' @param deltaAt Optional point mass location for the Normal family.
#' @param c0 Penalty constant.
#' @param scale Logical; whether to center and scale the spline basis.
#' @param pDegree Degrees of freedom for the natural spline basis.
#' @param aStart Initial value for the optimizer. A scalar is recycled to the
#'   number of columns of `Q`.
#' @param obs_weights Optional observation weights. For Poisson and Normal this
#'   should have length `nrow(P)` after construction; for Binomial it should
#'   have length equal to the number of rows of `X` (or of supplied `logP`).
#' @param ... Additional arguments passed to [stats::nlm()].
#'
#' @return A list containing the optimizer output, working matrices, covariance
#'   estimates, summary statistics, and internal objective helpers.
#' @export
deconv_fast <- function(tau, X, y, Q, P, n = 40,
                        family = c("Poisson", "Normal", "Binomial"),
                        ignoreZero = TRUE, deltaAt = NULL, c0 = 1,
                        scale = TRUE, pDegree = 5, aStart = 1.0,
                        obs_weights = NULL, compute_stats = TRUE, ...) {

  family <- match.arg(family)
  if (!is.null(deltaAt) && (family != "Normal")) {
    stop("deltaAt applies only to Normal.", call. = FALSE)
  }
  
  softmax <- function(z) {
    zmax <- max(z)
    ez <- exp(z - zmax)
    ez / sum(ez)
  }
  
  logsoftmax <- function(z) {
    zmax <- max(z)
    zz <- z - zmax
    zz - log(sum(exp(zz)))
  }
  
  logsumexp <- function(v) {
    vmax <- max(v)
    vmax + log(sum(exp(v - vmax)))
  }
  
  y_vec_make <- function(y, nrow_p) {
    if (length(y) == 1L && y == 1) {
      rep(1.0, nrow_p)
    } else {
      as.numeric(y)
    }
  }
  
  logP <- NULL

  if (missing(Q) && missing(P)) {
    if (missing(tau)) {
      stop("`tau` must be supplied when `P` and `Q` are not provided.",
        call. = FALSE
      )
    }

    if (missing(X)) {
      stop("`X` must be supplied when `P` and `Q` are not provided.",
        call. = FALSE
      )
    }

    if (family == "Poisson") {
      support_of_x <- if (ignoreZero) seq_len(n) else seq.int(0, n - 1)
      P <- outer(support_of_x, tau, function(x, lam) stats::dpois(x, lam))
      if (ignoreZero) {
        P <- sweep(P, 2, 1 - exp(-tau), "/")
      }

      if (missing(y)) {
        y <- tabulate(match(X, support_of_x), nbins = length(support_of_x))
      }

      Q1 <- splines::ns(tau, pDegree)
      if (scale) {
        Q1 <- scale(Q1, center = TRUE, scale = FALSE)
        Q1 <- sweep(Q1, 2, sqrt(colSums(Q1 * Q1)), "/")
      }
      Q <- cbind(1.0, Q1)
      Q <- sweep(Q, 2, sqrt(colSums(Q * Q)), "/")
      
      if (is.null(obs_weights)) {
        obs_weights <- rep(1.0, nrow(P))
      }

    } else if (family == "Normal") {
      r <- round(range(X), 1)
      x_bin <- seq(from = r[1], to = r[2], length.out = n)
      lo <- x_bin[-length(x_bin)]
      hi <- x_bin[-1]

      P <- outer(lo, tau, stats::pnorm) - outer(hi, tau, stats::pnorm)

      intervals <- findInterval(X, x_bin, all.inside = TRUE)
      y <- tabulate(
        pmax.int(1L, pmin.int(intervals, n - 1L)),
        nbins = n - 1L
      )

      Q1 <- splines::ns(tau, pDegree)
      if (scale) {
        Q1 <- scale(Q1, center = TRUE, scale = FALSE)
        Q1 <- sweep(Q1, 2, sqrt(colSums(Q1 * Q1)), "/")
      }
      if (!is.null(deltaAt)) {
        I0 <- as.numeric(abs(tau - deltaAt) < 1e-10)
        Q <- cbind(I0, Q1)
      } else {
        Q <- Q1
      }
      
      if (is.null(obs_weights)) {
        obs_weights <- rep(1.0, nrow(P))
      }

    } else {
      if (is.data.frame(X)) {
        X <- as.matrix(X)
      }
      if (!is.matrix(X) || ncol(X) != 2) {
        stop(
          "For family = 'Binomial', `X` must have two columns: (n_i, s_i).",
          call. = FALSE
        )
      }
      storage.mode(X) <- "double"

      logP <- vapply(
        tau,
        function(w) stats::dbinom(X[, 2], size = X[, 1], prob = w, log = TRUE),
        numeric(nrow(X))
      )
      logP <- as.matrix(logP)
      P <- matrix(NA_real_, nrow = nrow(logP), ncol = ncol(logP))
      y <- 1

      Q <- splines::ns(tau, pDegree)
      if (scale) {
        Q <- scale(Q, center = TRUE, scale = FALSE)
        Q <- sweep(Q, 2, sqrt(colSums(Q * Q)), "/")
      }
      
      if (is.null(obs_weights)) {
        obs_weights <- rep(1.0, nrow(logP))
      }
    }
  } else {
    if (!missing(X) || missing(y) || missing(P) || missing(Q)) {
      stop("If supplying `P` and `Q`, pass `P`, `Q`, and `y` together and omit `X`.",
        call. = FALSE
      )
    }
    
    if (family == "Binomial") {
      logP <- as.matrix(P)
      storage.mode(logP) <- "double"
      P <- matrix(NA_real_, nrow = nrow(logP), ncol = ncol(logP))
      if (is.null(obs_weights)) {
        obs_weights <- rep(1.0, nrow(logP))
      }
    } else if (is.null(obs_weights)) {
      obs_weights <- rep(1.0, nrow(P))
    }
  }
  
  Q <- as.matrix(Q)
  storage.mode(Q) <- "double"
  
  if (family != "Binomial") {
    P <- as.matrix(P)
    storage.mode(P) <- "double"
  }

  p <- ncol(Q)
  if (length(aStart) == 1L) {
    aStart <- rep(aStart, p)
  }
  if (length(aStart) != p) {
    stop("`aStart` has wrong length.", call. = FALSE)
  }
  
  obs_weights <- as.numeric(obs_weights)
  if (family == "Binomial") {
    if (length(obs_weights) != nrow(logP)) {
      stop(
        "`obs_weights` must have length equal to the number of Binomial observations.",
        call. = FALSE
      )
    }
  } else if (length(obs_weights) != nrow(P)) {
    stop("`obs_weights` must have length `nrow(P)`.", call. = FALSE)
  }
  if (anyNA(obs_weights) || any(obs_weights < 0)) {
    stop("`obs_weights` must be a non-missing, non-negative numeric vector.",
      call. = FALSE
    )
  }

  penalty_grad <- function(a) {
    aa <- sqrt(sum(a * a))
    if (aa == 0) {
      rep(0, length(a))
    } else {
      c0 * a / aa
    }
  }

  penalty_hessian <- function(a) {
    aa <- sqrt(sum(a * a))
    if (aa == 0) {
      matrix(0, nrow = length(a), ncol = length(a))
    } else {
      (c0 / aa) * (diag(length(a)) - (a %o% a) / (aa * aa))
    }
  }

  loglik <- function(a) {
    Qa <- as.vector(Q %*% a)
    g <- softmax(Qa)
    if (family == "Binomial") {
      logg <- logsoftmax(Qa)
      logf <- apply(logP, 1, function(row) logsumexp(row + logg))
      yv <- y_vec_make(y, nrow(logP))
      wy <- obs_weights * yv
      val <- -sum(wy * logf) + c0 * sqrt(sum(a * a))

      R <- exp(logP - logf)
      u <- as.vector(crossprod(R, wy))
      sy <- sum(wy)
    } else {
      f <- as.vector(P %*% g)
      f <- pmax(f, 1e-300)
      yv <- y_vec_make(y, nrow(P))
      wy <- obs_weights * yv

      val <- -sum(wy * log(f)) + c0 * sqrt(sum(a * a))

      u <- as.vector(t(P) %*% (wy / f))
      sy <- sum(wy)
    }

    v <- g * (u - sy)
    grad <- -as.vector(crossprod(Q, v)) + penalty_grad(a)

    attr(val, "gradient") <- grad
    val
  }

  opt <- stats::nlm(f = loglik, p = aStart, gradtol = 1e-10, ...)
  mle <- opt$estimate

  if (!is.logical(compute_stats) || length(compute_stats) != 1L ||
      is.na(compute_stats)) {
    stop("`compute_stats` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  if (!compute_stats) {
    g <- softmax(as.vector(Q %*% mle))

    return(list(
      mle = mle,
      Q = Q,
      P = if (family == "Binomial") NULL else P,
      logP = if (family == "Binomial") logP else NULL,
      S = NA_real_,
      cov = NULL,
      cov.g = NULL,
      stats = cbind(theta = tau, g = g),
      loglik = loglik,
      statsFunction = NULL
    ))
  }

  stats_function <- function(a) {
    g <- softmax(as.vector(Q %*% a))
    G <- cumsum(g)
    if (family == "Binomial") {
      logg <- logsoftmax(as.vector(Q %*% a))
      logf <- apply(logP, 1, function(row) logsumexp(row + logg))
      f <- exp(logf)
      Pt <- exp(logP - logf)
      y_hat <- obs_weights
      W <- sweep(t(Pt), 1, g, "*") - g %o% rep(1.0, nrow(logP))
    } else {
      f <- as.vector(P %*% g)
      f <- pmax(f, 1e-300)
      yv <- y_vec_make(y, nrow(P))
      y_hat <- sum(obs_weights * yv) * f
      Pt <- sweep(P, 1, f, "/")
      W <- sweep(t(Pt), 1, g, "*") - g %o% rep(1.0, nrow(P))
    }

    qw <- crossprod(Q, W)
    ywq <- (sweep(t(W), 1, y_hat, "*")) %*% Q
    I1 <- qw %*% ywq

    s_dot <- penalty_grad(a)
    s_dot_dot <- penalty_hessian(a)

    Rval <- sum(diag(s_dot_dot)) / max(1e-300, sum(diag(I1)))

    A <- I1 + s_dot_dot
    R <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(R)) {
      A <- A + diag(1e-8, nrow(A))
      R <- chol(A)
    }
    I2 <- chol2inv(R)

    bias <- as.vector(-I2 %*% s_dot)
    Cov <- I2 %*% (I1 %*% t(I2))

    Dq <- (diag(g) - (g %o% g)) %*% Q
    bias_g <- as.vector(Dq %*% bias)
    Cov_g <- Dq %*% Cov %*% t(Dq)
    se_g <- sqrt(pmax(0, diag(Cov_g)))

    CovG <- apply(Cov_g, 2, cumsum)
    CovG <- t(apply(CovG, 1, cumsum))
    se_G <- sqrt(pmax(0, diag(CovG)))

    mat <- cbind(
      theta = tau,
      g = g,
      SE.g = se_g,
      G = G,
      SE.G = se_G,
      Bias.g = bias_g
    )
    list(S = Rval, cov = Cov, cov.g = Cov_g, mat = mat)
  }

  stats <- stats_function(mle)

  if (family == "Poisson" && ignoreZero) {
    mat <- stats$mat
    tg <- mat[, "g"] / (1 - exp(-tau))
    stats$mat <- cbind(mat, tg = tg / sum(tg))
  }

  list(
    mle = mle,
    Q = Q,
    P = if (family == "Binomial") NULL else P,
    logP = if (family == "Binomial") logP else NULL,
    S = stats$S,
    cov = stats$cov,
    cov.g = stats$cov.g,
    stats = stats$mat,
    loglik = loglik,
    statsFunction = stats_function
  )
}
