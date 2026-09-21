## Row-wise Poisson SuSiE used by the MiSo initialization.
## F is row-normalized: Gamma shapes/rates need no factor dimension.
## The public interface always returns N x D summaries, including N = 1.


## Extract a sparse observation without converting its implicit zeros to a
## dense vector. Dense input keeps the original reference calculation.
.miso_susie_observation <- function(y, M) {
  sparse = .miso_is_sparse(y) || inherits(y, "sparseVector")
  if (.miso_is_sparse(y)) {
    y = .miso_prepare_counts(y)
    .miso_stop(nrow(y) == 1 && ncol(y) == M,
               "A sparse observation must be a 1 by M matrix.")
    features = rep.int(seq_len(M), diff(y@p))
    values = y@x
  } else if (inherits(y, "sparseVector")) {
    .miso_stop(length(y) == M, "y must have length M = ncol(F).")
    y = methods::as(y, "dsparseVector")
    .miso_stop(all(is.finite(y@x)) && all(y@x >= 0),
               "y must contain finite nonnegative counts.")
    positive = y@x > 0
    features = y@i[positive]
    values = y@x[positive]
  } else {
    .miso_stop(is.numeric(y) && length(y) == M,
               "y must contain one numeric count per column of F.")
    features = seq_len(M)
    values = as.numeric(y)
  }
  .miso_stop(all(is.finite(values)) && all(values >= 0),
             "y must contain finite nonnegative counts.")
  list(values = values, features = features, sparse = sparse)
}

.poisson_susie_elbo <- function(y, F, gamma, alpha, beta,
                                  alpha0, beta0, xi, eps = 1e-12,
                                  log_F = NULL) {
  D = nrow(gamma)
  K = ncol(gamma)
  if (is.null(log_F)) log_F = log(pmax(F, eps))
  expected_log_F = gamma %*% log_F
  expected_log_lambda = digamma(alpha) - log(beta)
  expected_lambda = alpha / beta

  data_bound = 0
  for (d in seq_len(D)) {
    data_bound = data_bound +
      sum(xi[d, ] * y *
            (expected_log_lambda[d] + expected_log_F[d, ] -
               log(pmax(xi[d, ], eps)))) -
      expected_lambda[d]
  }

  data_bound -
    sum(.miso_gamma_kl(alpha, beta, alpha0, beta0)) -
    sum(gamma * (log(pmax(gamma, eps)) - log(1 / K)))
}

## The row-wise driver normalizes F once, avoiding O(N K M) repeated
## normalization when only a few features per observation are nonzero.
.poisson_susie_fit_normalized <- function(y, F, D, alpha0, beta0, max_iters,
                                          tol, update_prior, seed, eps, keep_xi) {
  if (!is.null(seed)) set.seed(seed)
  K = nrow(F)
  M = ncol(F)
  observation = .miso_susie_observation(y, M)
  y = observation$values
  ## Normalize the full dictionary BEFORE selecting observed features. Its
  ## total rate remains one, including the omitted zero-count coordinates.
  work_F = F[, observation$features, drop = FALSE]
  log_F = log(pmax(work_F, eps))
  .miso_stop(length(alpha0) == D && length(beta0) == D &&
               all(alpha0 > 0) && all(beta0 > 0),
             "alpha0 and beta0 must contain D positive values.")

  gamma = matrix(rexp(D * K), D, K)
  gamma = .miso_normalize_rows(gamma, eps)
  xi = matrix(1 / D, D, length(y))
  alpha = alpha0 + sum(y) / D
  beta = beta0 + 1
  elbo = rep(NA_real_, max_iters)

  for (iteration in seq_len(max_iters)) {
    expected_log_lambda = digamma(alpha) - log(beta)
    log_xi = gamma %*% log_F + expected_log_lambda
    xi = .miso_softmax_columns(log_xi)

    alpha = alpha0 + as.vector(xi %*% y)
    beta = beta0 + 1

    ## Keep the per-slot products: batching them changes rounding in tied
    ## factors and can alter the discrete slot alignment during initialization.
    for (d in seq_len(D)) {
      score = as.vector(log_F %*% (xi[d, ] * y))
      gamma[d, ] = exp(score - max(score))
      gamma[d, ] = gamma[d, ] / sum(gamma[d, ])
    }

    if (update_prior) {
      mean_lambda = alpha / beta
      mean_log_lambda = digamma(alpha) - log(beta)
      alpha0 = .miso_gamma_shape_from_moments(mean_lambda, mean_log_lambda)
      beta0 = alpha0 / pmax(mean_lambda, eps)
    }

    elbo[iteration] = .poisson_susie_elbo(
      y, work_F, gamma, alpha, beta, alpha0, beta0, xi, eps, log_F
    )
    if (iteration > 1) {
      relative_change = abs(elbo[iteration] - elbo[iteration - 1]) /
        (abs(elbo[iteration - 1]) + 1)
      if (relative_change < tol) break
    }
  }

  ## Synchronize q(lambda), xi, and gamma with the final prior parameters.
  expected_log_lambda = digamma(alpha) - log(beta)
  xi_gamma = gamma
  log_xi = gamma %*% log_F + expected_log_lambda
  xi = .miso_softmax_columns(log_xi)
  alpha = alpha0 + as.vector(xi %*% y)
  beta = beta0 + 1
  for (d in seq_len(D)) {
    score = as.vector(log_F %*% (xi[d, ] * y))
    gamma[d, ] = exp(score - max(score))
    gamma[d, ] = gamma[d, ] / sum(gamma[d, ])
  }
  final_elbo = .poisson_susie_elbo(
    y, work_F, gamma, alpha, beta, alpha0, beta0, xi, eps, log_F
  )

  ## Full xi is a diagnostic output only. Reconstruct it from the parameters
  ## that produced the last xi, not the subsequently updated alpha and gamma.
  if (!keep_xi) {
    xi = NULL
  } else if (observation$sparse) {
    xi = .miso_softmax_columns(xi_gamma %*% log(pmax(F, eps)) + expected_log_lambda)
  }

  list(
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    alpha0 = alpha0,
    beta0 = beta0,
    xi = xi,
    elbo = elbo[seq_len(iteration)],
    final_elbo = final_elbo,
    n_iter = iteration
  )
}

#' Fit Poisson sum-of-single-effects models to observations
#'
#' Observations are fitted independently against a common row-normalized
#' dictionary. A vector is treated as one observation. Results retain matrix
#' and array dimensions even for a single observation, slot, or factor.
#'
#' @param Y A finite nonnegative numeric vector or matrix, or a Matrix sparse
#'   vector/matrix. Matrix rows are observations and columns are features.
#'   Fractional values are allowed and are not rounded.
#' @param F Nonnegative factor-by-feature dictionary, normalized internally.
#' @param D Positive integer number of factor slots per observation.
#' @param alpha0,beta0 Positive length-D loading prior shapes and rates.
#' @param max_iters Maximum iterations per observation.
#' @param tol Relative ELBO-change tolerance per observation.
#' @param update_prior Update the loading prior parameters.
#' @param seed Random seed; row i uses seed + i - 1.
#' @param eps Numerical floor.
#' @param keep_fits Retain individual fits, including their ELBO histories and
#'   D by M auxiliary allocation arrays. Defaults to FALSE to save memory.
#' @return A list containing `gamma` (N by D by K), `alpha` and `beta`
#'   (N by D), the normalized dictionary `F` (K by M), and `fits` (NULL or a
#'   list of row-level fits). The list can be cached and passed to [miso_init()].
#' @export
#' @examples
#' F <- rbind(c(4, 1, 1), c(1, 1, 4))
#' ps <- poisson_susie_fit(c(5, 1, 2), F, D = 2, max_iters = 5)
#' dim(ps$gamma)
poisson_susie_fit <- function(Y, F, D, alpha0 = rep(1, D),
                               beta0 = rep(1, D), max_iters = 100,
                               tol = 1e-6, update_prior = TRUE,
                               seed = 1, eps = 1e-12, keep_fits = FALSE) {
  if (inherits(Y, "sparseVector")) {
    observation = .miso_susie_observation(Y, ncol(F))
    Y = Matrix::sparseMatrix(i = rep.int(1L, length(observation$features)),
                             j = observation$features, x = observation$values,
                             dims = c(1L, ncol(F)))
  } else if (is.numeric(Y) && is.null(dim(Y))) {
    Y = matrix(Y, nrow = 1L)
  }
  Y = .miso_prepare_counts(Y)
  .miso_validate_dictionary(F, ncol(Y))
  .miso_positive_integer(D, "D")
  .miso_positive_integer(max_iters, "max_iters")
  .miso_stop(length(seed) == 1 && is.finite(seed), "seed must be finite.")
  F = .miso_normalize_rows(F, eps)
  sparse = .miso_is_sparse(Y)
  if (sparse) row_Y = methods::as(Y, "RsparseMatrix")
  N = nrow(Y)
  K = nrow(F)
  gamma = array(0, c(N, D, K))
  alpha = beta = matrix(0, N, D)
  fits = if (keep_fits) vector("list", N) else NULL

  for (i in seq_len(N)) {
    if (sparse) {
      first = row_Y@p[i] + 1L
      last = row_Y@p[i + 1L]
      index = if (first <= last) seq.int(first, last) else integer()
      y = Matrix::sparseVector(x = row_Y@x[index], i = row_Y@j[index] + 1L,
                               length = ncol(Y))
    } else {
      y = Y[i, ]
    }
    fit = .poisson_susie_fit_normalized(
      y, F, D, alpha0, beta0, max_iters, tol, update_prior,
      seed = seed + i - 1, eps = eps, keep_xi = keep_fits
    )
    gamma[i, , ] = fit$gamma
    alpha[i, ] = fit$alpha
    beta[i, ] = fit$beta
    if (keep_fits) fits[[i]] = fit
  }

  list(gamma = gamma, alpha = alpha, beta = beta, fits = fits, F = F)
}
