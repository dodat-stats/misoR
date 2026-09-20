## MiSo: variational empirical Bayes with a Dirichlet prior on pi.
##
## The implementation follows the row-normalized formulation in the paper.
## Local xi variables are computed in feature blocks and are never stored.


## Fit from a supplied dictionary or estimate one with optional Poisson NMF.
## init is a list containing gamma, alpha0, beta0, and optionally omega/F.
## init_control controls motif discovery and optional dictionary estimation.
#' Fit a mixture of sparse Poisson factor models
#'
#' Each outer iteration performs one Jacobi SuSiE update, alternating factor
#' selection and dictionary updates using the same sufficient statistics, then
#' mixture refinement and ELBO evaluation. Auxiliary allocations are computed
#' in feature blocks and are never retained. Updates are undamped.
#'
#' @inheritParams miso_init
#' @param init Optional initialization list from [miso_init()], or a list with
#'   `gamma` (S by D by K), `alpha0` and `beta0` (S by D), and optionally `omega`
#'   (N by S), `F`, and `phi0`. If `phi0` is omitted in the call, use `init$phi0`
#'   when available. D is inferred from supplied parameters when omitted.
#' @param max_iters Maximum outer iterations.
#' @param mf_iters Number of alternating gamma/F updates per outer iteration.
#'   Only one update is needed when either parameter is fixed.
#' @param mixture_max_iters Maximum responsibility/Dirichlet updates per outer
#'   iteration; these stop at relative Dirichlet-parameter change below 1e-10.
#' @param tol Relative ELBO-change tolerance. Zero disables early stopping.
#' @param min_iters Minimum outer iterations before convergence can be counted.
#' @param patience Number of consecutive small ELBO changes required to stop.
#' @param update_prior,update_gamma,update_F Whether to update loading priors,
#'   factor-selection probabilities, and the dictionary, respectively.
#' @param block_size Number of features in each allocation block.
#' @param verbose Print iteration diagnostics.
#' @param init_control Named list of initialization controls: `tau`,
#'   `min_fraction`, `gamma_floor`, `susie_max_iters`, `susie_tol`, and
#'   `nmf_max_iters`. The last applies only when F is absent. Must be empty
#'   when supplying `init`.
#' @param keep_initialization Retain starting parameters as `initialization`.
#'
#' @return A `miso_fit` list containing `F` (K by M), `gamma` (S by D by K),
#'   `alpha` (N by S by D), `beta`, `alpha0`, `beta0`, and `prior_mean`
#'   (each S by D), `omega` (N by S), `phi0`, `phi`, `pi`, and `allocation_mass`
#'   (each length S), and `z_hat` (length N). Also contains `component_elbo`
#'   (N by S), `elbo`, `final_elbo`, `converged`, `n_iter`, and `call`.
#'   The returned state is the last completed iteration, with
#'   `final_elbo == tail(elbo, 1)`. F rows sum to one, so expected factor
#'   loadings are on the count scale. `pi` is the posterior mean mixture weight;
#'   `allocation_mass` is the mean responsibility across observations.
#' @seealso [miso_init()], [summary.miso_fit()], [plot.miso_fit()],
#'   [predict.miso_fit()]
#' @export
#' @examples
#' set.seed(1)
#' Y <- matrix(rpois(60, 2), 10, 6)
#' F <- rbind(c(4, 4, 1, 1, 1, 1), c(1, 1, 1, 1, 4, 4))
#' fit <- miso_fit(Y, F, D = 2, max_iters = 5,
#'                 init_control = list(susie_max_iters = 5))
#' print(fit)
#' L <- predict(fit, type = "loadings")
#' stopifnot(all(dim(L) == c(10, 2)))
miso_fit <- function(Y, F = NULL, D = NULL, init = NULL, K = NULL,
                     phi0 = 0.1, max_iters = 100, mf_iters = 3,
                     mixture_max_iters = 30, tol = 1e-5,
                     min_iters = 5, patience = 2,
                     update_prior = TRUE, update_gamma = TRUE, update_F = TRUE,
                     block_size = 100, eps = 1e-12, verbose = FALSE,
                     seed = 1, init_control = list(), keep_initialization = FALSE) {
  Y = .miso_prepare_counts(Y)
  .miso_positive_integer(max_iters, "max_iters")
  .miso_positive_integer(mixture_max_iters, "mixture_max_iters")
  .miso_positive_integer(mf_iters, "mf_iters")
  .miso_stop(is.list(init_control) &&
               (length(init_control) == 0 ||
                  (!is.null(names(init_control)) &&
                   all(nzchar(names(init_control))) && !anyDuplicated(names(init_control)))),
             "init_control must be a named list without duplicate names.")
  allowed = c("tau", "min_fraction", "gamma_floor", "susie_max_iters", "susie_tol",
               "nmf_max_iters")
  .miso_stop(all(names(init_control) %in% allowed),
             paste("init_control accepts:", paste(allowed, collapse = ", ")))
  if (!is.null(init)) {
    .miso_stop(is.list(init) && all(c("gamma", "alpha0", "beta0") %in% names(init)),
               "init must contain gamma, alpha0, and beta0.")
    if (missing(phi0) && !is.null(init$phi0)) phi0 = init$phi0
    if (is.null(F)) F = init$F
    .miso_stop(!is.null(F), "Supply F or include F in init.")
    .miso_stop(length(init_control) == 0,
               "init_control is only used when init is NULL.")
    if (!is.null(D)) {
      .miso_positive_integer(D, "D")
      .miso_stop(length(dim(init$gamma)) == 3 && dim(init$gamma)[2] == D,
                 "D must agree with init$gamma.")
    }
  } else {
    .miso_stop(is.null(F) || is.null(init_control$nmf_max_iters),
               "nmf_max_iters only applies when F is NULL.")
    init = do.call(miso_init,
                   c(list(Y = Y, F = F, D = D, K = K, phi0 = phi0,
                          seed = seed, eps = eps), init_control))
    F = init$F
  }
  .miso_validate_dictionary(F, ncol(Y))
  if (!is.null(K)) {
    .miso_positive_integer(K, "K")
    .miso_stop(nrow(F) == K, "K must agree with nrow(F).")
  }
  fit = .miso_fit(Y, F, init$gamma, init$alpha0, init$beta0, init$omega,
                  phi0, max_iters, mf_iters, mixture_max_iters, tol,
                  min_iters, patience, update_prior, update_gamma, update_F,
                  block_size, eps, verbose)
  fit$call = match.call()
  if (keep_initialization) fit$initialization = init
  fit
}

.miso_update_dirichlet <- function(component_elbo, phi0, phi = NULL,
                                  max_iters = 30, tol = 1e-10) {
  N = nrow(component_elbo)
  S = ncol(component_elbo)
  if (!length(phi0) %in% c(1, S)) {
    stop("phi0 must contain one value or one per motif.", call. = FALSE)
  }
  phi0 = rep(phi0, length.out = S)
  if (is.null(phi)) phi = phi0 + N / S

  for (iteration in seq_len(max_iters)) {
    expected_log_pi = .miso_dirichlet_expected_log(phi)
    omega = .miso_softmax_rows(
      sweep(component_elbo, 2, expected_log_pi, "+")
    )
    phi_new = phi0 + colSums(omega)
    change = max(abs(phi_new - phi) / pmax(phi, 1))
    phi = phi_new
    if (change < tol) break
  }

  list(
    omega = omega,
    phi = phi,
    expected_log_pi = .miso_dirichlet_expected_log(phi),
    posterior_mean_pi = phi / sum(phi),
    allocation_mass = colMeans(omega)
  )
}

.miso_update_gamma <- function(F, C, gamma, eps = 1e-12) {
  S = dim(gamma)[1]
  D = dim(gamma)[2]
  log_F = log(pmax(F, eps))
  score = matrix(C, S * D, ncol(F)) %*% t(log_F)
  .miso_gamma_from_scores(score, gamma, eps)
}

## Full coordinate update, shared by dense counts and streamed scores.
.miso_gamma_from_scores <- function(score, gamma, eps = 1e-12) {
  S = dim(gamma)[1]
  D = dim(gamma)[2]
  K = dim(gamma)[3]
  score = matrix(score, S * D, K)
  gamma[] = .miso_softmax_rows(score)
  .miso_normalize_gamma(gamma, eps)
}

.miso_update_priors <- function(alpha, beta, omega, alpha0, beta0,
                               eps = 1e-12) {
  S = ncol(omega)
  D = ncol(alpha0)

  for (s in seq_len(S)) {
    weight = omega[, s]
    total_weight = sum(weight)
    if (total_weight <= eps) next

    for (d in seq_len(D)) {
      mean_lambda = sum(weight * alpha[, s, d] / beta[s, d]) / total_weight
      mean_log_lambda = sum(
        weight * (digamma(alpha[, s, d]) - log(beta[s, d]))
      ) / total_weight
      shape = .miso_gamma_shape_from_moments(mean_lambda, mean_log_lambda)
      alpha0[s, d] = shape
      beta0[s, d] = shape / pmax(mean_lambda, eps)
    }
  }
  list(alpha0 = alpha0, beta0 = beta0)
}

## Keep an unused factor's dictionary row unchanged.
.miso_update_F <- function(F, gamma, C, eps = 1e-12) {
  K = nrow(F)
  M = ncol(F)
  S = dim(gamma)[1]
  D = dim(gamma)[2]
  T = crossprod(matrix(gamma, S * D, K), matrix(C, S * D, M))
  active = rowSums(T) > 0
  F[active, ] = .miso_normalize_rows(T[active, , drop = FALSE], eps)
  F
}

.miso_elbo <- function(component_elbo, omega, phi, phi0, gamma,
                      eps = 1e-12) {
  S = ncol(omega)
  K = dim(gamma)[3]
  if (!length(phi0) %in% c(1, S)) {
    stop("phi0 must contain one value or one per motif.", call. = FALSE)
  }
  phi0 = rep(phi0, length.out = S)
  expected_log_pi = .miso_dirichlet_expected_log(phi)

  sum(omega * component_elbo) +
    sum(sweep(omega, 2, expected_log_pi, "*")) -
    sum(omega * log(pmax(omega, eps))) -
    .miso_dirichlet_kl(phi, phi0) -
    sum(gamma * (log(pmax(gamma, eps)) - log(1 / K)))
}

## SuSiE step: one Jacobi pass for both loading counts and MF statistics.
## beta is the stored posterior rate, not a rate reset from the latest prior.
.miso_susie_step <- function(prepared, F, gamma, alpha, beta, alpha0, beta0,
                             omega, update_gamma = TRUE, update_F = TRUE,
                             eps = 1e-12) {
  N = dim(alpha)[1]
  S = dim(gamma)[1]
  D = dim(gamma)[2]
  K = dim(gamma)[3]
  log_F = log(pmax(F, eps))
  C = if (update_F) array(0, c(S, D, ncol(F))) else NULL
  stream_scores = update_gamma && !update_F
  scores = if (stream_scores) array(0, c(S, D, K)) else NULL

  for (s in seq_len(S)) {
    pass = prepared$pass(
      prepared$blocks,
      digamma(matrix(alpha[, s, ], N, D)) -
        matrix(log(beta[s, ]), N, D, byrow = TRUE),
      matrix(gamma[s, , ], D, K) %*% log_F,
      omega = if (update_gamma || update_F) omega[, s] else NULL,
      log_F = if (stream_scores) log_F else NULL, eps = eps
    )
    ## Other motifs do not depend on this motif's alpha. Its entire feature
    ## pass is complete before replacing any of its shapes.
    alpha[, s, ] = sweep(pass$allocated, 2, alpha0[s, ], "+")
    if (update_F) C[s, , ] = pass$C
    if (stream_scores) scores[s, , ] = pass$score
  }
  list(alpha = alpha, beta = beta0 + 1, C = C, scores = scores)
}

## Mixture step's second pass: bounds at the updated parameters, no counts.
.miso_component_elbo <- function(prepared, F, gamma, alpha, beta, alpha0, beta0,
                                  eps = 1e-12) {
  N = dim(alpha)[1]
  S = dim(gamma)[1]
  D = dim(gamma)[2]
  K = dim(gamma)[3]
  log_F = log(pmax(F, eps))
  component_elbo = matrix(0, N, S)
  for (s in seq_len(S)) {
    alpha_s = matrix(alpha[, s, ], N, D)
    pass = prepared$pass(
      prepared$blocks,
      digamma(alpha_s) - matrix(log(beta[s, ]), N, D, byrow = TRUE),
      matrix(gamma[s, , ], D, K) %*% log_F,
      evaluate_bound = TRUE, compute_allocated = FALSE, eps = eps
    )
    component_elbo[, s] = pass$bound
    for (d in seq_len(D)) {
      component_elbo[, s] = component_elbo[, s] -
        alpha_s[, d] / beta[s, d] -
        .miso_gamma_kl(alpha_s[, d], beta[s, d], alpha0[s, d], beta0[s, d])
    }
  }
  component_elbo
}

.miso_fit <- function(Y, F, gamma_init, alpha0_init, beta0_init,
                     omega_init = NULL, phi0 = 0.1,
                     max_iters = 100, mf_iters = 3,
                     mixture_max_iters = 30, tol = 1e-5,
                     min_iters = 5, patience = 2,
                     update_prior = TRUE, update_gamma = TRUE, update_F = TRUE,
                     block_size = 100, eps = 1e-12,
                     verbose = FALSE) {
  Y = .miso_prepare_counts(Y)
  N = nrow(Y)
  S = dim(gamma_init)[1]
  D = dim(gamma_init)[2]
  K = nrow(F)
  if (is.null(omega_init)) omega_init = matrix(1 / S, N, S)
  .miso_validate_inputs(
    Y, F, gamma_init, alpha0_init, beta0_init, omega_init, phi0
  )

  F = .miso_normalize_rows(F, eps)
  gamma = .miso_normalize_gamma(gamma_init, eps)
  alpha0 = alpha0_init
  beta0 = beta0_init
  omega = .miso_normalize_rows(omega_init, eps)
  phi0 = rep(phi0, length.out = S)
  phi = phi0 + colSums(omega)
  D = dim(gamma)[2]
  alpha = array(0, c(N, S, D))
  totals = .miso_row_sums(Y)
  for (s in seq_len(S)) for (d in seq_len(D)) {
    alpha[, s, d] = alpha0[s, d] + totals / D
  }
  beta = beta0 + 1
  elbo = rep(NA_real_, max_iters)
  prepared = .miso_prepare_pass(Y, block_size)
  stable_iterations = 0
  converged = FALSE

  for (iteration in seq_len(max_iters)) {
    ## SuSiE: update alpha/beta and priors, using the current responsibilities.
    susie = .miso_susie_step(
      prepared, F, gamma, alpha, beta, alpha0, beta0, omega,
      update_gamma, update_F, eps
    )
    alpha = susie$alpha
    beta = susie$beta
    if (update_prior) {
      prior = .miso_update_priors(alpha, beta, omega, alpha0, beta0, eps)
      alpha0 = prior$alpha0
      beta0 = prior$beta0
    }

    ## MF: alternate full gamma/F updates with the same C. If either is fixed,
    ## a second repetition would be identical, so perform just one.
    for (mf in seq_len(if (update_gamma && update_F) mf_iters else 1L)) {
      if (update_gamma) {
        gamma = if (update_F) {
          .miso_update_gamma(F, susie$C, gamma, eps)
        } else {
          .miso_gamma_from_scores(susie$scores, gamma, eps)
        }
      }
      if (update_F) F = .miso_update_F(F, gamma, susie$C, eps)
    }
    rm(susie)

    ## Mixture: refresh bounds, update omega/phi, then record the full ELBO.
    component_elbo = .miso_component_elbo(
      prepared, F, gamma, alpha, beta, alpha0, beta0, eps
    )
    mixture = .miso_update_dirichlet(
      component_elbo, phi0, phi, max_iters = mixture_max_iters
    )
    omega = mixture$omega
    phi = mixture$phi
    elbo[iteration] = .miso_elbo(component_elbo, omega, phi, phi0, gamma, eps)

    if (verbose) {
      message(sprintf(
        "iteration %d: ELBO %.6f; allocation mass %s",
        iteration, elbo[iteration],
        paste(format(round(colMeans(omega), 4), nsmall = 4), collapse = ", ")
      ))
    }

    if (iteration > 1) {
      improvement = (elbo[iteration] - elbo[iteration - 1]) /
        (abs(elbo[iteration - 1]) + 1)
      stable_iterations = if (iteration >= min_iters && abs(improvement) < tol) {
        stable_iterations + 1
      } else {
        0
      }
      if (stable_iterations >= patience) {
        converged = TRUE
        break
      }
    }
  }

  n_iter = iteration
  elbo = elbo[seq_len(n_iter)]

  fit = list(
    omega = omega,
    z_hat = max.col(omega),
    phi0 = phi0,
    phi = phi,
    pi = phi / sum(phi),
    allocation_mass = colMeans(omega),
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    alpha0 = alpha0,
    beta0 = beta0,
    prior_mean = alpha0 / beta0,
    F = F,
    component_elbo = component_elbo,
    elbo = elbo,
    final_elbo = elbo[n_iter],
    converged = converged,
    n_iter = n_iter,
    call = match.call()
  )
  class(fit) = "miso_fit"
  fit
}

.miso_expected_loadings <- function(fit) {
  N = nrow(fit$omega)
  S = ncol(fit$omega)
  D = dim(fit$gamma)[2]
  K = dim(fit$gamma)[3]
  loading = matrix(0, N, K)

  for (s in seq_len(S)) {
    gamma_s = matrix(fit$gamma[s, , ], D, K)
    expected_lambda = matrix(fit$alpha[, s, ], N, D) /
      matrix(fit$beta[s, ], N, D, byrow = TRUE)
    loading = loading + fit$omega[, s] * (expected_lambda %*% gamma_s)
  }
  loading
}

#' Extract fitted mean counts or factor loadings
#'
#' Prediction currently describes the observations used for fitting.
#' The posterior mean loading is averaged over motif responsibilities and
#' factor selections. Means are L times F. Both outputs are dense matrices,
#' even when the input counts were sparse; no allocation array is constructed.
#'
#' @param object A fitted `miso_fit` object.
#' @param type `"mean"` for the N by M reconstructed rate matrix, or
#'   `"loadings"` for the N by K posterior mean factor loadings.
#' @param ... Reserved; supplying additional arguments, including `newdata`,
#'   raises an error. New-observation inference is not implemented.
#' @return A dense numeric matrix.
#' @export
predict.miso_fit <- function(object, type = c("mean", "loadings"), ...) {
  .miso_stop(length(list(...)) == 0,
             "Unused arguments: prediction is for fitted observations only.")
  type = match.arg(type)
  loading = .miso_expected_loadings(object)
  if (type == "loadings") loading else loading %*% object$F
}

#' Print a fitted MiSo model
#'
#' @param x A fitted `miso_fit` object.
#' @param ... Additional arguments; currently unused.
#' @return The fitted object invisibly.
#' @export
print.miso_fit <- function(x, ...) {
  cat("Dirichlet-pi MiSo fit\n")
  cat("  observations:", nrow(x$omega), "\n")
  cat("  motifs:", ncol(x$omega), "\n")
  cat("  iterations:", x$n_iter, "\n")
  cat("  converged:", x$converged, "\n")
  cat("  final ELBO:", format(x$final_elbo, digits = 8), "\n")
  cat("  allocation mass:",
      paste(format(round(x$allocation_mass, 4), nsmall = 4), collapse = ", "),
      "\n")
  invisible(x)
}
