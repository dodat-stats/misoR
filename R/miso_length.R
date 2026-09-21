## Document-length-adjusted MiSo: Y_ij ~ Poi(n_i (LF)_ij).

#' Fit MiSo with document-length-adjusted Poisson means
#'
#' Fit the Poisson working likelihood `Y[i,j] ~ Poisson(n[i] * (L %*% F)[i,j])`
#' with the same sparse motif structure as [miso_fit()]. Counts remain on their
#' original scale. Dictionary rows sum to one; loadings describe contributions
#' per unit of document length and are not constrained to sum to one. Both
#' fitters use the same SuSiE, factor, and mixture steps, with one update of
#' factor selections and the dictionary per outer iteration. Nonnegative
#' fractional inputs are accepted and optimize the generalized Poisson objective.
#'
#' @inheritParams miso_fit
#' @param D Positive integer number of factor slots per motif. Inferred from
#'   `init$gamma` when an initialization is supplied.
#' @param n Known positive document lengths, one per row of `Y`. Defaults to
#'   row totals over the supplied vocabulary. With this default, remove empty
#'   documents before fitting. Other positive lengths may be supplied explicitly;
#'   `rep(1, nrow(Y))` recovers the original likelihood and loading scale.
#' @param init Optional initialization with `gamma`, `alpha0`, `beta0`, and
#'   optionally `omega`, `F`, and `phi0`, as in [miso_fit()]. Loading priors must
#'   be on the per-unit-n scale. Reuse `fit$initialization` from a length-adjusted
#'   fit with `keep_initialization = TRUE`; its stored `n` must match. Do not
#'   directly reuse the count-scale priors returned by [miso_init()].
#'
#' @details
#' For integer counts, when `n` is the observed row total, this is a working likelihood, not the
#' exact sampling law conditional on that total. Write `r[i] = sum(L[i,])`
#' and `p[i,] = (L %*% F)[i,] / r[i]`. The likelihood equals a multinomial
#' likelihood for `p[i,]` times `Poisson(n[i]; n[i] * r[i])`. It therefore
#' encourages, but does not enforce, a total loading of one. Loading priors
#' mean that the Bayesian fit is not generally identical to a multinomial fit.
#'
#' Automatic initialization estimates the usual Poisson NMF dictionary if
#' necessary, fits count-scale Poisson SuSiE, and rescales each row's Gamma
#' posterior from loading `a` to `a / n[i]` by multiplying its rate by `n[i]`.
#' Motif discovery and pooled prior estimation then use these rescaled moments.
#' This reuses count-scale SuSiE for initialization only; subsequent MiSo
#' iterations optimize the length-adjusted objective.
#'
#' Each loading posterior has rate `beta0[s,d] + n[i]` at its coordinate update.
#' Subsequent empirical-Bayes updates hold that posterior fixed. Consequently,
#' returned `beta` need not equal returned `beta0 + n`, since the latter prior
#' may have been updated later in the same iteration.
#'
#' @return A `miso_fit` object with the fields described in [miso_fit()], except
#'   `beta` is N by S by D. Also includes `n` and `elbo_constant`, the omitted
#'   data-only term `sum(rowSums(Y) * log(n)) - sum(lgamma(Y + 1))`. Add this
#'   constant to `elbo` or `final_elbo` for the full working-likelihood bound.
#'   `predict(fit, "loadings")` returns per-unit-n loadings and `predict(fit)`
#'   returns mean counts, multiplying the reconstruction by `n` row-wise.
#'   Normalizing posterior mean loadings gives a descriptive composition, not
#'   the exact posterior mean of the normalized latent loadings.
#' @seealso [miso_fit()], [predict.miso_fit()], [plot.miso_fit()]
#' @export
#' @examples
#' Y <- matrix(c(12, 3, 2, 1, 3, 15, 2, 4, 6, 6, 5, 3), 3, 4, byrow = TRUE)
#' F <- rbind(c(5, 1, 1, 1), c(1, 5, 1, 1))
#' fit <- miso_fit_length(Y, F, D = 2, max_iters = 5,
#'                        init_control = list(susie_max_iters = 5))
#' L <- predict(fit, "loadings")
#' theta <- L / rowSums(L)
#' fitted_counts <- predict(fit)
miso_fit_length <- function(Y, F = NULL, D = NULL, init = NULL, K = NULL,
                            n = NULL, phi0 = 0.1, max_iters = 100,
                            mixture_max_iters = 30, tol = 1e-5,
                            min_iters = 5, patience = 2,
                            update_prior = TRUE, update_gamma = TRUE,
                            update_F = TRUE, block_size = 100, eps = 1e-12,
                            verbose = FALSE, seed = 1, init_control = list(),
                            keep_initialization = FALSE) {
  Y = .miso_prepare_counts(Y)
  if (is.null(n)) n = .miso_row_sums(Y)
  .miso_stop(is.numeric(n) && is.null(dim(n)) && length(n) == nrow(Y) &&
               all(is.finite(n)) && all(n > 0),
             paste("n must contain one positive finite document length per row;",
                   "remove empty documents when using row totals."))
  n = as.numeric(n)
  fit = .miso_fit_entry(
    Y, F, D, init, K, phi0, max_iters, mixture_max_iters, tol,
    min_iters, patience, update_prior, update_gamma, update_F, block_size,
    eps, verbose, seed, init_control, keep_initialization,
    phi0_missing = missing(phi0), n = n
  )
  fit$call = match.call()
  fit
}

## Transform row-level posteriors before pooling across documents, so that
## motif discovery and empirical-Bayes initialization use comparable rates.
.miso_init_length <- function(Y, n, F = NULL, D = NULL, K = NULL,
                              phi0 = 0.1, tau = 0.9, min_fraction = 0.05,
                              gamma_floor = 0.05, susie_max_iters = 100,
                              susie_tol = 1e-6, nmf_max_iters = 200,
                              seed = 1, eps = 1e-12) {
  .miso_stop(length(eps) == 1 && is.finite(eps) && eps > 0,
             "eps must be positive and finite.")
  .miso_positive_integer(D, "D")
  if (is.null(F)) {
    .miso_positive_integer(K, "K (required when F is NULL)")
    .miso_positive_integer(nmf_max_iters, "nmf_max_iters")
    F = .poisson_nmf_fit(Y, K, max_iters = nmf_max_iters, seed = seed)$F
  }
  .miso_validate_dictionary(F, ncol(Y))
  if (!is.null(K)) {
    .miso_positive_integer(K, "K")
    .miso_stop(nrow(F) == K, "K must agree with nrow(F).")
  }
  F = .miso_normalize_rows(F, eps)
  ps = poisson_susie_fit(Y, F, D, max_iters = susie_max_iters,
                         tol = susie_tol, seed = seed, eps = eps)
  ps$beta = sweep(ps$beta, 1, n, "*")
  initial = .miso_initialize_from_susie(
    ps, D, tau, min_fraction, gamma_floor, phi0, eps
  )
  initial$F = F
  initial$phi0 = rep(phi0, length.out = initial$S)
  initial$n = n
  initial$settings = list(D = D, K = nrow(F), tau = tau,
    min_fraction = min_fraction, gamma_floor = gamma_floor,
    susie_max_iters = susie_max_iters, susie_tol = susie_tol,
    nmf_max_iters = nmf_max_iters, seed = seed, loading_scale = "per_unit_n")
  initial
}
