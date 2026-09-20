#' Construct a reusable MiSo initialization
#'
#' Estimate a dictionary if needed, fit Poisson SuSiE independently to each
#' observation, discover motif supports, and initialize the mixture parameters.
#' Supply cached Poisson-SuSiE results to change support cutoffs without refitting
#' the observation-level models.
#'
#' @param Y A nonnegative numeric count matrix (observations by features), dense
#'   or a sparse matrix from the Matrix package.
#' @param F A nonnegative dictionary, factors by features. Rows are normalized
#'   internally. If omitted, use the cached dictionary or estimate one by NMF.
#' @param D Positive integer number of factor slots per motif. May be omitted
#'   when `poisson_susie` is supplied.
#' @param K Number of dictionary factors. Required when estimating `F`.
#' @param poisson_susie Optional result of [poisson_susie_fit()] for the same
#'   observations, in the same row order, and the same dictionary. Counts are
#'   not retained in that result, so the caller must ensure row correspondence.
#' @param phi0 Positive scalar or vector of Dirichlet prior parameters. A vector
#'   must have one entry per discovered motif.
#' @param tau Cumulative expected-loading fraction defining each initial support,
#'   in `(0, 1]`.
#' @param min_fraction Minimum initial motif fraction, in `(0, 1]`.
#' @param gamma_floor Initial factor-selection probability smoothing, in `[0, 1)`.
#'   This controls initialization, not update damping.
#' @param susie_max_iters,susie_tol Poisson-SuSiE iteration limit and tolerance;
#'   used only when cached results are absent.
#' @param nmf_max_iters Iteration limit for dictionary initialization.
#' @param seed Random seed for initialization.
#' @param eps Positive numerical floor.
#' @param keep_intermediates Retain the row-wise posterior summaries as
#'   `poisson_susie` in the returned list. Auxiliary allocation arrays are omitted.
#'
#' @return A list accepted by `miso_fit(init = ...)`: `F` (K by M), `gamma`
#'   (S by D by K), `alpha0` and `beta0` (S by D), `omega` (N by S), `phi0`,
#'   and `phi`. Also includes `S`, `cluster`, `pattern`, `support`, `anchors`,
#'   `expected_loading`, `slot_permutation`, `settings`, and `call`.
#' @export
#' @examples
#' set.seed(1)
#' Y <- matrix(rpois(60, 2), 10, 6)
#' F <- rbind(c(4, 4, 1, 1, 1, 1), c(1, 1, 1, 1, 4, 4))
#' ps <- poisson_susie_fit(Y, F, D = 2, max_iters = 5)
#' initial <- miso_init(Y, poisson_susie = ps, tau = 0.9)
#' fit <- miso_fit(Y, init = initial, max_iters = 5)
miso_init <- function(Y, F = NULL, D = NULL, K = NULL,
                      poisson_susie = NULL, phi0 = 0.1, tau = 0.9,
                      min_fraction = 0.05, gamma_floor = 0.05,
                      susie_max_iters = 100, susie_tol = 1e-6,
                      nmf_max_iters = 200, seed = 1, eps = 1e-12,
                      keep_intermediates = FALSE) {
  Y = .miso_prepare_counts(Y)
  .miso_stop(length(eps) == 1 && is.finite(eps) && eps > 0,
             "eps must be positive and finite.")
  cached = !is.null(poisson_susie)
  if (cached) {
    .miso_stop(is.list(poisson_susie) && length(dim(poisson_susie$gamma)) == 3,
               "poisson_susie must be a Poisson-SuSiE result.")
    if (is.null(D)) D = dim(poisson_susie$gamma)[2]
    if (is.null(F)) F = poisson_susie$F
    .miso_stop(!is.null(F), "Supply the dictionary used for poisson_susie.")
  }
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
  if (cached) {
    .miso_stop(identical(dim(poisson_susie$gamma), c(nrow(Y), as.integer(D), nrow(F))),
               "Cached Poisson-SuSiE dimensions must agree with Y, D, and F.")
    if (!is.null(poisson_susie$F)) {
      .miso_stop(isTRUE(all.equal(unname(F), unname(poisson_susie$F), tolerance = 1e-10)),
                 "F must match the dictionary used for poisson_susie.")
    }
  } else {
    poisson_susie = poisson_susie_fit(Y, F, D, max_iters = susie_max_iters,
                                     tol = susie_tol, seed = seed, eps = eps)
  }
  initial = .miso_initialize_from_susie(
    poisson_susie, D, tau, min_fraction, gamma_floor, phi0, eps
  )
  initial$F = F
  initial$phi0 = rep(phi0, length.out = initial$S)
  initial$settings = list(D = D, K = nrow(F), tau = tau,
    min_fraction = min_fraction, gamma_floor = gamma_floor,
    susie_max_iters = susie_max_iters, susie_tol = susie_tol,
    nmf_max_iters = nmf_max_iters, seed = seed, reused_poisson_susie = cached)
  initial$call = match.call()
  if (keep_intermediates) initial$poisson_susie = poisson_susie
  initial
}
