#' Construct a reusable MiSo initialization
#'
#' By default, initialize directly from Poisson NMF loadings: truncate each
#' loading profile, group identical supports, and optionally merge or split
#' groups to obtain S motifs. The earlier Poisson-SuSiE initializer remains
#' available with `method = "poisson_susie"` or cached `poisson_susie` results.
#'
#' @param Y A finite nonnegative numeric matrix (observations by features), dense
#'   or a sparse matrix from the Matrix package. Fractional values are allowed.
#' @param F A nonnegative dictionary, K by M. Rows are normalized internally;
#'   the NMF method requires each row to have positive total mass.
#'   If omitted, estimate it by Poisson NMF (or use the cached SuSiE dictionary).
#'   With method "nmf", supplying F alone fits loadings while holding F fixed.
#' @param D Positive integer number of factor slots per motif. May be omitted
#'   when `poisson_susie` is supplied.
#' @param K Number of factors. Required when estimating F.
#' @param poisson_susie Optional cached result of [poisson_susie_fit()] on the
#'   same observations and dictionary. Selects the legacy method when `method`
#'   is omitted. With n supplied, cached posteriors must be on the count scale.
#' @param phi0 Positive scalar or vector of Dirichlet prior parameters. A vector
#'   must have one entry per final motif.
#' @param tau Cumulative loading fraction defining each support, in `(0, 1]`.
#' @param min_fraction Minimum motif fraction for the legacy SuSiE method only.
#'   The direct NMF method does not consolidate rare patterns automatically.
#' @param gamma_floor Initial factor-selection probability smoothing, in `[0, 1)`.
#' @param susie_max_iters,susie_tol Iteration limit and tolerance for the legacy
#'   SuSiE method only, when cached results are absent.
#' @param nmf_max_iters Iteration limit for fitting NMF or fixed-F loadings.
#' @param seed Random seed for NMF and reproducible random group splitting.
#' @param eps Positive numerical floor.
#' @param keep_intermediates Retain row-wise posteriors for the legacy method.
#' @param S Target number of initial motifs for the NMF method, or NULL to keep
#'   all S0 distinct supports. Must be between 1 and nrow(Y). Smaller targets
#'   merge groups by largest cosine similarity of their mean loading vectors;
#'   larger targets randomly bisect the largest group until S is reached.
#' @param L Optional N by K NMF loadings, supplied together with F, on the count
#'   scale so that L times F approximates Y. Normalizing F rescales L to preserve
#'   their product. Supplying both skips NMF fitting. NMF method only.
#' @param method Initialization method: "nmf" (default) or "poisson_susie".
#' @param n Optional positive exposures, one per observation. NMF loadings or
#'   SuSiE posterior moments are divided by n before grouping and prior fitting.
#'   The returned initialization is then for [miso_fit_length()], with matching n.
#' @param unsupported_fraction For NMF initialization, combined prior mean of
#'   unsupported slots as a fraction of the group's mean total loading (default
#'   0.01). Each such slot has shape 1 and an equal share, floored at eps.
#'
#' @details
#' Supports use the shortest decreasing-loading prefix reaching tau, capped at
#' D. Ties use factor index. A zero-loading row has empty support. On merging,
#' centers are weighted by group size; the D largest center coordinates in the
#' union of original observation supports determine the new support. Full means
#' determine cosine similarity. Two zero centers have similarity 1; a zero and
#' a nonzero center have similarity 0. Pair ties are deterministic.
#'
#' Splits randomly partition the largest group into floor(size/2) and
#' ceiling(size/2) observations. Groups with identical supports remain distinct.
#' Supported slots use Gamma maximum-likelihood moment equations on loadings
#' floored at eps, with shapes bounded between 0.001 and 10000, including singleton
#' and constant groups. Gamma parameters are shape and rate. Unsupported slots
#' have uniform factor-selection probabilities and small positive prior means.
#' This is initialization only, not ARD thresholding.
#'
#' @return A list accepted by `miso_fit(init = ...)`, or `miso_fit_length` when n
#'   is supplied: F, gamma (S by D by K), a (N by S by D), b, omega
#'   (N by S), phi0, phi, S, cluster, pattern, support, anchors, expected_loading,
#'   settings and call. b is S by D or N by S by D with exposures.
#'   population_seed records initialization-only Gamma shape/rate estimates,
#'   not a prior specification or a population posterior. NMF initialization also
#'   returns S0, initial_cluster,
#'   centers and merge/split history. The legacy method returns slot_permutation.
#' @export
#' @examples
#' set.seed(1)
#' Y <- matrix(rpois(60, 2), 10, 6)
#' initial <- miso_init(Y, K = 3, D = 2, S = 2, nmf_max_iters = 10)
#' fit <- miso_fit(Y, init = initial, S = 2, max_iters = 5, tol = 0)
#' # For an existing NMF fit, supply both L (N by K) and F (K by M).
#' # initial <- miso_init(Y, L = L, F = F, D = 2, S = 3)
miso_init <- function(Y, F = NULL, D = NULL, K = NULL,
                      poisson_susie = NULL, phi0 = 0.1, tau = 0.9,
                      min_fraction = 0.05, gamma_floor = 0.05,
                      susie_max_iters = 100, susie_tol = 1e-6,
                      nmf_max_iters = 200, seed = 1, eps = 1e-12,
                      keep_intermediates = FALSE, S = NULL, L = NULL,
                      method = c("nmf", "poisson_susie"), n = NULL,
                      unsupported_fraction = 0.01) {
  Y = .miso_prepare_counts(Y)
  .miso_stop(length(eps) == 1 && is.finite(eps) && eps > 0,
             "eps must be positive and finite.")
  .miso_stop(length(seed) == 1 && is.finite(seed), "seed must be finite.")
  cached = !is.null(poisson_susie)
  if (missing(method) && cached) method = "poisson_susie"
  method = match.arg(method)
  if (!is.null(n)) {
    .miso_stop(is.numeric(n) && is.null(dim(n)) && length(n) == nrow(Y) &&
                 all(is.finite(n)) && all(n > 0),
               "n must contain one positive finite exposure per observation.")
    n = as.numeric(n)
  }
  if (method == "nmf") {
    .miso_stop(!cached, "poisson_susie requires method = 'poisson_susie'.")
    .miso_stop(missing(min_fraction) && missing(susie_max_iters) && missing(susie_tol),
               "min_fraction and susie_* controls require method = 'poisson_susie'.")
  } else {
    .miso_stop(is.null(L) && is.null(S) && missing(unsupported_fraction),
               "L, S and unsupported_fraction require method = 'nmf'.")
  }
  if (cached) {
    .miso_stop(is.list(poisson_susie) && length(dim(poisson_susie$gamma)) == 3,
               "poisson_susie must be a Poisson-SuSiE result.")
    if (is.null(D)) D = dim(poisson_susie$gamma)[2]
    if (is.null(F)) F = poisson_susie$F
    .miso_stop(!is.null(F), "Supply the dictionary used for poisson_susie.")
  }
  .miso_positive_integer(D, "D")
  if (!is.null(S)) {
    .miso_positive_integer(S, "S")
    .miso_stop(S <= nrow(Y), "S must not exceed the number of observations.")
  }
  supplied_L = !is.null(L)
  if (supplied_L) {
    .miso_stop(!is.null(F), "Supply F together with L.")
    .miso_stop(is.matrix(L) && is.numeric(L) &&
                 identical(dim(L), c(nrow(Y), nrow(F))) &&
                 all(is.finite(L)) && all(L >= 0),
               "L must be a finite nonnegative N by K matrix matching Y and F.")
  }
  if (is.null(F)) {
    .miso_positive_integer(K, "K (required when F is NULL)")
    .miso_positive_integer(nmf_max_iters, "nmf_max_iters")
    nmf = .poisson_nmf_fit(Y, K, max_iters = nmf_max_iters, seed = seed)
    F = nmf$F
    if (method == "nmf") L = nmf$L
  }
  .miso_validate_dictionary(F, ncol(Y))
  if (!is.null(K)) {
    .miso_positive_integer(K, "K")
    .miso_stop(nrow(F) == K, "K must agree with nrow(F).")
  }
  if (method == "nmf") {
    scale = rowSums(F)
    .miso_stop(all(is.finite(scale)) && all(scale > 0),
               "F must have finite positive row sums for NMF initialization.")
    F = F / scale
    if (supplied_L) L = sweep(L, 2, scale, "*")
    if (is.null(L)) {
      .miso_positive_integer(nmf_max_iters, "nmf_max_iters")
      L = .poisson_nmf_fit(Y, nrow(F), max_iters = nmf_max_iters,
                          seed = seed, F = F)$L
    }
    if (!supplied_L) L[.miso_row_sums(Y) == 0, ] = 0
    if (!is.null(n)) L = sweep(L, 1, n, "/")
    initial = .miso_initialize_from_loadings(
      L, D, S, tau, gamma_floor, phi0, seed, eps, unsupported_fraction)
    initial$settings = list(method = method, D = D, K = nrow(F), S = S,
      tau = tau, gamma_floor = gamma_floor, nmf_max_iters = nmf_max_iters,
      seed = seed, supplied_loadings = supplied_L,
      unsupported_fraction = unsupported_fraction)
  } else {
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
    if (!is.null(n)) poisson_susie$beta = sweep(poisson_susie$beta, 1, n, "*")
    initial = .miso_initialize_from_susie(
      poisson_susie, D, tau, min_fraction, gamma_floor, phi0, eps)
    initial$phi0 = rep(phi0, length.out = initial$S)
    initial$settings = list(method = method, D = D, K = nrow(F), tau = tau,
      min_fraction = min_fraction, gamma_floor = gamma_floor,
      susie_max_iters = susie_max_iters, susie_tol = susie_tol,
      nmf_max_iters = nmf_max_iters, seed = seed, reused_poisson_susie = cached)
    if (keep_intermediates) initial$poisson_susie = poisson_susie
  }
  initial$F = F
  if (!is.null(n)) initial$n = n
  initial$settings$loading_scale = if (is.null(n)) "counts" else "per_unit_n"
  initial$population_seed = list(shape=initial$shape_seed,rate=initial$rate_seed)
  N = nrow(Y); S = dim(initial$gamma)[1]; D = dim(initial$gamma)[2]
  initial$a = array(rep(as.vector(initial$shape_seed),each=N),c(N,S,D)) + .miso_row_sums(Y)/D
  initial$b = .miso_posterior_rates(initial$rate_seed,n)
  initial$shape_seed = initial$rate_seed = NULL
  initial$format_version = 2L
  initial$loading_units = if(is.null(n)) "counts" else "per_unit_n"
  initial$input_dimnames = dimnames(Y)
  initial$call = match.call()
  initial
}
