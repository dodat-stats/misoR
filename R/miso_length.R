#' Fit length-adjusted MiSO with Bayesian mixture parameters
#' @inheritParams miso_fit
#' @param n Positive exposures. Defaults to rowSums(Y), requiring nonempty rows.
#'   Counts remain unnormalized. For observed row totals the Poisson model is
#'   a working likelihood, not the exact conditional sampling law of the totals.
#' @details Local rates b are N by S by D and equal E(beta)+n at their coordinate
#'   update. Subsequent population updates hold them fixed. Unit exposures recover
#'   ordinary MiSO with matched initialization and priors. Input NMF loadings passed
#'   to miso_init remain on the count scale and are divided by n there.
#' @return A version-2 miso_fit object as in [miso_fit()], with n and per-exposure
#'   loadings. predict(fit) multiplies L times F by n. ELBO includes the full
#'   data-only constant; do not add elbo_constant a second time. Normalized posterior
#'   mean loadings are descriptive fractions, not E(L/sum(L)).
#' @export
miso_fit_length <- function(Y,F=NULL,D=NULL,init=NULL,K=NULL,n=NULL,phi0=0.1,
                     population_prior="row_totals",quadrature_control=list(),
                     max_iters=500,mixture_max_iters=30,tol=1e-8,min_iters=5,patience=2,
                     update_gamma=TRUE,update_F=TRUE,block_size=100,eps=1e-12,
                     verbose=FALSE,seed=1,init_control=list(),keep_initialization=FALSE,S,warm_up_iters=20L) {
  .miso_stop(!missing(S) && !is.null(S), "Supply S explicitly, including when init is supplied.")
  Y <- .miso_prepare_counts(Y)
  if(is.null(n)) n <- .miso_row_sums(Y)
  n <- as.numeric(n); .miso_check_exposure(n,nrow(Y))
  fit <- .miso_fit_entry(Y,F,D,init,K,phi0,population_prior,quadrature_control,
    max_iters,mixture_max_iters,tol,min_iters,patience,update_gamma,update_F,
    block_size,eps,verbose,seed,init_control,keep_initialization,n=n,S=S,
    phi0_missing=missing(phi0),prior_missing=missing(population_prior),
    warm_up_iters=warm_up_iters,warm_up_missing=missing(warm_up_iters))
  fit$call <- match.call(); fit
}
