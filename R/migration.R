#' Convert a legacy MiSO fit to the new field layout
#' @param object An old fitted object containing alpha, beta, alpha0, beta0.
#' @param Y Optional original data, required to restore the full likelihood constant
#'   for old ordinary fits. Must have matching dimensions and observation order.
#' @details This changes representation only: it does not turn an empirical-Bayes
#'   fit into a Bayesian fit. The converted object is marked legacy_eb and has no
#'   q_population. Use it as init in a new fit to perform Bayesian inference.
#'   Without Y, an ordinary legacy ELBO remains explicitly marked as omitting the
#'   data-only constant. Even after restoring that constant, legacy and Bayesian
#'   ELBOs have different population-prior terms and must not be mixed for model
#'   selection. Refit the converted state under the Bayesian model with explicit S.
#' @return A version-2 miso_fit object suitable for summaries, prediction and
#'   initialization, with a,b,alpha_mean,beta_mean,population_mean.
#' @export
miso_convert_fit <- function(object,Y=NULL) {
  if(identical(object$format_version,2L)) return(object)
  .miso_stop(is.list(object) && all(c("F","gamma","omega","alpha","beta","alpha0","beta0","elbo") %in% names(object)),
    "Not a supported legacy fit.")
  object[["a"]] <- object$alpha; object[["b"]] <- object$beta
  object$alpha_mean <- object$alpha0; object$beta_mean <- object$beta0
  object$population_mean <- object$alpha0/object$beta0
  object[c("alpha","beta","alpha0","beta0","prior_mean")] <- NULL
  object$format_version <- 2L; object$inference_model <- "legacy_eb"
  object$loading_units <- if(is.null(object[["n"]])) "counts" else "per_unit_n"
  constant <- object$elbo_constant
  if(!is.null(Y)) {
    Y <- .miso_prepare_counts(Y)
    .miso_stop(identical(dim(Y),c(nrow(object$omega),ncol(object$F))),"Y dimensions do not match fit.")
    values <- if(.miso_is_sparse(Y)) Y@x else Y
    constant <- sum(.miso_row_sums(Y)*log(if(is.null(object[["n"]])) 1 else object[["n"]]))-sum(lgamma(values+1))
    object$input_dimnames <- dimnames(Y)
  }
  object$elbo_convention <- if(is.null(constant)) "legacy_without_constant" else "full"
  if(!is.null(constant)) {object$elbo <- object$elbo+constant; object$elbo_constant <- constant}
  object$final_elbo <- tail(object$elbo,1)
  # Old initialization objects have a different contract and are intentionally omitted.
  object$initialization <- NULL
  object$elbo_stage <- rep("legacy",length(object$elbo))
  object$stop_reason <- "legacy_unknown"
  object$converted_from <- list(format_version=1L,inference_model="empirical_bayes")
  class(object) <- "miso_fit"; object
}
