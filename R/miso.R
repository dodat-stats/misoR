#' Fit MiSO with Bayesian mixture parameters and an estimated dictionary
#'
#' Joint population blocks q(alpha_sd,beta_sd) are updated by one-dimensional
#' quadrature. Local Gamma factors use shape a and rate b. The dictionary F
#' remains point-estimated by variational EM. All Gamma parameters use shape/rate.
#' @param Y Dense or sparse nonnegative N by M observations. Integer counts have
#'   a Poisson interpretation; fractional inputs use a generalized Poisson objective.
#' @param F Optional K by M initial dictionary. Rows are normalized. With supplied
#'   init, F must match its dictionary if both are given.
#' @param D,K,S Positive slot, factor, and motif counts. Inferred from init when
#'   supplied. K is required when estimating an initial dictionary. With S=NULL,
#'   automatic initialization retains its discovered supports, not an ELBO selection.
#' @param init Initialization from [miso_init()], a new-format fitted object, or
#'   a list with F and gamma (S by D by K). Optional omega is N by S; a and b must
#'   be supplied together. a is N by S by D; b is S by D for ordinary MiSO or
#'   N by S by D for length adjustment. Optional q_population contains joint
#'   posterior blocks in column-major (s,d) order. Prior incompatibility is an
#'   error; omit these blocks to initialize under a different prior. A fitted
#'   object continues inference on matching observations and feature order.
#' @param phi0 Positive Dirichlet prior scalar or S-vector; default 0.1. A supplied
#'   initialization's value is inherited when this argument is omitted.
#' @param population_prior A list with alpha and beta, each a named shape/rate
#'   pair, or "row_totals" (default). NULL also calibrates from row totals.
#'   An initialization's prior is inherited when this argument is omitted.
#' @param quadrature_control Named controls tol (1e-8), max_order (512), tail_drop
#'   (40), tail_expand (10). Refinement and tail-expansion discrepancies must
#'   both meet tolerance; these are numerical estimates, not certified error bounds.
#' @param max_iters Maximum main Jacobi iterations, excluding warm-up (default 500).
#'   With tol > 0, reaching this limit without convergence raises a warning.
#' @param mixture_max_iters Positive inner membership iteration limit.
#' @param warm_up_iters Nonnegative integer number of initial fixed-F Gauss-Seidel
#'   sweeps (default 20). If omitted when init is a fitted miso_fit object, use zero.
#'   An explicit value overrides this resume behavior. The full variational state
#'   passes to Jacobi unchanged. With update_F=FALSE, F stays fixed in both stages.
#'   tol=0 requests a fixed number of iterations and suppresses convergence warnings.
#' @param tol Relative ELBO change tolerance; zero disables early stopping.
#' @param min_iters,patience Minimum iterations and consecutive small changes.
#' @param update_gamma,update_F Whether to update factor selections and dictionary.
#' @param block_size Number of features per allocation block.
#' @param eps Positive numerical floor for probabilities and logarithms.
#' @param verbose Print iteration progress.
#' @param seed Initialization seed.
#' @param init_control Named controls passed to [miso_init()], used only without init.
#' @param keep_initialization Store the complete starting state in the fit.
#' @return A version-2 miso_fit object. Local variational parameters are a and b;
#'   alpha_mean and beta_mean are population posterior means. q_population holds
#'   joint moments, KL, normalizers, sufficient statistics and numerical diagnostics.
#'   population_mean is E(alpha/beta), not the ratio of means. Also contains F,
#'   gamma, omega, phi, phi0, pi, allocation_mass, z_hat, component_elbo, elbo,
#'   final_elbo, elbo_terms, elbo_constant, converged, n_iter, quadrature and
#'   quadrature_history. ELBO includes the data-only Poisson constant for both
#'   models. The final state is the last full coordinate sweep; coordinates need
#'   not be mutually optimal until convergence. Priors and calibration are retained.
#' @export
#' @examples
#' Y <- matrix(c(5,1,4,2,1,5,2,4),4,2,byrow=TRUE)
#' fit <- miso_fit(Y,F=diag(2),D=1,S=2,max_iters=3,tol=0,
#'                 init_control=list(nmf_max_iters=5))
#' predict(fit, "loadings")
miso_fit <- function(Y,F=NULL,D=NULL,init=NULL,K=NULL,phi0=0.1,
                     population_prior="row_totals",quadrature_control=list(),
                     max_iters=500,mixture_max_iters=30,tol=1e-8,min_iters=5,patience=2,
                     update_gamma=TRUE,update_F=TRUE,block_size=100,eps=1e-12,
                     verbose=FALSE,seed=1,init_control=list(),keep_initialization=FALSE,S,warm_up_iters=20L) {
  .miso_stop(!missing(S) && !is.null(S), "Supply S explicitly, including when init is supplied.")
  fit <- .miso_fit_entry(Y,F,D,init,K,phi0,population_prior,quadrature_control,
    max_iters,mixture_max_iters,tol,min_iters,patience,update_gamma,update_F,
    block_size,eps,verbose,seed,init_control,keep_initialization,S=S,
    phi0_missing=missing(phi0),prior_missing=missing(population_prior),
    warm_up_iters=warm_up_iters,warm_up_missing=missing(warm_up_iters))
  fit$call <- match.call(); fit
}

.miso_fit_entry <- function(Y,F,D,init,K,phi0,population_prior,quadrature_control,
    max_iters,mixture_max_iters,tol,min_iters,patience,update_gamma,update_F,
    block_size,eps,verbose,seed,init_control,keep_initialization,n=NULL,S=NULL,
    phi0_missing=FALSE,prior_missing=FALSE,warm_up_iters=20L,warm_up_missing=FALSE) {
  .miso_stop(is.numeric(warm_up_iters) && length(warm_up_iters)==1L &&
    is.finite(warm_up_iters) && warm_up_iters>=0 && warm_up_iters<.Machine$integer.max &&
    warm_up_iters==trunc(warm_up_iters), "warm_up_iters must be a nonnegative integer.")
  if(warm_up_missing && inherits(init,"miso_fit")) warm_up_iters <- 0L
  warm_up_iters <- as.integer(warm_up_iters)
  Y <- .miso_prepare_counts(Y); N <- nrow(Y)
  .miso_check_exposure(n,N)
  for (name in c("max_iters","mixture_max_iters","min_iters","patience","block_size"))
    .miso_positive_integer(get(name),name)
  .miso_stop(is.numeric(tol) && length(tol)==1L && is.finite(tol) && tol>=0,"Invalid tol.")
  .miso_stop(is.numeric(eps) && length(eps)==1L && is.finite(eps) && eps>0,"Invalid eps.")
  for (name in c("update_gamma","update_F","verbose","keep_initialization"))
    .miso_stop(is.logical(get(name)) && length(get(name))==1L && !is.na(get(name)),paste("Invalid",name))
  for(name in c("D","K","S")) if(!is.null(get(name))) .miso_positive_integer(get(name),name)
  if(!is.null(S)) .miso_stop(S<=N,"S must not exceed the number of observations.")
  .miso_stop(is.list(init_control) && (length(init_control)==0 ||
    (!is.null(names(init_control)) && all(nzchar(names(init_control))) &&
    !anyDuplicated(names(init_control)))),"init_control must be a uniquely named list.")
  allowed <- setdiff(names(formals(miso_init)),c("Y","F","D","K","phi0","seed","n","S"))
  .miso_stop(all(names(init_control) %in% allowed),"Unknown initialization control.")
  if(is.null(init)) {
    init <- do.call(miso_init,c(list(Y=Y,F=F,D=D,K=K,S=S,n=n,phi0=phi0,seed=seed),init_control))
  } else {
    .miso_stop(is.list(init) && !is.null(init$gamma),"init must contain gamma and F.")
    .miso_stop(length(init_control)==0,"init_control is only used without init.")
    .miso_stop(!any(c("alpha0","beta0","alpha","beta") %in% names(init)),
      "Legacy fields in init: convert an old fit with miso_convert_fit() or rebuild initialization.")
    if(!is.null(init$format_version)) .miso_stop(identical(init$format_version,2L),"Unsupported initialization format.")
    if(!is.null(init$loading_units)) .miso_stop(identical(init$loading_units,
      if(is.null(n)) "counts" else "per_unit_n"),"Initialization loading units do not match.")
    if(!is.null(init[["n"]])) .miso_stop(!is.null(n) && isTRUE(all.equal(as.numeric(init[["n"]]),as.numeric(n))),
      "Initialization exposures do not match.")
    if(!is.null(init$input_dimnames)) .miso_stop(identical(init$input_dimnames,dimnames(Y)),
      "Observation or feature names/order differ from initialization.")
    if(phi0_missing && !is.null(init$phi0)) phi0 <- init$phi0
    if(prior_missing && !is.null(init$population_prior)) population_prior <- init$population_prior
    if(!is.null(F) && !is.null(init$F)) .miso_stop(isTRUE(all.equal(
      .miso_normalize_rows(F,eps),.miso_normalize_rows(init$F,eps))),"F conflicts with init$F.")
    if(!is.null(F)) init$F <- F
  }
  .miso_stop(!is.null(init$F),"Supply F or include F in init.")
  F <- init$F; .miso_validate_dictionary(F,ncol(Y))
  .miso_stop(all(rowSums(F)>0),"F rows must have positive mass.")
  dimensions <- dim(init$gamma)
  .miso_stop(length(dimensions)==3L && all(dimensions>0) && dimensions[3]==nrow(F),"Invalid gamma dimensions.")
  for(name in c("S","D","K")) {
    position <- match(name,c("S","D","K"))
    if(!is.null(get(name))) .miso_stop(get(name)==dimensions[position],paste(name,"conflicts with init."))
  }
  S <- dimensions[1]; D <- dimensions[2]
  .miso_stop(S<=N,"S must not exceed the number of observations.")
  .miso_stop(all(is.finite(init$gamma) & init$gamma>=0) &&
    all(apply(init$gamma,c(1,2),sum)>0),"gamma must have positive slice totals.")
  .miso_stop(is.numeric(phi0) && length(phi0) %in% c(1L,S) && all(is.finite(phi0)&phi0>0),"Invalid phi0.")
  prior <- .miso_resolve_prior(population_prior,Y,D,n)
  control <- .miso_quadrature_control(quadrature_control)
  fit <- .miso_fit(Y,F,init,prior,control,phi0,max_iters,mixture_max_iters,tol,
                  min_iters,patience,update_gamma,update_F,block_size,eps,verbose,n,keep_initialization,warm_up_iters)
  if(!fit$converged && tol>0) warning(sprintf(
    "MiSO did not converge within max_iters = %d main iterations (%d warm-up sweeps). Resume with init = fit or increase max_iters.",
    max_iters,fit$warm_up_iters),call.=FALSE)
  fit$input_dimnames <- dimnames(Y); fit
}

.miso_fit <- function(Y,F,init,prior,control,phi0,max_iters,mixture_max_iters,tol,
                     min_iters,patience,update_gamma,update_F,block_size,eps,verbose,n,keep_initialization,warm_up_iters=0L) {
  N <- nrow(Y); S <- dim(init$gamma)[1]; D <- dim(init$gamma)[2]
  # Saved Bayesian posteriors may contain probabilities below eps. Preserve
  # already-normalized states on continuation instead of flooring them again.
  resuming <- inherits(init,"miso_fit") && identical(init$inference_model,"bayesian_mixture")
  if(!resuming || max(abs(rowSums(F)-1))>1e-8) F <- .miso_normalize_rows(F,eps)
  gamma <- init$gamma
  if(!resuming || max(abs(apply(gamma,c(1,2),sum)-1))>1e-8) gamma <- .miso_normalize_gamma(gamma,eps)
  omega <- init$omega
  if(is.null(omega)) omega <- matrix(1/S,N,S)
  .miso_stop(is.matrix(omega) && identical(dim(omega),c(N,S)) &&
    all(is.finite(omega)&omega>=0) && all(rowSums(omega)>0),"Invalid omega.")
  if(!resuming || max(abs(rowSums(omega)-1))>1e-8) omega <- .miso_normalize_rows(omega,eps)
  phi0 <- rep(phi0,length.out=S); phi <- phi0+colSums(omega)
  totals <- .miso_row_sums(Y)
  .miso_stop(is.null(init[["a"]])==is.null(init[["b"]]),"Supply a and b together.")
  a <- init[["a"]]; b <- init[["b"]]
  if(is.null(a)) {
    a <- array(prior$alpha[1]/prior$alpha[2],c(N,S,D)) + totals/D
    b <- .miso_posterior_rates(matrix(prior$beta[1]/prior$beta[2],S,D),n)
  }
  .miso_stop(is.numeric(a) && identical(dim(a),c(N,S,D)) && all(is.finite(a)&a>0),"Invalid a.")
  .miso_stop(is.numeric(b) && identical(dim(b),if(is.null(n)) c(S,D) else c(N,S,D)) &&
    all(is.finite(b)&b>0),"Invalid b dimensions or values for this model.")
  q_population <- init$q_population
  if(is.null(q_population)) q_population <- .miso_population_update(a,b,omega,prior,control) else {
    .miso_stop(is.list(q_population) && length(q_population)==S*D,"Invalid population blocks.")
    q_population <- lapply(q_population,function(h) {
      .miso_stop(is.list(h) && all(c("n","T","R","prior") %in% names(h)) &&
        isTRUE(all.equal(unname(h$prior),.miso_prior_vector(prior))),
        "Population block prior mismatch; omit q_population to initialize under a new prior.")
      .miso_gamma_pair(h$n,h$T,h$R,.miso_prior_vector(prior),control)
    })
  }
  starting <- NULL
  if(keep_initialization) {
    starting <- list(F=F,gamma=gamma,omega=omega,a=a,b=b,q_population=q_population,
      phi0=phi0,population_prior=prior,format_version=2L,
      loading_units=if(is.null(n)) "counts" else "per_unit_n",input_dimnames=dimnames(Y))
    if(!is.null(n)) starting$n <- n
    for(field in c("settings","population_seed","expected_loading","support","anchors",
                   "centers","history","S0","initial_cluster","cluster","pattern"))
      if(!is.null(init[[field]])) starting[[field]] <- init[[field]]
  }
  alpha_mean <- .miso_population_matrix(q_population,"mean_alpha",S,D)
  beta_mean <- .miso_population_matrix(q_population,"mean_beta",S,D)
  prepared <- .miso_prepare_pass(Y,block_size)
  values <- if(.miso_is_sparse(Y)) Y@x else Y
  constant <- sum(totals*log(if(is.null(n)) 1 else n))-sum(lgamma(values+1))
  initial_component <- .miso_component_elbo(prepared,F,gamma,a,b,q_population,eps,n)
  initial_elbo <- .miso_elbo_terms(initial_component,omega,phi,phi0,gamma,q_population,constant,eps)[["ELBO"]]
  rm(initial_component)
  total_limit <- as.double(warm_up_iters)+max_iters
  .miso_stop(total_limit<=.Machine$integer.max,"Combined iteration limit is too large.")
  trajectory <- matrix(NA_real_,total_limit,7L,dimnames=list(NULL,
    c("local_bound","assignment_term","KL_pi","KL_gamma","KL_population","constant","ELBO")))
  diagnostics <- matrix(NA_real_,total_limit,3L,dimnames=list(NULL,c("max_order","refinement_error","tail_check_error")))
  stable <- 0L; converged <- FALSE
  for(iteration in seq_len(total_limit)) {
    warming <- iteration<=warm_up_iters
    main_iteration <- iteration-warm_up_iters
    if(warming) {
      pass <- .miso_gs_sweep(prepared,F,gamma,a,b,alpha_mean,beta_mean,omega,n,eps,update_gamma)
      a <- pass$a; b <- pass$b; gamma <- pass$gamma
    } else {
      pass <- .miso_susie_step(prepared,F,gamma,a,b,alpha_mean,beta_mean,omega,
                              update_gamma,update_F,eps,n)
      a <- pass$a; b <- pass$b
      if(update_gamma) gamma <- if(update_F) .miso_update_gamma(F,pass$C,gamma,eps) else
        .miso_gamma_from_scores(pass$scores,gamma,eps)
    }
    q_population <- .miso_population_update(a,b,omega,prior,control)
    alpha_mean <- .miso_population_matrix(q_population,"mean_alpha",S,D)
    beta_mean <- .miso_population_matrix(q_population,"mean_beta",S,D)
    if(update_F && !warming) F <- .miso_update_F(F,gamma,pass$C,eps)
    rm(pass)
    component <- .miso_component_elbo(prepared,F,gamma,a,b,q_population,eps,n)
    mixture <- .miso_update_dirichlet(component,phi0,phi,mixture_max_iters)
    omega <- mixture$omega; phi <- mixture$phi
    terms <- .miso_elbo_terms(component,omega,phi,phi0,gamma,q_population,constant,eps)
    .miso_stop(all(is.finite(terms)),"Non-finite ELBO terms.")
    trajectory[iteration,] <- terms
    qd <- .miso_quadrature_diagnostics(q_population,S,D)
    diagnostics[iteration,] <- c(max(qd$order),max(qd$refinement_error),max(qd$tail_check_error))
    if(verbose) message(sprintf("%s iteration %d: ELBO %.6f; population KL %.4f",
      if(warming) "GS warm-up" else "Jacobi",if(warming) iteration else main_iteration,terms["ELBO"],terms["KL_population"]))
    if(main_iteration>1L) {
      relative <- abs(terms["ELBO"]-trajectory[iteration-1L,"ELBO"])/(1+abs(trajectory[iteration-1L,"ELBO"]-constant))
      stable <- if(main_iteration>=min_iters && relative<tol) stable+1L else 0L
      if(tol>0 && stable>=patience) {converged <- TRUE; break}
    }
  }
  trajectory <- as.data.frame(trajectory[seq_len(iteration),,drop=FALSE])
  fit <- list(F=F,gamma=gamma,omega=omega,a=a,b=b,alpha_mean=alpha_mean,beta_mean=beta_mean,
    q_population=q_population,population_prior=prior,
    population_mean=.miso_population_matrix(q_population,"mean_alpha_over_beta",S,D),
    phi0=phi0,phi=phi,pi=phi/sum(phi),allocation_mass=colMeans(omega),
    z_hat=max.col(omega,ties.method="first"),component_elbo=component,
    elbo=trajectory$ELBO,final_elbo=tail(trajectory$ELBO,1),elbo_terms=trajectory,
    elbo_constant=constant,elbo_convention="full",converged=converged,n_iter=iteration,
    n_iter_main=as.integer(main_iteration),warm_up_iters=warm_up_iters,
    elbo_stage=c(rep("warm_up",warm_up_iters),rep("main",main_iteration)),
    warm_up_elbo=head(trajectory$ELBO,warm_up_iters),initial_elbo=initial_elbo,
    stop_reason=if(converged) "converged" else if(tol==0) "fixed_iterations" else "max_iters",
    fit_control=list(max_iters=max_iters,tol=tol,min_iters=min_iters,patience=patience,
      warm_up_iters=warm_up_iters,update_F=update_F,update_gamma=update_gamma),
    quadrature=qd,quadrature_history=as.data.frame(diagnostics[seq_len(iteration),,drop=FALSE]),
    quadrature_control=control,format_version=2L,inference_model="bayesian_mixture",
    loading_units=if(is.null(n)) "counts" else "per_unit_n")
  if(!is.null(n)) fit[["n"]] <- n
  if(keep_initialization) fit$initialization <- starting
  class(fit) <- "miso_fit"; fit
}

.miso_component_elbo <- function(prepared,F,gamma,a,b,q_population,eps=1e-12,n=NULL) {
  N <- dim(a)[1]; S <- dim(a)[2]; D <- dim(a)[3]; K <- nrow(F)
  component <- matrix(0,N,S); log_F <- log(pmax(F,eps))
  for(s in seq_len(S)) {
    shape <- matrix(a[,s,],N,D); rate <- .miso_rate_matrix(b,s,N,D)
    elog <- digamma(shape)-log(rate); mu <- shape/rate
    pass <- prepared$pass(prepared$blocks,elog,matrix(gamma[s,,],D,K)%*%log_F,
      evaluate_bound=TRUE,compute_allocated=FALSE,eps=eps)
    component[,s] <- pass$bound
    for(d in seq_len(D)) {
      h <- q_population[[s+(d-1L)*S]]
      entropy <- shape[,d]-log(rate[,d])+lgamma(shape[,d])+(1-shape[,d])*digamma(shape[,d])
      logprior <- h$mean_alpha_log_beta-h$mean_lgamma_alpha+
        (h$mean_alpha-1)*elog[,d]-h$mean_beta*mu[,d]
      component[,s] <- component[,s]-(if(is.null(n)) 1 else n)*mu[,d]+logprior+entropy
    }
  }
  component
}
.miso_elbo_terms <- function(component,omega,phi,phi0,gamma,q_population,constant,eps) {
  local <- sum(omega*component)
  assignment <- sum(sweep(omega,2,.miso_dirichlet_expected_log(phi),"*"))-sum(omega*log(pmax(omega,eps)))
  klpi <- .miso_dirichlet_kl(phi,phi0)
  klgamma <- sum(gamma*(log(pmax(gamma,eps))+log(dim(gamma)[3])))
  klpop <- sum(vapply(q_population,function(h) h$KL,0.0))
  c(local_bound=local,assignment_term=assignment,KL_pi=klpi,KL_gamma=klgamma,
    KL_population=klpop,constant=constant,ELBO=local+assignment-klpi-klgamma-klpop+constant)
}

.miso_susie_step <- function(prepared, F, gamma, a, b, alpha_mean, beta_mean,
                             omega, update_gamma = TRUE, update_F = TRUE,
                             eps = 1e-12, n = NULL) {
  N = dim(a)[1]
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
      digamma(matrix(a[, s, ], N, D)) -
        log(.miso_rate_matrix(b, s, N, D)),
      matrix(gamma[s, , ], D, K) %*% log_F,
      omega = if (update_gamma || update_F) omega[, s] else NULL,
      log_F = if (stream_scores) log_F else NULL, eps = eps
    )
    ## Other motifs do not depend on this motif's a. Its entire feature
    ## pass is complete before replacing any of its shapes.
    a[, s, ] = sweep(pass$allocated, 2, alpha_mean[s, ], "+")
    if (update_F) C[s, , ] = pass$C
    if (stream_scores) scores[s, , ] = pass$score
  }
  list(a = a, b = .miso_posterior_rates(beta_mean, n), C = C, scores = scores)
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

.miso_expected_loadings <- function(fit) {
  N = nrow(fit$omega)
  S = ncol(fit$omega)
  D = dim(fit$gamma)[2]
  K = dim(fit$gamma)[3]
  loading = matrix(0, N, K)

  for (s in seq_len(S)) {
    gamma_s = matrix(fit$gamma[s, , ], D, K)
    expected_lambda = matrix(fit[["a"]][, s, ], N, D) /
      .miso_rate_matrix(fit[["b"]], s, N, D)
    loading = loading + fit$omega[, s] * (expected_lambda %*% gamma_s)
  }
  loading
}

#' Extract fitted mean counts or factor loadings
#'
#' Prediction currently describes the observations used for fitting.
#' The posterior mean loading is averaged over motif responsibilities and
#' factor selections. Means are L times F for [miso_fit()], and
#' `n * (L %*% F)` for [miso_fit_length()]. Length-adjusted loadings are per
#' unit of `n`; they are not constrained to sum to one. Both outputs are dense matrices,
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
  if (type == "loadings") return(loading)
  mean = loading %*% object$F
  if (!is.null(object[["n"]])) mean = object[["n"]] * mean
  mean
}

#' Print a fitted MiSo model
#'
#' @param x A fitted `miso_fit` object.
#' @param ... Additional arguments; currently unused.
#' @return The fitted object invisibly.
#' @export
print.miso_fit <- function(x, ...) {
  cat("MiSO fit (",x$inference_model,")\n",sep="")
  if (!is.null(x[["n"]])) cat("  document-length-adjusted Poisson working likelihood\n")
  cat("  observations:", nrow(x$omega), "\n")
  cat("  motifs:", ncol(x$omega), "\n")
  cat("  iterations:", x$n_iter, "\n")
  if(!is.null(x$n_iter_main)) cat("  warm-up:",x$warm_up_iters,"; main:",x$n_iter_main,"\n")
  cat("  converged:", x$converged, "\n")
  if(!is.null(x$stop_reason)) cat("  stop reason:",x$stop_reason,"\n")
  cat("  final ELBO:", format(x$final_elbo, digits = 8), "\n")
  cat("  allocation mass:",
      paste(format(round(x$allocation_mass, 4), nsmall = 4), collapse = ", "),
      "\n")
  invisible(x)
}
