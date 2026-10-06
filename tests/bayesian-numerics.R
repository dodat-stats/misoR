library(misoR)
local({
  ns <- asNamespace('misoR')
  close <- function(a,b,tol=1e-8) {
    msg <- all.equal(a,b,tolerance=tol,check.attributes=FALSE)
    if(!isTRUE(msg)) stop(paste(msg,collapse='\n'))
  }
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
  prior <- list(alpha=c(shape=2,rate=2),beta=c(shape=2,rate=3))
  vec <- unname(c(prior$alpha,prior$beta))
  # Independent integration in alpha (not log alpha), with analytic beta moments.
  Nw <- 3.2; T <- 12; R <- 2.5
  h <- ns$.miso_gamma_pair(Nw,T,R,vec)
  logkernel <- function(x) dgamma(x,2,2,log=TRUE)-Nw*lgamma(x)+(x-1)*R+
    lgamma(2+Nw*x)+2*log(3)-lgamma(2)-(2+Nw*x)*log(3+T)
  shift <- optimize(logkernel,c(.001,20),maximum=TRUE)$objective
  f <- function(x) exp(logkernel(x)-shift)
  z <- integrate(f,0,Inf,rel.tol=1e-10)$value
  E <- function(fun) integrate(function(x) fun(x)*f(x),0,Inf,rel.tol=1e-9)$value/z
  close(h$mean_alpha,E(identity))
  close(h$mean_beta,E(function(x) (2+Nw*x)/(3+T)))
  close(h$mean_alpha_over_beta,E(function(x) x*(3+T)/(1+Nw*x)))
  close(h$mean_alpha_log_beta,E(function(x) x*(digamma(2+Nw*x)-log(3+T))))
  close(h$mean_lgamma_alpha,E(lgamma)); close(h$logZ,log(z)+shift)
  empty <- ns$.miso_gamma_pair(0,0,0,vec)
  close(empty$mean_alpha,1); close(empty$mean_beta,2/3)
  close(empty$mean_alpha_over_beta,3); close(empty$KL,0)
  rejects(ns$.miso_gamma_pair(2,5,0,vec,list(max_order=100)))
  rejects(ns$.miso_gamma_pair(2,5,0,vec,list(unknown=1)))
  # Calibration identifies total scale, independent of F or initialization.
  Y <- rbind(c(.25,1.1,0),c(2.7,0,.35),c(0,0,0),c(3,1,2))
  pr <- miso_population_prior(Y,2)
  close(pr$alpha['shape']/pr$alpha['rate']*pr$beta['rate']/(pr$beta['shape']-1),mean(rowSums(Y))/2)
  close(miso_population_prior(Y,4)$beta['rate'],pr$beta['rate']/2)
  lengthpr <- miso_population_prior(Y[-3,],2,n=rowSums(Y[-3,]))
  close(lengthpr$calibration$mean_slot,.5)
  rejects(miso_population_prior(Y*0,2))
  # One slot: allocation equals row totals, so local updates are analytic.
  initial <- list(F=matrix(c(.2,.3,.5),1),gamma=array(1,c(1,1,1)))
  for(adjusted in c(FALSE,TRUE)) for(sparse in c(FALSE,TRUE)) {
    n <- if(adjusted) c(1.3,4.2,2.1,3) else NULL
    args <- list(Y=if(sparse) Matrix::Matrix(Y,sparse=TRUE) else Y,init=initial,
      population_prior=prior,max_iters=1,tol=0,update_F=FALSE,keep_initialization=TRUE)
    if(adjusted) args$n <- n
    args$S <- dim(args$init$gamma)[1]; args$warm_up_iters <- 0L
    fit <- do.call(if(adjusted) miso_fit_length else miso_fit,args)
    h0 <- fit$initialization$q_population[[1]]
    close(as.vector(fit$a),h0$mean_alpha+rowSums(Y))
    close(as.vector(fit$b),h0$mean_beta+(if(is.null(n)) 1 else n))
    a <- as.vector(fit$a); b <- as.vector(fit$b); h1 <- fit$q_population[[1]]
    mu <- a/b; elog <- digamma(a)-log(b)
    entropy <- a-log(b)+lgamma(a)+(1-a)*digamma(a)
    logprior <- h1$mean_alpha_log_beta-h1$mean_lgamma_alpha+(h1$mean_alpha-1)*elog-h1$mean_beta*mu
    exposure <- if(is.null(n)) rep(1,nrow(Y)) else n
    likelihood <- rowSums(sweep(Y,2,log(as.vector(initial$F)),'*'))+rowSums(Y)*elog-exposure*mu
    constant <- sum(rowSums(Y)*log(exposure))-sum(lgamma(Y+1))
    close(fit$component_elbo[,1],likelihood+logprior+entropy)
    close(fit$final_elbo,sum(likelihood+logprior+entropy)-h1$KL+constant)
    close(predict(fit),outer(exposure*mu,as.vector(initial$F)))
  }
  # Several motifs/slots: full sweeps, sparse/dense and fixed-F streamed scores.
  set.seed(19); Y <- matrix(rpois(30,2),6,5); Y[1,] <- 0
  ini <- miso_init(Y,K=3,D=2,S=3,nmf_max_iters=4)
  for(fixed in c(FALSE,TRUE)) {
    a <- miso_fit(S = dim(ini$gamma)[1], warm_up_iters = 0L, Y,init=ini,population_prior=prior,update_F=!fixed,max_iters=12,tol=0)
    b <- miso_fit(S = dim(ini$gamma)[1], warm_up_iters = 0L, Matrix::Matrix(Y,sparse=TRUE),init=ini,population_prior=prior,
      update_F=!fixed,max_iters=12,tol=0,block_size=1)
    fields <- c('F','gamma','omega','a','b','alpha_mean','beta_mean','population_mean','elbo')
    close(a[fields],b[fields]); stopifnot(min(diff(a$elbo)) > -1e-7)
    if(fixed) close(a$F,ini$F)
    stopifnot(max(a$quadrature$tail_check_error)<1e-8,max(a$quadrature$refinement_error)<1e-8)
    # Unit exposure is identical, with the rate matrix expanded across rows.
    il <- ini; il$n <- rep(1,6); il$loading_units <- 'per_unit_n'
    il$b <- array(rep(as.vector(ini$b),each=6),dim(ini$a))
    l <- miso_fit_length(S = dim(il$gamma)[1], warm_up_iters = 0L, Y,init=il,n=il$n,population_prior=prior,update_F=!fixed,max_iters=12,tol=0)
    close(l[setdiff(fields,'b')],a[setdiff(fields,'b')]); close(predict(l),predict(a))
    # Exact warm restart follows the same coordinate trajectory.
    first <- miso_fit(S = dim(ini$gamma)[1], warm_up_iters = 0L, Y,init=ini,population_prior=prior,update_F=!fixed,max_iters=5,tol=0)
    rest <- miso_fit(S = dim(first$gamma)[1], warm_up_iters = 0L, Y,init=first,update_F=!fixed,max_iters=7,tol=0)
    close(tail(a$elbo,7),rest$elbo)
    close(a[setdiff(fields,'elbo')],rest[setdiff(fields,'elbo')])
  }
  # Exposure unit transformation also transforms the beta hyperprior rate.
  n <- seq_len(6)+1; il <- miso_init(Y,K=3,D=2,S=2,n=n,nmf_max_iters=4)
  f <- miso_fit_length(S = dim(il$gamma)[1], warm_up_iters = 0L, Y,init=il,n=n,population_prior=prior,max_iters=10,tol=0)
  scale <- 7; scaled <- il; scaled$n <- n*scale; scaled$b <- il$b*scale
  p <- prior; p$beta['rate'] <- p$beta['rate']/scale
  g <- miso_fit_length(S = dim(scaled$gamma)[1], warm_up_iters = 0L, Y,init=scaled,n=n*scale,population_prior=p,max_iters=10,tol=0)
  close(f$elbo,g$elbo); close(predict(f),predict(g)); close(predict(f,'loadings'),scale*predict(g,'loadings'))
  # Invalid state, prior mismatch, and feature order must not pass silently.
  bad <- ini; bad$b <- NULL; rejects(miso_fit(tol = 0, S = dim(bad$gamma)[1], warm_up_iters = 0L, Y,init=bad))
  rejects(miso_fit(tol = 0, S = dim(f$gamma)[1], warm_up_iters = 0L, Y,init=f)); rejects(miso_fit_length(tol = 0, S = dim(f$gamma)[1], warm_up_iters = 0L, Y,init=f,n=n+1))
  rejects(miso_fit_length(tol = 0, S = dim(f$gamma)[1], warm_up_iters = 0L, Y,init=f,n=n,population_prior=list(alpha=c(shape=2,rate=4),beta=c(shape=2,rate=3))))
  rejects(miso_fit_length(S = dim(il$gamma)[1], warm_up_iters = 0L, Y,init=il,n=n,quadrature_control=list(tol=NA)))
  rejects(miso_fit_length(tol = 0, S = dim(il$gamma)[1], warm_up_iters = 0L, Y,init=il,n=n,update_F=NA))
  for(value in c(-1,NA,Inf)) {badY<-Y;badY[1,1]<-value;rejects(miso_fit(tol = 0, S = dim(ini$gamma)[1], warm_up_iters = 0L, badY,init=ini))}
  # No fitting code stores the old ambiguous field names.
  stopifnot(!any(c('alpha','beta','alpha0','beta0','prior_mean') %in% names(f)))
})
cat('Bayesian population integration, analytic bounds, calibration, coordinate updates, and invariance checks passed.\n')
