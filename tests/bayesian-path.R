library(misoR)
local({
  ns <- asNamespace('misoR')
  close <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-8,check.attributes=FALSE)))
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
  set.seed(6); Y <- matrix(rpois(60,3),10,6)
  for(adjusted in c(FALSE,TRUE)) {
    n <- if(adjusted) rowSums(Y) else NULL
    ini <- miso_init(Y,K=3,D=2,S=3,n=n,nmf_max_iters=5)
    # Distinct slot distributions make the alignment optimum unique.
    for(s in 1:3) {ini$gamma[s,1,] <- c(.8,.1,.1); ini$gamma[s,2,] <- c(.1,.8,.1)}
    args <- list(Y=Y,init=ini,max_iters=8,tol=0,update_gamma=FALSE)
    if(adjusted) args$n <- n
    args$S <- dim(args$init$gamma)[1]; args$warm_up_iters <- 0L
    f <- do.call(if(adjusted) miso_fit_length else miso_fit,args)
    merged <- ns$.miso_path_merge(f,Y,block_size=2)
    u <- merged$record$first; v <- merged$record$second
    perm <- merged$record$slot_order[[1]]; g <- merged$init
    close(rowSums(g$omega),rep(1,nrow(Y)))
    close(g$omega[,u],f$omega[,u]+f$omega[,v])
    close(g$phi[u],f$phi[u]+f$phi[v]-f$phi0[1])
    close(g$population_prior,f$population_prior); close(g$F,f$F)
    N <- nrow(Y); D <- dim(f$a)[3]; w <- colSums(f$omega)
    for(d in seq_len(D)) {
      mu1 <- f$a[,u,d]/ns$.miso_rate_matrix(f$b,u,N,D)[,d]
      mu2 <- f$a[,v,perm[d]]/ns$.miso_rate_matrix(f$b,v,N,D)[,perm[d]]
      l1 <- digamma(f$a[,u,d])-log(ns$.miso_rate_matrix(f$b,u,N,D)[,d])
      l2 <- digamma(f$a[,v,perm[d]])-log(ns$.miso_rate_matrix(f$b,v,N,D)[,perm[d]])
      T <- sum(f$omega[,u]*mu1+f$omega[,v]*mu2)
      R <- sum(f$omega[,u]*l1+f$omega[,v]*l2)
      h <- g$q_population[[u+(d-1)*2]]
      close(h$n,w[u]+w[v]); close(h$T,T); close(h$R,R)
      reference <- ns$.miso_gamma_pair(w[u]+w[v],T,R,unname(c(f$population_prior$alpha,f$population_prior$beta)))
      close(h$mean_alpha,reference$mean_alpha); close(h$mean_beta,reference$mean_beta)
      close(g$gamma[u,d,],(w[u]*f$gamma[u,d,]+w[v]*f$gamma[v,perm[d],])/(w[u]+w[v]))
    }
    # A permutation of the second component's slots must be undone by alignment.
    changed <- f
    changed$gamma[v,,] <- f$gamma[v,2:1,]
    changed$a[,v,] <- f$a[,v,2:1]
    if(adjusted) changed$b[,v,] <- f$b[,v,2:1] else changed$b[v,] <- f$b[v,2:1]
    changed$population_mean[v,] <- f$population_mean[v,2:1]
    changed$q_population[c(v,v+3)] <- f$q_population[c(v+3,v)]
    same <- ns$.miso_path_merge(changed,Y)$init
    close(same$gamma,g$gamma); close(same$a,g$a); close(same$b,g$b)
    # Independent, scalar merged allocation pass.
    ell <- (f$omega[,u]*(digamma(matrix(f$a[,u,],N,D))-log(ns$.miso_rate_matrix(f$b,u,N,D)))+
      f$omega[,v]*(digamma(matrix(f$a[,v,],N,D)[,perm])-log(ns$.miso_rate_matrix(f$b,v,N,D)[,perm])))/g$omega[,u]
    allocated <- matrix(0,N,D)
    for(i in seq_len(N)) for(m in seq_len(ncol(Y))) {
      score <- ell[i,]+as.vector(matrix(g$gamma[u,,],D,3)%*%log(pmax(f$F[,m],1e-12)))
      p <- exp(score-max(score)); p <- p/sum(p)
      allocated[i,] <- allocated[i,]+Y[i,m]*p
    }
    for(d in seq_len(D)) close(g$a[,u,d],allocated[,d]+g$q_population[[u+(d-1)*2]]$mean_alpha)
    path <- miso_fit_path(Y,3:1,init=ini,model=if(adjusted) 'length' else 'count',n=n,max_iters=6,tol=0)
    stopifnot(nrow(path$merges)==2,all(is.finite(path$elbo)))
    for(fit in path$fits) {
      close(fit$population_prior,f$population_prior)
      stopifnot(fit$elbo_convention=='full',min(diff(fit$elbo)) > -1e-7)
    }
    tmp <- tempfile(fileext='.pdf'); grDevices::pdf(tmp)
    plot(path); grDevices::dev.off(); unlink(tmp)
  }
  # Empty responsibilities do not produce NaN merged states.
  f$omega[,] <- 0; f$omega[,3] <- 1
  z <- ns$.miso_path_merge(f,Y)$init
  stopifnot(all(is.finite(z$a)),all(is.finite(z$b)))
  rejects(miso_fit_path(Y,c(3,1),K=3,D=2))
  rejects(miso_fit_path(Y,3:1,K=3,D=2,phi0=c(.1,.2,.1)))
})
cat('Merge alignment, pooled statistics, allocation initialization, and Bayesian paths passed.\n')
