library(misoR)
local({
  ns <- asNamespace('misoR')
  close <- function(x,y) stopifnot(isTRUE(all.equal(x,y,tolerance=1e-8,check.attributes=FALSE)))
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
  set.seed(812);Y<-matrix(rpois(6*7,3),6,7);Y[1,]<-0;Y[,1]<-0
  prior<-list(alpha=c(shape=2,rate=2),beta=c(shape=2,rate=4))
  fields<-c('F','gamma','omega','a','b','phi','alpha_mean','beta_mean','population_mean')
  for(adjusted in c(FALSE,TRUE)) {
    n<-if(adjusted)seq_len(nrow(Y))+1 else NULL
    ini<-miso_init(Y,K=3,D=2,S=2,n=n,nmf_max_iters=5)
    fitter<-if(adjusted)miso_fit_length else miso_fit
    args<-list(Y=Y,S=2,init=ini,population_prior=prior,max_iters=2,tol=0,warm_up_iters=1L,keep_initialization=TRUE)
    if(adjusted)args$n<-n
    fit<-do.call(fitter,args)
    # Independent scalar allocations for one GS sweep (not the blocked kernel).
    ref<-fit$initialization;N<-nrow(Y);D<-2;K<-3;S<-2
    h<-ref$q_population
    for(s in 1:S)for(d in 1:D) {
      counts<-matrix(0,N,ncol(Y));rates<-ns$.miso_rate_matrix(ref$b,s,N,D)
      for(i in 1:N)for(m in seq_len(ncol(Y))) {
        eta<-vapply(1:D,function(j)digamma(ref$a[i,s,j])-log(rates[i,j])+sum(ref$gamma[s,j,]*log(ref$F[,m])),0.0)
        prob<-exp(eta-max(eta));prob<-prob/sum(prob)
        counts[i,m]<-Y[i,m]*prob[d]
      }
      ref$a[,s,d]<-h[[s+(d-1)*S]]$mean_alpha+rowSums(counts)
      rate<-h[[s+(d-1)*S]]$mean_beta
      if(adjusted)ref$b[,s,d]<-rate+n else ref$b[s,d]<-rate+1
      score<-as.vector(log(ref$F)%*%colSums(counts*ref$omega[,s]))
      p<-exp(score-max(score));p<-p/sum(p);p<-pmax(p,1e-12)
      ref$gamma[s,d,]<-p/sum(p)
    }
    control<-ns$.miso_quadrature_control()
    ref$q_population<-ns$.miso_population_update(ref$a,ref$b,ref$omega,prior,control)
    comp<-ns$.miso_component_elbo(ns$.miso_prepare_pass(Y,100),ref$F,ref$gamma,ref$a,ref$b,ref$q_population,1e-12,n)
    mix<-ns$.miso_update_dirichlet(comp,ref$phi0,ref$phi0+colSums(ref$omega),30L)
    ref$omega<-mix$omega;ref$phi<-mix$phi
    rargs<-args;rargs$init<-ref;rargs$warm_up_iters<-0L
    expected<-do.call(fitter,rargs)
    close(fit[fields],expected[fields]);close(tail(fit$elbo,2),expected$elbo)
    stopifnot(fit$n_iter==3,fit$n_iter_main==2,fit$warm_up_iters==1,
      identical(fit$elbo_stage,c('warm_up','main','main')),fit$stop_reason=='fixed_iterations',
      all(diff(c(fit$initial_elbo,fit$elbo))> -1e-7))
    sparse<-args;sparse$Y<-Matrix::Matrix(Y,sparse=TRUE)
    close(do.call(fitter,sparse)[c(fields,'elbo')],fit[c(fields,'elbo')])
    # A split run preserves all state and automatically skips a second warm-up.
    first<-args;first$max_iters<-1;first<-do.call(fitter,first)
    rest<-args;rest$init<-first;rest$max_iters<-1;rest$warm_up_iters<-NULL
    rest<-do.call(fitter,rest)
    close(rest[fields],fit[fields]);close(rest$elbo,tail(fit$elbo,1));stopifnot(rest$warm_up_iters==0)
    # Tiny valid fitted probabilities are retained exactly on continuation.
    tiny<-first;tiny$F[,1]<-1e-16;tiny$F<-tiny$F/rowSums(tiny$F)
    tiny$omega[1,]<-c(1,1e-20);tiny$omega[1,]<-tiny$omega[1,]/sum(tiny$omega[1,])
    tiny$gamma[1,1,]<-c(1,1e-20,1e-20);tiny$gamma[1,1,]<-tiny$gamma[1,1,]/sum(tiny$gamma[1,1,])
    keep<-args;keep$init<-tiny;keep$warm_up_iters<-NULL;keep$max_iters<-1
    kept<-do.call(fitter,keep)$initialization
    stopifnot(identical(kept$F,tiny$F),identical(kept$omega,tiny$omega),identical(kept$gamma,tiny$gamma))
    forced<-args;forced$init<-first;forced$max_iters<-1
    stopifnot(do.call(fitter,forced)$warm_up_iters==1)
    fixed<-args;fixed$update_F<-FALSE;fixed$update_gamma<-FALSE
    f<-do.call(fitter,fixed);close(f$F,fit$initialization$F);close(f$gamma,fit$initialization$gamma)
    fresh<-args;fresh$warm_up_iters<-NULL;fresh$max_iters<-1
    stopifnot(do.call(fitter,fresh)$warm_up_iters==20)
    for(v in list(-1,1.5,NA_real_,Inf,c(1,2),NULL)) {
      bad<-args;bad['warm_up_iters']<-list(v);rejects(do.call(fitter,bad))
    }
    bad<-args;bad$S<-NULL;rejects(do.call(fitter,bad))
    bad$S<-3;rejects(do.call(fitter,bad))
    # Warn for a cap hit; tol=0 deliberately requests fixed work and stays quiet.
    warnings<-character();short<-args;short$max_iters<-1;short$tol<-1e-8
    f<-withCallingHandlers(do.call(fitter,short),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    stopifnot(length(warnings)==1,grepl('max_iters = 1',warnings),f$stop_reason=='max_iters',!f$converged)
    warnings<-character();short$tol<-0
    withCallingHandlers(do.call(fitter,short),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    stopifnot(length(warnings)==0)
    path<-miso_fit_path(Y,2:1,init=ini,model=if(adjusted)'length'else'count',n=n,
      population_prior=prior,max_iters=2,tol=0,warm_up_iters=3L)
    stopifnot(identical(path$summary$warm_up_iters,c(3L,0L)),all(path$summary$main_iterations==2))
    tmp<-tempfile(fileext='.pdf');grDevices::pdf(tmp);plot(fit);grDevices::dev.off();unlink(tmp)
    sm<-summary(fit);stopifnot(sm$warm_up_iters==1,sm$main_iterations==2,sm$quadrature$problematic_blocks==0)
  }
  stopifnot(formals(miso_fit)$max_iters==500,formals(miso_fit_length)$max_iters==500,
    formals(miso_fit)$tol==1e-8,formals(miso_fit_length)$tol==1e-8)
})
cat('GS scalar reference, stages, resume, sparse/count/length, fixed parameters, warnings and paths passed.\n')
