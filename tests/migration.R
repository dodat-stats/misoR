library(misoR)
local({
  file <- if(file.exists('fixtures/legacy-fits.rds')) 'fixtures/legacy-fits.rds' else
    'development/misoR-bayes/tests/fixtures/legacy-fits.rds'
  fixture <- readRDS(file)
  close <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-8,check.attributes=FALSE)))
  for(model in c('count','length')) {
    old <- fixture[[model]]; new <- miso_convert_fit(old,fixture$Y)
    close(predict(new),fixture[[paste0(model,'_prediction')]])
    close(new$elbo,old$elbo+new$elbo_constant)
    close(new$population_mean,old$alpha0/old$beta0)
    stopifnot(new$inference_model=='legacy_eb',is.null(new$q_population),
      all(new$elbo_stage=='legacy'),new$stop_reason=='legacy_unknown',
      all(is.na(summary(new)$slots$alpha_sd)),all(is.na(summary(new)$slots$beta_sd)))
    close(miso_convert_fit(new),new)
    args <- list(Y=fixture$Y,init=new,max_iters=3,tol=0)
    if(model=='length')args$n<-new$n
    args$S <- dim(args$init$gamma)[1]; args$warm_up_iters <- 0L
    fitted <- do.call(if(model=='length') miso_fit_length else miso_fit,args)
    stopifnot(fitted$inference_model=='bayesian_mixture',length(fitted$q_population)>0)
  }
  partial <- miso_convert_fit(fixture$count)
  stopifnot(partial$elbo_convention=='legacy_without_constant')
  close(partial$elbo,fixture$count$elbo)
})
cat('Legacy conversion preserves predictions, identifies inference type, and restores ELBO constants.\n')
