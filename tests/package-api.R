library(misoR)
local({
  close <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-8)))
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
  set.seed(41)
  Y <- matrix(rpois(20 * 12, 1), 20, 12)
  F <- matrix(rexp(3 * 12), 3, 12)
  ps <- poisson_susie_fit(Y, F, D = 2, max_iters = 5, seed = 4)
  fresh <- miso_init(Y, F, D = 2, method = "poisson_susie", susie_max_iters = 5, seed = 4,
                    phi0 = 0.2, keep_intermediates = TRUE)
  cached <- miso_init(Y, poisson_susie = ps, phi0 = 0.2)
  fields <- c("F", "gamma", "a", "b", "omega", "phi0", "support")
  close(fresh[fields], cached[fields])
  close(fresh$poisson_susie, ps)
  stopifnot(cached$settings$reused_poisson_susie)
  fit <- miso_fit(S = dim(cached$gamma)[1], warm_up_iters = 0L, Y, init = cached, max_iters = 4, tol = 0)
  close(fit$phi0, cached$phi0)
  override <- miso_fit(tol = 0, S = dim(cached$gamma)[1], warm_up_iters = 0L, Y, init = cached, phi0 = 0.3, max_iters = 1)
  close(override$phi0, rep(0.3, ncol(fit$omega)))
  sparse <- miso_fit(S = dim(cached$gamma)[1], warm_up_iters = 0L, Matrix::Matrix(Y, sparse = TRUE), init = cached, max_iters = 4, tol = 0)
  close(predict(fit), predict(sparse))
  close(predict(fit), predict(fit, type = "loadings") %*% fit$F)
  rejects(miso_init(Y[-1, ], poisson_susie = ps))
  rejects(miso_init(Y, poisson_susie = ps, D = 3))
  rejects(miso_init(Y, F = F + 1, poisson_susie = ps))
  rejects(predict(fit, newdata = Y))
  sm <- summary(fit)
  stopifnot(inherits(sm, "summary.miso_fit"), sum(sm$motifs$n_assigned) == nrow(Y),
            sm$final_elbo == tail(fit$elbo, 1), sm$elbo_decreases == 0,
            nrow(sm$slots) == dim(fit$gamma)[1] * dim(fit$gamma)[2])
  for (j in seq_len(nrow(sm$slots))) {
    slot <- sm$slots[j, ]
    close(slot$probability, fit$gamma[slot$motif, slot$slot, slot$factor])
    close(slot$population_mean, fit$population_mean[slot$motif, slot$slot])
  }
  invisible(capture.output(print(sm)))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  plot(fit)
  order_clusters <- rev(sort(unique(fit$z_hat)))
  data <- plot(fit, type = "loadings", factor_order = 3:1,
               cluster_order = order_clusters, sort_by = 1)
  close(rowSums(data$values), rep(1, nrow(Y)))
  close(data$loadings, predict(fit, "loadings"))
  stopifnot(identical(sort(data$observation_order), seq_len(nrow(Y))),
            identical(data$factor_order, 3:1), identical(data$cluster_order, order_clusters))
  counts <- plot(fit, type = "loadings", normalize = FALSE)
  close(counts$values, counts$loadings[counts$observation_order, , drop = FALSE])
  rejects(plot(fit, type = "loadings", factor_order = c(1, 1, 3)))
  rejects(plot(fit, type = "loadings", cluster_order = 0))
  rejects(plot(fit, type = "loadings", normalize = NA))
  rejects(plot(fit, type = "loadings", sort_by = "unknown"))
  rejects(plot(fit, type = "loadings", sort_by = NA_real_))
  # Different dominant factors per group; unequal totals distinguish fractions
  # from raw loadings. Ties retain observation order despite reordered colors.
  L <- rbind(c(90,10,0), c(1,9,0), c(1,4,0), c(1,9,0),
             c(2,0,8), c(8,0,2), c(2,0,8), c(1,0,0))
  fixture <- structure(list(F=diag(3), z_hat=rep(1:2,each=4),
    omega=cbind(rep(c(1,0),each=4),rep(c(0,1),each=4)),
    a=array(0,c(8,2,3)), b=matrix(1,2,3), gamma=array(0,c(2,3,3))),
    class="miso_fit")
  for(s in 1:2) { fixture$a[,s,] <- L; fixture$gamma[s,,] <- diag(3) }
  ordered <- plot(fixture, type="loadings", sort_by="dominant",
                  factor_order=3:1, cluster_order=2:1)
  stopifnot(identical(ordered$observation_order,c(8L,6L,5L,7L,2L,4L,3L,1L)))
  absolute <- plot(fixture, type="loadings", sort_by="dominant",
                   normalize=FALSE, cluster_order=2:1)
  close(absolute$observation_order,ordered$observation_order)
  close(absolute$values,L[ordered$observation_order,])
  numeric_sort <- plot(fixture,type="loadings",sort_by=2)
  stopifnot(identical(numeric_sort$observation_order,c(1L,3L,2L,4L,5L,6L,7L,8L)))
  # Equal mean fractions choose factor 1, independently of display order.
  fixture$a[,1,] <- fixture$a[,2,] <- rbind(c(8,2,0),c(2,8,0),
    c(2,8,0),c(8,2,0),matrix(0,4,3))
  tied <- plot(fixture,type="loadings",sort_by="dominant",factor_order=3:1)
  stopifnot(identical(tied$observation_order,c(1L,4L,2L,3L,5L,6L,7L,8L)))
  single <- miso_fit(tol = 0, S = 1L, warm_up_iters = 0L, matrix(0, 1, 1), F = matrix(1), D = 1, max_iters = 2,
    population_prior=list(alpha=c(shape=2,rate=2),beta=c(shape=2,rate=2)))
  plot(single)
  data <- plot(single, type = "loadings")
  stopifnot(identical(dim(data$loadings), c(1L, 1L)), nrow(summary(single)$slots) == 1)
})
cat("Initialization, summary, prediction, and plot API tests passed.\n")
