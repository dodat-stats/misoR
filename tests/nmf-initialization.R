library(misoR)
local({
  close <- function(a, b, tolerance = 1e-8) {
    stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)))
  }
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
  # A and B differ strongly in magnitude, but have almost the same direction.
  # C is closer to B in Euclidean distance. Cosine must merge A and B.
  L <- rbind(c(90, 10, 0), c(90, 10, 0), c(90, 10, 0),
             c(8, 2, 0), c(8, 2, 0), c(0, 0, 10), c(0, 0, 10))
  F <- diag(3)
  Y <- L %*% F
  initial <- miso_init(Y, L = L, F = F, D = 2, tau = .85)
  stopifnot(initial$S0 == 3, initial$S == 3,
            identical(initial$cluster, c(1L, 1L, 1L, 2L, 2L, 3L, 3L)),
            identical(initial$support, list(1L, 1:2, 3L)))
  merged <- miso_init(Y, L = L, F = F, D = 2, tau = .85, S = 2)
  stopifnot(length(merged$history) == 1,
            identical(merged$history[[1]]$groups, c(1L, 2L)),
            identical(merged$cluster, c(1L, 1L, 1L, 1L, 1L, 2L, 2L)))
  close(merged$centers[1, ], colMeans(L[1:5, ]))
  close(merged$phi, .1 + c(5, 2))
  close(merged$history[[1]]$cosine, sum(c(90, 10) * c(8, 2)) /
          sqrt(sum(c(90, 10)^2) * sum(c(8, 2)^2)))
  # A union larger than D is truncated by the pooled mean, with factor-index ties.
  pair <- rbind(c(8, 2, 0), c(8, 0, 2))
  truncated <- miso_init(pair, L = pair, F = F, D = 2, S = 1, tau = .95)
  stopifnot(truncated$S0 == 2, identical(truncated$support[[1]], 1:2))
  close(truncated$centers[1, ], c(8, 1, 1))
  close(truncated$gamma[1, 1, ], c(1 - .05 + .05/3, .05/3, .05/3))
  # All original observations remain exactly once after sequential splits.
  split <- miso_init(Y, L = L, F = F, D = 2, tau = .85, S = 6, seed = 19)
  again <- miso_init(Y, L = L, F = F, D = 2, tau = .85, S = 6, seed = 19)
  close(split$cluster, again$cluster)
  stopifnot(split$S == 6, all(tabulate(split$cluster, nbins = 6) > 0),
            all(rowSums(split$omega) == 1), length(split$history) == 3,
            length(unique(vapply(split$support, paste, "", collapse = ","))) < 6)
  singletons <- miso_init(Y, L = L, F = F, D = 2, tau = .85, S = nrow(Y))
  stopifnot(all(colSums(singletons$omega) == 1),
            all(is.finite(singletons$population_seed$shape)), all(is.finite(singletons$population_seed$rate)),
            all(singletons$population_seed$shape <= 1e4), all(singletons$population_seed$rate > 0))
  # Gamma point-estimate MLE equations on a nondegenerate group.
  varying <- rbind(c(2, 0), c(4, 0), c(8, 0))
  g <- miso_init(varying, L = varying, F = diag(2), D = 1)
  close(as.numeric(g$population_seed$shape/g$population_seed$rate), mean(varying[, 1]))
  close(as.numeric(log(g$population_seed$shape) - digamma(g$population_seed$shape)),
        log(mean(varying[, 1])) - mean(log(varying[, 1])), tolerance = 1e-7)
  # Unsupported-slot means jointly use 1% of the group's total loading.
  soft <- miso_init(varying, L = varying, F = diag(2), D = 3)
  close(sum((soft$population_seed$shape/soft$population_seed$rate)[1, 2:3]), .01 * mean(rowSums(varying)))
  close(soft$gamma[1, 2, ], c(.5, .5))
  # Rescaling a dictionary and its loadings must preserve the initialization.
  scale <- c(2, 5, 3)
  scaled <- miso_init(Y, L = sweep(L, 2, scale, "/"),
                     F = diag(scale), D = 2, tau = .85, S = 2)
  fields <- c("F", "expected_loading", "centers", "gamma", "population_seed", "omega")
  close(scaled[fields], merged[fields])
  close(scaled$expected_loading %*% scaled$F, Y)
  # Fixed-F Poisson NMF recovers loadings for a disjoint dictionary, dense/sparse.
  for (counts in list(Y, Matrix::Matrix(Y, sparse = TRUE))) {
    fixed <- miso_init(counts, F = F, D = 2, tau = .85, S = 2, nmf_max_iters = 10)
    close(fixed$expected_loading, L)
    close(fixed$F, F)
    close(fixed$omega, merged$omega)
  }
  # tau=1, empty observations, and singleton dimensions have finite priors.
  zeros <- miso_init(matrix(0, 3, 2), L = matrix(0, 3, 2), F = diag(2), D = 2, S = 3)
  stopifnot(all(lengths(zeros$support) == 0), all(is.finite(zeros$population_seed$rate)))
  full <- miso_init(matrix(c(.1, .2, .3), 1), L = matrix(c(.1, .2, .3), 1),
                   F = F, D = 2, tau = 1)
  stopifnot(identical(full$support[[1]], 2:3))
  one <- miso_init(matrix(0, 1, 1), F = matrix(1), D = 1, S = 1)
  stopifnot(one$S0 == 1, all(is.finite(one$population_seed$rate)))
  # Exposure conversion precedes centers and Gamma prior estimation.
  n <- rowSums(Y)
  adjusted <- miso_init(Y, L = L, F = F, D = 2, tau = .85, S = 2, n = n)
  direct <- miso_init(Y, L = L/n, F = F, D = 2, tau = .85, S = 2)
  close(adjusted[fields], direct[fields])
  close(adjusted$n, n)
  # Both fitters accept target S and can reuse exactly the same initialization.
  controls <- list(L = L, tau = .85)
  rawfit <- miso_fit(warm_up_iters = 0L, Y, F = F, D = 2, S = 2, init_control = controls,
                    keep_initialization = TRUE, max_iters = 6, tol = 0)
  lengthfit <- miso_fit_length(warm_up_iters = 0L, Y, F = F, D = 2, S = 2, init_control = controls,
                              keep_initialization = TRUE, max_iters = 6, tol = 0)
  close(rawfit$initialization[fields], merged[fields])
  close(lengthfit$initialization[fields], adjusted[fields])
  for (fit in list(rawfit, lengthfit)) {
    stopifnot(all(is.finite(fit$elbo)),
              all(diff(fit$elbo) >= -1e-8 * (1 + abs(head(fit$elbo, -1)))))
  }
  reuse <- miso_fit_length(warm_up_iters = 0L, Y, init = adjusted, S = 2, max_iters = 6, tol = 0)
  close(reuse$elbo, lengthfit$elbo)
  # Fully automatic defaults use NMF; explicit legacy initialization still works.
  auto <- miso_fit(tol = 0, warm_up_iters = 0L, Y, K = 3, D = 2, S = 2, max_iters = 2,
                  init_control = list(nmf_max_iters = 10), keep_initialization = TRUE)
  stopifnot(auto$initialization$settings$method == "nmf", ncol(auto$omega) == 2)
  rejects(miso_init(Y, L = L, D = 2))
  rejects(miso_init(Y, L = L[-1, ], F = F, D = 2))
  rejects(miso_init(Y, L = L, F = matrix(0, 3, 3), D = 2))
  for (bad in list(0, 1.5, nrow(Y) + 1, NA_real_, Inf, c(1, 2))) {
    rejects(miso_init(Y, L = L, F = F, D = 2, S = bad))
  }
  rejects(miso_fit(tol = 0, warm_up_iters = 0L, Y, init = merged, S = 3))
  rejects(miso_fit(tol = 0, S = dim(adjusted$gamma)[1], warm_up_iters = 0L, Y, init = adjusted))
  rejects(miso_init(Y, F = F, D = 2, min_fraction = .05))
  rejects(miso_init(Y, F = F, D = 2, method = "poisson_susie", S = 2))
  rejects(miso_init(Y, L = L, F = F, D = 2, unsupported_fraction = 0))
})
cat("Direct NMF initialization tests passed.\n")
