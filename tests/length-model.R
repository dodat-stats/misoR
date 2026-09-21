library(misoR)

local({
  close <- function(a, b, tolerance = 1e-8) {
    result <- all.equal(a, b, tolerance = tolerance, check.attributes = TRUE)
    if (!isTRUE(result)) stop(paste(result, collapse = "\n"))
  }
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

  # An exactly conjugate one-topic model: posterior and full ELBO must equal
  # the analytic Gamma-Poisson posterior and integrated evidence, respectively.
  Y <- rbind(c(1, 2, 0), c(10, 4, 6), c(0, 0, 0))
  n <- c(3, 20, 7)
  F <- matrix(c(.2, .3, .5), 1)
  a0 <- 2.4; b0 <- 3.1
  initial <- list(F = F, gamma = array(1, c(1, 1, 1)),
                  alpha0 = matrix(a0), beta0 = matrix(b0))
  fit <- miso_fit_length(Y, n = n, init = initial, max_iters = 1,
                         update_prior = FALSE, update_F = FALSE)
  a <- a0 + rowSums(Y); b <- b0 + n
  close(as.vector(fit$alpha), a)
  close(as.vector(fit$beta), b)
  close(as.vector(predict(fit, "loadings")), a / b)
  close(predict(fit), outer(n * a / b, as.vector(F)))
  evidence <- sum(Y * log(outer(n, as.vector(F))) - lgamma(Y + 1)) +
    sum(a0 * log(b0) - lgamma(a0) + lgamma(a) - a * log(b))
  close(fit$final_elbo + fit$elbo_constant, evidence)
  stopifnot(identical(dim(fit$beta), c(3L, 1L, 1L)))

  # Independent scalar allocation and entropy-form bound for unequal lengths,
  # multiple motifs, slots, and factors. Fix global parameters to isolate the
  # coordinate update and its mixture evidence.
  Y <- rbind(c(6, 2, 1, 0), c(1, 7, 3, 2), c(2, 0, 1, 3))
  n <- rowSums(Y)
  F <- rbind(c(.5, .2, .2, .1), c(.1, .5, .1, .3))
  gamma <- array(c(.7, .2, .3, .8, .3, .8, .7, .2), c(2, 2, 2))
  initial <- list(F = F, gamma = gamma, alpha0 = matrix(c(1, 2, 3, 4), 2),
                  beta0 = matrix(c(2, 1, 5, 3), 2))
  alpha <- beta <- array(0, c(3, 2, 2))
  softmax <- function(x) { p <- exp(x - max(x)); p / sum(p) }
  kl <- function(a, b, a0, b0) {
    entropy <- a - log(b) + lgamma(a) + (1 - a) * digamma(a)
    log_prior <- a0 * log(b0) - lgamma(a0) + (a0 - 1) *
      (digamma(a) - log(b)) - b0 * a / b
    -entropy - log_prior
  }
  for (i in 1:3) for (s in 1:2) {
    beta[i, s, ] <- initial$beta0[s, ] + n[i]
    old_alpha <- initial$alpha0[s, ] + sum(Y[i, ]) / 2
    alpha[i, s, ] <- initial$alpha0[s, ]
    for (j in 1:4) {
      score <- vapply(1:2, function(d) digamma(old_alpha[d]) - log(beta[i, s, d]) +
                        sum(gamma[s, d, ] * log(F[, j])), 0.0)
      alpha[i, s, ] <- alpha[i, s, ] + Y[i, j] * softmax(score)
    }
  }
  bound <- matrix(0, 3, 2)
  for (i in 1:3) for (s in 1:2) {
    for (j in 1:4) {
      score <- vapply(1:2, function(d) digamma(alpha[i, s, d]) - log(beta[i, s, d]) +
                        sum(gamma[s, d, ] * log(F[, j])), 0.0)
      p <- softmax(score)
      bound[i, s] <- bound[i, s] + Y[i, j] * sum(p * (score - log(p)))
    }
    for (d in 1:2) bound[i, s] <- bound[i, s] - n[i] * alpha[i, s, d] / beta[i, s, d] -
      kl(alpha[i, s, d], beta[i, s, d], initial$alpha0[s, d], initial$beta0[s, d])
  }
  fit <- miso_fit_length(Y, init = initial, max_iters = 1,
                         update_prior = FALSE, update_gamma = FALSE, update_F = FALSE)
  close(fit$alpha, alpha); close(fit$beta, beta); close(fit$component_elbo, bound)

  # Independent EB moments/root solve, evaluated at the stored posterior.
  learned <- miso_fit_length(Y, init = initial, max_iters = 1,
                             update_gamma = FALSE, update_F = FALSE)
  close(learned$alpha, alpha); close(learned$beta, beta)
  for (s in 1:2) for (d in 1:2) {
    mu <- mean(alpha[, s, d] / beta[, s, d])
    logmu <- mean(digamma(alpha[, s, d]) - log(beta[, s, d]))
    shape <- uniroot(function(x) log(x) - digamma(x) - log(mu) + logmu,
                     c(1e-3, 1e4), tol = 1e-12)$root
    close(learned$alpha0[s, d], shape)
    close(learned$beta0[s, d], shape / mu)
  }
  stopifnot(max(abs(learned$beta[, 1, 1] - learned$beta0[1, 1] - n)) > 1e-4)

  # Unit lengths recover ordinary MiSo; a common change of units together
  # with transformed priors preserves count predictions and all updates.
  for (prior in c(FALSE, TRUE)) for (dictionary in c(FALSE, TRUE)) {
    original <- miso_fit(Y, init = initial, max_iters = 4, tol = 0,
                         update_prior = prior, update_F = dictionary)
    for (scale in c(1, 7)) {
      scaled_init <- initial
      scaled_init$beta0 <- scale * initial$beta0
      scaled <- miso_fit_length(Y, n = rep(scale, 3), init = scaled_init,
                                max_iters = 4, tol = 0,
                                update_prior = prior, update_F = dictionary)
      fields <- c("alpha", "F", "gamma", "omega", "phi", "alpha0")
      close(scaled[fields], original[fields])
      close(scaled$beta0, scale * original$beta0)
      for (i in 1:3) close(scaled$beta[i, , ], scale * original$beta)
      close(predict(scaled), predict(original))
      close(scale * predict(scaled, "loadings"), predict(original, "loadings"))
      close(scaled$elbo + sum(Y) * log(scale), original$elbo)
    }
  }

  # Sparse/dense, block size, streamed factor scores, and ELBO ascent with
  # unequal lengths, including initialization and all parameter updates.
  for (dictionary in c(FALSE, TRUE)) {
    dense <- miso_fit_length(Y, init = initial, max_iters = 15, tol = 0,
                             update_F = dictionary)
    for (width in c(1L, 3L)) {
      sparse <- miso_fit_length(Matrix::Matrix(Y, sparse = TRUE), init = initial,
                                max_iters = 15, tol = 0, block_size = width,
                                update_F = dictionary)
      fields <- c("alpha", "beta", "alpha0", "beta0", "F", "gamma", "omega", "elbo",
                   "elbo_constant", "n")
      close(dense[fields], sparse[fields])
      close(predict(dense), predict(sparse))
    }
    stopifnot(all(diff(dense$elbo) >= -1e-8 * (1 + abs(head(dense$elbo, -1)))))
  }
  automatic <- miso_fit_length(Y, F, D = 2, max_iters = 3, tol = 0, phi0 = .3,
                               init_control = list(susie_max_iters = 4),
                               keep_initialization = TRUE)
  resumed <- miso_fit_length(Y, init = automatic$initialization, max_iters = 3, tol = 0)
  fields <- c("alpha", "beta", "F", "gamma", "omega", "phi0", "elbo")
  close(automatic[fields], resumed[fields])
  sparse_auto <- miso_fit_length(Matrix::Matrix(Y, sparse = TRUE), F, D = 2,
                                 max_iters = 3, tol = 0, phi0 = .3,
                                 init_control = list(susie_max_iters = 4))
  close(automatic[fields], sparse_auto[fields])
  # Initialization moments must be rescaled before motif-prior estimation.
  ps <- poisson_susie_fit(Y, F, D = 2, max_iters = 4)
  ps$beta <- n * ps$beta
  reference_init <- miso_init(Y, poisson_susie = ps, phi0 = .3)
  init_fields <- c("gamma", "alpha0", "beta0", "omega", "expected_loading")
  close(automatic$initialization[init_fields], reference_init[init_fields])
  nmf_fit <- miso_fit_length(Y, K = 2, D = 2, max_iters = 2,
                             init_control = list(nmf_max_iters = 3, susie_max_iters = 3))
  stopifnot(all(is.finite(predict(nmf_fit))),
            identical(dim(nmf_fit$F), c(2L, 4L)),
            inherits(summary(automatic), "summary.miso_fit"))
  invisible(capture.output(print(automatic)))
  for (input in list(matrix(3, 1, 1), matrix(c(1, 8), 2, 1))) {
    singleton <- miso_fit_length(input, F = matrix(1), D = 1, max_iters = 2)
    stopifnot(identical(dim(singleton$beta), c(nrow(input), 1L, 1L)),
              identical(dim(predict(singleton)), dim(input)))
  }
  rejects(miso_fit_length(rbind(Y, 0), F, D = 2))
  for (bad in list(0, -1, NA_real_, Inf, "3", c(1, 2), c(1, 2, 0), matrix(1, 3))) {
    rejects(miso_fit_length(Y, n = bad, init = initial))
  }
  rejects(miso_fit_length(Y, n = n + 1, init = automatic$initialization))
  rejects(miso_fit(Y, init = automatic$initialization))
  rejects(predict(automatic, newdata = Y))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  plotted <- plot(automatic, type = "loadings", normalize = FALSE)
  close(plotted$loadings, predict(automatic, "loadings"))
})
cat("Length-adjusted MiSo tests passed: conjugacy, evidence, updates, scaling, sparse counts, and API.\n")
