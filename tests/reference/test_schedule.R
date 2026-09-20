## Independent scalar-loop reference for the two-pass Jacobi schedule.
local({
  close = function(x, y, tolerance = 1e-7) {
    result = all.equal(x, y, tolerance = tolerance)
    if (!isTRUE(result)) stop(paste(result, collapse = "\n"))
  }
  probability = function(score) {
    value = exp(score - max(score))
    value / sum(value)
  }
  normalize = function(value) {
    value = pmax(value, 1e-12)
    value / sum(value)
  }
  gamma_kl = function(a, b, a0, b0) {
    entropy = a - log(b) + lgamma(a) + (1 - a) * digamma(a)
    log_prior = a0 * log(b0) - lgamma(a0) +
      (a0 - 1) * (digamma(a) - log(b)) - b0 * a / b
    -entropy - log_prior
  }
  reference_fit = function(Y, F, initial, iterations, mf_iters,
                           update_prior = TRUE, update_gamma = TRUE, update_F = TRUE) {
    N = nrow(Y); M = ncol(Y); K = nrow(F)
    gamma = initial$gamma; S = dim(gamma)[1]; D = dim(gamma)[2]
    for (k in seq_len(K)) F[k, ] = normalize(F[k, ])
    for (s in seq_len(S)) for (d in seq_len(D)) gamma[s, d, ] = normalize(gamma[s, d, ])
    omega = initial$omega
    for (i in seq_len(N)) omega[i, ] = normalize(omega[i, ])
    alpha0 = initial$alpha0; beta0 = initial$beta0
    alpha = array(0, c(N, S, D)); beta = beta0 + 1
    for (i in seq_len(N)) for (s in seq_len(S)) for (d in seq_len(D)) {
      alpha[i, s, d] = alpha0[s, d] + sum(Y[i, ]) / D
    }
    phi0 = rep(0.1, S); phi = phi0 + colSums(omega)
    history = numeric(iterations)
    for (iteration in seq_len(iterations)) {
      xi = array(0, c(N, S, D, M))
      for (i in seq_len(N)) for (s in seq_len(S)) for (m in seq_len(M)) {
        score = vapply(seq_len(D), function(d) digamma(alpha[i, s, d]) -
          log(beta[s, d]) + sum(gamma[s, d, ] * log(pmax(F[, m], 1e-12))), 0.0)
        xi[i, s, , m] = probability(score)
      }
      C = array(0, c(S, D, M))
      for (s in seq_len(S)) for (d in seq_len(D)) {
        for (i in seq_len(N)) alpha[i, s, d] = alpha0[s, d] + sum(Y[i, ] * xi[i, s, d, ])
        for (m in seq_len(M)) C[s, d, m] = sum(omega[, s] * Y[, m] * xi[, s, d, m])
      }
      beta = beta0 + 1
      if (update_prior) for (s in seq_len(S)) for (d in seq_len(D)) {
        if (sum(omega[, s]) <= 1e-12) next
        mean_lambda = sum(omega[, s] * alpha[, s, d] / beta[s, d]) / sum(omega[, s])
        mean_log = sum(omega[, s] * (digamma(alpha[, s, d]) - log(beta[s, d]))) / sum(omega[, s])
        delta = max(log(mean_lambda) - mean_log, 1e-12)
        equation = function(a) log(a) - digamma(a) - delta
        shape = if (equation(1e-3) <= 0) 1e-3 else if (equation(1e4) >= 0) 1e4 else
          uniroot(equation, c(1e-3, 1e4), tol = 1e-12)$root
        alpha0[s, d] = shape
        beta0[s, d] = shape / max(mean_lambda, 1e-12)
      }
      for (r in seq_len(mf_iters)) {
        if (update_gamma) for (s in seq_len(S)) for (d in seq_len(D)) {
          score = vapply(seq_len(K), function(k) sum(C[s, d, ] * log(pmax(F[k, ], 1e-12))), 0.0)
          gamma[s, d, ] = normalize(probability(score))
        }
        if (update_F) for (k in seq_len(K)) {
          T = numeric(M)
          for (m in seq_len(M)) for (s in seq_len(S)) for (d in seq_len(D)) {
            T[m] = T[m] + gamma[s, d, k] * C[s, d, m]
          }
          if (sum(T) > 0) F[k, ] = normalize(T)
        }
      }
      # Explicit entropy form, independent of the production log-sum-exp bound.
      L = matrix(0, N, S)
      for (i in seq_len(N)) for (s in seq_len(S)) {
        for (m in seq_len(M)) {
          score = vapply(seq_len(D), function(d) digamma(alpha[i, s, d]) -
            log(beta[s, d]) + sum(gamma[s, d, ] * log(pmax(F[, m], 1e-12))), 0.0)
          x = probability(score); positive = x > 0
          L[i, s] = L[i, s] + Y[i, m] * sum(x[positive] * (score[positive] - log(x[positive])))
        }
        for (d in seq_len(D)) L[i, s] = L[i, s] - alpha[i, s, d] / beta[s, d] -
          gamma_kl(alpha[i, s, d], beta[s, d], alpha0[s, d], beta0[s, d])
      }
      for (r in seq_len(30)) {
        for (i in seq_len(N)) omega[i, ] = probability(L[i, ] + digamma(phi) - digamma(sum(phi)))
        new_phi = phi0 + colSums(omega)
        change = max(abs(new_phi - phi) / pmax(phi, 1))
        phi = new_phi
        if (change < 1e-10) break
      }
      expected_log_pi = digamma(phi) - digamma(sum(phi))
      log_B = function(a) sum(lgamma(a)) - lgamma(sum(a))
      pi_kl = log_B(phi0) - log_B(phi) + sum((phi - phi0) * expected_log_pi)
      value = -pi_kl - sum(gamma * log(pmax(gamma, 1e-12) * K))
      for (i in seq_len(N)) for (s in seq_len(S)) {
        value = value + omega[i, s] * (L[i, s] + expected_log_pi[s] - log(max(omega[i, s], 1e-12)))
      }
      history[iteration] = value
    }
    list(alpha = alpha, beta = beta, gamma = gamma, F = F, omega = omega, phi = phi,
         alpha0 = alpha0, beta0 = beta0, component_elbo = L, elbo = history)
  }

  set.seed(20260920)
  N = 6L; S = 2L; D = 2L; K = 3L; M = 7L
  Y = matrix(rpois(N * M, 3), N, M); Y[1, ] = 0; Y[, M] = 0
  F = .miso_normalize_rows(matrix(runif(K * M, .1, 1), K, M))
  initial = list(gamma = .miso_normalize_gamma(array(runif(S * D * K), c(S, D, K))),
                 alpha0 = matrix(c(1, 4, 2, 3), S, D),
                 beta0 = matrix(c(1, 2, 3, 1), S, D),
                 omega = matrix(rep(c(.9, .1), each = N), N, S))
  inputs = list(Y)
  if (requireNamespace("Matrix", quietly = TRUE)) inputs[[2]] = Matrix::Matrix(Y, sparse = TRUE)
  for (mf_iters in c(1L, 3L)) for (iterations in c(1L, 4L)) {
    expected = reference_fit(Y, F, initial, iterations, mf_iters)
    for (input in inputs) for (width in c(1L, 3L, M)) {
      fit = miso_fit(input, F, init = initial, max_iters = iterations,
                     mf_iters = mf_iters, tol = 0, block_size = width)
      close(fit[names(expected)], expected)
      stopifnot(identical(fit$final_elbo, tail(fit$elbo, 1)),
                max(abs(fit$beta - fit$beta0 - 1)) > .01,
                all(diff(fit$elbo) >= -1e-8))
    }
  }
  # Exercise all optional fixed-parameter combinations against the same equations.
  for (prior in c(FALSE, TRUE)) for (gamma in c(FALSE, TRUE)) for (dictionary in c(FALSE, TRUE)) {
    expected = reference_fit(Y, F, initial, 3, 2, prior, gamma, dictionary)
    for (input in inputs) {
      fit = miso_fit(input, F, init = initial, max_iters = 3, tol = 0,
                     mf_iters = 2, update_prior = prior, update_gamma = gamma, update_F = dictionary)
      close(fit[names(expected)], expected)
    }
  }
  one = miso_fit(Y, F, init = initial, max_iters = 1, mf_iters = 1)
  three = miso_fit(Y, F, init = initial, max_iters = 1, mf_iters = 3)
  stopifnot(max(abs(one$F - three$F)) > 1e-5)

  # Count the passes, including early stopping: no hidden final refresh and
  # no additional xi passes for extra MF repetitions. Bound passes return no counts.
  for (input in inputs) for (learn_F in c(FALSE, TRUE)) {
    exported = .miso_test_clone()
    internal = environment(exported$miso_fit)
    kernel_name = if (inherits(input, "sparseMatrix")) ".miso_sparse_pass" else ".miso_dense_pass"
    kernel = internal[[kernel_name]]; phases = character()
    internal[[kernel_name]] = function(...) {
      args = list(...); result = kernel(...)
      bound = isTRUE(args$evaluate_bound)
      phases <<- c(phases, if (bound) "mixture" else "susie")
      if (bound) stopifnot(is.null(result$allocated), is.null(result$C), is.null(result$score))
      if (!learn_F) stopifnot(is.null(result$C))
      result
    }
    fit = exported$miso_fit(input, F, init = initial, max_iters = 10, mf_iters = 4,
                             update_F = learn_F, tol = 1e6, min_iters = 2, patience = 1)
    stopifnot(fit$converged, fit$n_iter == 2L,
              identical(phases, rep(c(rep("susie", S), rep("mixture", S)), fit$n_iter)),
              identical(fit$final_elbo, tail(fit$elbo, 1)))
  }
  # Multiple starting states: monotonicity to numerical tolerance with full updates.
  for (seed in 1:8) {
    set.seed(seed)
    initial$gamma = .miso_normalize_gamma(array(rexp(S * D * K), c(S, D, K)))
    fit = miso_fit(Y, F, init = initial, max_iters = 20, tol = 0)
    stopifnot(all(diff(fit$elbo) >= -1e-8 * (1 + abs(head(fit$elbo, -1)))))
  }
})
cat("Two-pass schedule tests passed: explicit equations, pass counts, and ELBO ascent.\n")
