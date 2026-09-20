
set.seed(21)

## The hand-constructed row-wise fit checks the expected-loading initializer,
## including support recovery and slot alignment.
N = 12
D = 2
K = 4
gamma_ps = array(1e-4, c(N, D, K))
alpha_ps = matrix(5, N, D)
beta_ps = matrix(1, N, D)
for (i in seq_len(N)) {
  active = if (i <= N / 2) 1:2 else 3:4
  active = if (i %% 2 == 0) rev(active) else active
  gamma_ps[i, 1, active[1]] = 1
  gamma_ps[i, 2, active[2]] = 1
}
for (i in seq_len(N)) {
  for (d in seq_len(D)) gamma_ps[i, d, ] = gamma_ps[i, d, ] /
    sum(gamma_ps[i, d, ])
}

initialization = .miso_initialize_from_susie(
  list(gamma = gamma_ps, alpha = alpha_ps, beta = beta_ps),
  D = D, tau = 0.9, min_fraction = 0.2, phi0 = 0.1
)
stopifnot(
  initialization$S == 2,
  all(sort(lengths(split(initialization$cluster, initialization$cluster))) ==
        c(6, 6)),
  all(vapply(initialization$support, length, 0L) == D),
  all(is.finite(initialization$alpha0)),
  all(is.finite(initialization$beta0)),
  max(abs(initialization$phi -
            (0.1 + colSums(initialization$omega)))) < 1e-10
)

## The stochastic test exercises every variational update on a small data set.
M = 25
F = .miso_normalize_rows(matrix(rexp(K * M), K, M))
true_motif = rep(1:2, each = N / 2)
true_loading = matrix(0, N, K)
for (i in seq_len(N)) {
  active = if (true_motif[i] == 1) 1:2 else 3:4
  true_loading[i, active] = rgamma(D, shape = 12, rate = 1)
}
Y = matrix(rpois(N * M, as.vector(true_loading %*% F)), N, M)

nmf = .poisson_nmf_fit(Y, K, max_iters = 10, seed = 21)
stopifnot(
  all(abs(rowSums(nmf$F) - 1) < 1e-8),
  tail(nmf$log_likelihood, 1) >= nmf$log_likelihood[1] - 1e-8
)

fit = .miso_fit(
  Y, F,
  gamma_init = initialization$gamma,
  alpha0_init = initialization$alpha0,
  beta0_init = initialization$beta0,
  omega_init = initialization$omega,
  phi0 = 0.1,
  max_iters = 15,
  min_iters = 3,
  mf_iters = 3,
  update_F = TRUE
)

stopifnot(
  is.finite(fit$final_elbo),
  all(abs(rowSums(fit$omega) - 1) < 1e-8),
  all(abs(apply(fit$gamma, c(1, 2), sum) - 1) < 1e-8),
  all(abs(rowSums(fit$F) - 1) < 1e-8),
  max(abs(fit$phi - (fit$phi0 + colSums(fit$omega)))) < 1e-8,
  all(dim(.miso_expected_loadings(fit)) == c(N, K)),
  all(dim(predict(fit)) == c(N, M))
)

## Compact initialization agrees with optional diagnostics, including singleton dimensions.
local({
  check_initialization = function(Y, F, D) {
    compact = poisson_susie_fit(Y, F, D, max_iters = 5, seed = 17)
    diagnostic = poisson_susie_fit(
      Y, F, D, max_iters = 5, seed = 17, keep_fits = TRUE
    )
    N = nrow(Y)
    K = nrow(F)
    stopifnot(
      identical(dim(compact$alpha), c(N, as.integer(D))),
      identical(dim(compact$beta), c(N, as.integer(D))),
      is.null(compact$fits), length(diagnostic$fits) == N,
      identical(compact$gamma, diagnostic$gamma),
      identical(compact$alpha, diagnostic$alpha),
      identical(compact$beta, diagnostic$beta),
      object.size(compact) < object.size(diagnostic),
      isTRUE(all.equal(.miso_initialize_from_susie(compact, D),
                       .miso_initialize_from_susie(diagnostic, D), tolerance = 1e-12))
    )
  }
  check_initialization(Y, F, D)
  check_initialization(Y[1, , drop = FALSE], F[1, , drop = FALSE], 1L)
  check_initialization(matrix(0, 3, 1), matrix(1, 1, 1), 1L)

  lean = miso_fit(Y, K = K, D = D, max_iters = 1,
                   init_control = list(nmf_max_iters = 5, susie_max_iters = 5))
  full = miso_fit(Y, K = K, D = D, max_iters = 1, keep_initialization = TRUE,
                   init_control = list(nmf_max_iters = 5, susie_max_iters = 5))
  lean$call = full$call = NULL
  stopifnot(is.null(lean$initialization), !is.null(full$initialization),
            isTRUE(all.equal(lean, full[names(lean)], check.attributes = FALSE)),
            object.size(lean) < object.size(full))
})
