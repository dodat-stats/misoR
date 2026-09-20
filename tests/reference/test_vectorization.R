## Mathematical references for exact vectorization changes.
local({
  close = function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-10)))
  set.seed(731)
  for (N in c(1L, 7L)) for (K in c(1L, 5L)) {
    scores = matrix(rnorm(N * K, sd = 1000), N, K,
                    dimnames = list(paste0("r", seq_len(N)), paste0("k", seq_len(K))))
    weight = exp(scores - apply(scores, 1, max))
    stopifnot(identical(.miso_softmax_rows(scores), weight / rowSums(weight)))
  }
  stopifnot(identical(.miso_softmax_rows(matrix(numeric(), 0, 5)),
                      matrix(numeric(), 0, 5)))
  for (S in c(1L, 3L)) for (D in c(1L, 5L)) for (K in c(1L, 4L)) {
    M = 11L
    gamma = array(rexp(S * D * K), c(S, D, K),
                   dimnames = list(paste0("s", seq_len(S)),
                                   paste0("d", seq_len(D)), paste0("k", seq_len(K))))
    gamma[1] = 0
    reference = pmax(gamma, 1e-12)
    for (s in seq_len(S)) reference[s, , ] = .miso_normalize_rows(matrix(reference[s, , ], D))
    close(.miso_normalize_gamma(gamma), reference)
    gamma = reference
    F = .miso_normalize_rows(matrix(rexp(K * M), K, M,
          dimnames = list(paste0("factor", seq_len(K)), paste0("feature", seq_len(M)))))
    C = array(rexp(S * D * M), c(S, D, M))
    C[, , M] = 0
    reference = gamma
    for (s in seq_len(S)) for (d in seq_len(D)) {
      score = as.vector(log(F) %*% C[s, d, ])
      probability = exp(score - max(score))
      reference[s, d, ] = probability / sum(probability)
    }
    close(.miso_update_gamma(F, C, gamma), .miso_normalize_gamma(reference))
    allocated = matrix(0, K, M)
    for (s in seq_len(S)) for (d in seq_len(D)) {
      allocated = allocated + gamma[s, d, ] %o% C[s, d, ]
    }
    dimnames(allocated) = dimnames(F)
    close(.miso_update_F(F, gamma, C), .miso_normalize_rows(allocated))
    close(.miso_update_F(F, gamma, C * 0), F)
    ## A factor with no expected counts must retain its normalized row.
    if (K > 1) {
      gamma[, , K] = 0
      updated_F = .miso_update_F(F, gamma, C)
      close(updated_F[K, ], F[K, ])
    }
  }
  ## ELBO's cached-log-dictionary option, including a
  ## sparse all-zero observation represented by no observed features.
  for (M in c(0L, 1L, 13L)) for (D in c(1L, 5L)) {
    K = 4L
    y = rpois(M, 2)
    F = matrix(runif(K*M, 0.01, 1), K, M)
    gamma = .miso_normalize_rows(matrix(rexp(D*K), D, K))
    xi = .miso_softmax_columns(matrix(rnorm(D*M), D, M))
    alpha = runif(D, 1, 4); beta = runif(D, 1, 4)
    alpha0 = rep(1, D); beta0 = rep(2, D)
    expected_log_F = gamma %*% log(F)
    reference = 0
    for (d in seq_len(D)) {
      reference = reference + sum(xi[d, ]*y*(digamma(alpha[d])-log(beta[d])+
                          expected_log_F[d, ]-log(pmax(xi[d, ], 1e-12))))-alpha[d]/beta[d]
    }
    reference = reference - sum(.miso_gamma_kl(alpha, beta, alpha0, beta0)) -
      sum(gamma*(log(pmax(gamma, 1e-12))-log(1/K)))
    close(.poisson_susie_elbo(y,F,gamma,alpha,beta,alpha0,beta0,xi), reference)
    close(.poisson_susie_elbo(y,F,gamma,alpha,beta,alpha0,beta0,xi,log_F=log(F)), reference)
  }
})
cat("Vectorization reference tests passed.\n")
