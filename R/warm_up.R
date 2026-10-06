# Fixed-F GS sweep. Recompute normalized allocations over all slots before
# each slot update; population moments and memberships stay fixed for the sweep.
.miso_gs_sweep <- function(prepared, F, gamma, a, b, alpha_mean, beta_mean,
                          omega, n = NULL, eps = 1e-12,
                          update_gamma = TRUE) {
  N <- dim(a)[1]; S <- dim(a)[2]; D <- dim(a)[3]; K <- nrow(F)
  log_F <- log(pmax(F, eps))
  next_b <- .miso_posterior_rates(beta_mean, n)
  for(s in seq_len(S)) for(d in seq_len(D)) {
    pass <- prepared$pass(prepared$blocks,
      digamma(matrix(a[,s,],N,D)) - log(.miso_rate_matrix(b,s,N,D)),
      matrix(gamma[s,,],D,K) %*% log_F,
      omega = omega[,s], log_F = log_F, eps = eps)
    a[,s,d] <- alpha_mean[s,d] + pass$allocated[,d]
    if(is.null(n)) b[s,d] <- next_b[s,d] else b[,s,d] <- next_b[,s,d]
    if(update_gamma) {
      p <- as.numeric(.miso_softmax_rows(matrix(pass$score[d,],1,K)))
      p <- pmax(p,eps); gamma[s,d,] <- p/sum(p)
    }
  }
  list(a=a,b=b,gamma=gamma)
}
