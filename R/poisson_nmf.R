## Small Poisson NMF implementation for producing the initialization dictionary.


## Evaluate LF only at observed entries, with O(nnz(Y)) workspace rather
## than constructing either N x M means or nnz(Y) x K gathered factors.
.miso_nmf_observed_mean <- function(L, F, rows, columns) {
  value = numeric(length(rows))
  for (k in seq_len(ncol(L))) value = value + L[rows, k] * F[k, columns]
  value
}

.poisson_nmf_fit <- function(Y, K, max_iters = 200, tol = 1e-6,
                            seed = 1, eps = 1e-10, F = NULL) {
  Y = .miso_prepare_counts(Y)
  sparse = .miso_is_sparse(Y)
  .miso_stop(K >= 1 && K == as.integer(K), "K must be a positive integer.")
  set.seed(seed)
  N = nrow(Y)
  M = ncol(Y)
  L = matrix(rexp(N * K), N, K)
  update_F = is.null(F)
  F = if (update_F) {
    .miso_normalize_rows(matrix(rexp(K * M), K, M), eps)
  } else {
    .miso_validate_dictionary(F, M)
    .miso_stop(nrow(F) == K, "K must agree with nrow(F).")
    .miso_stop(all(is.finite(rowSums(F))) && all(rowSums(F) > 0),
               "Fixed F must have finite positive row sums.")
    F / rowSums(F)
  }
  objective = rep(NA_real_, max_iters)
  if (sparse) {
    rows = Y@i + 1L
    columns = rep.int(seq_len(M), diff(Y@p))
    log_factorial = sum(lgamma(Y@x + 1))
  } else {
    log_factorial = sum(lgamma(Y + 1))
  }

  for (iteration in seq_len(max_iters)) {
    if (sparse) {
      ratio = Y
      ratio@x = Y@x / pmax(.miso_nmf_observed_mean(L, F, rows, columns), eps)
    } else {
      ratio = Y / pmax(L %*% F, eps)
    }
    L = L * as.matrix(ratio %*% t(F)) /
      matrix(rowSums(F), N, K, byrow = TRUE)
    L = pmax(L, eps)

    if (update_F) {
      if (sparse) {
        ratio@x = Y@x / pmax(.miso_nmf_observed_mean(L, F, rows, columns), eps)
      } else {
        ratio = Y / pmax(L %*% F, eps)
      }
      F = F * as.matrix(t(L) %*% ratio) /
        matrix(colSums(L), K, M)
      F = pmax(F, eps)
      row_scale = rowSums(F)
      F = F / row_scale
      L = sweep(L, 2, row_scale, "*")
    }

    log_term = if (sparse) {
      sum(Y@x * log(pmax(.miso_nmf_observed_mean(L, F, rows, columns), eps)))
    } else {
      sum(Y * log(pmax(L %*% F, eps)))
    }
    ## Sum the actual rate analytically; eps protects log/division only.
    objective[iteration] = log_term -
      sum(colSums(L) * rowSums(F)) - log_factorial
    if (iteration > 1) {
      relative_change = abs(objective[iteration] - objective[iteration - 1]) /
        (abs(objective[iteration - 1]) + 1)
      if (relative_change < tol) break
    }
  }

  list(
    L = L,
    F = F,
    log_likelihood = objective[seq_len(iteration)],
    n_iter = iteration
  )
}
