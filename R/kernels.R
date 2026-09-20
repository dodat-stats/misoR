## Numerical helpers shared by the MiSo implementation.

.miso_stop <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

.miso_normalize_rows <- function(x, eps = 1e-12) {
  x = pmax(x, eps)
  x / rowSums(x)
}

.miso_normalize_gamma <- function(gamma, eps = 1e-12) {
  .miso_stop(length(dim(gamma)) == 3, "gamma must be an S by D by K array.")
  ## Flatten (s, d) into rows; retaining the original array preserves dimnames.
  gamma[] = .miso_normalize_rows(
    matrix(gamma, nrow = dim(gamma)[1] * dim(gamma)[2]), eps
  )
  gamma
}

.miso_softmax_rows <- function(log_weight) {
  ## D (or the number of motifs) is small; reduce whole columns in C rather
  ## than invoking an R function separately for every observation/nonzero.
  if (nrow(log_weight) == 0) return(log_weight)
  row_max = log_weight[, 1]
  if (ncol(log_weight) > 1) {
    for (k in 2:ncol(log_weight)) row_max = pmax(row_max, log_weight[, k])
  }
  weight = exp(log_weight - row_max)
  weight / rowSums(weight)
}

## Stable log sum exp over columns, without normalized allocations or entropy.
.miso_log_sum_exp_rows <- function(score) {
  if (nrow(score) == 0) return(numeric())
  maximum = score[, 1]
  if (ncol(score) > 1) {
    for (d in 2:ncol(score)) maximum = pmax(maximum, score[, d])
  }
  shift = maximum
  shift[!is.finite(shift)] = 0
  shift + log(rowSums(exp(score - shift)))
}

.miso_gamma_kl <- function(shape, rate, prior_shape, prior_rate) {
  (shape - prior_shape) * digamma(shape) +
    prior_shape * (log(rate) - log(prior_rate)) -
    lgamma(shape) + lgamma(prior_shape) -
    (rate - prior_rate) * shape / rate
}

.miso_dirichlet_expected_log <- function(parameter) {
  digamma(parameter) - digamma(sum(parameter))
}

.miso_dirichlet_kl <- function(parameter, prior) {
  lgamma(sum(parameter)) - sum(lgamma(parameter)) -
    lgamma(sum(prior)) + sum(lgamma(prior)) +
    sum((parameter - prior) * .miso_dirichlet_expected_log(parameter))
}

.miso_gamma_shape_from_moments <- function(mean, mean_log,
                                          min_shape = 1e-3,
                                          max_shape = 1e4,
                                          eps = 1e-12) {
  delta = pmax(log(pmax(mean, eps)) - mean_log, eps)
  shape = ifelse(delta > 0.5, 1 / (2 * delta), 1 / delta)
  shape = pmin(pmax(shape, min_shape), max_shape)

  for (iteration in seq_len(12)) {
    objective = log(shape) - digamma(shape) - delta
    derivative = 1 / shape - trigamma(shape)
    candidate = shape - objective / derivative
    shape = pmin(pmax(candidate, min_shape), max_shape)
  }
  shape
}

.miso_feature_blocks <- function(M, block_size) {
  .miso_stop(length(block_size) == 1 && is.finite(block_size) &&
               block_size >= 1 && block_size == as.integer(block_size),
             "block_size must be a positive integer.")
  split(seq_len(M), ceiling(seq_len(M) / block_size))
}

.miso_is_sparse <- function(Y) inherits(Y, "sparseMatrix")

## Validate stored values only: is.finite(Y) or Y >= 0 can materialize a
## mostly-TRUE dense logical matrix when Y is sparse.
.miso_prepare_counts <- function(Y) {
  if (.miso_is_sparse(Y)) {
    .miso_stop(requireNamespace("Matrix", quietly = TRUE),
               "Install the Matrix package to use sparse counts.")
    values = if ("x" %in% methods::slotNames(Y)) Y@x else 1
    .miso_stop(all(is.finite(values)) && all(values >= 0),
               "Y must contain finite nonnegative counts.")
    Y = methods::as(methods::as(methods::as(Y, "dMatrix"),
                               "generalMatrix"), "CsparseMatrix")
    Y = Matrix::drop0(Y)
    .miso_stop(all(is.finite(Y@x)), "Y must contain finite counts.")
  } else {
    .miso_stop(is.matrix(Y) && is.numeric(Y),
               "Y must be a numeric matrix or a Matrix sparse matrix.")
    .miso_stop(all(is.finite(Y)) && all(Y >= 0),
               "Y must contain finite nonnegative counts.")
  }
  .miso_stop(all(dim(Y) > 0), "Y must have at least one row and one column.")
  Y
}

.miso_row_sums <- function(Y) {
  if (.miso_is_sparse(Y)) as.numeric(Matrix::rowSums(Y)) else rowSums(Y)
}

## Each block contains only positive observations, with global row/column
## indices. Empty feature blocks are omitted. Y must be canonical dgCMatrix.
.miso_sparse_blocks <- function(Y, block_size) {
  blocks = lapply(.miso_feature_blocks(ncol(Y), block_size), function(features) {
    first = Y@p[features[1]] + 1L
    last = Y@p[tail(features, 1) + 1L]
    if (first > last) return(NULL)
    index = seq.int(first, last)
    list(i = Y@i[index] + 1L,
         j = rep.int(features, diff(Y@p[c(features, tail(features, 1) + 1L)])),
         x = Y@x[index])
  })
  Filter(Negate(is.null), blocks)
}

.miso_softmax_columns <- function(log_weight) {
  if (ncol(log_weight) == 0) return(log_weight)
  t(.miso_softmax_rows(t(log_weight)))
}

.miso_positive_integer <- function(value, name) {
  .miso_stop(length(value) == 1 && is.finite(value) && value >= 1 &&
               value == as.integer(value), paste(name, "must be a positive integer."))
}

.miso_validate_dictionary <- function(F, M) {
  .miso_stop(is.matrix(F) && is.numeric(F) && nrow(F) >= 1 && ncol(F) == M,
             "F must be a numeric K by M matrix with K >= 1.")
  .miso_stop(all(is.finite(F)) && all(F >= 0), "F must be finite and nonnegative.")
}

.miso_validate_inputs <- function(Y, F, gamma, alpha0, beta0,
                                 omega, phi0) {
  .miso_prepare_counts(Y)
  .miso_stop(is.matrix(F) && ncol(F) == ncol(Y),
             "F must be a K by M matrix with M = ncol(Y).")
  .miso_stop(all(is.finite(F)) && all(F >= 0),
             "F must be finite and nonnegative.")

  dimensions = dim(gamma)
  .miso_stop(length(dimensions) == 3 && dimensions[3] == nrow(F),
             "gamma must be an S by D by K array with K = nrow(F).")
  .miso_stop(all(dim(alpha0) == dimensions[1:2]) &&
               all(dim(beta0) == dimensions[1:2]),
             "alpha0 and beta0 must be S by D matrices.")
  .miso_stop(all(is.finite(alpha0)) && all(alpha0 > 0) &&
               all(is.finite(beta0)) && all(beta0 > 0),
             "Gamma prior parameters must be finite and positive.")
  .miso_stop(all(dim(omega) == c(nrow(Y), dimensions[1])) &&
               all(is.finite(omega)) && all(omega >= 0),
             "omega must be a nonnegative N by S matrix.")
  .miso_stop(length(phi0) %in% c(1, dimensions[1]) &&
               all(is.finite(phi0)) && all(phi0 > 0),
             "phi0 must contain one positive value or one per motif.")
  invisible(TRUE)
}

.miso_xi <- function(expected_log_lambda, expected_log_F, block,
                      log_normalizer = FALSE) {
  N = nrow(expected_log_lambda)
  D = ncol(expected_log_lambda)
  B = length(block)
  log_xi = array(0, c(N, D, B))

  for (d in seq_len(D)) {
    log_xi[, d, ] = expected_log_lambda[, d] +
      matrix(expected_log_F[d, block], N, B, byrow = TRUE)
  }

  log_max = matrix(log_xi[, 1, ], N, B)
  if (D > 1) {
    for (d in 2:D) {
      log_max = pmax(log_max, matrix(log_xi[, d, ], N, B))
    }
  }
  xi = array(0, c(N, D, B))
  for (d in seq_len(D)) {
    xi[, d, ] = exp(matrix(log_xi[, d, ], N, B) - log_max)
  }
  denominator = matrix(xi[, 1, ], N, B)
  if (D > 1) {
    for (d in 2:D) denominator = denominator + matrix(xi[, d, ], N, B)
  }
  if (log_normalizer) return(log_max + log(denominator))
  for (d in seq_len(D)) {
    xi[, d, ] = matrix(xi[, d, ], N, B) / denominator
  }
  xi
}

## Dense and sparse passes share an output contract; only accumulation differs.
## Passing log_F requests streamed gamma scores instead of feature counts C.
.miso_prepare_pass <- function(Y, block_size) {
  Y = .miso_prepare_counts(Y)
  if (.miso_is_sparse(Y)) {
    list(blocks = .miso_sparse_blocks(Y, block_size), pass = .miso_sparse_pass)
  } else {
    list(blocks = list(Y = Y, features = .miso_feature_blocks(ncol(Y), block_size)),
         pass = .miso_dense_pass)
  }
}

.miso_dense_pass <- function(blocks, expected_log_lambda, expected_log_F,
                             omega = NULL, evaluate_bound = FALSE,
                             compute_allocated = TRUE, eps = 1e-12, log_F = NULL) {
  N = nrow(expected_log_lambda)
  D = ncol(expected_log_lambda)
  M = ncol(expected_log_F)
  allocated = if (compute_allocated) matrix(0, N, D) else NULL
  C = if (!is.null(omega) && is.null(log_F)) matrix(0, D, M) else NULL
  score_sum = if (!is.null(log_F)) matrix(0, D, nrow(log_F)) else NULL
  bound = if (evaluate_bound) numeric(N) else NULL
  for (block in blocks$features) {
    Y_block = blocks$Y[, block, drop = FALSE]
    if (evaluate_bound) {
      bound = bound + rowSums(Y_block * .miso_xi(
        expected_log_lambda, expected_log_F, block, log_normalizer = TRUE))
    }
    if (!compute_allocated && is.null(omega)) next
    xi = .miso_xi(expected_log_lambda, expected_log_F, block)
    counts = if (!is.null(omega)) matrix(0, D, length(block)) else NULL
    for (d in seq_len(D)) {
      xi_d = matrix(xi[, d, ], N, length(block))
      xi_y = xi_d * Y_block
      if (compute_allocated) allocated[, d] = allocated[, d] + rowSums(xi_y)
      if (!is.null(omega)) {
        counts[d, ] = colSums((omega * Y_block) * xi_d)
      }
    }
    if (!is.null(omega)) {
      if (is.null(log_F)) C[, block] = counts else {
        score_sum = score_sum + counts %*% t(log_F[, block, drop = FALSE])
      }
    }
  }
  result = list(allocated = allocated, C = C, bound = bound)
  if (!is.null(log_F)) result$score = score_sum
  result
}

## Exact sparse MiSo kernels. Counts are processed at nonzero entries;
## dictionary, posterior, and global sufficient-statistic arrays stay dense.

## One motif, one frozen local state. No N x D x M allocation tensor is made.
.miso_sparse_pass <- function(blocks, expected_log_lambda, expected_log_F,
                              omega = NULL, evaluate_bound = FALSE,
                              compute_allocated = TRUE, eps = 1e-12, log_F = NULL) {
  N = nrow(expected_log_lambda)
  D = ncol(expected_log_lambda)
  M = ncol(expected_log_F)
  allocated = if (compute_allocated) matrix(0, N, D) else NULL
  C = if (!is.null(omega) && is.null(log_F)) matrix(0, D, M) else NULL
  score_sum = if (!is.null(log_F)) matrix(0, D, nrow(log_F)) else NULL
  bound = if (evaluate_bound) numeric(N) else NULL
  ## Transpose the small dictionary once, instead of transposing a gathered
  ## D x nnz(block) matrix for every block.
  expected_log_F_by_feature = t(expected_log_F)

  for (block in blocks) {
    score = expected_log_lambda[block$i, , drop = FALSE] +
      expected_log_F_by_feature[block$j, , drop = FALSE]
    if (evaluate_bound) {
      contribution = block$x * .miso_log_sum_exp_rows(score)
      by_row = rowsum(matrix(contribution, ncol = 1), block$i, reorder = FALSE)
      index = as.integer(rownames(by_row))
      bound[index] = bound[index] + by_row[, 1]
    }
    if (!compute_allocated && is.null(omega)) next
    xi = .miso_softmax_rows(score)
    xi_y = xi * block$x
    if (compute_allocated) {
      by_row = rowsum(xi_y, block$i, reorder = FALSE)
      index = as.integer(rownames(by_row))
      allocated[index, ] = allocated[index, , drop = FALSE] + by_row
    }
    if (!is.null(omega)) {
      by_column = rowsum(xi_y * omega[block$i], block$j, reorder = FALSE)
      index = as.integer(rownames(by_column))
      if (is.null(log_F)) C[, index] = t(by_column) else {
        score_sum = score_sum + t(by_column) %*% t(log_F[, index, drop = FALSE])
      }
    }
  }
  result = list(allocated = allocated, C = C, bound = bound)
  if (!is.null(log_F)) result$score = score_sum
  result
}

