## Run from miso-main/: Rscript --vanilla tests/run.R
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Sparse tests require the Matrix package.")
}

local({
  close = function(x, y, tolerance = 1e-8) {
    comparison = all.equal(x, y, tolerance = tolerance, check.attributes = TRUE)
    if (!isTRUE(comparison)) stop(paste(comparison, collapse = "\n"))
  }
  rejects = function(expression) {
    stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
  }
  compare_pipeline = function(Y) {
    N = nrow(Y)
    M = ncol(Y)
    K = 3L
    S = 2L
    D = 2L
    set.seed(46)
    F = .miso_normalize_rows(matrix(rexp(K * M), K, M))
    gamma = .miso_normalize_gamma(array(rexp(S * D * K), c(S, D, K)))
    alpha0 = matrix(c(1, 3, 2, 1), S, D)
    beta0 = matrix(c(2, 1, 1, 4), S, D)
    omega = .miso_normalize_rows(matrix(rexp(N * S), N, S))
    sparse = methods::as(Matrix::Matrix(Y, sparse = TRUE), "generalMatrix")
    sparse = methods::as(sparse, "CsparseMatrix")
    original = sparse

    for (block_size in unique(c(1L, 4L, M + 2L))) {
      dense = .miso_fit(Y, F, gamma, alpha0, beta0, omega, max_iters = 2,
                        tol = 0, block_size = block_size)
      sparse_fit = .miso_fit(sparse, F, gamma, alpha0, beta0, omega, max_iters = 2,
                             tol = 0, block_size = block_size)
      dense$call = sparse_fit$call = NULL
      close(sparse_fit, dense)
    }

    dense_fit = .miso_fit(Y, F, gamma, alpha0, beta0, omega, max_iters = 5,
                         tol = 0, update_F = TRUE)
    sparse_fit = .miso_fit(sparse, F, gamma, alpha0, beta0, omega, max_iters = 5,
                          tol = 0, update_F = TRUE)
    dense_fit$call = sparse_fit$call = NULL
    close(sparse_fit, dense_fit, tolerance = 1e-7)
    close(predict(sparse_fit), predict(dense_fit), tolerance = 1e-7)
    stopifnot(identical(sparse, original))

    ## The rate is summed without clipping, even when eps is deliberately large.
    for (eps in c(1e-10, 0.02)) {
      dense_nmf = .poisson_nmf_fit(Y, K, max_iters = 5, tol = 0, eps = eps)
      sparse_nmf = .poisson_nmf_fit(sparse, K, max_iters = 5, tol = 0, eps = eps)
      close(sparse_nmf, dense_nmf)
      mean_Y = sparse_nmf$L %*% sparse_nmf$F
      close(tail(sparse_nmf$log_likelihood, 1),
            sum(Y * log(pmax(mean_Y, eps)) - mean_Y - lgamma(Y + 1)))
    }
    dense_ps = poisson_susie_fit(Y, F, D, max_iters = 5, tol = 0, keep_fits = TRUE)
    sparse_ps = poisson_susie_fit(sparse, F, D, max_iters = 5, tol = 0,
                                   keep_fits = TRUE)
    close(sparse_ps, dense_ps)
    close(.miso_initialize_from_susie(sparse_ps, D),
          .miso_initialize_from_susie(dense_ps, D))
    dense_auto = miso_fit(Y, K = K, D = D, max_iters = 2,
                           init_control = list(nmf_max_iters = 5, susie_max_iters = 5))
    sparse_auto = miso_fit(sparse, K = K, D = D, max_iters = 2,
                            init_control = list(nmf_max_iters = 5, susie_max_iters = 5))
    dense_auto$call = sparse_auto$call = NULL
    close(sparse_auto, dense_auto)

    ## A sparse single-row matrix and vector use the same public single-row API.
    row = Y[1, ]
    vector = Matrix::sparseVector(x = row[row > 0], i = which(row > 0), length = M)
    close(poisson_susie_fit(vector, F, D, max_iters = 5, seed = 4),
          poisson_susie_fit(row, F, D, max_iters = 5, seed = 4))
    close(poisson_susie_fit(sparse[1, , drop = FALSE], F, D, max_iters = 5, seed = 4),
          poisson_susie_fit(row, F, D, max_iters = 5, seed = 4))
    stopifnot(is.null(poisson_susie_fit(vector, F, D, max_iters = 2,
                                       keep_fits = FALSE)$fits))
    invisible(NULL)
  }

  set.seed(92)
  Y = matrix(rpois(9 * 13, 0.6), 9, 13)
  Y[1, ] = 0
  Y[, 13] = 0
  Y[2, 2] = 6
  compare_pipeline(Y)
  compare_pipeline(matrix(c(0, 3, 0, 2), nrow = 1))
  compare_pipeline(matrix(c(0, 3, 0, 2), ncol = 1))
  compare_pipeline(matrix(4, 1, 1))
  compare_pipeline(matrix(0, 1, 1))
  compare_pipeline(matrix(0, 4, 5))

  ## Canonicalization supports CSR, triplets with duplicate indices and stored
  ## zeros, symmetric storage, logical/pattern matrices, and unit diagonals.
  csc = Matrix::sparseMatrix(i = c(1, 3), j = c(2, 4), x = c(2, 5), dims = c(4, 4))
  triplet = methods::new("dgTMatrix", i = c(0L, 0L, 2L), j = c(1L, 1L, 3L),
                         x = c(1, 2, 0), Dim = c(4L, 4L))
  symmetric = Matrix::sparseMatrix(i = c(1, 2), j = c(1, 3), x = c(2, 4),
                                   symmetric = TRUE, dims = c(4, 4))
  pattern = Matrix::sparseMatrix(i = c(1, 3), j = c(2, 4), dims = c(4, 4))
  triangular = methods::new("dtCMatrix", Dim = c(3L, 3L), p = integer(4), diag = "U")
  for (input in list(csc, methods::as(csc, "RsparseMatrix"), triplet,
                     symmetric, pattern, csc > 0, triangular)) {
    canonical = .miso_prepare_counts(input)
    stopifnot(inherits(canonical, "dgCMatrix"), all(canonical@x > 0))
    close(as.matrix(canonical), as.matrix(input) * 1)
    F = matrix(1 / ncol(input), 1, ncol(input))
    gamma = array(1, c(1, 1, 1))
    prior = matrix(1, 1, 1)
    omega = matrix(1, nrow(input), 1)
    sparse_fit = .miso_fit(input, F, gamma, prior, prior, omega, max_iters = 2)
    dense_fit = .miso_fit(as.matrix(input) * 1, F, gamma, prior, prior, omega, max_iters = 2)
    sparse_fit$call = dense_fit$call = NULL
    close(sparse_fit, dense_fit)
  }
  for (bad in c(-1, NA_real_, Inf)) {
    invalid = Matrix::sparseMatrix(i = 1, j = 1, x = bad, dims = c(2, 3))
    rejects(.poisson_nmf_fit(invalid, 1, max_iters = 1))
    rejects(poisson_susie_fit(invalid, matrix(1/3, 1, 3), 1, max_iters = 1))
    rejects(.miso_fit(invalid, matrix(1/3, 1, 3), array(1, c(1,1,1)),
                     matrix(1), matrix(1), max_iters = 1))
  }

  ## Allocation guard: this input would occupy 48 MB if made dense. Sparse
  ## local inference and NMF should never allocate an object near that size.
  if (isTRUE(capabilities("profmem"))) {
    large = Matrix::sparseMatrix(i = seq_len(2000), j = seq_len(2000),
                                 x = 1, dims = c(2000, 3000))
    F = matrix(1 / 3000, 2, 3000)
    gamma = array(0.5, c(2, 2, 2))
    prior = matrix(1, 2, 2)
    omega = matrix(0.5, 2000, 2)
    path = tempfile()
    on.exit({ Rprofmem(NULL); unlink(path) })
    Rprofmem(path)
    fit = .miso_fit(large, F, gamma, prior, prior, omega,
                    max_iters = 1, block_size = 100)
    nmf = .poisson_nmf_fit(large, 2, max_iters = 1)
    ps = poisson_susie_fit(large[1:3, , drop = FALSE], F, 2, max_iters = 1)
    Rprofmem(NULL)
    allocations = suppressWarnings(as.numeric(sub(" .*", "", readLines(path))))
    stopifnot(max(allocations, na.rm = TRUE) < 2000 * 3000 * 8 / 4)
  }
})
cat("Sparse tests passed: dense equivalence, edge cases, formats, and allocation guard.\n")
