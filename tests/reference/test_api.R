## Public initialization, fitting, validation, and prediction contracts.
local({
  close = function(x, y, tolerance = 1e-8) {
    comparison = all.equal(x, y, tolerance = tolerance)
    if (!isTRUE(comparison)) stop(paste(comparison, collapse = "\n"))
  }
  rejects = function(expression) {
    stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
  }
  set.seed(121)
  Y = matrix(rpois(12 * 15, 1), 12, 15); Y[1, ] = 0
  F = .miso_normalize_rows(matrix(rexp(4 * 15), 4, 15))
  inputs = list(Y)
  if (requireNamespace("Matrix", quietly = TRUE)) inputs[[2]] = Matrix::Matrix(Y, sparse = TRUE)
  for (counts in inputs) {
    initial = .miso_initialize(counts, F, D = 2, susie_max_iters = 5, seed = 23)
    automatic = miso_fit(counts, F, D = 2, seed = 23,
                         init_control = list(susie_max_iters = 5),
                         keep_initialization = TRUE, max_iters = 4, tol = 0,
                         update_F = TRUE, mf_iters = 2)
    close(automatic$initialization[names(initial)], initial)
    automatic$initialization = NULL
    supplied = miso_fit(counts, init = initial, max_iters = 4, tol = 0,
                        update_F = TRUE, mf_iters = 2)
    supplied$call = automatic$call = NULL
    close(automatic, supplied)
    close(predict(supplied, "loadings"), .miso_expected_loadings(supplied))
    close(predict(supplied), predict(supplied, "loadings") %*% supplied$F)
    explicit = miso_fit(counts, init = initial, max_iters = 4, tol = 0,
                        update_prior = TRUE, update_F = TRUE, mf_iters = 3)
    defaults = miso_fit(counts, init = initial, max_iters = 4, tol = 0)
    explicit$call = defaults$call = NULL
    close(defaults, explicit)
    stopifnot(identical(defaults$final_elbo, tail(defaults$elbo, 1)))
    learned = miso_fit(counts, K = 3, D = 2, seed = 23,
                       init_control = list(nmf_max_iters = 5, susie_max_iters = 5),
                       keep_initialization = TRUE, max_iters = 1)
    close(learned$initialization$F,
          .poisson_nmf_fit(counts, 3, max_iters = 5, seed = 23)$F)

  }
  single = poisson_susie_fit(Y[2, ], F, 2, max_iters = 3)
  close(single, poisson_susie_fit(Y[2, , drop = FALSE], F, 2, max_iters = 3))
  stopifnot(identical(dim(single$alpha), c(1L, 2L)),
            identical(dim(single$gamma), c(1L, 2L, 4L)), is.null(single$fits))
  for (input in list(matrix(0, 1, 1), matrix(3, 1, 1))) {
    fit = miso_fit(input, matrix(1, 1, 1), 1,
                   init_control = list(susie_max_iters = 2), max_iters = 2)
    stopifnot(is.finite(fit$final_elbo), identical(dim(fit$alpha), c(1L, 1L, 1L)))
  }
  rejects(miso_fit(Y, F)); rejects(miso_fit(Y, D = 2))
  rejects(miso_fit(Y, F, D = 2, max_iters = 0))
  for (bad in list(0, -1, 1.5, NA_real_, Inf, numeric())) {
    rejects(miso_fit(Y, init = initial, mf_iters = bad))
  }
  rejects(miso_fit(Y, init = initial, n_inner = 3))
  rejects(miso_fit(Y, init = initial, gamma_step = 0.5))
  rejects(miso_fit(Y, init = initial, F_step = 0.5))
  rejects(miso_fit(Y, F, D = 2, init_control = list(typo = 1)))
  rejects(miso_fit(Y, F, D = 2, init_control = list(nmf_max_iters = 3)))
  rejects(miso_fit(Y, init = list())); rejects(miso_fit(Y, init = initial, D = 3))
  rejects(poisson_susie_fit(Y, F, D = 0)); rejects(poisson_susie_fit(Y, F[, 1:2], D = 2))
  stopifnot(setequal(getNamespaceExports("misoR"),
                     c("miso_fit", "miso_init", "poisson_susie_fit")))
})
cat("Public API tests passed.\n")
