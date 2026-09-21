library(misoR)

local({
  close <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-8)))
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

  # One factor gives an analytic Gamma update even for fractional observations.
  # Fractional row totals detect any accidental integer conversion or rounding.
  Y <- rbind(c(0.25, 1.1, 0), c(2.7, 0, 0.35), c(0, 0, 0))
  F <- matrix(c(.2, .3, .5), 1)
  a0 <- 2.4; b0 <- 3.1
  initial <- list(F = F, gamma = array(1, c(1, 1, 1)),
                  alpha0 = matrix(a0), beta0 = matrix(b0))
  for (adjusted in c(FALSE, TRUE)) for (sparse in c(FALSE, TRUE)) {
    n <- if (adjusted) c(1.3, 4.2, 2.1) else rep(1, nrow(Y))
    input <- if (sparse) Matrix::Matrix(Y, sparse = TRUE) else Y
    args <- list(Y = input, init = initial, max_iters = 1,
                 update_prior = FALSE, update_F = FALSE)
    if (adjusted) args$n <- n
    fit <- do.call(if (adjusted) miso_fit_length else miso_fit, args)
    a <- a0 + rowSums(Y); b <- b0 + n
    close(as.vector(fit$alpha), a)
    close(as.vector(predict(fit, "loadings")), a / b)
    close(predict(fit), outer(n * a / b, as.vector(F)))
    constant <- sum(rowSums(Y) * log(n)) - sum(lgamma(Y + 1))
    reference <- sum(Y * log(outer(n, as.vector(F))) - lgamma(Y + 1)) +
      sum(a0 * log(b0) - lgamma(a0) + lgamma(a) - a * log(b))
    close(fit$final_elbo + constant, reference)
    if (adjusted) close(fit$elbo_constant, constant)
  }

  # Exercise fractional data through NMF, row-wise SuSiE initialization, and
  # every outer update, for both public fitters and both storage formats.
  Y <- rbind(c(.25, 1.1, 0, 2.3), c(2.7, 0, .35, 1.25),
             c(.9, .4, 1.8, 0), c(0, 2.2, .7, 1.1))
  fields <- c("alpha", "beta", "alpha0", "beta0", "gamma", "F", "omega", "elbo")
  for (fitter in list(miso_fit, miso_fit_length)) {
    args <- list(K = 2, D = 2, max_iters = 12, tol = 0, seed = 4,
                 init_control = list(nmf_max_iters = 5, susie_max_iters = 5))
    dense <- do.call(fitter, c(list(Y = Y), args))
    sparse <- do.call(fitter, c(list(Y = Matrix::Matrix(Y, sparse = TRUE),
                                    block_size = 1), args))
    close(dense[fields], sparse[fields])
    stopifnot(all(is.finite(predict(dense))),
              all(diff(dense$elbo) >= -1e-8 * (1 + abs(head(dense$elbo, -1)))))
    for (value in c(-.1, NA_real_, Inf)) {
      invalid <- Y; invalid[1, 1] <- value
      rejects(do.call(fitter, c(list(Y = invalid), args)))
      rejects(do.call(fitter, c(list(Y = Matrix::Matrix(invalid, sparse = TRUE)), args)))
    }
  }
})
cat("Nonnegative real-valued input tests passed: analytic updates, objective, initialization, and sparse/dense fits.\n")
