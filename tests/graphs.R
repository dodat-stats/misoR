library(misoR)
local({
  close <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10)))
  rejects <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

  # Two slots selecting factor 1 must both contribute to its loading and edges.
  f <- list(gamma = array(0, c(1, 3, 2)), omega = matrix(1),
            a = array(c(2, 3, 7), c(1, 1, 3)), b = matrix(1, 1, 3))
  f$gamma[1, 1:2, 1] <- 1
  f$gamma[1, 3, 2] <- 1
  g <- miso_fit_graphs(f)
  close(g$node_loadings, matrix(c(5, 7), 1))
  close(g$aggregate_loadings, c(5, 7))
  close(g$aggregate, matrix(c(0, 35, 35, 0), 2))
  stopifnot(identical(g$vertices[[1]], 1:2))
  # The known-selection graph is invariant to combining the duplicate slots.
  combined <- list(gamma = array(diag(2), c(1, 2, 2)), omega = matrix(1),
                   a = array(c(5, 7), c(1, 1, 2)), b = matrix(1, 1, 2))
  close(miso_fit_graphs(combined), g)

  # Independent scalar reference: uncertainty, soft membership, and per-row rates.
  set.seed(54)
  N <- 4L; S <- 2L; D <- 3L; K <- 4L
  u <- list(gamma = array(runif(S*D*K), c(S, D, K)),
            omega = matrix(runif(N*S), N, S),
            a = array(runif(N*S*D, 1, 4), c(N, S, D)),
            b = array(runif(N*S*D, .5, 2), c(N, S, D)), n = 1:N)
  u$gamma <- sweep(u$gamma, c(1, 2), apply(u$gamma, c(1, 2), sum), "/")
  u$omega <- u$omega / rowSums(u$omega)
  ug <- miso_fit_graphs(u)
  E <- array(0, c(K, K, S)); A <- matrix(0, S, K)
  for (s in 1:S) for (i in 1:N) for (d in 1:D) for (k in 1:K) {
    w <- u$omega[i, s] / sum(u$omega[, s])
    md <- u$a[i, s, d] / u$b[i, s, d]
    A[s, k] <- A[s, k] + w*md*u$gamma[s, d, k]
    for (e in 1:D) for (l in 1:K) {
      if (d == e || k == l || !all(c(k, l) %in% ug$vertices[[s]])) next
      E[k, l, s] <- E[k, l, s] + w*md*u$a[i, s, e]/u$b[i, s, e] *
        u$gamma[s, d, k]*u$gamma[s, e, l]
    }
  }
  close(ug$node_loadings, A)
  close(ug$edge_weights, E)
  close(ug$aggregate, E[, , 1]*mean(u$omega[, 1]) + E[, , 2]*mean(u$omega[, 2]))
  close(ug$aggregate_loadings, colSums(A*colMeans(u$omega)))
  stopifnot(ug$units == "loading per unit of n")
  # An uncertain single slot cannot yield a between-factor edge.
  one <- list(gamma = array(c(.6, .4), c(1, 1, 2)), omega = matrix(1),
              a = array(5, c(1, 1, 1)), b = matrix(1))
  close(miso_fit_graphs(one)$aggregate, matrix(0, 2, 2))
  # Empty motifs and singleton shapes are valid.
  u$omega[,] <- 0; u$omega[, 1] <- 1
  zero <- miso_fit_graphs(u)
  stopifnot(all(zero$edge_weights[, , 2] == 0), all(zero$node_loadings[2, ] == 0))
  bad <- f; bad$b[1] <- 0; rejects(miso_fit_graphs(bad))
  bad <- f; bad$gamma[] <- NA; rejects(miso_fit_graphs(bad))
  bad <- f; bad$omega[] <- .5; rejects(miso_fit_graphs(bad))

  Y <- matrix(c(5, 1, 4, 2, 1, 5, 2, 4), 4, 2, byrow = TRUE)
  fit <- miso_fit(tol = 0, S = 2L, warm_up_iters = 0L, Y, F = diag(2), D = 2, max_iters = 5)
  length_fit <- miso_fit_length(tol = 0, S = 2L, warm_up_iters = 0L, Y, F = diag(2), D = 2, max_iters = 5)
  stopifnot(miso_fit_graphs(fit)$units == "counts",
            miso_fit_graphs(length_fit)$units == "loading per unit of n")
  for (x in list(fit, length_fit))
    close(miso_fit_graphs(x)$aggregate_loadings, colMeans(predict(x, "loadings")))
  single <- miso_fit(tol = 0, S = 1L, warm_up_iters = 0L, matrix(0, 1, 1), F = matrix(1), D = 1, max_iters = 2,
    population_prior=list(alpha=c(shape=2,rate=2),beta=c(shape=2,rate=2)))
  before <- serialize(fit, NULL)
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = 9, height = 4)
  on.exit({ grDevices::dev.off(); unlink(file) })
  old <- graphics::par(c("mfrow", "mar", "oma"))
  rng <- .Random.seed
  layouts <- c("circle", "hub")
  if (requireNamespace("igraph", quietly = TRUE)) layouts <- c(layouts, "fr", "kk", "tree")
  for (layout in layouts) for (type in c("subgraphs", "aggregate")) {
    for (x in list(fit, length_fit, single)) for (size in c("constant", "loading")) {
      expected <- miso_fit_graphs(x)
      close(plot_miso_graphs(x, type, layout = layout, node_size = size), expected)
      close(plot(x, type = type, layout = layout, node_size = size), expected)
      close(plot_miso_graphs(expected, type, layout = layout, node_size = size), expected)
      close(graphics::par(c("mfrow", "mar", "oma")), old)
    }
  }
  triangle <- list(gamma = array(diag(3), c(1, 3, 3)), omega = matrix(1),
    a = array(c(2, 3, 7), c(1, 1, 3)), b = matrix(1, 1, 3))
  for (layout in layouts) for (type in c("subgraphs", "aggregate"))
    close(plot_miso_graphs(triangle, type, layout = layout), miso_fit_graphs(triangle))
  close(plot(fit, type = "subgraphs", colors = c("red", "blue")), miso_fit_graphs(fit))
  # Composition must advance through the caller's layout without resetting it.
  graphics::layout(matrix(1:4, 2, 2, byrow = TRUE))
  graphics::par(mar = c(1, 1, 2, 1))
  titles <- c("First motif", "Second motif")
  close(plot_miso_graphs(ug, main = titles, arrange = FALSE, show_legend = FALSE), ug)
  stopifnot(identical(unname(graphics::par("mfg")), c(1L, 2L, 2L, 2L)))
  plot_miso_graphs(ug, "aggregate", arrange = FALSE, show_legend = FALSE)
  stopifnot(identical(unname(graphics::par("mfg")), c(2L, 1L, 2L, 2L)))
  graphics::plot.new()
  stopifnot(identical(unname(graphics::par("mfg")), c(2L, 2L, 2L, 2L)))
  graphics::par(mfrow = c(1, 1))
  plot_miso_graphs(zero, node_size = "loading")
  legacy <- g[c("vertices", "edge_weights", "aggregate", "mass")]
  close(plot_miso_graphs(legacy), legacy)
  rejects(plot_miso_graphs(legacy, node_size = "loading"))
  close(.Random.seed, rng)
  stopifnot(identical(serialize(fit, NULL), before))
  # No random state is introduced if the caller had none.
  rm(".Random.seed", envir = .GlobalEnv)
  plot_miso_graphs(triangle, layout = tail(layouts, 1))
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  assign(".Random.seed", rng, envir = .GlobalEnv)
})
cat("Graph moments, duplicate slots, length adjustment, and plotting API checks passed.\n")
