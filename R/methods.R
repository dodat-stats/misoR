#' Summarize a MiSo fit
#'
#' @param object A fitted `miso_fit` object.
#' @param ... Additional arguments; currently unused.
#' @return A `summary.miso_fit` list with dimensions, convergence diagnostics,
#'   a `motifs` table, and a `slots` table. `motifs` contains posterior mixture
#'   weights (`pi`), mean responsibilities (`allocation_mass`), hard-assignment
#'   sizes, and unordered MAP supports. `slots` contains the most probable factor,
#'   its selection probability, and the loading prior shape, rate, and mean.
#' @export
#' @examples
#' Y <- matrix(c(5, 1, 4, 2, 1, 5, 2, 4), 4, 2, byrow = TRUE)
#' fit <- miso_fit(Y, F = diag(2), D = 1, max_iters = 5)
#' summary(fit)
summary.miso_fit <- function(object, ...) {
  S = ncol(object$omega)
  D = dim(object$gamma)[2]
  K = nrow(object$F)
  probabilities = matrix(object$gamma, S * D, K)
  selected = max.col(probabilities, ties.method = "first")
  slots = data.frame(motif = rep(seq_len(S), D), slot = rep(seq_len(D), each = S),
    factor = selected, probability = probabilities[cbind(seq_len(S * D), selected)],
    prior_shape = as.vector(object$alpha0), prior_rate = as.vector(object$beta0),
    prior_mean = as.vector(object$prior_mean))
  slots = slots[order(slots$motif, slots$slot), ]
  rownames(slots) = NULL
  supports = vapply(seq_len(S), function(s)
    paste(sort(unique(slots$factor[slots$motif == s])), collapse = "+"), "")
  increments = diff(object$elbo)
  result = list(
    dimensions = c(N = nrow(object$omega), M = ncol(object$F), K = K, S = S, D = D),
    converged = object$converged, iterations = object$n_iter,
    final_elbo = object$final_elbo,
    min_elbo_increment = if (length(increments)) min(increments) else NA_real_,
    elbo_decreases = sum(increments < -1e-8 * (1 + abs(head(object$elbo, -1L)))),
    motifs = data.frame(motif = seq_len(S), pi = object$pi,
      allocation_mass = object$allocation_mass,
      n_assigned = tabulate(object$z_hat, nbins = S), support = supports),
    slots = slots)
  class(result) = "summary.miso_fit"
  result
}

#' @rdname summary.miso_fit
#' @param x A summary returned by `summary()`.
#' @param digits Number of significant digits to print.
#' @export
print.summary.miso_fit <- function(x, digits = 4, ...) {
  cat("MiSo fit:", paste(names(x$dimensions), x$dimensions, sep = "=", collapse = ", "), "\n")
  cat("Converged:", x$converged, "after", x$iterations, "iterations\n")
  cat("Final ELBO:", format(x$final_elbo, digits = max(digits, 8)), "\n")
  cat("ELBO decreases:", x$elbo_decreases, "\n\nMotifs:\n")
  print(x$motifs, digits = digits, row.names = FALSE)
  cat("\nFactor slots and loading priors:\n")
  print(x$slots, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Plot ELBO convergence or learned loadings
#'
#' Loading bars use posterior means averaged over mixture uncertainty.
#' Hard cluster labels determine grouping only. No truth labels are needed.
#' Empty hard-assignment clusters are omitted; fitted components are not merged.
#'
#' @param x A fitted `miso_fit` object.
#' @param type Either `"elbo"` or `"loadings"`.
#' @param normalize Show within-observation fractions rather than absolute
#'   expected loadings (counts for [miso_fit()], per unit of `n` for
#'   [miso_fit_length()]). Used for loading plots only.
#' @param cluster_order Optional permutation of occupied fitted cluster indices.
#' @param factor_order Optional permutation of all fitted factor indices.
#' @param sort_by Optional fitted factor index. Within each cluster, order bars
#'   by the loading fraction of that factor. Otherwise preserve observation order.
#' @param col Optional color (ELBO) or vector with one color per factor in the
#'   original factor order (loadings).
#' @param main Optional plot title.
#' @param ... Additional graphical arguments to `plot()` or `barplot()`.
#' @return Invisibly, for loading plots, a list containing the N by K posterior
#'   mean `loadings` in original order, displayed `values`, `observation_order`,
#'   `cluster_order`, `factor_order`, and `bar_midpoints`. For ELBO plots, returns
#'   the fitted object invisibly.
#' @export
#' @examples
#' Y <- matrix(c(5, 1, 4, 2, 1, 5, 2, 4), 4, 2, byrow = TRUE)
#' fit <- miso_fit(Y, F = diag(2), D = 1, max_iters = 5)
#' plot(fit)
#' plot(fit, type = "loadings")
plot.miso_fit <- function(x, type = c("elbo", "loadings"), normalize = TRUE,
                          cluster_order = NULL, factor_order = NULL,
                          sort_by = NULL, col = NULL, main = NULL, ...) {
  type = match.arg(type)
  if (type == "elbo") {
    if (is.null(main)) main = "MiSo ELBO"
    if (is.null(col)) col = "steelblue"
    graphics::plot(seq_along(x$elbo), x$elbo, type = "l", col = col,
                   xlab = "Outer iteration", ylab = "ELBO", main = main, ...)
    return(invisible(x))
  }
  .miso_stop(is.logical(normalize) && length(normalize) == 1 && !is.na(normalize),
             "normalize must be TRUE or FALSE.")
  L = .miso_expected_loadings(x)
  K = ncol(L)
  occupied = sort(unique(x$z_hat))
  if (is.null(cluster_order)) cluster_order = occupied
  if (is.null(factor_order)) factor_order = seq_len(K)
  .miso_stop(length(cluster_order) == length(occupied) &&
               setequal(cluster_order, occupied),
             "cluster_order must permute the occupied cluster indices.")
  .miso_stop(length(factor_order) == K && setequal(factor_order, seq_len(K)),
             "factor_order must permute all factor indices.")
  if (!is.null(sort_by)) {
    .miso_stop(length(sort_by) == 1 && sort_by %in% seq_len(K),
               "sort_by must be a fitted factor index.")
  }
  fractions = L / pmax(rowSums(L), .Machine$double.xmin)
  rows = unlist(lapply(cluster_order, function(s) {
    i = which(x$z_hat == s)
    if (is.null(sort_by)) i else i[order(fractions[i, sort_by], i)]
  }), use.names = FALSE)
  sizes = tabulate(match(x$z_hat, cluster_order), nbins = length(cluster_order))
  starts = c(1L, head(cumsum(sizes), -1L) + 1L)
  ends = cumsum(sizes)
  gaps = numeric(length(rows))
  gaps[starts[-1L]] = max(1, nrow(L) / 125)
  if (is.null(col)) col = grDevices::hcl.colors(K, "Dark 3")
  .miso_stop(length(col) == K, "col must contain one color per factor.")
  values = if (normalize) fractions else L
  displayed = values[rows, factor_order, drop = FALSE]
  if (is.null(main)) main = if (normalize) "Learned loading fractions" else "Learned loadings"
  labels = rownames(x$F)
  if (is.null(labels)) labels = paste("Factor", seq_len(K))
  old_mar = graphics::par(mar = pmax(graphics::par("mar"), c(4.5, 4, 2.5, 1)))
  on.exit(graphics::par(old_mar))
  loading_label = if (is.null(x[["n"]])) "Expected loading (counts)" else
    "Expected loading (per unit of n)"
  midpoints = as.vector(graphics::barplot(
    t(displayed), space = gaps, border = NA, col = col[factor_order],
    axes = FALSE, axisnames = FALSE, xaxs = "i", yaxs = "i",
    ylim = c(0, 1.25 * max(1e-12, rowSums(displayed))),
    main = main, ylab = if (normalize) "Loading fraction" else loading_label, ...))
  graphics::axis(2, at = if (normalize) seq(0, 1, 0.25) else NULL, las = 1)
  graphics::axis(1, at = (midpoints[starts] + midpoints[ends]) / 2,
    labels = paste0("C", cluster_order, "\nn=", sizes), tick = FALSE, line = -0.3, cex.axis = 0.8)
  graphics::legend("top", legend = labels[factor_order], fill = col[factor_order],
    border = NA, ncol = min(K, 6L), bty = "n", cex = 0.8)
  graphics::box(bty = "l")
  invisible(list(loadings = L, values = displayed, observation_order = rows,
    cluster_order = cluster_order, factor_order = factor_order, bar_midpoints = midpoints))
}
