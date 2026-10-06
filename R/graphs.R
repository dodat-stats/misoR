#' Summarize posterior factor co-loading graphs
#'
#' Compute motif graphs without changing the fit. Duplicate slots contribute
#' additively to factor loadings and co-loading edges; only vertex labels are
#' deduplicated. All slot probabilities enter the moment calculations.
#'
#' @param fit A [miso_fit()] or [miso_fit_length()] result, or a list containing
#'   compatible `gamma`, `a`, `b`, and `omega` arrays.
#' @details For motif s, let `G[d,k] = gamma[s,d,k]` and
#'   `m[i,d] = a[i,s,d] / b[s,d]` (or `/ b[i,s,d]` for length adjustment).
#'   Responsibility-weighted mean factor loadings are the weighted average of
#'   `m %*% G`. For distinct slots d,e, `Q[d,e]` is the weighted average of
#'   `m[i,d]*m[i,e]`, with `Q[d,d]=0`. Off-diagonal entries of `t(G) %*% Q %*% G`
#'   give expected co-loading. Same-slot products are excluded because one slot
#'   cannot select two different factors. These are moments, not covariances.
#'
#'   Display vertices are the union of the slots' MAP factors. Edge weights are
#'   zeroed outside those vertices and on the diagonal. Loading summaries keep
#'   all K factors, including probability mass outside the displayed vertices.
#'   Aggregate edges and loadings average motifs using mean responsibilities,
#'   not the Dirichlet posterior mean. A zero-mass motif has zero summaries.
#'   Length-adjusted fits retain per-unit-n loading units; edges have squared units.
#' @return A `miso_graphs` list with `vertices` (one factor-index vector per motif),
#'   `edge_weights` (K by K by S), `aggregate` (K by K), `mass` (length S),
#'   `node_loadings` (S by K), `aggregate_loadings` (length K), and `units`.
#' @seealso [plot_miso_graphs()]
#' @export
miso_fit_graphs <- function(fit) {
  dims <- dim(fit$gamma)
  .miso_stop(length(dims) == 3L && all(dims > 0), "gamma must be S by D by K.")
  S <- dims[1]; D <- dims[2]; K <- dims[3]
  N <- nrow(fit$omega)
  .miso_stop(length(N) == 1L && N > 0 && ncol(fit$omega) == S &&
    identical(dim(fit[["a"]]), c(N, S, D)) &&
    (identical(dim(fit[["b"]]), c(S, D)) || identical(dim(fit[["b"]]), c(N, S, D))),
    "Incompatible a, b, omega, or gamma dimensions.")
  for (field in c("gamma", "omega", "a", "b"))
    .miso_stop(is.numeric(fit[[field]]) && all(is.finite(fit[[field]])) &&
      all(fit[[field]] >= 0), paste(field, "must be finite and nonnegative."))
  .miso_stop(all(fit[["b"]] > 0), "b must be positive.")
  .miso_stop(max(abs(apply(fit$gamma, c(1, 2), sum) - 1)) < 1e-8 &&
    max(abs(rowSums(fit$omega) - 1)) < 1e-8, "gamma and omega must be normalized.")
  selected <- matrix(apply(fit$gamma, c(1, 2), which.max), S, D)
  vertices <- lapply(seq_len(S), function(s) sort(unique(selected[s, ])))
  totals <- colSums(fit$omega)
  edges <- array(0, c(K, K, S))
  loadings <- matrix(0, S, K)
  for (s in seq_len(S)) {
    if (totals[s] == 0) next
    G <- matrix(fit$gamma[s, , ], D, K)
    m <- matrix(fit[["a"]][, s, ], N, D) / .miso_rate_matrix(fit[["b"]], s, N, D)
    weights <- fit$omega[, s] / totals[s]
    loadings[s, ] <- colSums((m %*% G) * weights)
    Q <- crossprod(m, m * weights)
    diag(Q) <- 0
    E <- t(G) %*% Q %*% G
    diag(E) <- 0
    v <- vertices[[s]]
    edges[v, v, s] <- E[v, v]
  }
  mass <- totals / N
  structure(list(vertices = vertices, edge_weights = edges,
    aggregate = apply(sweep(edges, 3, mass, "*"), c(1, 2), sum), mass = mass,
    node_loadings = loadings, aggregate_loadings = colSums(loadings * mass),
    units = if (is.null(fit[["n"]])) "counts" else "loading per unit of n"), class = "miso_graphs")
}

# Vertex coordinates only: layout never changes the displayed edges or weights.
.miso_graph_layout <- function(weights, layout, root = NULL, seed = 1) {
  n <- nrow(weights)
  circle <- function(n) {
    if (n == 1) return(matrix(c(0, 0), 1, 2))
    if (n == 2) return(rbind(c(-1, 0), c(1, 0)))
    angle <- base::pi / 2 + 2 * base::pi * (seq_len(n) - 1) / n
    cbind(cos(angle), sin(angle))
  }
  xy <- circle(n)
  if (layout == "circle" || n <= 2 || !any(weights > 0)) return(xy)
  if (is.null(root)) root <- which.max(rowSums(weights))
  if (layout == "hub") {
    xy[root, ] <- 0
    xy[-root, ] <- circle(n - 1)
    return(xy)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop('Install igraph for this layout: install.packages("igraph")')
  }
  # Keep stochastic layouts reproducible without affecting simulation/fitting.
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv)
    else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(seed)
  g <- igraph::graph_from_adjacency_matrix(weights, mode = "upper",
                                          weighted = TRUE, diag = FALSE)
  strength <- igraph::E(g)$weight / mean(igraph::E(g)$weight)
  xy <- switch(layout,
    fr = igraph::layout_with_fr(g, weights = strength, niter = 1000),
    kk = igraph::layout_with_kk(g, weights = 1 / strength),
    tree = igraph::layout_as_tree(g, root = root, mode = "all"))
  # Fit the existing plotting window, including disconnected/flat layouts.
  xy <- sweep(xy, 2, (apply(xy, 2, min) + apply(xy, 2, max)) / 2)
  span <- apply(abs(xy), 2, max)
  sweep(xy, 2, ifelse(span > 0, span, 1), "/")
}

#' Plot MiSO motif subgraphs or the aggregate factor graph
#'
#' @param graphs A MiSO fit (recommended), a [miso_fit_graphs()] result, or a
#'   legacy graph list with `vertices`, `edge_weights`, `aggregate`, and `mass`.
#' @param type Plot each motif (`"subgraphs"`) or their weighted `"aggregate"`.
#' @param colors,labels Vectors of K colors and factor labels in fitted factor order.
#' @param edge_max Common denominator for edge widths. Defaults to the maximum
#'   across motif and aggregate edges, so separate calls share a scale.
#' @param layout Vertex layout: `"circle"`, `"hub"`, `"fr"`, `"kk"`, or `"tree"`.
#'   The last three require the optional igraph package. Layout never changes edges.
#' @param ncol Maximum number of subgraph columns.
#' @param main Optional overall title. For subgraphs, a vector of S titles labels
#'   the individual panels; with `arrange = FALSE`, use NULL or S panel titles.
#' @param root Optional factor index for hub or tree layouts.
#' @param seed Random seed for layouts; the caller's random state is preserved.
#' @param node_size `"constant"` for equal circles (default), or `"loading"` for
#'   areas proportional to posterior mean factor loadings, summed over all slots.
#'   Loading areas use a common scale across all motifs and the aggregate.
#' @param arrange Set FALSE to draw into the caller's existing base-graphics
#'   layout, consuming one panel per motif or one aggregate panel. In this mode
#'   margins and layout are controlled by the caller and are not reset afterward.
#' @param show_legend Show the co-loading annotation or edge-width legend.
#' @return Invisibly, the graph summary used for plotting. Graphics parameters
#'   are restored on exit when `arrange = TRUE`. No fit parameters are changed.
#' @details Duplicate slots add their expected loadings and their contributions
#'   to edges with other factors. They do not create duplicate vertices or self-edges.
#'   Factor uncertainty is preserved; slots are not collapsed by their MAP labels.
#'   The percentage above each motif is its mean posterior responsibility.
#' @seealso [miso_fit_graphs()], [plot.miso_fit()]
#' @export
#' @examples
#' Y <- matrix(c(5, 1, 4, 2, 1, 5, 2, 4), 4, 2, byrow = TRUE)
#' fit <- miso_fit(Y, F = diag(2), D = 2, S = 2, max_iters = 5, tol = 0)
#' plot_miso_graphs(fit, "subgraphs", node_size = "loading")
#' plot(fit, type = "aggregate")
plot_miso_graphs <- function(graphs, type = c("subgraphs", "aggregate"),
                             colors = NULL, labels = NULL, edge_max = NULL,
                             layout = c("circle", "hub", "fr", "kk", "tree"), ncol = 4,
                             main = NULL, root = NULL, seed = 1,
                             node_size = c("constant", "loading"),
                             arrange = TRUE, show_legend = TRUE) {
  if (!is.null(graphs$gamma)) graphs <- miso_fit_graphs(graphs)
  node_size <- match.arg(node_size)
  type <- match.arg(type)
  layout <- match.arg(layout)
  K <- nrow(graphs$aggregate)
  S <- length(graphs$vertices)
  stopifnot(is.logical(arrange), length(arrange) == 1, !is.na(arrange),
            is.logical(show_legend), length(show_legend) == 1, !is.na(show_legend))
  if (!is.null(main)) stopifnot(is.character(main), !anyNA(main),
    length(main) %in% if (type == "subgraphs" && arrange) c(1L, S) else
      if (type == "subgraphs") S else 1L)
  if (is.null(colors)) colors <- grDevices::hcl.colors(K, "Dark 3")
  if (is.null(labels)) labels <- paste0("F", seq_len(K))
  if (is.null(edge_max)) edge_max <- max(graphs$edge_weights, graphs$aggregate)
  units <- if (is.null(graphs$units)) "counts" else graphs$units
  edge_units <- if (units == "counts") "counts squared" else paste0("(", units, ") squared")
  if (node_size == "loading") {
    .miso_stop(!is.null(graphs$node_loadings) && !is.null(graphs$aggregate_loadings),
               "Loading-sized nodes require a fit or a miso_fit_graphs() summary.")
    node_max <- max(graphs$node_loadings, graphs$aggregate_loadings, 1e-12)
  }
  stopifnot(length(colors) == K, length(labels) == K,
            length(edge_max) == 1, is.finite(edge_max), edge_max >= 0,
            length(ncol) == 1, ncol >= 1, ncol == as.integer(ncol),
            length(seed) == 1, is.finite(seed))
  if (!is.null(root)) {
    stopifnot(length(root) == 1, is.finite(root), root == as.integer(root),
              root >= 1, root <= K)
  }
  draw <- function(weights, v, xy, title, loading = NULL) {
    graphics::plot(NA, xlim = c(-1.35, 1.35), ylim = c(-1.6, 1.35),
         asp = 1, axes = FALSE, xlab = "", ylab = "", main = title)
    pairs <- which(upper.tri(weights) & weights > 0, arr.ind = TRUE)
    for (p in seq_len(nrow(pairs))) {
      j <- pairs[p, 1]; k <- pairs[p, 2]
      graphics::segments(xy[j, 1], xy[j, 2], xy[k, 1], xy[k, 2],
               lwd = 10 * weights[j, k] / max(edge_max, 1e-12), col = "grey45")
    }
    radius <- if (node_size == "constant") rep(.16, length(v)) else
      .22 * sqrt(loading / node_max)
    graphics::symbols(xy[, 1], xy[, 2], circles = radius,
            inches = FALSE, add = TRUE, bg = colors[v], fg = "white")
    graphics::text(xy[, 1], xy[, 2], labels = labels[v], col = "white", font = 2)
  }
  if (arrange) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
  }
  if (type == "subgraphs") {
    columns <- min(ncol, S)
    panel_titles <- !is.null(main) && length(main) == S
    if (arrange) graphics::par(mfrow = c(ceiling(S / columns), columns),
        mar = c(0.5, 0.5, 2, 0.5),
        oma = c(if (show_legend) 1 else 0, 0, if (is.null(main) || panel_titles) 0 else 2, 0))
    for (s in seq_len(S)) {
      v <- graphs$vertices[[s]]
      weights <- matrix(graphs$edge_weights[v, v, s], length(v), length(v))
      local_root <- if (!is.null(root) && root %in% v) match(root, v) else NULL
      xy <- .miso_graph_layout(weights, layout, local_root, seed)
      draw(weights, v, xy,
           if (panel_titles) main[s] else sprintf("C%d (%.1f%%)", s, 100 * graphs$mass[s]),
           if (node_size == "loading") graphs$node_loadings[s, v] else NULL)
    }
    if (show_legend) graphics::mtext(paste0("Edges: co-loading (", edge_units, ")",
                 if (node_size == "loading") "; node area: mean loading" else ""),
          side = 1, outer = TRUE, cex = 0.8)
    if (!is.null(main) && !panel_titles) graphics::mtext(main, side = 3, outer = TRUE, font = 2)
  } else {
    if (arrange) graphics::par(mfrow = c(1, 1), mar = c(1, 1, 3, 1), oma = c(0, 0, 0, 0))
    xy <- .miso_graph_layout(graphs$aggregate, layout, root, seed)
    if (is.null(main)) main <- "Aggregated factor graph"
    draw(graphs$aggregate, seq_len(K), xy, main, graphs$aggregate_loadings)
    if (show_legend && edge_max > 0) {
      graphics::legend("bottom", legend = signif(edge_max * c(0.25, 0.5, 1), 3),
             title = paste0("Co-loading (", edge_units, ")"), col = "grey45",
             lwd = c(2.5, 5, 10), ncol = 3, bty = "n", cex = 0.85)
    }
  }
  invisible(graphs)
}
