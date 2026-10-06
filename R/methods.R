#' Summarize a MiSO fit
#' @param object A miso_fit object in format version 2. Convert legacy objects
#'   explicitly with [misoR::miso_convert_fit()].
#' @param ... Reserved.
#' @return Dimensions, convergence, motif and slot summaries, and worst quadrature
#'   discrepancies. alpha_mean/beta_mean are posterior means; their SDs come from
#'   the joint population block. population_mean is E(alpha/beta). Converted legacy
#'   fits have NA SDs and no population posterior until refitted.
#' @export
summary.miso_fit <- function(object,...) {
  .miso_stop(identical(object$format_version,2L),"Convert legacy fits with miso_convert_fit().")
  S <- ncol(object$omega); D <- dim(object$gamma)[2]; K <- nrow(object$F)
  probabilities <- matrix(object$gamma,S*D,K)
  selected <- max.col(probabilities,ties.method="first")
  sd_field <- function(field) if(is.null(object$q_population)) rep(NA_real_,S*D) else
    sqrt(pmax(0,vapply(object$q_population,function(h) h[[field]],0.0)))
  slots <- data.frame(motif=rep(seq_len(S),D),slot=rep(seq_len(D),each=S),
    factor=selected,probability=probabilities[cbind(seq_len(S*D),selected)],
    alpha_mean=as.vector(object$alpha_mean),alpha_sd=sd_field("var_alpha"),
    beta_mean=as.vector(object$beta_mean),beta_sd=sd_field("var_beta"),
    population_mean=as.vector(object$population_mean))
  slots <- slots[order(slots$motif,slots$slot),]; rownames(slots) <- NULL
  supports <- vapply(seq_len(S),function(s)
    paste(sort(unique(slots$factor[slots$motif==s])),collapse="+"),"")
  increments <- diff(object$elbo); q <- object$quadrature
  result <- list(dimensions=c(N=nrow(object$omega),M=ncol(object$F),K=K,S=S,D=D),
    inference_model=object$inference_model,converged=object$converged,iterations=object$n_iter,
    warm_up_iters=object$warm_up_iters,main_iterations=object$n_iter_main,stop_reason=object$stop_reason,
    final_elbo=object$final_elbo,elbo_convention=object$elbo_convention,
    min_elbo_increment=if(length(increments)) min(increments) else NA_real_,
    elbo_decreases=sum(increments < -1e-8*(1+abs(head(object$elbo,-1L)-if(is.null(object$elbo_constant)) 0 else object$elbo_constant))),
    motifs=data.frame(motif=seq_len(S),pi=object$pi,allocation_mass=object$allocation_mass,
      n_assigned=tabulate(object$z_hat,nbins=S),support=supports),slots=slots,
    quadrature=if(is.null(q)) NULL else list(max_order=max(q$order),
      max_refinement_error=max(q$refinement_error),max_tail_check_error=max(q$tail_check_error),
      problematic_blocks=sum(!q$converged)))
  class(result) <- "summary.miso_fit"; result
}
#' @rdname summary.miso_fit
#' @param x A summary object.
#' @param digits Display precision.
#' @export
print.summary.miso_fit <- function(x,digits=4,...) {
  cat("MiSO:",paste(names(x$dimensions),x$dimensions,sep="=",collapse=", "),"\n")
  cat("Inference:",x$inference_model,"; converged:",x$converged,"after",x$iterations,"iterations\n")
  if(!is.null(x$warm_up_iters)) cat("Stages:",x$warm_up_iters,"warm-up +",x$main_iterations,"main; stop:",x$stop_reason,"\n")
  cat("Final ELBO:",format(x$final_elbo,digits=max(digits,8)),"(",x$elbo_convention,")\n")
  cat("ELBO decreases:",x$elbo_decreases,"\n")
  if(!is.null(x$quadrature)) {
    cat("Quadrature: max order",x$quadrature$max_order,"; max refinement discrepancy",
      format(x$quadrature$max_refinement_error,digits=3),"; max tail discrepancy",
      format(x$quadrature$max_tail_check_error,digits=3),"\n")
  }
  print(x$motifs,digits=digits,row.names=FALSE)
  cat("\nPopulation posterior and slot summaries:\n")
  print(x$slots,digits=digits,row.names=FALSE); invisible(x)
}

#' Plot ELBO convergence, learned loadings, or factor graphs
#'
#' Loading bars use posterior means averaged over mixture uncertainty.
#' Hard cluster labels determine grouping only. No truth labels are needed.
#' Empty hard-assignment clusters are omitted; fitted components are not merged.
#'
#' @param x A fitted `miso_fit` object.
#' @param type `"elbo"`, `"loadings"`, `"subgraphs"`, or `"aggregate"`.
#'   Graph types delegate to [plot_miso_graphs()] and accept its options via `...`.
#' @param normalize Show within-observation fractions rather than absolute
#'   expected loadings (counts for [miso_fit()], per unit of `n` for
#'   [miso_fit_length()]). Used for loading plots only.
#' @param cluster_order Optional permutation of occupied fitted cluster indices.
#' @param factor_order Optional permutation of all fitted factor indices.
#' @param sort_by Optional fitted factor index. Within each cluster, order bars
#'   by the loading fraction of that factor. Otherwise preserve observation order.
#' @param col Optional color (ELBO) or vector with one color per factor in the
#'   original factor order (loadings or graphs). For graphs, overrides `colors`.
#' @param main Optional plot title.
#' @param ... Additional arguments to `plot()`, `barplot()`, or [plot_miso_graphs()].
#' @return Invisibly, for loading plots, a list containing the N by K posterior
#'   mean `loadings` in original order, displayed `values`, `observation_order`,
#'   `cluster_order`, `factor_order`, and `bar_midpoints`. For ELBO plots, returns
#'   the fitted object invisibly. Graph plots return the [miso_fit_graphs()] summary.
#' @export
#' @examples
#' Y <- matrix(c(5, 1, 4, 2, 1, 5, 2, 4), 4, 2, byrow = TRUE)
#' fit <- miso_fit(Y, F = diag(2), D = 1, S = 2, max_iters = 5, tol = 0)
#' plot(fit)
#' plot(fit, type = "loadings")
plot.miso_fit <- function(x, type = c("elbo", "loadings", "subgraphs", "aggregate"), normalize = TRUE,
                          cluster_order = NULL, factor_order = NULL,
                          sort_by = NULL, col = NULL, main = NULL, ...) {
  type = match.arg(type)
  if (type %in% c("subgraphs", "aggregate")) {
    args = list(...)
    if (!is.null(col)) args$colors = col
    return(invisible(do.call(plot_miso_graphs, c(list(graphs = x, type = type, main = main), args))))
  }
  if (type == "elbo") {
    if (is.null(main)) main = "MiSo ELBO"
    if (is.null(col)) col = "steelblue"
    graphics::plot(seq_along(x$elbo), x$elbo, type = "l", col = col,
                   xlab = "Outer iteration", ylab = "ELBO", main = main, ...)
    if(!is.null(x$warm_up_iters) && x$warm_up_iters>0L) {
      graphics::abline(v=x$warm_up_iters+.5,lty=2,col="grey50")
      graphics::legend("bottomright","End of fixed-F GS warm-up",lty=2,col="grey50",bty="n",cex=.8)
    }
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
