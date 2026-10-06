#' Fit a descending path of MiSO motif counts
#'
#' Fit once at `S[1]`, then repeatedly merge two fitted motifs and refit through
#' the remaining values of `S`. This is a warm-start heuristic, not independent fits or a guarantee
#' of global optima. No motif count is selected automatically.
#'
#' @inheritParams miso_fit
#' @param S A nonempty integer vector forming a consecutive descending sequence,
#'   e.g. `6:1` or `6:3`. A scalar fits one model and still returns a path.
#'   Ascending, repeated, skipped, or nonpositive counts are rejected.
#' @param D,K Positive slot and factor counts, inferred from init when supplied.
#' @param model `"count"` for [miso_fit()] or `"length"` for [miso_fit_length()].
#' @param n Exposures for `model = "length"`; NULL uses row totals, as in
#'   [miso_fit_length()]. All must be positive. Must be NULL for `model = "count"`.
#' @param init Optional initialization for the first fit, with exactly `S[1]`
#'   motifs. For the length model its loading parameters must already be on the per-unit-n
#'   scale. Reuse [miso_init()] or a fit's retained initialization.
#' @param seed Random seed set once at the beginning of the path and passed to
#'   the first fitter's initialization.
#' @param ... Named fitting controls shared across the path, e.g. `max_iters`,
#'   `tol`, `min_iters`, `patience`, `update_F`, `verbose`, `phi0`, and
#'   `keep_initialization`. `init_control` and `warm_up_iters` apply only to the
#'   first fit. Post-merge fits use zero warm-up sweeps and continue with Jacobi.
#'   Unknown, unnamed, or duplicated controls are rejected. The path requires
#'   a symmetric Dirichlet prior: `phi0` must be scalar (or inherited from a
#'   symmetric initialization). Its per-motif value stays fixed across S.
#'
#' @details The initial model is fitted and successive models merge one pair of
#'   motifs. Pair scores use mass-weighted cosine distance of population profiles
#'   based on E(alpha/beta). Slots are aligned by maximum selection-probability
#'   overlap. Responsibilities and local sufficient statistics N,T,R are pooled;
#'   the joint population posterior is recomputed and one allocation pass initializes
#'   the merged local Gamma factors. The dictionary and other components are retained.
#'   The same numerical population prior is used throughout the path. All ELBOs
#'   include the data-only constant. Merges are approximate initializations; bounds
#'   need not be monotone as S changes. Compare converged fits, preferably from
#'   several starts; the path is not a guarantee of consistent model selection.
#' @return A `miso_fit_path` list containing:
#'   \item{fits}{A descending list of `miso_fit` objects named `S6`, `S5`, etc.
#'     Use names to select a motif count, e.g. `path$fits[["S3"]]`. Each fit
#'     retains its entire ELBO trajectory and supports prediction and plotting.}
#'   \item{elbo}{Named vector of final ELBOs in the same order as `fits`.}
#'   \item{summary}{Data frame with `S`, `ELBO`, `iterations`, and `converged`.}
#'   \item{merges}{One row per merge: `from_S`, `to_S`, `first`, `second`,
#'     `cosine`, `cost`, and a list column `slot_order` mapping retained slots
#'     to slots in the second motif. Indices refer to the fit at `from_S`.}
#'   \item{model, S, seed, n, call}{Path settings, resolved exposures
#'     (NULL for count fits), and the original call.}
#' @seealso [miso_fit()], [miso_fit_length()], [plot.miso_fit_path()]
#' @export
#' @examples
#' Y <- matrix(c(5, 1, 4, 2, 1, 5, 2, 4), 4, 2, byrow = TRUE)
#' path <- miso_fit_path(Y, S = 3:1, K = 2, D = 2,
#'                       max_iters = 5, tol = 0, seed = 1)
#' path$summary
#' plot(path)
#' plot(path$fits[["S2"]], type = "subgraphs")
#' length_path <- miso_fit_path(Y, S = 3:1, K = 2, D = 2,
#'   model = "length", n = rowSums(Y), max_iters = 5, tol = 0)
miso_fit_path <- function(Y, S, F = NULL, D = NULL, K = NULL,
                          init = NULL, model = c("count", "length"), n = NULL,
                          seed = 1, ...) {
  call <- match.call()
  model <- match.arg(model)
  .miso_stop(is.numeric(S) && is.null(dim(S)) && length(S) > 0 &&
    all(is.finite(S)) && all(S >= 1 & S <= .Machine$integer.max) &&
    all(S == trunc(S)) && all(diff(S) == -1),
    "S must be a consecutive descending vector of positive integers, e.g. 6:1.")
  S <- as.integer(S)
  Y <- .miso_prepare_counts(Y)
  .miso_stop(S[1] <= nrow(Y), "S[1] must not exceed the number of observations.")
  .miso_stop(length(seed) == 1 && is.numeric(seed) && is.finite(seed) &&
    abs(seed) <= .Machine$integer.max && seed == as.integer(seed), "seed must be an integer.")
  controls <- list(...)
  reserved <- c("Y", "F", "D", "K", "init", "S", "n", "seed")
  .miso_stop(length(controls) == 0 || (!is.null(names(controls)) &&
    all(nzchar(names(controls))) && !anyDuplicated(names(controls)) &&
    all(names(controls) %in% setdiff(names(formals(miso_fit)), reserved))),
    "Supply named, unique fitting controls supported by miso_fit().")
  if ("phi0" %in% names(controls)) {
    .miso_stop(is.numeric(controls$phi0) && length(controls$phi0) == 1 &&
      is.finite(controls$phi0) && controls$phi0 > 0, "phi0 must be a positive scalar for a path.")
  } else if (!is.null(init$phi0)) {
    prior <- init$phi0
    .miso_stop(is.numeric(prior) && length(prior) %in% c(1L, S[1]) &&
      all(is.finite(prior)) && all(prior > 0) && all(prior == prior[1]),
      "The initialization must have a symmetric phi0 for a path.")
  }
  if (model == "count") {
    .miso_stop(is.null(n), "n requires model = 'length'.")
  } else {
    if (is.null(n)) n <- .miso_row_sums(Y)
    .miso_stop(is.numeric(n) && is.null(dim(n)) && length(n) == nrow(Y) &&
      all(is.finite(n)) && all(n > 0), "n must contain one positive exposure per observation.")
    n <- as.numeric(n)
  }
  set.seed(seed)
  fits <- stats::setNames(vector("list", length(S)), paste0("S", S))
  merges <- data.frame(from_S = integer(), to_S = integer(), first = integer(),
    second = integer(), cosine = numeric(), cost = numeric())
  merges$slot_order <- I(list())
  fitter <- if (model == "count") miso_fit else miso_fit_length
  args <- c(list(Y = Y, F = F, D = D, K = K, init = init, S = S[1], seed = seed), controls)
  if (model == "length") args$n <- n
  for (j in seq_along(S)) {
    if (j > 1L) {
      merged <- .miso_path_merge(fits[[j-1L]],Y,
        block_size=if(is.null(controls$block_size)) 100L else controls$block_size)
      merges <- rbind(merges, merged$record)
      # Use the preceding learned dictionary; initialization controls apply once.
      args$warm_up_iters <- 0L
      args$F <- NULL
      args$init_control <- NULL
      args$init <- merged$init
      args$S <- S[j]
      args$population_prior <- fits[[j-1L]]$population_prior
    }
    if (isTRUE(controls$verbose)) message("MiSO path: S = ", S[j])
    fits[[j]] <- do.call(fitter, args)
  }
  elbo <- vapply(fits, function(f) f$final_elbo, 0.0)
  summary <- data.frame(S = S, ELBO = unname(elbo),
    iterations = vapply(fits, function(f) f$n_iter, 0L),
    warm_up_iters = vapply(fits, function(f) f$warm_up_iters, 0L),
    main_iterations = vapply(fits, function(f) f$n_iter_main, 0L),
    converged = vapply(fits, function(f) f$converged, TRUE), row.names = NULL)
  if (any(!summary$converged) && (is.null(controls$tol) || controls$tol != 0))
    warning("Non-converged fits at S = ", paste(summary$S[!summary$converged], collapse = ", "),
            ". Inspect path$summary and the per-fit ELBO histories.", call. = FALSE)
  structure(list(fits = fits, elbo = elbo, summary = summary, merges = merges,
    model = model, S = S, seed = seed,
    n = n, call = call), class = "miso_fit_path")
}

# Merge aligned components by pooling local sufficient statistics, then perform
# one local initialization pass. This is a starting state, not an ELBO ascent step.
.miso_path_merge <- function(f,Y,block_size=100L,eps=1e-12) {
  S <- dim(f$gamma)[1]; D <- dim(f$gamma)[2]; K <- dim(f$gamma)[3]; N <- nrow(f$omega)
  .miso_stop(S>=2 && identical(f$inference_model,"bayesian_mixture"),"Merge requires a Bayesian fit with at least two motifs.")
  means <- f$population_mean
  profile <- t(matrix(vapply(seq_len(S),function(s)
    colSums(matrix(f$gamma[s,,],D,K)*means[s,]),numeric(K)),K,S))
  norm <- sqrt(rowSums(profile^2)); unit <- profile/pmax(norm,eps)
  cosine <- pmin(pmax(tcrossprod(unit),-1),1)
  cosine[outer(norm==0,norm==0,"&")] <- 1
  mass <- colSums(f$omega)
  cost <- (1-cosine)*outer(mass,mass)/pmax(outer(mass,mass,"+"),eps)
  cost[!upper.tri(cost)] <- Inf
  pair <- arrayInd(which.min(cost),dim(cost)); u <- pair[1]; v <- pair[2]
  gu <- matrix(f$gamma[u,,],D,K); gv <- matrix(f$gamma[v,,],D,K)
  perm <- .miso_anchor_assignment(gu %*% t(gv))
  weight <- mass[c(u,v)]; weight <- if(sum(weight)>0) weight/sum(weight) else c(.5,.5)
  merged_gamma <- weight[1]*gu+weight[2]*gv[perm,,drop=FALSE]
  wu <- f$omega[,u]; wv <- f$omega[,v]; wg <- wu+wv
  au <- matrix(f[["a"]][,u,],N,D); av <- matrix(f[["a"]][,v,],N,D)[,perm,drop=FALSE]
  bu <- .miso_rate_matrix(f[["b"]],u,N,D); bv <- .miso_rate_matrix(f[["b"]],v,N,D)[,perm,drop=FALSE]
  lu <- digamma(au)-log(bu); lv <- digamma(av)-log(bv)
  T <- colSums(wu*(au/bu)+wv*(av/bv))
  R <- colSums(wu*lu+wv*lv)
  control <- f$quadrature_control
  blocks <- lapply(seq_len(D),function(d)
    .miso_gamma_pair(sum(wg),T[d],R[d],.miso_prior_vector(f$population_prior),control))
  elog <- (wu*lu+wv*lv)/pmax(wg,.Machine$double.xmin)
  for(d in seq_len(D)) elog[wg==0,d] <- weight[1]*lu[wg==0,d]+weight[2]*lv[wg==0,d]
  prepared <- .miso_prepare_pass(.miso_prepare_counts(Y),block_size)
  pass <- prepared$pass(prepared$blocks,elog,merged_gamma %*% log(pmax(f$F,eps)),eps=eps)
  keep <- setdiff(seq_len(S),v); newS <- S-1L
  init <- list(F=f$F,gamma=f$gamma[keep,,,drop=FALSE],a=f[["a"]][,keep,,drop=FALSE],
    b=if(length(dim(f[["b"]]))==3L) f[["b"]][,keep,,drop=FALSE] else f[["b"]][keep,,drop=FALSE],
    omega=f$omega[,keep,drop=FALSE],phi0=rep(f$phi0[1],newS),population_prior=f$population_prior,
    format_version=2L,loading_units=f$loading_units,input_dimnames=f$input_dimnames)
  init$gamma[u,,] <- merged_gamma; init$omega[,u] <- wg
  init$phi <- init$phi0+colSums(init$omega)
  init$q_population <- vector("list",newS*D)
  for(d in seq_len(D)) for(s in seq_len(newS))
    init$q_population[[s+(d-1L)*newS]] <- if(s==u) blocks[[d]] else f$q_population[[keep[s]+(d-1L)*S]]
  am <- vapply(blocks,function(h) h$mean_alpha,0.0)
  bm <- vapply(blocks,function(h) h$mean_beta,0.0)
  init[["a"]][,u,] <- sweep(pass$allocated,2,am,"+")
  if(is.null(f[["n"]])) init[["b"]][u,] <- bm+1 else {
    init[["n"]] <- f[["n"]]; init[["b"]][,u,] <- matrix(bm,N,D,byrow=TRUE)+f[["n"]]
  }
  record <- data.frame(from_S=S,to_S=newS,first=u,second=v,cosine=cosine[u,v],cost=cost[u,v])
  record$slot_order <- I(list(perm))
  list(init=init,record=record,pooled=list(N=sum(wg),T=T,R=R))
}

#' Print a MiSO fitting path
#' @param x A [miso_fit_path()] result.
#' @param ... Additional arguments passed to the summary table's print method.
#' @return The path, invisibly.
#' @export
print.miso_fit_path <- function(x, ...) {
  cat("MiSO", x$model, "path: S =", x$S[1], "to", tail(x$S, 1), "\n")
  do.call(print, c(list(x = x$summary), utils::modifyList(list(row.names = FALSE), list(...))))
  invisible(x)
}

#' Plot the ELBO along a MiSO fitting path
#' @param x A [miso_fit_path()] result.
#' @param ... Additional arguments to [graphics::plot.default()], including
#'   `xlab`, `ylab`, `main`, `type`, and `col`. The x axis is motif count S.
#' @return The path's summary table, invisibly. No S is selected automatically.
#' @export
plot.miso_fit_path <- function(x, ...) {
  tab <- x$summary[order(x$summary$S), , drop = FALSE]
  options <- utils::modifyList(list(type = "b", pch = 16, xlab = "Number of motifs, S",
    ylab = "Final ELBO", main = "MiSO motif-count path", xaxt = "n"), list(...))
  do.call(graphics::plot, c(list(x = tab$S, y = tab$ELBO), options))
  if (identical(options$xaxt, "n") && !identical(options$axes, FALSE))
    graphics::axis(1, at = tab$S)
  invisible(x$summary)
}
