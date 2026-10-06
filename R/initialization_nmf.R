## Direct NMF initialization: exact supports, optional cosine merges or splits.
.miso_initialize_from_loadings <- function(L, D, S = NULL, tau = 0.9,
                                          gamma_floor = 0.05, phi0 = 0.1,
                                          seed = 1, eps = 1e-12,
                                          unsupported_fraction = 0.01) {
  .miso_positive_integer(D, "D")
  .miso_stop(is.matrix(L) && is.numeric(L) && all(dim(L) > 0) &&
               all(is.finite(L)) && all(L >= 0),
             "L must be a finite nonnegative N by K matrix.")
  .miso_stop(length(tau) == 1 && is.finite(tau) && tau > 0 && tau <= 1,
             "tau must be in (0, 1].")
  .miso_stop(length(gamma_floor) == 1 && is.finite(gamma_floor) &&
               gamma_floor >= 0 && gamma_floor < 1,
             "gamma_floor must be in [0, 1).")
  .miso_stop(length(unsupported_fraction) == 1 && is.finite(unsupported_fraction) &&
               unsupported_fraction > 0 && unsupported_fraction <= 1,
             "unsupported_fraction must be in (0, 1].")
  N = nrow(L)
  K = ncol(L)
  if (!is.null(S)) {
    .miso_positive_integer(S, "S")
    .miso_stop(S <= N, "S must not exceed the number of observations.")
  }
  totals = rowSums(L)
  .miso_stop(all(is.finite(totals)), "L row totals must be finite.")
  pattern = lapply(seq_len(N), function(i) {
    if (totals[i] == 0) return(integer())
    ranked = order(-L[i, ], seq_len(K))
    fraction = cumsum(L[i, ranked]) / totals[i]
    fraction[K] = 1 # Avoid rounding below one when tau = 1.
    r = which(fraction >= tau)[1]
    sort(ranked[seq_len(min(r, D))])
  })
  keys = vapply(pattern, paste, "", collapse = ",")
  raw_cluster = match(keys, unique(keys))
  S0 = max(raw_cluster)
  members = lapply(seq_len(S0), function(s) which(raw_cluster == s))
  if (is.null(S)) S = S0
  centers = do.call(rbind, lapply(members, function(ix) colMeans(L[ix, , drop = FALSE])))
  sizes = lengths(members)
  # Preserve full original support unions: a previously truncated factor may
  # become important after a subsequent merge.
  unions = lapply(members, function(ix) pattern[[ix[1]]])
  history = list()
  if (S < S0) {
    unit_centers = function(x) {
      scale = apply(x, 1, max)
      x = x / ifelse(scale > 0, scale, 1)
      norm = sqrt(rowSums(x^2))
      x / ifelse(norm > 0, norm, 1)
    }
    unit = unit_centers(centers)
    similarity = tcrossprod(unit)
    zero = which(rowSums(abs(centers)) == 0)
    similarity[zero, zero] = 1
    diag(similarity) = -Inf
    while (length(members) > S) {
      score = similarity
      score[!upper.tri(score)] = -Inf
      pair = arrayInd(which.max(score), dim(score))
      a = pair[1]; b = pair[2]
      history[[length(history) + 1L]] = list(
        action = "merge", groups = c(a, b), sizes = sizes[c(a, b)],
        cosine = similarity[a, b])
      weight = sizes[a] / (sizes[a] + sizes[b])
      centers[a, ] = weight * centers[a, ] + (1 - weight) * centers[b, ]
      sizes[a] = sizes[a] + sizes[b]
      members[[a]] = c(members[[a]], members[[b]])
      unions[[a]] = sort(unique(c(unions[[a]], unions[[b]])))
      members = members[-b]; unions = unions[-b]; sizes = sizes[-b]
      centers = centers[-b, , drop = FALSE]
      similarity = similarity[-b, -b, drop = FALSE]
      unit = unit_centers(centers)
      values = as.numeric(unit %*% unit[a, ])
      if (all(centers[a, ] == 0)) values[rowSums(abs(centers)) == 0] = 1
      similarity[a, ] = similarity[, a] = values
      diag(similarity) = -Inf
    }
  } else if (S > S0) {
    set.seed(seed)
    while (length(members) < S) {
      a = which.max(lengths(members))
      ix = members[[a]]
      shuffled = ix[sample.int(length(ix))]
      cut = floor(length(ix) / 2)
      members[[a]] = shuffled[seq_len(cut)]
      members[[length(members) + 1L]] = shuffled[seq.int(cut + 1L, length(ix))]
      unions[[length(unions) + 1L]] = unions[[a]]
      history[[length(history) + 1L]] = list(
        action = "split", group = a, sizes = c(cut, length(ix) - cut))
    }
    centers = do.call(rbind, lapply(members, function(ix) colMeans(L[ix, , drop = FALSE])))
  }
  .miso_stop(length(phi0) %in% c(1L, S) && all(is.finite(phi0)) && all(phi0 > 0),
             "phi0 must contain one positive value or one per final motif.")
  support = lapply(seq_len(S), function(s) {
    candidates = unions[[s]]
    ranked = candidates[order(-centers[s, candidates], candidates)]
    sort(head(ranked, D))
  })
  anchors = lapply(seq_len(S), function(s) {
    ix = support[[s]]
    ix[order(-centers[s, ix], ix)]
  })
  cluster = integer(N)
  for (s in seq_len(S)) cluster[members[[s]]] = s
  omega = matrix(0, N, S)
  omega[cbind(seq_len(N), cluster)] = 1
  gamma = array(1 / K, c(S, D, K))
  shape_seed = rate_seed = matrix(1, S, D)
  for (s in seq_len(S)) {
    for (d in seq_along(anchors[[s]])) {
      k = anchors[[s]][d]
      gamma[s, d, ] = gamma_floor / K
      gamma[s, d, k] = 1 - gamma_floor + gamma_floor / K
      # Gamma MLE from point estimates, bounded to avoid infinite shape for
      # singleton/constant groups. The same shape bounds are used by EB.
      values = pmax(L[members[[s]], k], eps)
      mean_value = mean(values)
      shape_seed[s, d] = .miso_gamma_shape_from_moments(mean_value, mean(log(values)))
      rate_seed[s, d] = shape_seed[s, d] / mean_value
    }
    unsupported = D - length(anchors[[s]])
    if (unsupported > 0) {
      mean_value = max(unsupported_fraction * mean(totals[members[[s]]]) /
                         unsupported, eps)
      rate_seed[s, seq.int(length(anchors[[s]]) + 1L, D)] = 1 / mean_value
    }
  }
  list(S = S, S0 = S0, F = NULL, omega = omega,
       phi0 = rep(phi0, length.out = S), phi = rep(phi0, length.out = S) + colSums(omega),
       gamma = gamma, shape_seed = shape_seed, rate_seed = rate_seed,
       cluster = cluster, initial_cluster = raw_cluster, pattern = pattern,
       support = support, anchors = anchors, centers = centers,
       expected_loading = L, history = history)
}
