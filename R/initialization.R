## Expected-loading initialization in Section 3.3 of the paper.


.miso_cosine <- function(x, y, eps = 1e-12) {
  sum(x * y) / pmax(sqrt(sum(x^2) * sum(y^2)), eps)
}

.miso_nearest_profile <- function(profile, candidates) {
  which.max(apply(candidates, 1, .miso_cosine, y = profile))
}

.miso_anchor_assignment <- function(score) {
  A = nrow(score)
  D = ncol(score)
  if (A == 0) return(integer())

  best_value = -Inf
  best_assignment = integer(A)
  search = function(d, available, assignment, value) {
    if (d > A) {
      if (value > best_value) {
        best_value <<- value
        best_assignment <<- assignment
      }
      return(invisible(NULL))
    }
    for (slot in available) {
      search(
        d + 1, available[available != slot],
        c(assignment, slot), value + score[d, slot]
      )
    }
    invisible(NULL)
  }
  search(1, seq_len(D), integer(), 0)
  best_assignment
}

.miso_consolidate_patterns <- function(pattern, loading, min_fraction) {
  N = nrow(loading)
  keys = vapply(pattern, function(x) paste(sort(x), collapse = ","), "")
  raw_keys = unique(keys)
  raw_cluster = match(keys, raw_keys)
  raw_members = lapply(seq_along(raw_keys), function(s) which(raw_cluster == s))
  raw_support = lapply(raw_members, function(index) pattern[[index[1]]])
  raw_profile = do.call(rbind, lapply(
    raw_members, function(index) colMeans(loading[index, , drop = FALSE])
  ))
  raw_size = lengths(raw_members)
  cutoff = min_fraction * N
  retained = which(raw_size >= cutoff)

  if (length(retained) >= 2) {
    raw_destination = rep(NA_integer_, length(raw_members))
    raw_destination[retained] = seq_along(retained)
    for (s in setdiff(seq_along(raw_members), retained)) {
      raw_destination[s] = .miso_nearest_profile(
        raw_profile[s, ], raw_profile[retained, , drop = FALSE]
      )
    }
    cluster = raw_destination[raw_cluster]
    support = raw_support[retained]
  } else {
    members = raw_members
    support = raw_support
    while (length(members) > 1 && any(lengths(members) < cutoff)) {
      profiles = do.call(rbind, lapply(
        members, function(index) colMeans(loading[index, , drop = FALSE])
      ))
      smallest = which.min(lengths(members))
      neighbors = setdiff(seq_along(members), smallest)
      nearest = neighbors[.miso_nearest_profile(
        profiles[smallest, ], profiles[neighbors, , drop = FALSE]
      )]
      members[[nearest]] = c(members[[nearest]], members[[smallest]])
      members = members[-smallest]
      support = support[-smallest]
    }
    cluster = integer(N)
    for (s in seq_along(members)) cluster[members[[s]]] = s
  }

  list(cluster = cluster, support = support)
}

.miso_initialize_from_susie <- function(ps, D, tau = 0.9,
                                       min_fraction = 0.05,
                                       gamma_floor = 0.05,
                                       phi0 = 0.1, eps = 1e-12) {
  gamma_ps = ps$gamma
  alpha_ps = ps$alpha
  beta_ps = ps$beta
  dimensions = dim(gamma_ps)
  .miso_stop(length(dimensions) == 3 && dimensions[2] == D,
             "ps$gamma must be an N by D by K array.")
  .miso_stop(identical(dim(alpha_ps), dimensions[1:2]) &&
               identical(dim(beta_ps), dimensions[1:2]),
             "ps$alpha and ps$beta must both be N by D matrices.")
  .miso_stop(all(is.finite(gamma_ps)) && all(gamma_ps >= 0) &&
               all(is.finite(alpha_ps)) && all(alpha_ps > 0) &&
               all(is.finite(beta_ps)) && all(beta_ps > 0),
             "Poisson SuSiE parameters must be finite and valid.")
  .miso_stop(tau > 0 && tau <= 1, "tau must be in (0, 1].")
  .miso_stop(min_fraction > 0 && min_fraction <= 1,
             "min_fraction must be in (0, 1].")
  .miso_stop(gamma_floor >= 0 && gamma_floor < 1,
             "gamma_floor must be in [0, 1).")

  N = dimensions[1]
  K = dimensions[3]
  loading = matrix(0, N, K)
  posterior_mean = posterior_mean_log = matrix(0, N, D)
  for (i in seq_len(N)) {
    gamma_i = matrix(gamma_ps[i, , ], D, K)
    alpha_i = matrix(alpha_ps[i, ], D, K)
    beta_i = matrix(beta_ps[i, ], D, K)
    loading[i, ] = colSums(
      gamma_i * alpha_i / beta_i
    )
    posterior_mean[i, ] = rowSums(gamma_i * alpha_i / beta_i)
    posterior_mean_log[i, ] = rowSums(
      gamma_i * (digamma(alpha_i) - log(beta_i))
    )
  }

  pattern = vector("list", N)
  for (i in seq_len(N)) {
    order_i = order(loading[i, ], decreasing = TRUE)
    r = which(cumsum(loading[i, order_i]) / sum(loading[i, ]) >= tau)[1]
    pattern[[i]] = sort(order_i[seq_len(min(r, D))])
  }

  consolidated = .miso_consolidate_patterns(
    pattern, loading, min_fraction
  )
  cluster = consolidated$cluster
  support = consolidated$support
  S = max(cluster)
  omega = matrix(0, N, S)
  omega[cbind(seq_len(N), cluster)] = 1
  .miso_stop(length(phi0) %in% c(1, S) &&
               all(is.finite(phi0)) && all(phi0 > 0),
             "phi0 must contain one positive value or one per motif.")
  phi0 = rep(phi0, length.out = S)
  phi = phi0 + colSums(omega)

  cluster_loading = do.call(rbind, lapply(seq_len(S), function(s) {
    colMeans(loading[cluster == s, , drop = FALSE])
  }))
  anchors = vector("list", S)
  gamma = array(1 / K, c(S, D, K))
  for (s in seq_len(S)) {
    anchors[[s]] = support[[s]][
      order(cluster_loading[s, support[[s]]], decreasing = TRUE)
    ]
    for (d in seq_along(anchors[[s]])) {
      gamma[s, d, ] = gamma_floor / K
      gamma[s, d, anchors[[s]][d]] =
        1 - gamma_floor + gamma_floor / K
    }
  }

  permutation = matrix(0L, N, D)
  for (i in seq_len(N)) {
    s = cluster[i]
    A = length(anchors[[s]])
    anchor_score = matrix(0, A, D)
    if (A > 0) {
      for (d in seq_len(A)) {
        k = anchors[[s]][d]
        anchor_score[d, ] = gamma_ps[i, , k] *
          alpha_ps[i, ] /
          beta_ps[i, ]
      }
    }
    anchored_slots = .miso_anchor_assignment(anchor_score)
    remaining_slots = setdiff(seq_len(D), anchored_slots)
    slot_loading = posterior_mean[i, ]
    remaining_slots = remaining_slots[
      order(slot_loading[remaining_slots], decreasing = TRUE)
    ]
    permutation[i, ] = c(anchored_slots, remaining_slots)
  }

  alpha0 = beta0 = matrix(0, S, D)
  for (s in seq_len(S)) {
    members = which(cluster == s)
    for (d in seq_len(D)) {
      aligned = cbind(members, permutation[members, d])
      mean_lambda = mean(posterior_mean[aligned])
      mean_log_lambda = mean(posterior_mean_log[aligned])
      alpha0[s, d] = .miso_gamma_shape_from_moments(
        mean_lambda, mean_log_lambda
      )
      beta0[s, d] = alpha0[s, d] / pmax(mean_lambda, eps)
    }
  }

  list(
    S = S,
    omega = omega,
    phi = phi,
    gamma = gamma,
    alpha0 = alpha0,
    beta0 = beta0,
    cluster = cluster,
    pattern = pattern,
    support = support,
    anchors = anchors,
    expected_loading = loading,
    slot_permutation = permutation
  )
}

.miso_initialize <- function(Y, F, D, tau = 0.9,
                            min_fraction = 0.05,
                            gamma_floor = 0.05, phi0 = 0.1,
                            susie_max_iters = 100,
                            susie_tol = 1e-6, seed = 1, eps = 1e-12,
                            keep_intermediates = FALSE) {
  F = .miso_normalize_rows(F, eps)
  ps = poisson_susie_fit(
    Y, F, D, max_iters = susie_max_iters, tol = susie_tol,
    seed = seed, eps = eps, keep_fits = keep_intermediates
  )
  initialization = .miso_initialize_from_susie(
    ps, D, tau, min_fraction, gamma_floor, phi0, eps
  )
  initialization$F = F
  if (keep_intermediates) initialization$poisson_susie = ps
  initialization
}

