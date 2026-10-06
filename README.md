# misoR

MiSO models nonnegative observations as mixtures of sparse Poisson factor models.
Version 0.2.0 uses Bayesian mixture parameters (`pi`, population `alpha` and `beta`,
and factor selections `gamma`) and point-estimates the dictionary `F` by variational EM.
Both ordinary and length-adjusted models share the same blocked inference engine.

```r
library(misoR)
fit <- miso_fit(Y, K = 12, D = 5, S = 6)
summary(fit)
L <- predict(fit, "loadings")
fitted_counts <- predict(fit)
plot(fit, type = "loadings")
plot(fit, type = "subgraphs")

# Exposures default to row totals; omit empty observations with this default.
fit_length <- miso_fit_length(Y, K = 12, D = 5, S = 6, n = rowSums(Y))
```

`Y` can be dense or a sparse Matrix. Fractional nonnegative values are accepted as
a generalized Poisson objective; the count likelihood interpretation requires integers.
Rows of `F` sum to one. Length-model loadings are per unit of exposure and are not
constrained to sum to one. With observed row totals as exposures this is a working
Poisson likelihood, not the exact conditional multinomial model.

## Population prior

The default calibration uses `Y` directly, without an NMF fit. If `t` is the mean
row total (ordinary model), or mean row total divided by exposure (length model),
the prior targets an average slot loading of `t/D`. Under shape/rate parameterization,
`alpha ~ Gamma(2,2)` and `beta ~ Gamma(2,t/D)` independently achieve this target.
For exposures equal to row totals, `t=1`. Shape 2 is the hyperprior shape, distinct
from the random population loading shape `alpha`, whose prior mean is 1.

```r
prior <- miso_population_prior(Y, D = 5)   # retains calibration diagnostics
fit <- miso_fit(Y, K = 12, D = 5, S = 6, population_prior = prior)
# Or specify positive shape/rate pairs; beta's shape must exceed 1.
prior <- list(alpha = c(shape = 2, rate = 2),
              beta = c(shape = 2, rate = 10))
```

Use the same prior across starts and candidate S. Across D the calibration rule
stays fixed, keeping the expected *total* loading constant; the beta prior rate
therefore scales as 1/D. All-zero data require an explicit prior. Calibration uses
data, so describe it as a data-calibrated prior rather than a prespecified prior.

## Initialization and continuing a fit

```r
initial <- miso_init(Y, K = 12, D = 5, S = 6)
fit <- miso_fit(Y, init = initial, S = 6, keep_initialization = TRUE)
continued <- miso_fit(Y, init = fit, S = 6, max_iters = 500)

# Reuse count-scale NMF factors and loadings without refitting NMF:
initial <- miso_init(Y, F = F_nmf, L = L_nmf, D = 5, S = 6)
# Holding a supplied dictionary fixed is explicit:
fit <- miso_fit(Y, F = F_nmf, D = 5, S = 6, update_F = FALSE)
```

A custom initialization minimally supplies `F` and `gamma` (S by D by K).
It may also supply `omega` (N by S), `a` and `b` together, and `q_population`.
Ordinary `b` is S by D; length-model `b` is N by S by D. Local `a` is always
N by S by D. Exposures, data order, dimensions and loading units must match.
Named rows/columns are checked; ordering of unnamed data remains the caller's
responsibility. A fitted object's prior is inherited unless explicitly overridden;
incompatible `q_population` blocks are rejected. Remove those blocks to use a fit
as a new initialization under a changed prior. This is no longer an exact continuation.

NMF estimates determine starting values only. `initial$population_seed` records
initialization-only Gamma fits; it is neither the prior specification nor the
population posterior. New global posterior blocks are computed during fitting.

## Returned quantities and ELBO

- `a`, `b`: local Gamma variational parameters.
- `alpha_mean`, `beta_mean`: posterior means of population parameters.
- `q_population`: joint posterior moments, covariance, KL, normalizer and sufficient statistics.
- `population_mean`: E(alpha/beta), not the ratio of posterior means.
- `gamma`: selection probabilities; `omega`: observation responsibilities; `pi`: mean mixture weights.
- `elbo_terms`: local bound, assignment term, KL_pi, KL_gamma, KL_population, constant and ELBO.
- `quadrature`: per-block convergence, order, refinement discrepancy, bounds and tail sensitivity.
- `quadrature_history`: maximum numerical discrepancies in each fitted iteration.

The ELBO **includes** the data-only Poisson constant in both models. Do not add
`elbo_constant` again. `final_elbo` is the final entry of `elbo`. The state is the
last complete coordinate sweep; local rates need not equal the *latest* population
mean rate plus exposure until convergence. Numerical quadrature discrepancies
are estimates, not certified bounds. Integration failures raise an informative error.

```r
path <- miso_fit_path(Y, S = 8:4, K = 12, D = 5,
  max_iters = 2000, tol = 1e-8, min_iters = 30, patience = 8)
plot(path)
```

The path aligns slots, pools population-update sufficient statistics, and refits
following each merge. It holds the prior fixed. It is a warm-start heuristic;
ELBO need not rise as S decreases and does not guarantee correct model selection.

## Migration from version 0.0.1

This is an intentional field-layout change; ambiguous old aliases are not retained.
Use `miso_convert_fit(old_fit, Y)` to read an old fitted object in the new layout.
The converted object remains marked `legacy_eb` until refitted. Prediction is unchanged.
With Y, its old ELBO is converted to the full-constant convention. Without Y, an old
ordinary fit remains explicitly marked `legacy_without_constant`. Even with the
constant restored, legacy and Bayesian objectives have different population-prior
terms and must not be mixed for model selection; refit under the new model first.

Old `alpha`, `beta` become `a`, `b`; old `alpha0`, `beta0` become `alpha_mean`,
`beta_mean`; `prior_mean` becomes `population_mean`. Population posterior SDs are
unavailable for converted EB fits. Rebuild old standalone initialization lists with
`miso_init()`; do not pass `alpha0`/`beta0` to the new fitter. The old `update_prior`
switch is removed: Bayesian population blocks are always updated. The exported
row-wise `poisson_susie_fit()` retains its existing empirical-Bayes API as an
initialization utility.

The prior empirical-Bayes equation tests are retained with the complete old source
snapshot. The new package tests cover joint quadrature against independent integration,
analytic local updates and ELBOs, sparse/dense and exposure invariance, exact warm
continuation, merge statistics, migration, initialization and plotting.

## Warm-up, convergence, and explicit motif counts

`S` is required in `miso_fit()` and `miso_fit_length()`, including when `init`
is supplied. It must match that state. `miso_init(S = NULL)` may still discover
supports, but callers must explicitly choose the motif count for fitting.

Fresh fits default to `warm_up_iters = 20L`: fixed-dictionary Gauss-Seidel sweeps,
followed by ordinary Jacobi updates with dictionary learning. The full variational
state is retained at the transition. `warm_up_iters = 0L` restores immediate Jacobi.
Supplying a fitted `miso_fit` as `init` automatically skips warm-up when the argument
is omitted; an explicit integer overrides this. Paths warm up the first fresh fit
only and use zero warm-up after merges. `update_F = FALSE` keeps the dictionary
fixed throughout both stages; `update_gamma = FALSE` also applies to warm-up.

`max_iters = 500` caps the main stage **in addition to** warm-up. The default
`tol = 1e-8` avoids the premature stopping observed with the older 1e-5 default. If that cap is reached without convergence, the fit is
returned with `converged = FALSE`, `stop_reason = "max_iters"`, and a warning.
Use `tol = 0` for intentional fixed-iteration runs (no convergence warning).
For ELBO comparison use converged multiple starts and tighter controls, e.g.
`tol = 1e-8, min_iters = 30, patience = 8, max_iters = 2000`.

`n_iter` counts the combined trajectory; `warm_up_iters` and `n_iter_main` give
stage counts. `elbo_stage` marks the stages; `warm_up_elbo` is the warm-up prefix.
The ELBO plot marks its end. A resumed call returns that call's trajectory only.
`summary()` reports stage counts, stopping reason, and quadrature discrepancies.

Legacy conversion preserves predictions and labels the result `legacy_eb`.
It does not create a Bayesian posterior: refit with explicit `S` to obtain one.
Retained legacy uncertainty is not reported as Bayesian population uncertainty.
