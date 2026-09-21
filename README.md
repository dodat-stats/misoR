# misoR

MiSo fits a mixture of sparse Poisson factor models to a count matrix with observations in rows and features in columns. The package supports dense matrices and `Matrix` sparse matrices. The implementation is pure R.

## Install locally

From this directory in a terminal:

```sh
R CMD INSTALL .
```

`Matrix` must be installed. In R, use `install.packages("Matrix")` if needed.

## Fit and inspect

```r
library(misoR)

# Learn a dictionary with K factors, using D slots per motif.
fit <- miso_fit(Y, K = 6, D = 2)

# Alternatively, initialize with a supplied K-by-M dictionary.
fit <- miso_fit(Y, F = F, D = 2)
# F is refined by default; set update_F = FALSE to keep it fixed.

print(fit)
summary(fit)
L <- predict(fit, type = "loadings")  # N by K, posterior mean loadings
mu <- predict(fit, type = "mean")    # N by M, equal to L %*% fit$F
plot(fit, type = "elbo")
plot(fit, type = "loadings")
plot(fit, type = "loadings", normalize = FALSE)
```

Prediction currently describes fitted observations. Both prediction outputs are dense, including for sparse count inputs. The loading plot uses posterior means averaged over mixture uncertainty; hard cluster assignments only group the bars. `cluster_order`, `factor_order`, and `sort_by` control the display without changing the fit. No simulation truth or component merging is used.

For a runnable simulation, see `inst/examples/small-fit.R` or run
`source(system.file("examples", "small-fit.R", package = "misoR"))`.

## Fit document counts with document-length adjustment

Use raw document-by-word counts, after selecting the vocabulary. Remove empty
documents before fitting with the default `n = rowSums(Y)`.

```r
fit <- miso_fit_length(Y, K = 6, D = 2, keep_initialization = TRUE)
# Or supply an initial dictionary: miso_fit_length(Y, F = F, D = 2)
L <- predict(fit, "loadings")  # Contributions per unit of document length
theta <- L / rowSums(L)        # Descriptive topic proportions
counts <- predict(fit)        # fit$n * (L %*% fit$F)
plot(fit, type = "loadings")

# Reuse the initialization on the same counts and lengths.
fit2 <- miso_fit_length(Y, init = fit$initialization)
```

This fits the Poisson working likelihood `Y[i,j] ~ Poisson(n[i] * (L %*% F)[i,j])`.
The factor rows sum to one; the loading rows are unconstrained. For observed
document lengths, the likelihood decomposes into a multinomial composition
likelihood and a term encouraging total loadings of one. It is not the exact
conditional sampling law given the observed total.

Automatic initialization rescales row-wise Poisson-SuSiE posteriors before
estimating shared loading priors. Do not directly pass count-scale priors from
`miso_init()` to this function. Custom `init` priors must be per unit of `n`.
The returned `beta` has dimensions N by S by D, while `alpha0` and `beta0`
remain S by D. The retained posterior rates can differ from `beta0 + n`
after an empirical-Bayes update. `elbo` omits data-only constants; add
`fit$elbo_constant` to recover the full bound. Normalized posterior mean
loadings are descriptive proportions, not exact posterior means of proportions.
An explicit positive vector `n` is also accepted; unit lengths recover the
original likelihood. See `?miso_fit_length` for details.

## Reuse initialization

```r
initial <- miso_init(Y, F = F, D = 2, seed = 1)
initial$support
fit <- miso_fit(Y, init = initial)

# Cache Poisson SuSiE when exploring support cutoffs.
ps <- poisson_susie_fit(Y, F, D = 2)
initial90 <- miso_init(Y, poisson_susie = ps, tau = 0.90)
initial95 <- miso_init(Y, poisson_susie = ps, tau = 0.95)
fit90 <- miso_fit(Y, init = initial90)
fit95 <- miso_fit(Y, init = initial95)
```

Cached results retain the normalized dictionary. They must refer to the same observations in the same row order. A reusable initialization stores its Dirichlet prior; `miso_fit()` uses that prior unless `phi0` is explicitly supplied.

## Algorithm and returned parameters

Both `miso_fit()` and `miso_fit_length()` use the same optimization loop. Each iteration follows SuSiE → factor → mixture → ELBO. The SuSiE step uses one Jacobi allocation pass to accumulate loading counts and `C`, then updates the loading posterior, gamma, and loading priors. The factor step updates F once using `C` and the updated gamma. The mixture step evaluates bounds in a second feature pass, then refines responsibilities and Dirichlet parameters. No damping is used. Loading priors, gamma, and F are updated by default.

Length adjustment changes initialization, posterior rates, and the Poisson-rate term, while retaining the same update schedule. Both fitters accept finite nonnegative real-valued matrices, including fractional entries, without rounding. For fractional data, the optimized criterion is the generalized Poisson objective.

The fit is a documented list. Principal fields are:

| Field | Dimensions / meaning |
|---|---|
| `F` | K by M, row-normalized factor profiles |
| `gamma` | S by D by K, factor-selection probabilities |
| `alpha` | N by S by D, posterior loading shapes |
| `beta` | S by D, posterior loading rates |
| `alpha0`, `beta0` | S by D, loading prior shapes and rates |
| `omega` | N by S, posterior motif responsibilities |
| `phi0`, `phi`, `pi` | Length S, Dirichlet prior, posterior, and mean |
| `allocation_mass` | Length S, mean responsibilities |
| `z_hat` | Length N, hard cluster assignments |
| `elbo`, `final_elbo` | Convergence history and last recorded value |

Allocation variables xi are never stored by MiSo. See `?miso_fit`, `?miso_init`, `?poisson_susie_fit`, and `?plot.miso_fit` for complete argument and return-value documentation.

## Development and publishing

See [DEVELOPMENT.md](DEVELOPMENT.md) for documentation and package checks, and [PUBLISHING.md](PUBLISHING.md) for private GitHub setup. This initial version carries a private-use license notice; no open-source license has been selected.
