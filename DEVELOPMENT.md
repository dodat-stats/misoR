# Developing misoR

Version 0.2.0 retains joint Bayesian population updates in `R/population_posterior.R`.
The core fitter shares blocked kernels between ordinary and length models. Do not
replace joint moments such as E(alpha log beta) or E(alpha/beta) by products or ratios
of separate posterior means.

Documentation is generated from roxygen comments. Package tests use base R.

```sh
Rscript -e 'roxygen2::roxygenise(".")'
R CMD build .
R CMD check --no-manual misoR_0.2.0.tar.gz
```

The prior implementation is preserved outside the package in a source snapshot;
new Bayesian equations replace its EB-specific equation tests. Initializer and graph
regressions continue to run. An additional repository-level validation compares the
package against the independent standalone research implementation over 35 iterations.

Fixed-F GS warm-up uses the same prepared data and population blocks as the main
Jacobi engine. Warm-up is not repeated automatically on resumes or after merges.
Regression tests with warm_up_iters=0 retain their original equation checks.
