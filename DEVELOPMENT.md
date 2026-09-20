# Developing misoR

This directory is a standalone package. It needs no files from `../code`, the manuscript, or the simulation results. The inference equations were copied from the tested two-pass implementation; further package development belongs in `R/`. The older research scripts continue using their existing source files.

Documentation is written in roxygen comments in `R/`; generated `man/*.Rd` and `NAMESPACE` are committed so installation needs no documentation tools.

```sh
# From misoR/, regenerate documentation when editing the public API.
Rscript -e 'roxygen2::roxygenise(".")'

# From the parent directory, build and check an isolated source package.
R CMD build misoR
R CMD check --no-manual misoR_0.0.1.tar.gz
R CMD INSTALL misoR_0.0.1.tar.gz
```

`--no-manual` skips PDF-manual compilation and its LaTeX dependency. It still checks the help pages and runs examples and tests. To compile a PDF manual as well, omit that option when LaTeX is available.

Tests use base R without testthat. The numerical reference suite checks the independent update equations, ELBO ascent, exact feature-pass counts, dense/sparse agreement, and edge cases. The package API suite checks reusable initialization, dictionary consistency, S3 dispatch, summaries, loading plots, and singleton dimensions.
