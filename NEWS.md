# misoR 0.2.0

- Loading plots accept `sort_by = "dominant"` to order observations within each
  motif by decreasing fraction of its largest-mean factor.
- Bayesian mixture inference uses joint population posteriors for alpha and beta;
  the normalized dictionary F remains point-estimated.
- New `warm_up_iters` control: 20 fixed-F Gauss-Seidel sweeps for fresh fits,
  followed by Jacobi. Fitted-state resumes skip warm-up unless explicitly asked;
  descending paths warm up only their first fit.
- `S` is now required in both fitters, including when supplying `init`. Existing
  calls must add their desired motif count, matching the supplied state.
- The default `tol` is tightened to 1e-8; the former 1e-5 default stopped too
  early on the BBC and cancer validation data.
- The default `max_iters` is 500 main iterations, excluding warm-up. Fits warn if
  this cap is reached without convergence. `tol = 0` intentionally runs a fixed
  number of iterations and does not warn.
- Continuing a normalized Bayesian fit preserves its saved probabilities exactly,
  including values below the numerical floor; they are not floored a second time.
- Fits retain combined ELBO histories, stage labels, separate iteration counts,
  stopping reason, and fitting controls. ELBO plots mark the warm-up transition.
- Legacy conversion remains representation-only: predictions are preserved,
  inference is labeled legacy_eb, and legacy objectives are not comparable to
  Bayesian objectives for model selection. Refit with explicit S.
