library(misoR)

# Equation-level checks access private helpers, but production namespaces stay
# locked. Pass-count instrumentation uses isolated copies of the functions.
test_environment <- new.env(parent = asNamespace("misoR"))
test_environment$.miso_test_clone <- function() {
  ns <- asNamespace("misoR")
  copy <- new.env(parent = ns)
  for (name in ls(ns, all.names = TRUE)) {
    value <- get(name, envir = ns, inherits = FALSE)
    if (is.function(value) && identical(environment(value), ns)) {
      environment(value) <- copy
      assign(name, value, envir = copy)
    }
  }
  copy
}
for (file in c("test_core.R", "test_vectorization.R", "test_api.R",
               "test_lse_streaming.R", "test_schedule.R", "test_sparse.R")) {
  sys.source(file.path("reference", file), envir = test_environment)
}
cat("All numerical regression tests passed.\n")
