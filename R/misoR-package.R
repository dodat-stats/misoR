#' Mixture of sparse Poisson factor models
#'
#' Fit MiSo using [miso_fit()], create reusable starting parameters with
#' [miso_init()], and fit row-wise Poisson SuSiE with [poisson_susie_fit()].
#' Use [miso_fit_length()] for document-length-adjusted Poisson means.
#' Observations are rows and features are columns throughout the package.
#'
#' @importFrom stats predict rexp rgamma
#' @importFrom utils head tail
#' @keywords internal
"_PACKAGE"
