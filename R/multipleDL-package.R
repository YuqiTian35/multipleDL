#' The 'multipleDL' package.
#' @title Address Detection Limits by Cumulative Probability Models
#'
#' @description The package allows fitting regression models on continuous/ordinal response data subject to detection limits (DLs)
#' based on cumulative probability models (CPMs). Both single and multiple DLs can be handled. Conditional quantiles and CDFs
#' (cumulative distribution functions) can obtained from fitted models.
#'
#'
#'
#' @docType package
#' @name multipleDL-package
#' @useDynLib multipleDL, .registration = TRUE
#' @import methods stats Rcpp rstantools
#' @importFrom rstan sampling
#' @importFrom RcppParallel CxxFlags RcppParallelLibs
#' @importFrom SparseM t as.matrix solve
#'
#' @references
#' Stan Development Team (2020). RSroxygen2::roxygenize()tan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#' Harrell, F. (2020). rms: Regression modeling strategies. R package version 6.1.0. https://CRAN.R-project.org/package=rms
#' Tian et al. "Addressing detection limits by semiparametric cumulative probability models." (2022) (to be submitted)
#'
NULL
