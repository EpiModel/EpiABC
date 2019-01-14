
#' EpiABC: Approximate Bayesian Computation for EpiModel on High-Performance Computing Clusters
#'
#' Implementation of Approximate Bayesian Computation with Sequential
#' Monte Carlo Methods (ABC-SMC) optimized for use with EpiModel R package epidemic
#' simulations running on High-Performance Computing (HPC) cluster environments such
#' as MPI and Slurm.
#'
#' \tabular{ll}{
#'    Package: \tab EpiABC\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.0.0\cr
#'    Date: \tab 2019-01-14\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @name EpiABC-package
#' @aliases EpiABC-package EpiABC
#' @docType package
#'
#' @keywords package
#'
#' @import MASS mnormt parallel pls abc lhs tensorA
#' @importFrom stats as.formula cov cov.wt lm optimize rnorm runif sd uniroot var
#' @importFrom utils flush.console read.table write.table
#'
NULL
