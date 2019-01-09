
#' EasyABC: performing efficient approximate Bayesian computation sampling
#' schemes using R
#' 
#' The package EasyABC enables to launch a series of simulations of a computer
#' code from the R platform, and to retrieve the simulation outputs in an
#' appropriate format for post-processing treatments. Four sequential sampling
#' schemes, three coupled-to-MCMC schemes and a Simulated Annealing scheme are
#' implemented. EasyABC further enables to launch the simulations in parallel
#' on multiple cores of a computer.
#' 
#' \tabular{ll}{ Package: EasyABC Type: Package Version: 1.5 Date: 2015-06-30
#' License: GPL-3 LazyLoad: yes }
#' 
#' @name EasyABC-package
#' @aliases EasyABC-package EasyABC
#' @docType package
#' @author Franck Jabot, Thierry Faure, Nicolas Dumoulin, Carlo Albert
#' @seealso \code{\link{ABC_rejection}}, \code{\link{ABC_sequential}},
#' \code{\link{ABC_mcmc}}, \code{\link{SABC}}, \code{\link{binary_model}},
#' \code{\link{binary_model_cluster}}
#' @keywords package abc
#' 
#' @import MASS mnormt parallel pls abc lhs tensorA
#' @importFrom stats as.formula cov cov.wt lm optimize rnorm runif sd uniroot var
#' @importFrom utils flush.console read.table write.table
#' 
NULL
