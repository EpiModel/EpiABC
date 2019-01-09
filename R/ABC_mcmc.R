## FUNCTION ABC_mcmc: ABC coupled to MCMC (Marjoram et al. 2003, Wegmann et al.
## 2009)


#' Coupled to MCMC schemes for ABC
#' 
#' This function implements three different algorithms to perform coupled to
#' MCMC ABC.
#' 
#' See the package's vignette for details on ABC-MCMC.
#' 
#' @param method a character string indicating the ABC-MCMC algorithm to be
#' used. Possible values are \code{"Marjoram_original"}, \code{"Marjoram"} and
#' \code{"Wegmann"}.  Note that the method \code{"Marjoram_original"} cannot be
#' used with multiple cores.
#' @param model a \code{R} function implementing the model to be simulated. It
#' must take as arguments a vector of model parameter values and it must return
#' a vector of summary statistics. When using the option \code{use_seed=TRUE},
#' \code{model} must take as arguments a vector containing a seed value and the
#' model parameter values.  A tutorial is provided in the package's vignette to
#' dynamically link a binary code to a \code{R} function.  Users may
#' alternatively wish to wrap their binary executables using the provided
#' functions \code{\link{binary_model}} and \code{\link{binary_model_cluster}}.
#' The use of these functions is associated with slightly different constraints
#' on the design of the binary code (see \code{\link{binary_model}} and
#' \code{\link{binary_model_cluster}}).
#' @param prior a list of prior information. Each element of the list
#' corresponds to a model parameter. The list element must be a vector whose
#' first argument determines the type of prior distribution: possible values
#' are \code{"unif"} for a uniform distribution on a segment, \code{"normal"}
#' for a normal distribution, \code{"lognormal"} for a lognormal distribution
#' or \code{"exponential"} for an exponential distribution.  The following
#' arguments of the list elements contain the characteritiscs of the prior
#' distribution chosen: for \code{"unif"}, two numbers must be given: the
#' minimum and maximum values of the uniform distribution; for \code{"normal"},
#' two numbers must be given: the mean and standard deviation of the normal
#' distribution; for \code{"lognormal"}, two numbers must be given: the mean
#' and standard deviation on the log scale of the lognormal distribution; for
#' \code{"exponential"}, one number must be given: the rate of the exponential
#' distribution. User-defined prior distributions can also be provided. See the
#' vignette for additional information on this topic.
#' @param summary_stat_target a vector containing the targeted (observed)
#' summary statistics.
#' @param prior_test a string expressing the constraints between model
#' parameters.  This expression will be evaluated as a logical expression, you
#' can use all the logical operators including \code{"<"}, \code{">"}, \ldots{}
#' Each parameter should be designated with \code{"X1"}, \code{"X2"}, \ldots{}
#' in the same order as in the prior definition.  If not provided, no
#' constraint will be applied.
#' @param n_rec a positive integer equal to the desired number of sampled
#' points along the MCMC.
#' @param n_between_sampling a positive integer equal to the desired spacing
#' between sampled points along the MCMC.
#' @param n_cluster a positive integer. If larger than 1 (the default value),
#' \code{ABC_mcmc} will launch \code{model} simulations in parallel on
#' \code{n_cluster} cores of the computer.
#' @param use_seed logical. If \code{FALSE} (default), \code{ABC_mcmc} provides
#' as input to the function \code{model} a vector containing the model
#' parameters used for the simulation.  If \code{TRUE}, \code{ABC_mcmc}
#' provides as input to the function \code{model} a vector containing an
#' integer seed value and the model parameters used for the simulation.  In
#' this last case, the seed value should be used by \code{model} to initialize
#' its pseudo-random number generators (if \code{model} is stochastic).
#' @param verbose logical. \code{FALSE} by default. If \code{TRUE},
#' \code{ABC_mcmc} writes in the current directory intermediary results at the
#' end of each step of the algorithm in the file "output_mcmc".  These outputs
#' have a matrix format, in wich each raw is a different simulation, the first
#' columns are the parameters used for this simulation, the following columns
#' are the summary statistics of this simulation, and the last column is the
#' distance between the simulation and the data.
#' @param dist_weights a vector containing the weights to apply to the distance
#' between the computed and the targeted statistics.  These weights can be used
#' to give more importance to a summary statistisc for example. The weights
#' will be normalized before applying them. If not provided, no weights will be
#' applied.
#' @param \dots Additional arguments can be passed depending on the choosen
#' method (see below)
#' @return The returned value is a list containing the following components:
#' \item{param}{ The model parameters used in the \code{model} simulations.}
#' \item{stats}{ The summary statistics obtained at the end of the \code{model}
#' simulations.} \item{dist}{ The distance of the simulations to the data.}
#' \item{stats_normalization}{ The standard deviation of the summary statistics
#' across the \code{model} simulations of the initial step.  These values are
#' used to normalize the summary statistics before the computation of the
#' Euclidean distance between simulations and data.  If \code{method} is
#' \code{"Marjoram_original"}, this is equal to \code{tab_normalization}.  If
#' \code{method} is \code{"Wegmann"}, this is not provided.} \item{epsilon}{
#' The final maximal distance between simulations and data in the retained
#' sample of particles.} \item{nsim}{ The number of \code{model} simulations
#' performed.} \item{n_between_sampling}{ The spacing between two sampled
#' points in the MCMC.} \item{computime}{ The computing time to perform the
#' simulations.} \item{min_stats}{ The minimal values of each summary
#' statistics during the calibration step, given when \code{method} is
#' \code{"Wegmann"}.} \item{max_stats}{ The maximal values of each summary
#' statistics during the calibration step, given when \code{method} is
#' \code{"Wegmann"}.} \item{lambda}{ The lambda values of the Box-Cox
#' transformation, given when \code{method} is \code{"Wegmann"}.}
#' \item{geometric_mean}{ The geometric means, given when \code{method} is
#' \code{"Wegmann"}.} \item{boxcox_mean}{ The means of Box-Cox transforms,
#' given when \code{method} is \code{"Wegmann"}.} \item{boxcox_sd}{ The
#' standard deviations of Box-Cox transforms, given when \code{method} is
#' \code{"Wegmann"}.} \item{pls_transform}{ The matrix of PLS transformation,
#' given when \code{method} is \code{"Wegmann"}.} \item{numcomp}{ The number of
#' used components for the PLS transformation, given when \code{method} is
#' \code{"Wegmann"}.}
#' @section Additional parameters: Depending on the choosen method, you can
#' specify the following arguments: \describe{ \item{dist_max}{ a positive
#' number, used when \code{method} is \code{"Marjoram_original"}.  This is the
#' tolerance threshold used during the MCMC.  If not provided by the user, it
#' is automatically computed as half the distance between the first simulation
#' and the target summary statistics and a warning is printed.}
#' \item{tab_normalization}{ a vector of the same length as
#' \code{summary_stat_target}, used when \code{method} is
#' \code{"Marjoram_original"}.  Each element contains a positive number by
#' which each summary statistics must be divided before the computation of the
#' Euclidean distance between simulations and data.  If not provided by the
#' user, the simulated summary statistics are divided by the target summary
#' statistics and a warning is printed.} \item{proposal_range}{ a vector of the
#' same length as the number of model parameters, used when \code{method} is
#' \code{"Marjoram_original"}.  Each element contains a positive number
#' defining the range of MCMC jumps for each model parameter.  If not provided
#' by the user, a default value is used for each parameter and a warning is
#' printed. The default value is 1/50 of the prior range for uniform
#' distributions, 1/20 of the standard deviation of the prior distribution for
#' normal distributions, 1/20 * exp ( sigma * sigma } for lognormal
#' distributions where sigma is the standard deviation of the prior
#' distribution in the log scale, and 1/20 of the inverse of the rate for
#' exponential distributions.  \item{n_calibration}{ a positive integer, used
#' when \code{method} is \code{"Marjoram"} or \code{"Wegmann"}.  This is the
#' number of simulations performed during the calibration step.  Default value
#' is 10000.} \item{tolerance_quantile}{ a positive number between 0 and 1
#' (strictly), used when \code{method} is \code{"Marjoram"} or
#' \code{"Wegmann"}.  This is the percentage of simulations retained during the
#' calibration step to determine the tolerance threshold to be used during the
#' MCMC.  Default value is 0.01.} \item{proposal_phi}{ a positive number, used
#' when \code{method} is \code{"Marjoram"} or \code{"Wegmann"}.  This is a
#' scaling factor defining the range of MCMC jumps.  Default value is 1.}
#' \item{numcomp}{ a positive integer, used when \code{method} is
#' \code{"Wegmann"}.  This is the number of components to be used for PLS
#' transformations.  Default value is 0 which encodes that this number is equal
#' to the number of summary statistics.} \item{seed_count}{ a positive integer,
#' the initial seed value provided to the function \code{model} (if
#' \code{use_seed=TRUE}). This value is incremented by 1 at each call of the
#' function \code{model}.} \item{progress_bar}{ logical, \code{FALSE} by
#' default. If \code{TRUE}, \code{ABC_mcmc} will output a bar of progression
#' with the estimated remaining computing time. Option not available with
#' multiple cores.  } \item{max_pick}{ a positive number, the max number of
#' fails when moving particle inside the prior. Enabled only if inside_prior is
#' to \code{TRUE}.  \code{10000} by default.  } }
#' @author Franck Jabot, Thierry Faure and Nicolas Dumoulin
#' @seealso \code{\link{binary_model}}, \code{\link{binary_model_cluster}},
#' \code{\link{ABC_rejection}}, \code{\link{ABC_emulation}},
#' \code{\link{ABC_sequential}}
#' @references Marjoram, P., Molitor, J., Plagnol, V. and Tavar\'e, S. (2003)
#' Markov chain Monte Carlo without likelihoods. \emph{PNAS}, \bold{100},
#' 15324--15328.
#' 
#' Wegmann, D., Leuenberger, C. and Excoffier, L. (2009) Efficient approximate
#' Bayesian computation coupled with Markov chain Monte Carlo without
#' likelihood. \emph{Genetics}, \bold{182}, 1207-1218.
#' @keywords abc model inference mcmc
#' @examples
#' 
#'  \dontrun{ 
#'     ##### EXAMPLE 1 #####
#'     #####################
#' 
#'     ## the model has two parameters and outputs two summary statistics.
#'     ## defining a simple toy model:
#'     toy_model<-function(x){ c( x[1] + x[2] + rnorm(1,0,0.1) , x[1] * x[2] + rnorm(1,0,0.1) ) }
#' 
#'     ## define prior information
#'     toy_prior=list(c("unif",0,1),c("normal",1,2))
#'     # a uniform prior distribution between 0 and 1 for parameter 1, and a normal distribution
#'     # of mean 1 and standard deviation of 2 for parameter 2.
#' 
#'     ## define the targeted summary statistics
#'     sum_stat_obs=c(1.5,0.5)
#' 
#'     ## to perform the Marjoram et al. (2003)'s method:
#'     ##
#'     ABC_Marjoram_original<-ABC_mcmc(method="Marjoram_original", model=toy_model, prior=toy_prior,
#'       summary_stat_target=sum_stat_obs)
#'     ABC_Marjoram_original
#' 
#'     ## artificial example to perform the Marjoram et al. (2003)'s method, with modifications
#'     # drawn from Wegmann et al. (2009) without Box-Cox and PLS transformations.
#'     ##
#'     ABC_Marjoram<-ABC_mcmc(method="Marjoram", model=toy_model, prior=toy_prior,
#'       summary_stat_target=sum_stat_obs)
#'     ABC_Marjoram
#' 
#' 
#'     ## artificial example to perform the Wegmann et al. (2009)'s method.
#'     ##
#'     ABC_Wegmann<-ABC_mcmc(method="Wegmann", model=toy_model, prior=toy_prior,
#'       summary_stat_target=sum_stat_obs)
#'     ABC_Wegmann
#' 
#' 
#'     ##### EXAMPLE 2 #####
#'     #####################
#' 
#'     ## this time, the model is a C++ function packed into a R function -- this time,
#'     # the option 'use_seed' must be turned to TRUE.
#' 
#'     ## define prior information
#'     trait_prior=list(c("unif",3,5),c("unif",-2.3,1.6),c("unif",-25,125),c("unif",-0.7,3.2))
#'     trait_prior
#' 
#'     ## define the targeted summary statistics
#'     sum_stat_obs=c(100,2.5,20,30000)
#' 
#' 
#'     ## artificial example to perform the Marjoram et al. (2003)'s method.
#'     ##
#'     n=10
#'     ABC_Marjoram_original<-ABC_mcmc(method="Marjoram_original", model=trait_model,
#'     prior=trait_prior, summary_stat_target=sum_stat_obs, n_rec=n, use_seed=TRUE)
#'     ABC_Marjoram_original
#' 
#'     ## artificial example to perform the Marjoram et al. (2003)'s method, with modifications
#'     # drawn from Wegmann et al. (2009) without Box-Cox and PLS transformations.
#'     ##
#'     n=10
#'     n_calib=10
#'     tol_quant=0.2 
#'     ABC_Marjoram<-ABC_mcmc(method="Marjoram", model=trait_model, prior=trait_prior,
#'       summary_stat_target=sum_stat_obs, n_rec=n, n_calibration=n_calib,
#'       tolerance_quantile=tol_quant, use_seed=TRUE)
#'     ABC_Marjoram
#' 
#' 
#'     ## artificial example to perform the Wegmann et al. (2009)'s method.
#'     ##
#'     n=10
#'     n_calib=10
#'     tol_quant=0.2 
#'     ABC_Wegmann<-ABC_mcmc(method="Wegmann", model=trait_model, prior=trait_prior,
#'         summary_stat_target=sum_stat_obs, n_rec=n, n_calibration=n_calib,
#'         tolerance_quantile=tol_quant, use_seed=TRUE)
#'     ABC_Wegmann
#'  }%dontrun
#' 
#' @export ABC_mcmc
ABC_mcmc <- function(method, model, prior, summary_stat_target, prior_test = NULL, 
    n_rec = 100, n_between_sampling = 10, n_cluster = 1, use_seed = FALSE, verbose = FALSE, 
    dist_weights=NULL, ...) {
    ## checking errors in the inputs
    if (missing(method)) 
        stop("'method' is missing")
    if (missing(model)) 
        stop("'model' is missing")
    if (missing(prior)) 
        stop("'prior' is missing")
    if (!is.null(prior_test)) 
        .check_prior_test(length(prior), prior_test)
    data = .wrap_constants_in_model(prior, model, use_seed)
    prior = data$new_prior
    model = data$new_model
    prior = .process_prior(prior)
    if (missing(summary_stat_target)) 
        stop("'summary_stat_target' is missing")
    if (!any(method == c("Marjoram_original", "Marjoram", "Wegmann"))) {
        stop("Method must be Marjoram_original, Marjoram or wegmann")
    }
    if (!is.vector(summary_stat_target)) 
        stop("'summary_stat_target' has to be a vector.")
    if (!is.vector(n_cluster)) 
        stop("'n_cluster' has to be a number.")
    if (length(n_cluster) > 1) 
        stop("'n_cluster' has to be a number.")
    if (n_cluster < 1) 
        stop("'n_cluster' has to be a positive number.")
    n_cluster = floor(n_cluster)
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean")
    if (!is.logical(verbose)) 
        stop("'verbose' has to be boolean")
    if (!is.null(dist_weights) && length(dist_weights)!=length(summary_stat_target)) {
        stop("'dist_weights' has to be the same length than 'summary_stat_target'")
    }
    mcmc = NULL
    if (n_cluster == 1) {
        mcmc = .ABC_mcmc_internal(method, model, prior, prior_test, n_rec, n_between_sampling, 
            summary_stat_target, use_seed, verbose, dist_weights=dist_weights, ...)
    } else {
        if (use_seed == FALSE) {
            stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly - see the package's vignette for more details.")
        }
        mcmc = .ABC_mcmc_cluster(method, model, prior, prior_test, n_rec, n_between_sampling, 
            summary_stat_target, n_cluster, use_seed, verbose, dist_weights=dist_weights, ...)
    }
    mcmc
} 
