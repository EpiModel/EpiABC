## FUNCTION ABC_rejection: brute-force ABC (Pritchard et al. 1999)


#' Rejection sampling scheme for ABC
#' 
#' This function launches a series of \code{nb_simul} model simulations with
#' model parameters drawn in the prior distribution specified in
#' \code{prior_matrix}.
#' 
#' 
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
#' @param nb_simul a positive integer equal to the desired number of
#' simulations of the model.
#' @param prior_test a string expressing the constraints between model
#' parameters.  This expression will be evaluated as a logical expression, you
#' can use all the logical operators including \code{"<"}, \code{">"}, \ldots{}
#' Each parameter should be designated with \code{"X1"}, \code{"X2"}, \ldots{}
#' in the same order as in the prior definition.  If not provided, no
#' constraint will be applied.
#' @param summary_stat_target a vector containing the targeted (observed)
#' summary statistics.  If not provided, \code{ABC_rejection} only launches the
#' simulations and outputs the simulation results.
#' @param tol tolerance, a strictly positive number (between 0 and 1)
#' indicating the proportion of simulations retained nearest the targeted
#' summary statistics.
#' @param use_seed logical. If \code{FALSE} (default), \code{ABC_rejection}
#' provides as input to the function \code{model} a vector containing the model
#' parameters used for the simulation.  If \code{TRUE}, \code{ABC_rejection}
#' provides as input to the function \code{model} a vector containing an
#' integer seed value and the model parameters used for the simulation.  In
#' this last case, the seed value should be used by \code{model} to initialize
#' its pseudo-random number generators (if \code{model} is stochastic).
#' @param seed_count a positive integer, the initial seed value provided to the
#' function \code{model} (if \code{use_seed=TRUE}). This value is incremented
#' by 1 at each call of the function \code{model}.
#' @param n_cluster a positive integer. If larger than 1 (the default value),
#' \code{ABC_rejection} will launch \code{model} simulations in parallel on
#' \code{n_cluster} cores of the computer.
#' @param verbose logical. \code{FALSE} by default. If \code{TRUE},
#' \code{ABC_rejection} writes in the current directory intermediary results at
#' the end of each step of the algorithm in the file "output".  These outputs
#' have a matrix format, in wich each raw is a different simulation, the first
#' columns are the parameters used for this simulation, and the last columns
#' are the summary statistics of this simulation.
#' @param progress_bar logical, \code{FALSE} by default. If \code{TRUE},
#' \code{ABC_rejection} will output a bar of progression with the estimated
#' remaining computing time. Option not available with multiple cores.
#' @return The returned value is a list containing the following components:
#' \item{param}{ The model parameters used in the \code{model} simulations.}
#' \item{stats}{ The summary statistics obtained at the end of the \code{model}
#' simulations.} \item{weights}{ The weights of the different \code{model}
#' simulations. In the standard rejection scheme, all \code{model} simulations
#' have the same weights.} \item{stats_normalization}{ The standard deviation
#' of the summary statistics across the \code{model} simulations.} \item{nsim}{
#' The number of \code{model} simulations performed.} \item{nrec}{ The number
#' of retained simulations (if targeted summary statistics are provided).}
#' \item{computime}{ The computing time to perform the simulations.}
#' @author Franck Jabot, Thierry Faure and Nicolas Dumoulin
#' @seealso \code{\link{binary_model}}, \code{\link{binary_model_cluster}},
#' \code{\link{ABC_sequential}}, \code{\link{ABC_mcmc}}
#' @references Pritchard, J.K., and M.T. Seielstad and A. Perez-Lezaun and M.W.
#' Feldman (1999) Population growth of human Y chromosomes: a study of Y
#' chromosome microsatellites. \emph{Molecular Biology and Evolution},
#' \bold{16}, 1791--1798.
#' @keywords abc model inference
#' @examples
#' 
#'     ##### EXAMPLE 1 #####
#'     #####################
#'     set.seed(1)
#' 
#'     ## artificial example to show how to use the 'ABC_rejection' function.
#'     ## defining a simple toy model:
#'     toy_model<-function(x){ 2 * x + 5 + rnorm(1,0,0.1) }
#' 
#'     ## define prior information
#'     toy_prior=list(c("unif",0,1)) # a uniform prior distribution between 0 and 1
#' 
#'     ## only launching simulations with parameters drawn in the prior distributions
#'     set.seed(1)
#'     n=10
#'     ABC_sim<-ABC_rejection(model=toy_model, prior=toy_prior, nb_simul=n)
#'     ABC_sim
#' 
#'     ## launching simulations with parameters drawn in the prior distributions
#'     # and performing the rejection step
#'     sum_stat_obs=6.5
#'     tolerance=0.2
#'     ABC_rej<-ABC_rejection(model=toy_model, prior=toy_prior, nb_simul=n,
#'       summary_stat_target=sum_stat_obs, tol=tolerance)
#' 
#'     ## NB: see the package's vignette to see how to pipeline 'ABC_rejection' with the function
#'     # 'abc' of the package 'abc' to perform other rejection schemes.
#'  \dontrun{ 
#' 
#'     ##### EXAMPLE 2 #####
#'     #####################
#' 
#'     ## this time, the model has two parameters and outputs two summary statistics.
#'     ## defining a simple toy model:
#'     toy_model2<-function(x){ c( x[1] + x[2] + rnorm(1,0,0.1) , x[1] * x[2] + rnorm(1,0,0.1) ) }
#' 
#'     ## define prior information
#'     toy_prior2=list(c("unif",0,1),c("normal",1,2))
#'     # a uniform prior distribution between 0 and 1 for parameter 1, and a normal distribution
#'     # of mean 1 and standard deviation of 2 for parameter 2.
#' 
#'     ## only launching simulations with parameters drawn in the prior distributions
#'     set.seed(1)
#'     n=10
#'     ABC_sim<-ABC_rejection(model=toy_model2, prior=toy_prior2, nb_simul=n)
#'     ABC_sim
#' 
#'     ## launching simulations with parameters drawn in the prior distributions
#'     # and performing the rejection step
#'     sum_stat_obs2=c(1.5,0.5)
#'     tolerance=0.2
#'     ABC_rej<-ABC_rejection(model=toy_model2, prior=toy_prior2, nb_simul=n,
#'       summary_stat_target=sum_stat_obs2, tol=tolerance)
#' 
#'     ## NB: see the package's vignette to see how to pipeline 'ABC_rejection' with the function
#'     # 'abc' of the package 'abc' to perform other rejection schemes.
#' 
#' 
#'     ##### EXAMPLE 3 #####
#'     #####################
#' 
#'     ## this time, the model is a C++ function packed into a R function -- this time, the option
#'     # 'use_seed' must be turned to TRUE.
#'     n=10
#'     trait_prior=list(c("unif",3,5),c("unif",-2.3,1.6),c("unif",-25,125),c("unif",-0.7,3.2))
#'     trait_prior
#' 
#'     ## only launching simulations with parameters drawn in the prior distributions
#'     ABC_sim<-ABC_rejection(model=trait_model, prior=trait_prior, nb_simul=n, use_seed=TRUE)
#'     ABC_sim
#' 
#'     ## launching simulations with parameters drawn in the prior distributions and performing
#'     # the rejection step
#'     sum_stat_obs=c(100,2.5,20,30000)
#'     tolerance=0.2
#'     ABC_rej<-ABC_rejection(model=trait_model, prior=trait_prior, nb_simul=n,
#'       summary_stat_target=sum_stat_obs, tol=tolerance, use_seed=TRUE)
#' 
#'     ## NB: see the package's vignette to see how to pipeline 'ABC_rejection' with the function
#'     # 'abc' of the package 'abc' to perform other rejection schemes.
#' 
#' 
#'     ##### EXAMPLE 4 - Parallel implementations #####
#'     ################################################
#' 
#'     ## NB: the option use_seed must be turned to TRUE.
#' 
#'     ## For models already running with the option use_seed=TRUE, simply change
#'     # the value of n_cluster:
#'     sum_stat_obs=c(100,2.5,20,30000)
#'     ABC_simb<-ABC_rejection(model=trait_model, prior=trait_prior, nb_simul=n,
#'       use_seed=TRUE, n_cluster=2)
#' 
#'     ## For other models, change the value of n_cluster and modify the model so that the first
#'     # parameter becomes a seed information value:
#'     toy_model_parallel<-function(x){ 
#' 	set.seed(x[1])
#' 	2 * x[2] + 5 + rnorm(1,0,0.1) }
#'     sum_stat_obs=6.5
#' 
#'     ABC_simb<-ABC_rejection(model=toy_model_parallel, prior=toy_prior, nb_simul=n,
#'       use_seed=TRUE, n_cluster=2)
#'  }%dontrun
#' 
#' @export ABC_rejection
ABC_rejection <- function(model, prior, nb_simul, prior_test = NULL, summary_stat_target = NULL, 
    tol = NULL, use_seed = FALSE, seed_count = 0, n_cluster = 1, verbose = FALSE, 
    progress_bar = FALSE) {
    ## checking errors in the inputs
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
    if (missing(nb_simul)) 
        stop("'nb_simul' is missing")
    if (nb_simul < 1) 
        stop("'nb_simul' must be a number larger than 1")
    if ((!is.null(summary_stat_target)) && (!is.vector(summary_stat_target))) {
        stop("'summary_stat_target' has to be a number")
    }
    if ((!is.null(tol)) && (!is.vector(tol))) {
        stop("'tol' has to be a number")
    }
    if (!is.logical(use_seed)) 
        stop("'use_seed' has to be boolean")
    if (!is.vector(seed_count)) 
        stop("'seed_count' has to be a number")
    if (length(seed_count) > 1) 
        stop("'seed_count' has to be a number")
    if (seed_count < 0) 
        stop("'seed_count' has to be a positive number")
    if (!is.vector(n_cluster)) 
        stop("'n_cluster' has to be a number.")
    if (length(n_cluster) > 1) 
        stop("'n_cluster' has to be a number.")
    if (n_cluster < 1) 
        stop("'n_cluster' has to be a positive number.")
    n_cluster = floor(n_cluster)
    if (!is.logical(verbose)) 
        stop("'verbose' has to be boolean")
    if (!is.logical(progress_bar)) 
        stop("'progress_bar' has to be boolean")
    nb_simul = floor(nb_simul)
    seed_count = floor(seed_count)
    rejection = NULL
    if ((!is.null(summary_stat_target)) && (is.null(tol))) {
        stop("'tol' is missing")
    }
    if (n_cluster == 1) {
        rejection = .ABC_rejection(model, prior, prior_test, nb_simul, use_seed, 
            seed_count, verbose, progress_bar)
    } else {
        if (use_seed == FALSE) {
            stop("For parallel implementations, you must specify the option 'use_seed=TRUE' and modify your model accordingly - see the package's vignette for more details.")
        }
        rejection = .ABC_rejection_cluster(model, prior, prior_test, nb_simul, seed_count, 
            n_cluster, verbose)
    }
    res = NULL
    if (is.null(summary_stat_target)) {
        res = list(param = rejection$param, stats = rejection$stats, weights = rejection$weights, 
            stats_normalization = rejection$stats_normalization, nsim = rejection$nsim, 
            computime = rejection$computime)
    } else {
        options(warn = -1)
        rej = abc(summary_stat_target, rejection$param, rejection$stats, tol, method = "rejection")
        options(warn = 0)
        nr = dim(rej$unadj.values)[1]
        res = list(param = rej$unadj.values, stats = rej$ss, weights = array(1/nr, 
            nr), stats_normalization = rejection$stats_normalization, nsim = rejection$nsim, 
            nrec = nr, computime = rejection$computime)
    }
    res
} 
