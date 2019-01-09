\name{ABC_sequential}
\alias{ABC_sequential}
\title{
	Sequential sampling schemes for ABC}
\description{
	This function implements four different algorithms to perform sequential sampling schemes for ABC.
	Sequential sampling schemes consist in sampling initially model parameters in the prior distribution, just like in a standard rejection-based ABC, in order to obtain a rough posterior distribution of parameter values, and in subsequently sampling close to this rough posterior distribution to refine it.
	Sequential sampling schemes have been shown to be more efficient than standard rejection-based procedures.
}
\usage{
	ABC_sequential(method, model, prior, nb_simul, summary_stat_target,
          prior_test=NULL, n_cluster = 1, use_seed = FALSE, verbose = FALSE,
          dist_weights=NULL, ...)
}

\arguments{
  \item{method}{
	a character string indicating the sequential algorithm to be used. Possible values are \code{"Beaumont"}, \code{"Drovandi"}, \code{"Delmoral"}, \code{"Lenormand"} and \code{"Emulation"}.}
  \item{model}{
    a \code{R} function implementing the model to be simulated. It must take as arguments a vector of model parameter values and it must return a vector of summary statistics. When using the option \code{use_seed=TRUE}, \code{model} must take as arguments a vector containing a seed value and the model parameter values.
    A tutorial is provided in the package's vignette to dynamically link a binary code to a \code{R} function.  Users may alternatively wish to wrap their binary executables using the provided functions \code{\link{binary_model}} and \code{\link{binary_model_cluster}}. The use of these functions is associated with slightly different constraints on the design of the binary code (see \code{\link{binary_model}} and \code{\link{binary_model_cluster}}).}
  \item{prior}{
    a list of prior information. Each element of the list corresponds to a model parameter. The list element must be a vector whose first argument determines the type of prior distribution: possible values are \code{"unif"} for a uniform distribution on a segment, \code{"normal"} for a normal distribution, \code{"lognormal"} for a lognormal distribution or \code{"exponential"} for an exponential distribution.
    The following arguments of the list elements contain the characteritiscs of the prior distribution chosen: for \code{"unif"}, two numbers must be given: the minimum and maximum values of the uniform distribution; for \code{"normal"}, two numbers must be given: the mean and standard deviation of the normal distribution; for \code{"lognormal"}, two numbers must be given: the mean and standard deviation on the log scale of the lognormal distribution; for \code{"exponential"}, one number must be given: the rate of the exponential distribution.
    Note that when using the method "Lenormand", solely uniform prior distributions are supported. User-defined prior distributions can also be provided. See the vignette for additional information on this topic.}
  \item{nb_simul}{
    	a positive integer equal to the desired number of simulations of the model below the tolerance threshold when \code{method} is \code{"Beaumont"}, \code{"Drovandi"} and \code{"Delmoral"}.
	When \code{method} is \code{"Lenormand"}, the number of simulations below the tolerance threshold is equal to \code{nb_simul * alpha}. See the package's vignette and Lenormand et al. (2012) for details.}
  \item{summary_stat_target}{
	a vector containing the targeted (observed) summary statistics.}
  \item{prior_test}{
        a string expressing the constraints between model parameters.
        This expression will be evaluated as a logical expression, you can use all the logical operators including \code{"<"}, \code{">"}, \ldots
        Each parameter should be designated with \code{"X1"}, \code{"X2"}, \ldots in the same order as in the prior definition.
        If not provided, no constraint will be applied.}
  \item{n_cluster}{
    a positive integer. If larger than 1 (the default value), \code{ABC_sequential} will launch \code{model} simulations in parallel on \code{n_cluster} cores of the computer.}
  \item{use_seed}{
    logical. If \code{FALSE} (default), \code{ABC_sequential} provides as input to the function \code{model} a vector containing the model parameters used for the simulation.
    If \code{TRUE}, \code{ABC_sequential} provides as input to the function \code{model} a vector containing an integer seed value and the model parameters used for the simulation.
    In this last case, the seed value should be used by \code{model} to initialize its pseudo-random number generators (if \code{model} is stochastic).}
  \item{verbose}{
    logical. \code{FALSE} by default. If \code{TRUE}, \code{ABC_sequential} writes in the current directory intermediary results at the end of each step of the algorithm various files.
    The file "n_simul_tot_step_iteration" (where iteration is the step number) contains the total number of simulations performed since the beginning of the algorithm at the end of the step "iteration".
    The file "R_step_iteration" (when using the method "Drovandi") is the parameter R used during the step "iteration" (see Drovandi and Pettitt 2011 for details).
    The file "p_acc_iteration" (when using the method "Lenormand") is the parameter p_acc computed at the end of the step "iteration" (see Lenormand et al. 2012 for details).
    The file "tolerance_step_iteration" (when using the method "Drovandi", "Delmoral" or "Lenormand") is the tolerance computed at the end of the step "iteration".
    The file "output_step_iteration" gives the simulations kept after each iteration and has a matrix format, in wich each row is a different simulation, the first column is the weight of the simulation, the following columns are the parameters used for this simulation, and the last columns are the summary statistics of this simulation.
    The file "model_step_iteration" gives the simulations performed at each iteration and has a matrix format, in which each row is a different simulation, the first column is the weight of the simulation, the following columns are the parameters used for this simulation, and the last columns are the summary statistics of this simulation.
    All these informations are further stored in a list (with the same formats) and are accessible from R - see \code{intermediary} in the value section below.}
  \item{dist_weights}{
    a vector containing the weights to apply to the distance between the computed and the targeted statistics.
    These weights can be used to give more importance to a summary statistisc for example. The weights will
    be normalized before applying them. If not provided, no weights will be applied.
  }
  \item{\dots}{Additional arguments can be passed depending on the choosen method (see below)}
}

\section{Additional parameters}{
  Depending on the choosen method, you can specify the following arguments:
  \describe{
    \item{seed_count}{
	  a positive integer, the initial seed value provided to the function \code{model} (if \code{use_seed=TRUE}). This value is incremented by 1 at each call of the function \code{model}.}
    \item{inside_prior}{
	  logical used when \code{method} is \code{"Beaumont"}, \code{"Lenormand"} or \code{"Emulation"}. \code{TRUE} by default.
	  If \code{FALSE}, parameter sampling is not restricted to the initial ranges of the prior distribution during the subsequent algorithm steps.}
    \item{tolerance_tab}{
	  a vector containing the sequence of tolerance thresholds when \code{method} is \code{"Beaumont"}, or the targeted final tolerance threshold when \code{method} is \code{"Drovandi"}.}
    \item{alpha}{
	  a positive number between 0 and 1 (strictly) used when \code{method} is \code{"Drovandi"}, \code{"Delmoral"}, \code{"Lenormand"} or \code{"Emulation"}.
	  \code{alpha} is the proportion of particles rejected at each step in the algorithm \code{"Drovandi"}.
	  This is the proportion of particles kept at each step in the algorithms \code{"Delmoral"}, \code{"Lenormand"} and \code{"Emulation"}.
	  Default values are 0.5 when \code{method} is \code{"Drovandi"}, \code{"Lenormand"} or \code{"Emulation"} and 0.9 for \code{"Delmoral"}.
	  See the package's vignette for details.}
    \item{c}{
	  a positive number between 0 and 1 (strictly) used when \code{method} is \code{"Drovandi"}.
	  This is the expected proportion of particles which are going to be duplicated at each step.
	  Default value is 0.01.
	  See the package's vignette and Drovandi and Pettitt (2011) for details.}
    \item{first_tolerance_level_auto}{
	  logical used when \code{method} is \code{"Drovandi"}.
	  Default value is \code{TRUE}. In this case, the first tolerance threshold is determined by the algorithm, by taking the 1-\code{alpha} quantile of the distances between the simulated and targeted summary statistics.
	  If \code{FALSE}, the initial tolerance threshold for the first step has to be provided as the first element of the vector \code{tolerance_tab}. In this case, the targeted final tolerance threshold is the second element of \code{tolerance_tab}.}
    \item{M}{
	  a positive integer used when \code{method} is \code{"Delmoral"}.
	  This is the number of \code{model} simulations performed for each parameter set.
	  Default value is 1. See the package's vignette and Del Moral et al. (2012) for details.}
    \item{nb_threshold}{
	  a positive integer used when \code{method} is \code{"Delmoral"}. Default value is 0.5*\code{nb_simul}.
	  This is the minimal effective sample size below which a resampling step is launched. See the package's vignette and Del Moral et al. (2012) for details.}
    \item{tolerance_target}{
	  a positive number used when \code{method} is \code{"Delmoral"}.
	  This is the targeted final tolerance threshold.}
    \item{p_acc_min}{
	  a positive number between 0 and 1 (strictly) used when \code{method} is \code{"Lenormand"} or \code{"Emulation"}.
	  This is the stopping criterion of the algorithm: a small number ensures a better convergence of the algorithm, but at a cost in computing time.
	  Default value is 0.05. See the package's vignette and Lenormand et al. (2012) for details.}
    \item{n_step_emulation}{
          a positive integer, the number of times the emulation is repeated. \code{9} by default.
          }
    \item{emulator_span}{
          a positive number, the number of design points selected for the local regression.
          \code{50} by default.
          }
    \item{progress_bar}{
      logical, \code{FALSE} by default. If \code{TRUE}, \code{ABC_sequential} will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.
          }
    \item{max_pick}{
          a positive number, the max number of fails when moving particle inside the prior. Enabled only if inside_prior is to \code{TRUE}.
          \code{10000} by default.
          }
  }
}
\details{
  See the package's vignette for details on the four algorithms.
}
\value{
  The returned value is a list containing the following components: 
  \item{param}{
  	The model parameters used in the \code{model} simulations.}
  \item{stats}{
  	The summary statistics obtained at the end of the \code{model} simulations.}
  \item{weights}{
  	The weights of the different \code{model} simulations.}
  \item{stats_normalization}{
  	The standard deviation of the summary statistics across the \code{model} simulations of the initial step.
	These values are used to normalize the summary statistics before the computation of the Euclidean distance between simulations and data.}
  \item{epsilon}{
  	The final maximal distance between simulations and data in the retained sample of particles.}
  \item{nsim}{
  	The number of \code{model} simulations performed.}
  \item{computime}{
  	The computing time to perform the simulations.}
  \item{intermediary}{
  	Intermediary results stored when the option \code{"verbose=TRUE"} is chosen. Each element of this list corresponds to a different step. See the argument \code{verbose} above for more details on the information stored.}
}
\references{
	Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009) Adaptive approximate Bayesian computation. \emph{Biometrika},\bold{96}, 983--990.

	Del Moral, P., Doucet, A., and Jasra, A. (2012) An adaptive sequential Monte Carlo method for approximate Bayesian computation. \emph{Statistics and Computing}, \bold{22}, 1009--1020.

	Drovandi, C. C. and Pettitt, A. N. (2011) Estimation of parameters for macroparasite population evolution using approximate Bayesian computation. \emph{Biometrics}, \bold{67}, 225--233.

	Lenormand, M., Jabot, F., Deffuant G. (2012) Adaptive approximate Bayesian computation for complex models. http://arxiv.org/pdf/1111.1308.pdf
	
	Jabot, F., Lagarrigues G., Courbaud B., Dumoulin N. (2015). A comparison of emulation methods for Approximate Bayesian Computation. To be published.
}
\author{Franck Jabot, Thierry Faure and Nicolas Dumoulin}
\seealso{
  \code{\link{binary_model}}, \code{\link{binary_model_cluster}}, \code{\link{ABC_rejection}}, \code{\link{ABC_emulation}}, \code{\link{ABC_mcmc}}
}
\examples{
 \dontrun{ 
    ##### EXAMPLE 1 #####
    #####################
    set.seed(1)

    ## artificial example to show how to use the 'ABC_sequential' function.
    ## defining a simple toy model:
    toy_model<-function(x){ 2 * x + 5 + rnorm(1,0,0.1) }

    ## define prior information
    toy_prior=list(c("unif",0,1)) # a uniform prior distribution between 0 and 1

    ## define the targeted summary statistics
    sum_stat_obs=6.5

    ## to perform the Beaumont et al. (2009)'s method:
    ##
    tolerance=c(1.5,0.5)
    ABC_Beaumont<-ABC_sequential(method="Beaumont", model=toy_model, prior=toy_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, tolerance_tab=tolerance)
    ABC_Beaumont

    ## to perform the Drovandi and Pettitt (2011)'s method:
    ##
    tolerance=0.5
    c_drov=0.7
    ABC_Drovandi<-ABC_sequential(method="Drovandi", model=toy_model, prior=toy_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, tolerance_tab=tolerance, c=c_drov)
    ABC_Drovandi

    ## to perform the Del Moral et al. (2012)'s method:
    ##
    alpha_delmo=0.5
    tolerance=0.5
    ABC_Delmoral<-ABC_sequential(method="Delmoral", model=toy_model, prior=toy_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, alpha=alpha_delmo, tolerance_target=tolerance)
    ABC_Delmoral

    ## to perform the Lenormand et al. (2012)'s method:
    ##
    pacc=0.4
    ABC_Lenormand<-ABC_sequential(method="Lenormand", model=toy_model, prior=toy_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, p_acc_min=pacc)
    ABC_Lenormand


    ##### EXAMPLE 2 #####
    #####################

    ## this time, the model has two parameters and outputs two summary statistics.
    ## defining a simple toy model:
    toy_model2<-function(x){ c( x[1] + x[2] + rnorm(1,0,0.1) , x[1] * x[2] + rnorm(1,0,0.1) ) }

    ## define prior information
    toy_prior2=list(c("unif",0,1),c("normal",1,2))
    # a uniform prior distribution between 0 and 1 for parameter 1, and a normal distribution
    # of mean 1 and standard deviation of 2 for parameter 2.

    ## define the targeted summary statistics
    sum_stat_obs2=c(1.5,0.5)

    ## to perform the Beaumont et al. (2009)'s method:
    ##
    tolerance=c(1.5,0.5)
    ABC_Beaumont<-ABC_sequential(method="Beaumont", model=toy_model2, prior=toy_prior2,
    nb_simul=20, summary_stat_target=sum_stat_obs2, tolerance_tab=tolerance)
    ABC_Beaumont

    ## to perform the Drovandi and Pettitt (2011)'s method:
    ##
    tolerance=0.5
    c_drov=0.7
    ABC_Drovandi<-ABC_sequential(method="Drovandi", model=toy_model2, prior=toy_prior2,
    nb_simul=20, summary_stat_target=sum_stat_obs2, tolerance_tab=tolerance, c=c_drov)
    ABC_Drovandi

    ## to perform the Del Moral et al. (2012)'s method:
    ##
    alpha_delmo=0.5
    tolerance=0.5
    ABC_Delmoral<-ABC_sequential(method="Delmoral", model=toy_model2, prior=toy_prior2,
    nb_simul=20, summary_stat_target=sum_stat_obs2, alpha=alpha_delmo, tolerance_target=tolerance)
    ABC_Delmoral

    ## to perform the Lenormand et al. (2012)'s method:
    ##
    pacc=0.4
    # Only uniform priors are supported for the method "Lenormand" (since it performs a Latin
    # Hypercube sampling at the beginning):
    toy_prior2=list(c("unif",0,1),c("unif",0.5,1.5))
    # a uniform prior distribution between 0 and 1 for parameter 1, and a normal distribution of
    # mean 1 and standard deviation of 1 for parameter 2.
    ABC_Lenormand<-ABC_sequential(method="Lenormand", model=toy_model2, prior=toy_prior2,
    nb_simul=20, summary_stat_target=sum_stat_obs2, p_acc_min=pacc)
    ABC_Lenormand

    ##### EXAMPLE 3 #####
    #####################

    ## this time, the model is a C++ function packed into a R function -- this time, the option
    # 'use_seed' must be turned to TRUE.
    n=10
    ## define prior information
    trait_prior=list(c("unif",3,5),c("unif",-2.3,1.6),c("unif",-25,125),c("unif",-0.7,3.2))
    trait_prior

    ## define the targeted summary statistics
    sum_stat_obs=c(100,2.5,20,30000)

    ## to perform the Beaumont et al. (2009)'s method:
    ##
    tolerance=c(8,5)
    ABC_Beaumont<-ABC_sequential(method="Beaumont", model=trait_model, prior=trait_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, tolerance_tab=tolerance, use_seed=TRUE)
    ABC_Beaumont

    ## to perform the Drovandi and Pettitt (2011)'s method:
    ##
    tolerance=3
    c_drov=0.7
    ABC_Drovandi<-ABC_sequential(method="Drovandi", model=trait_model, prior=trait_prior,
      nb_simul=20, summary_stat_target=sum_stat_obs, tolerance_tab=tolerance, c=c_drov,
      use_seed=TRUE)
    ABC_Drovandi

    ## to perform the Del Moral et al. (2012)'s method:
    ##
    alpha_delmo=0.5
    tolerance=3
    ABC_Delmoral<-ABC_sequential(method="Delmoral", model=trait_model, prior=trait_prior,
      nb_simul=20, summary_stat_target=sum_stat_obs, alpha=alpha_delmo,
      tolerance_target=tolerance, use_seed=TRUE)
    ABC_Delmoral

    ## to perform the Lenormand et al. (2012)'s method:
    ##
    pacc=0.4
    ABC_Lenormand<-ABC_sequential(method="Lenormand", model=trait_model, prior=trait_prior,
      nb_simul=20, summary_stat_target=sum_stat_obs, p_acc_min=pacc, use_seed=TRUE)
    ABC_Lenormand


    ##### EXAMPLE 4 - Parallel implementations #####
    ################################################

    ## NB: the option use_seed must be turned to TRUE.

    ## For models already running with the option use_seed=TRUE, simply change
    # the value of n_cluster:
    sum_stat_obs=c(100,2.5,20,30000)
    ABC_Lenormand<-ABC_sequential(method="Lenormand", model=trait_model, prior=trait_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, p_acc_min=pacc, use_seed=TRUE, n_cluster=2)
    ABC_Lenormand

    ## For other models, change the value of n_cluster and modify the model so that the
    # first parameter becomes a seed information value:
    toy_model_parallel<-function(x){ 
	set.seed(x[1])
	2 * x[2] + 5 + rnorm(1,0,0.1) }
    sum_stat_obs=6.5

    ABC_Lenormand<-ABC_sequential(method="Lenormand", model=toy_model_parallel, prior=toy_prior,
    nb_simul=20, summary_stat_target=sum_stat_obs, p_acc_min=pacc, use_seed=TRUE, n_cluster=2)
    ABC_Lenormand
 }%dontrun
}
\keyword{abc}
\keyword{model}
\keyword{inference}
\keyword{sequential_monte_carlo}
\keyword{population_monte_carlo}
