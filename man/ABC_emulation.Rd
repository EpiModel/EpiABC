\name{ABC_emulation}
\alias{ABC_emulation}
\title{
  Rejection sampling scheme for ABC using an emulator}
\description{
  This function launches a series of \code{nb_design_pts} model simulations with model parameters drawn in the prior distribution specified in \code{prior_matrix}, build an emulator with these computed design points and then launches a series of \code{nb_simul} emulator simulations.
}
\usage{
  ABC_emulation(model, prior,  nb_design_pts, nb_simul, prior_test=NULL,
  summary_stat_target=NULL, emulator_span = 50, tol=NULL, use_seed=FALSE,
  seed_count=0, n_cluster=1, verbose=FALSE, progress_bar=FALSE)
}
\arguments{
  \item{model}{
    a \code{R} function implementing the model to be simulated. It must take as arguments a vector of model parameter values and it must return a vector of summary statistics. When using the option \code{use_seed=TRUE}, \code{model} must take as arguments a vector containing a seed value and the model parameter values.
    A tutorial is provided in the package's vignette to dynamically link a binary code to a \code{R} function.  Users may alternatively wish to wrap their binary executables using the provided functions \code{\link{binary_model}} and \code{\link{binary_model_cluster}}. The use of these functions is associated with slightly different constraints on the design of the binary code (see \code{\link{binary_model}} and \code{\link{binary_model_cluster}}).}
  \item{prior}{
    a list of prior information. Each element of the list corresponds to a model parameter. The list element must be a vector whose first argument determines the type of prior distribution: possible values are \code{"unif"} for a uniform distribution on a segment, \code{"normal"} for a normal distribution, \code{"lognormal"} for a lognormal distribution or \code{"exponential"} for an exponential distribution.
    The following arguments of the list elements contain the characteritiscs of the prior distribution chosen: for \code{"unif"}, two numbers must be given: the minimum and maximum values of the uniform distribution; for \code{"normal"}, two numbers must be given: the mean and standard deviation of the normal distribution; for \code{"lognormal"}, two numbers must be given: the mean and standard deviation on the log scale of the lognormal distribution; for \code{"exponential"}, one number must be given: the rate of the exponential distribution. User-defined prior distributions can also be provided. See the vignette for additional information on this topic.}
  \item{nb_design_pts}{
    a positive integer equal to the desired number of simulations of the model used to build the emulator.}
  \item{nb_simul}{
    a positive integer equal to the desired number of simulations of the emulator.}
  \item{prior_test}{
        a string expressing the constraints between model parameters.
        This expression will be evaluated as a logical expression, you can use all the logical operators including \code{"<"}, \code{">"}, \ldots
        Each parameter should be designated with \code{"X1"}, \code{"X2"}, \ldots in the same order as in the prior definition.
        If not provided, no constraint will be applied.}
  \item{summary_stat_target}{
    a vector containing the targeted (observed) summary statistics.
    If not provided, \code{ABC_rejection} only launches the simulations and outputs the simulation results.}
  \item{emulator_span}{
    a positive number, the number of design points selected for the local regression.
    \code{50} by default.}
  \item{tol}{
    tolerance, a strictly positive number (between 0 and 1) indicating the proportion of simulations retained nearest the targeted summary statistics.}
  \item{use_seed}{
    logical. If \code{FALSE} (default), \code{ABC_rejection} provides as input to the function \code{model} a vector containing the model parameters used for the simulation.
    If \code{TRUE}, \code{ABC_rejection} provides as input to the function \code{model} a vector containing an integer seed value and the model parameters used for the simulation.
    In this last case, the seed value should be used by \code{model} to initialize its pseudo-random number generators (if \code{model} is stochastic).}
  \item{seed_count}{
    a positive integer, the initial seed value provided to the function \code{model} (if \code{use_seed=TRUE}). This value is incremented by 1 at each call of the function \code{model}.}
  \item{n_cluster}{
    a positive integer. If larger than 1 (the default value), \code{ABC_rejection} will launch \code{model} simulations in parallel on \code{n_cluster} cores of the computer.}
  \item{verbose}{
    logical. \code{FALSE} by default. If \code{TRUE}, \code{ABC_rejection} writes in the current directory intermediary results at the end of each step of the algorithm in the file "output".
    These outputs have a matrix format, in wich each raw is a different simulation, the first columns are the parameters used for this simulation, and the last columns are the summary statistics of this simulation.}
  \item{progress_bar}{
    logical, \code{FALSE} by default. If \code{TRUE}, \code{ABC_rejection} will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}
}
%\details{}
\value{
  The returned value is a list containing the following components: 
  \item{param}{
  	The model parameters used in the \code{model} simulations.}
  \item{stats}{
  	The summary statistics obtained at the end of the \code{model} simulations.}
  \item{weights}{
  	The weights of the different \code{model} simulations. In the standard rejection scheme, all \code{model} simulations have the same weights.}
  \item{stats_normalization}{
  	The standard deviation of the summary statistics across the \code{model} simulations.}
  \item{nsim}{
  	The number of \code{model} simulations performed.}
  \item{nrec}{
  	The number of retained simulations (if targeted summary statistics are provided).}
  \item{computime}{
  	The computing time to perform the simulations.}
}
\references{
  Jabot, F., Lagarrigues G., Courbaud B., Dumoulin N. (2015). A comparison of emulation methods for Approximate Bayesian Computation. To be published.
}
\author{Franck Jabot, Thierry Faure and Nicolas Dumoulin}
\seealso{
  \code{\link{binary_model}}, \code{\link{binary_model_cluster}}, \code{\link{ABC_sequential}}, \code{\link{ABC_mcmc}}
}
\examples{
 \dontrun{ 

    ##### EXAMPLE 1 #####
    #####################

    ## the model is a C++ function packed into a R function -- the option 'use_seed'
    ##  must be turned to TRUE.
    trait_prior=list(c("unif",3,5),c("unif",-2.3,1.6),c("unif",-25,125),c("unif",-0.7,3.2))
    trait_prior

    ## only launching simulations with parameters drawn in the prior distributions
    ABC_emul = ABC_emulation(model=trait_model, prior=trait_prior,
      nb_design_pts=10, nb_simul=300, use_seed=TRUE, progress=TRUE)
    ABC_emul

    ## launching simulations with parameters drawn in the prior distributions and performing
    # the rejection step
    sum_stat_obs=c(100,2.5,20,30000)
    ABC_emul = ABC_emulation(model=trait_model, prior=trait_prior, tol=0.2, nb_design_pts=10,
      nb_simul=100, summary_stat_target=sum_stat_obs, use_seed=TRUE, progress=TRUE)
    ABC_emul

 }%dontrun
}
\keyword{abc}
\keyword{model}
\keyword{inference}
\keyword{emulation}
