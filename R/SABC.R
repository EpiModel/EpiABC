#========================================================
#
#  Wrapper function for SABC algorithms
#
#========================================================

## Wrapper function for SABC algorithms. Runs function SABC.noninf() or 
## SABC.inf() depending on value of "method" argument. See documentation for 
## details of those functions.
## Furthermore, the distance functino f.dist is determined.
##
## ARGUMENTS:
## method:  Which algorithm to use. Available algorithms are "noninformative" 
##          and "informative", for the case of a non-informative or an 
##          informative prior.
## other arguments are passed to SABC.noninf() or SABC.inf()
##
## VALUES:
## Output of function SABC.noninf() or SABC.inf()



#' Simulated Annealing approach to Approximate Bayesian Computation (SABC)
#' 
#' Algorithms for the Simulated Annealing approach to Approximate Bayesian
#' Computation (SABC).
#' 
#' SABC defines a class of algorithms for particle ABC that are inspired by
#' Simulated Annealing. Unlike other algorithms, this class is not based on
#' importance sampling, and hence does not suffer from a loss of effective
#' sample size due to re-sampling. The approach is presented in detail in
#' Albert, Kuensch, and Scheidegger (2014; see references).
#' 
#' This package implements two versions of SABC algorithms, for the cases of a
#' non-informative or an informative prior. These are described in detail in
#' the paper. The algorithms can be selected using the \code{method} argument:
#' \code{method=noninformative} or \code{method=informative}. In the
#' informative case, the algorithm corrects for the bias caused by an over- or
#' under-representation of the prior.
#' 
#' The argument \code{adaptjump} allows a choice of whether to adapt the
#' covariance of the jump distribution. Default is TRUE.
#' 
#' Furthermore, the package allows for three different ways of using the data.
#' If \code{y} is not provided, the algorithm expects \code{r.model} to return
#' a scalar measuring the distance between a random sample from the likelihood
#' and the data. If \code{y} is provided and \code{summarystats = FALSE}, the
#' algorithm expects \code{r.model} to return a random sample from the
#' likelihood and uses the relative sum of squares to measure the distances
#' between \code{y} and random likelihood samples. If \code{summarystats =
#' TRUE} the algorithm calculates summary statistics semi-automatically, as
#' described in detail in the paper by Fearnhead et al. (2012; see references).
#' The summary statistics are calculated by means of a linear regression
#' applied to a sample from the prior and the image of \code{f.summarystats} of
#' an associated sample from the likelihood.
#' 
#' @aliases SABC SABC.noninf SABC.inf
#' @param r.model Function that returns either a random sample from the
#' likelihood or a scalar distance between such a sample and the data. The
#' first argument must be the parameter vector.
#' @param r.prior Function that returns a random sample from the prior.
#' @param d.prior Function that returns the density of the prior distribution.
#' @param n.sample Size of the ensemble.
#' @param eps.init Initial tolerance or temperature.
#' @param iter.max Total number of simulations from the likelihood.
#' @param v Tuning parameter that governs the annealing speed. Defaults to 1.2,
#' for the \code{noninformative} algorithm and 0.4, for the \code{informative}
#' algorithm.
#' @param beta Tuning parameter that governs the mixing in parameter space.
#' Defaults to 0.8.
#' @param delta Tuning parameter for the resampling steps. Defaults to 0.1.
#' @param resample Number of accepted particle updates after which a resampling
#' step is performed. Defaults to 5*\code{n.sample}.
#' @param verbose Shows the iteration progress each \code{verbose} simulations
#' from the likelihood. NULL for no output. Defaults to \code{verbose =
#' n.sample}.
#' @param adaptjump Whether to adapt covariance of jump distribution. Default
#' is TRUE.
#' @param method Argument to select algorithm. Accepts \code{noninformative} or
#' \code{informative}.
#' @param summarystats Whether summary statistics shall be calculated (semi-)
#' automatically. Defaults to FALSE.
#' @param y Data vector. Needs to be provided if either \code{summarystats =
#' TRUE} or if \code{r.model} returns a sample from the likelihood.
#' @param f.summarystats If \code{summarystats = TRUE} this function is needed
#' for the calculation of the summary statistics. Defaults to
#' \code{f.summarystats(x)=(x,x^2,x^3)}, where the powers are to be understood
#' element-wise.
#' @param ... further arguments passed to \code{r.model}
#' @return Returns a list with the following components: \item{E}{Matrix with
#' ensemble of samples.} \item{P}{Matrix with prior ensemble of samples.}
#' \item{eps}{Value of tolerance (temperature) at final iteration.}
#' \item{ESS}{Effective sample size, due to final bias correction
#' (\code{informative} algorithm only).}
#' @author Carlo Albert <carlo.albert@@eawag.ch>, Andreas Scheidegger, Tobia
#' Fasciati. Package initially compiled by Lukas M. Weber.
#' @references C. Albert, H. R. Kuensch and A. Scheidegger, Statistics and
#' Computing 0960-3174 (2014), arXiv:1208.2157, \emph{A Simulated Annealing
#' Approach to Approximate Bayes Computations}.
#' 
#' P. Fearnhead and D. Prangle, J. R. Statist. Soc. B 74 (2012),
#' \emph{Constructing summary statistics for approximate Bayesian computation:
#' semi-automatic approximate Bayesian computation}.
#' @examples
#' 
#'  \dontrun{ 
#' ## Example for "noninformative" case
#' # Prior is uniform on [-10,10]
#' d.prior <- function(par)
#'     dunif(par,-10,10)
#' r.prior <- function()
#'     runif(1,-10,10)
#' 
#' # Model is the sum of two normal distributions. Return distance to observation 0:
#' f.dist <- function(par)
#'     return( abs(rnorm( 1 , par , ifelse(runif(1)<0.5,1,0.1 ) )))
#' 
#' # Run algorithm ("noninformative" case)
#' res <- SABC(f.dist,r.prior,d.prior,n.sample=500,eps.init=2,iter.max=50000)
#'  }%dontrun
#' 
#'  \dontrun{%
#' # Histogram of results
#' hist(res$E[,1],breaks=200)
#'  }
#' 
#' @export SABC
SABC <- function(r.model, r.prior, d.prior, 
                 n.sample, eps.init, iter.max, 
                 v=ifelse(method=="informative",0.4,1.2), 
                 beta=0.8,  
                 delta=0.1, resample=5*n.sample, 
                 verbose=n.sample, 
                 method="noninformative", adaptjump=TRUE, 
                 summarystats=FALSE, y=NULL, f.summarystats=NULL, ...)
{
  f.dist <- r.model
  
  if( !summarystats & !is.null(y) )
  {
    if(!all(y !=0))
      f.dist <- function(par,...)
      {
        x <- r.model(par,...)
        return(sum((x-y)^2))
      }
    else
      f.dist <- function(par,...)
      {
        x <- r.model(par,...)
        return(sum((x-y)/y)^2)
      }
  }
  
  if( !summarystats & is.null(y) & length( r.model( r.prior(),... ) )>1 )
  {
      stop("r.model needs to return distance to data or data vector y needs to be provided!")
  }
  
  if(summarystats & is.null(y)) stop("data vector y needs to be provided!")
  if(summarystats & is.null(f.summarystats) )
    f.summarystats <- function(x) return(c(x,x^2,x^3))
  
  if (method == "noninformative" | method == "noninf" | method == "non") {
    SABC.noninf(f.dist=f.dist, d.prior=d.prior, r.prior=r.prior, 
                n.sample=n.sample, eps.init=eps.init, iter.max=iter.max, 
                v=v, beta=beta, delta=delta, resample=resample, 
                verbose=verbose, adaptjump=adaptjump, 
                summarystats=summarystats, y=y, f.summarystats=f.summarystats, ...)
  }
  else if (method == "informative" | method == "inf") {
    SABC.inf(f.dist=f.dist, d.prior=d.prior, r.prior=r.prior, 
             n.sample=n.sample, eps.init=eps.init, iter.max=iter.max, 
             v=v, beta=beta, delta=delta, resample=resample, 
             verbose=verbose, adaptjump=adaptjump, 
             summarystats=summarystats, y=y, f.summarystats=f.summarystats, ...)
  }
  else {
    print("error in SABC(): method name not recognized")
  }
}
