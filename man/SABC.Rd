\name{SABC}

\alias{SABC}
\alias{SABC.noninf}
\alias{SABC.inf}

\title{Simulated Annealing approach to Approximate Bayesian Computation (SABC)}

\description{
Algorithms for the Simulated Annealing approach to Approximate Bayesian Computation (SABC).
}

\usage{
SABC(r.model, r.prior, d.prior, n.sample, eps.init, iter.max,
     v=ifelse(method=="informative",0.4,1.2), beta=0.8,
     delta=0.1, resample=5*n.sample, verbose=n.sample,
     method="noninformative", adaptjump=TRUE,
     summarystats=FALSE, y=NULL, f.summarystats=NULL, ...)
}

\arguments{
  \item{r.model}{Function that returns either a random sample from the likelihood or a scalar distance between such a sample and the data. The first argument must be the parameter vector.}
  \item{r.prior}{Function that returns a random sample from the prior.}
  \item{d.prior}{Function that returns the density of the prior distribution.}
  \item{n.sample}{Size of the ensemble.}
  \item{eps.init}{Initial tolerance or temperature.}
  \item{iter.max}{Total number of simulations from the likelihood.}
  \item{v}{Tuning parameter that governs the annealing speed. Defaults to 1.2, for the \code{noninformative} algorithm and 0.4, for the \code{informative} algorithm.}
  \item{beta}{Tuning parameter that governs the mixing in parameter space. Defaults to 0.8.}
  \item{delta}{Tuning parameter for the resampling steps. Defaults to 0.1.}
  \item{resample}{Number of accepted particle updates after which a resampling step is performed. Defaults to 5*\code{n.sample}.}
  \item{verbose}{Shows the iteration progress each \code{verbose} simulations from the likelihood. NULL for no output. Defaults to \code{verbose = n.sample}.}
  \item{adaptjump}{Whether to adapt covariance of jump distribution. Default is TRUE.}
  \item{method}{Argument to select algorithm. Accepts \code{noninformative} or \code{informative}.}
  \item{summarystats}{Whether summary statistics shall be calculated (semi-) automatically. Defaults to FALSE.}
  \item{y}{Data vector. Needs to be provided if either \code{summarystats = TRUE} or if \code{r.model} returns a sample from the likelihood.}
  \item{f.summarystats}{If \code{summarystats = TRUE} this function is needed for the calculation of the summary statistics. Defaults to \code{f.summarystats(x)=(x,x^2,x^3)}, where the powers are to be understood element-wise.}
  \item{...}{further arguments passed to \code{r.model}}
}

\details{
SABC defines a class of algorithms for particle ABC that are inspired by Simulated Annealing. Unlike other algorithms, this class is not based on importance sampling, and hence does not suffer from a loss of effective sample size due to re-sampling. The approach is presented in detail in Albert, Kuensch, and Scheidegger (2014; see references).

This package implements two versions of SABC algorithms, for the cases of a non-informative or an informative prior. These are described in detail in the paper. The algorithms can be selected using the \code{method} argument: \code{method=noninformative} or \code{method=informative}.
In the informative case, the algorithm corrects for the bias caused by an over- or under-representation of the prior.

The argument \code{adaptjump} allows a choice of whether to adapt the covariance of the jump distribution. Default is TRUE.

Furthermore, the package allows for three different ways of using the data.
If \code{y} is not provided, the algorithm expects \code{r.model} to return a scalar measuring the distance between a random sample from the likelihood and the data.
If \code{y} is provided and \code{summarystats = FALSE}, the algorithm expects \code{r.model} to return a random sample from the likelihood and uses the relative sum of squares to measure the distances between \code{y} and random likelihood samples.
If \code{summarystats = TRUE} the algorithm calculates summary statistics semi-automatically, as described in detail in the paper by Fearnhead et al. (2012; see references).
The summary statistics are calculated by means of a linear regression applied to a sample from the prior and the image of \code{f.summarystats} of an associated sample from the likelihood.
}

\value{
Returns a list with the following components:
  \item{E}{Matrix with ensemble of samples.}
  \item{P}{Matrix with prior ensemble of samples.}
  \item{eps}{Value of tolerance (temperature) at final iteration.}
  \item{ESS}{Effective sample size, due to final bias correction (\code{informative} algorithm only).}
}

\references{
C. Albert, H. R. Kuensch and A. Scheidegger, Statistics and Computing 0960-3174 (2014),  arXiv:1208.2157, \emph{A Simulated Annealing Approach to Approximate Bayes Computations}.


P. Fearnhead and D. Prangle, J. R. Statist. Soc. B  74 (2012), \emph{Constructing summary statistics for approximate Bayesian computation: semi-automatic approximate Bayesian computation}.}

\author{
Carlo Albert <carlo.albert@eawag.ch>, Andreas Scheidegger, Tobia Fasciati. Package initially compiled by Lukas M. Weber.
}

\examples{
 \dontrun{ 
## Example for "noninformative" case
# Prior is uniform on [-10,10]
d.prior <- function(par)
    dunif(par,-10,10)
r.prior <- function()
    runif(1,-10,10)

# Model is the sum of two normal distributions. Return distance to observation 0:
f.dist <- function(par)
    return( abs(rnorm( 1 , par , ifelse(runif(1)<0.5,1,0.1 ) )))

# Run algorithm ("noninformative" case)
res <- SABC(f.dist,r.prior,d.prior,n.sample=500,eps.init=2,iter.max=50000)
 }%dontrun

 \dontrun{%
# Histogram of results
hist(res$E[,1],breaks=200)
 }
}
