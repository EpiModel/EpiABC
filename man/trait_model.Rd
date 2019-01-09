\name{trait_model}
\alias{trait_model}
\title{
  Trait model for EasyABC tutorial
}
\description{
  A stochastic individual-based model to demonstrate how the EasyABC functions can be used
}
\details{
  This model is drawn from Jabot (2010) A stochastic dispersal-limited trait-based model of community dynamics. Journal of Theoretical Biology, 262, 650-661.
  It represents the stochastic dynamics of an ecological community where each species is represented by a set of traits (i.e. characteristics) which determine its competitive ability.
  The model requires four parameters: an immigration rate \code{I}, and three additional parameters (\code{h}, \code{A} and \eqn{\sigma}{sigma}) describing the way traits determine species competitive ability.
  The model additionnally requires two fixed variables: the total number of individuals in the local community \code{J} and the number of traits used \code{n_t}.
  The model outputs four summary statistics: the species richness of the community \code{S}, its Shannon's index \code{H}, the mean of the trait value among individuals \code{MTV} and the skewness of the trait value distribution \code{STV}.
  The simulation code \code{trait_model} reads sequentially \code{J}, \code{I}, \code{A}, \code{n_t}, \code{h} and \eqn{\sigma}{sigma}.
}
\keyword{sample}
\examples{
  trait_model(input=c(1,500,1,1,1,1))
}
