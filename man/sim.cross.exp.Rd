\name{sim.cross.exp}
\alias{sim.cross.exp}
\title{Simulate a QTL experiment}
\description{
  Simulates data for a QTL backcross experiment using a model in which QTLs act additively.
}
\usage{
sim.cross.exp(map, model = NULL, n.ind = 100, type = c("bc"), error.prob = 0, 
              missing.prob = 0, keep.qtlgeno = TRUE, keep.errorind = TRUE, 
              m = 0, p = 0, map.function = c("haldane", "kosambi", "c-f", "morgan"),
              dist, params, logNorm = FALSE)
}
\arguments{
 \item{map}{A list whose components are vectors containing the marker
 locations on each of the chromosomes.}
 \item{model}{A matrix where each row corresponds to a
 different QTL, and gives the chromosome number, cM position and
 effects of the QTL.}
 \item{n.ind}{Number of individuals to simulate.}
 \item{type}{Indicates type of experimental cross to simulate. Only backcross supported.}
 \item{error.prob}{The genotyping error rate.}
 \item{missing.prob}{The rate of missing genotypes.}
 \item{keep.qtlgeno}{If TRUE, genotypes for the simulated QTLs will be
   included in the output.}
 \item{keep.errorind}{If TRUE, and if \code{error.prob > 0}, the
   identity of genotyping errors will be included in the output.}
 \item{m}{Interference parameter; a non-negative integer.  0 corresponds
   to no interference.}
 \item{p}{Probability that a chiasma comes from the no-interference
   mechanism} 
 \item{map.function}{Indicates whether to use the Haldane, Kosambi, 
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{dist}{Parametric distribution from which to generate phenotypes}
 \item{params}{Vector of parameters for parametric distributions}
 \item{logNorm}{If TRUE, continuous observations in point-mass mixture generated from
  log normal distribution. If FALSE, continuous observations are from truncated normal.}
}
\details{
  This function extends \code{\link[qtl]{sim.cross}} to simulate phenotypes from several
  additional parametric distributions including point-mass mixtures. The 
  available distributions and their specification in \code{dist} are
  normal (\code{"norm"}), log normal (\code{"lognorm"}), gamma (\code{"gamma"}), 
  truncated normal (\code{"truncNorm"}), t (\code{"t"}), Cauchy (\code{"cauchy"}),
  logistic (\code{"logistic"}), exponential (\code{"expon"}), Uniform (\code{"uniform"},
  Poisson (\code{"pois"}), and point-mass mixture (\code{"pointmass"}). 
  
  Parameters for the distribution are specified as a vector in the \code{param} 
  argument as follows: normal (\code{mean, sd}), log normal (\code{mean, sd}),
  gamma (\code{shape, scale}), truncated normal (\code{mean, sd, floor}) where
  \code{mean} and \code{sd} are the parameters of the normal distribution prior 
  to truncation and \code{floor} is the lower value at which the distribution
  is truncated, t (\code{df}), Cauchy (\code{location, scale}), logistic (\code{location, scale}),
  exponential (\code{rate}), uniform (\code{min, max}), Poisson (\code{lambda}), 
  point-mass mixture (\code{proportion, mean, sd}) where proportion is the proportion 
  of observations in the point-mass and mean and sd are parameters for generating 
  the continuous observations from a truncated normal or log normal distribution. 
  Point-mass values are always simulated at 0.
 
}
\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for
  details. 

  If \code{keep.qtlgeno} is TRUE, the cross object will contain a
  component \code{qtlgeno} which is a matrix containing the QTL
  genotypes (with complete data and no errors), coded as in the genotype
  data.

  If \code{keep.errorind} is TRUE and errors were simulated, each
  component of \code{geno} will each contain a matrix \code{errors},
  with 1's indicating simulated genotyping errors.
}
\author{Sandra L. Taylor, \email{sltaylor@ucdavis.edu} }
\seealso{ \code{\link[qtl]{sim.cross}}}
\examples{
\dontshow{library(qtl)}
# Simulate a genetic mapp
Map <- sim.map(c(100,100), n.mar=c(22), include.x=FALSE, eq.spacing=TRUE)
num.ind <- 200

# simulate 200 backcross individuals with one qtl affecting a point-mass mixture trait
# define a qtl at 30 cm on chromosome 1 with an additive effect of 0.6
# on the mean of the continuous component and -0.1 on the proportion of 
# observations in the point-mass at 0
QTL.info <- rbind(c(1,30,0.6,-0.01))
sim <- sim.cross.exp(map=Map, model=QTL.info, n.ind=num.ind, type="bc",
          keep.qtlgeno=TRUE, map.function="haldane", dist="pointmass", params=c(0.5,1,0.6))

}  
\keyword{datagen}