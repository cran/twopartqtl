\name{pmscan}
\alias{pmscan}
\title{Interval mapping with covariates for point-mass mixtures}
\description{
  Conducts interval mapping with covariates for experimental crosses for which the phenotype has a point-mass mixture distribution. Identifies QTLs by testing the joint null hypothesis of no differences in the proportion of observations in the point-mass and no difference the means of the non-point mass observations.
}
\usage{
pmscan(cross, chr, pheno.col = 1, pm.value = 0, maxit = 4000, tol = 1e-04, 
       addcovar, verbose = FALSE, imp.method = c("imp", "argmax"), 
       error.prob = 1e-04, map.function = c("haldane", "kosambi", "c-v", 
       "morgan"), n.perm, use.log = FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
  \item{chr}{Optional vector indicating the chromosomes for which LOD
   scores should be calculated.}
  \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype} 
  \item{pm.value}{Value of the point-mass observations}
  \item{maxit}{Maximum number of iterations}
  \item{tol}{Tolerance value for determining convergence }
  \item{addcovar}{A matrix of additive covariates with dimensions \emph{n} X \emph{number of covariates}}
  \item{verbose}{In the case \code{n.perm} is specified, displays
   information about the progress of the permutation tests.}
  \item{imp.method}{Method used to impute any missing marker genotype
   data. See \code{\link[qtl]{fill.geno}} for details.}
  \item{error.prob}{Genotyping error probability assumed when imputing
   the missing marker genotype data.}
  \item{map.function}{Map function used when imputing the missing marker
   genotype data.}
  \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
  \item{use.log}{If TRUE phenotype values not in the point-mass are log transformed. }
}
\details{
  The multipoint genotype probabilities are first calculated using 
  \code{\link[qtl]{calc.genoprob}}.

  The method is currently implemented only for experimental crosses consisting of 
  two genotypes (e.g., backcross, recombinant inbred lines).
 
  Only autosomal chromosomes are evaluated.
  
  Individuals with any missing phenotypes or covariates are dropped.
  
  Only additive covariates are modeled. Covariates must be numeric matrices.   
}
\value{
   The function returns an object of the same form as the function
  \code{\link[qtl]{scanone}}.

  If \code{n.perm} is missing, the function returns the scan results as
  a data.frame with three columns: chromosome, position, LOD score.
  Attributes indicate the names and positions of the chosen marker
  covariates.

  If \code{n.perm} > 0, the function results the results of a
  permutation test: a vector giving the genome-wide maximum LOD score in
  each of the permutations. 

}
\references{ 
  Taylor, S.L. and K.S. Pollard 20XX. Composite interval mapping to identify 
  quantitative trait loci for point-mass mixture phenotypes.  \emph{Genetics Research}, \bold{XX}, xxx--xxx
  }
\author{Sandra L. Taylor, \email{sltaylor@ucdavis.edu}}
\seealso{ \code{\link[qtl]{scanone}}}
\examples{
\dontshow{library(qtl)}
# Simulate backcross experiment
Map <- sim.map(c(100,100), n.mar=c(22), include.x=FALSE, eq.spacing=TRUE)
num.ind <- 200
QTL.info <- rbind(c(1,30,0.6,-0.02))
sim <- sim.cross.exp(map=Map, model=QTL.info, n.ind=num.ind, type="bc",
    keep.qtlgeno=TRUE, map.function="haldane", dist="pointmass", params=c(0.5,6,0.6))
sim <- calc.genoprob(sim, step=1)
Covs <- matrix(rnorm(400), nrow=200)

out <- pmscan(sim, chr=c(1:2), pheno.col=1, pm.value=0, addcovar=Covs)

# Permutation tests
\dontrun{out.perm <- pmscan(sim, chr=c(1:2), pheno.col=1, pm.value=0, addcovar=Covs, n.perm=1000)}
}
\keyword{models}
