\name{twopartCIM}
\alias{twopartCIM}
\title{Composite interval mapping for point-mass mixtures}
\description{
  Conducts composite interval mapping for phenotypes distributed as a point-mass
  mixture following the same procedure as QTL Cartographer. Significant covariate
  markers are chosen through a forward selection procedure up to a user-specified number.
  Selected markers are used as covariates in interval mapping procedure. Covariate
  marker covariates within a fixed, user-specified window of the location being
  tested are dropped. 
}
\usage{
twopartCIM(cross, pheno.col = 1, n.marcovar = 3, window = 10, pm.value = 0, 
          threshold = 1, maxit = 4000, tol = 1e-04, verbose = FALSE, 
          imp.method = c("imp", "argmax"), error.prob = 1e-04, 
          map.function = c("haldane", "kosambi", "c-v", "morgan"), use.log = FALSE, 
          n.perm)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
  \item{pheno.col}{Column number in the phenotype matrix to be used as the phenotype. Only one column can be analyzed at a time.} 
  \item{n.marcovar}{Number of marker covariates to use}
  \item{window}{Window size, in cM}
  \item{pm.value}{Value of the point-mass observations}
  \item{threshold}{Significance threshold (p-value) for retaining covariates}
  \item{maxit}{Maximum number of iterations}
  \item{tol}{Tolerance value for determining convergence }
  \item{verbose}{In the case \code{n.perm} is specified, displays
   information about the progress of the permutation tests.}
  \item{imp.method}{Method used to impute any missing marker genotype
   data. See \code{\link[qtl]{fill.geno}} for details.}
  \item{error.prob}{Genotyping error probability assumed when imputing
   the missing marker genotype data.}
  \item{map.function}{Map function used when imputing the missing marker
   genotype data.}
 \item{use.log}{If TRUE phenotype values not in the point-mass are log transformed. }
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
}
\details{
  Missing marker genotype data are first imputed via \code{\link[qtl]{fill.geno}}
  according to the specified \code{imp.method}.

  Covariate markers are then identified through a forward selection process as described
  in Taylor and Pollard (2009). The most significant marker is retained at each step
  until \code{n.marcovar} markers are identified. If \code{threshold} is less
  than 1, a marker will be used as a covariate only if it is significant (i.e.,
  p-value < threshold). Significant markers up to \code{n.marcovar} markers will
  be retained.
  
  Selected markers are then used as covariates in interval mapping procedure. Covariate
  marker covariates within a fixed, user-specified window of the location being
  tested are dropped. The number of covariate markers is adjusted if necessary
  to be at least twice the number of continuous observations.     
}
\value{
    The function returns an object of the same form as the function
  \code{\link[qtl]{scanone}}:

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
\seealso{ \code{\link[qtl]{cim}}}
\examples{
\dontshow{library(qtl)}
# Simulate backcross experiment
Map <- sim.map(c(100), n.mar=c(22), include.x=FALSE, eq.spacing=TRUE)
num.ind <- 200
QTL.info <- rbind(c(1,30,0.6,-0.02))
sim <- sim.cross.exp(map=Map, model=QTL.info, n.ind=num.ind, type="bc",
    keep.qtlgeno=TRUE, map.function="haldane", dist="pointmass", params=c(0.5,6,0.6))
sim <- calc.genoprob(sim, step=1)

out <- twopartCIM(sim, pheno.col=1, n.marcovar=3, pm.value=0)

# Permutation tests
\dontrun{out.perm <- twopartCIM(sim, pheno.col=1, n.marcovar=3, pm.value=0, n.perm=1000)}
}
\keyword{models}
