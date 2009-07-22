\name{twopartqtl-internal}
\alias{logitreg}
\alias{PMforwardSelect}
\alias{rpointmass}
\alias{rtruncnorm}
\alias{sim.cross.bc.exp}
\title{Internal twopartqtl functions}

\description{
  Internal qtl functions. These are generally not to be called by the
  user.
}

\usage{
logitreg <- function(x, y, wt = rep(1, length(y)), intercept = T, start = rep(0, p))
PMforwardSelect(geno, pheno, threshold, n.marcovar, pm.value=0)
rpointmass(num.ind, mix, mean, sd,logNorm)
rtruncnorm(num.ind, mean, sd, floor,logNorm)
sim.cross.bc.exp(map, model, n.ind, error.prob, missing.prob,
                    keep.errorind, m, p, map.function, dist, params, logNorm)
}


\author{Sandra L. Taylor, \email{sltaylor@ucdavis.edu} }

\keyword{internal}
