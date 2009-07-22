################################################################################
#
# twopartCIM.R
#
# copyright (c) 2009, Sandra L. Taylor
# 
# first written Aug, 2009
# Modified from cim.r in R/qtl written by Karl W Broman
#
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the twopartqtl package
# Contains: twopartCIM, PMforwardSelect, logitreg
#
################################################################################

################################################################################
# Composite interval mapping for point mass mixture data for experimental
# crosses with 2 genotypes
################################################################################

twopartCIM <- function(cross, pheno.col = 1, n.marcovar =3, window = 10, pm.value=0, 
              threshold=1, maxit=4000, tol=1e-4, verbose=FALSE, imp.method = c("imp",
               "argmax"), error.prob = 1e-4, map.function = c("haldane", "kosambi", "c-v", "morgan"), 
               use.log=FALSE, n.perm){

  imp.method <- match.arg(imp.method)
  map.function <- match.arg(map.function)
  type <- class(cross)[1]

  if (any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")
  
   if (use.log){
        for (i in 1:num.ind){
          if(!is.na(cross$pheno[i,pheno.col])){
            if (cross$pheno[i,pheno.col]!=pm.value) cross$pheno[i,pheno.col] <- log(cross$pheno[i,pheno.col])
          }
        }
   }

   pheno <-cross$pheno[,pheno.col]

  if (min(pheno, na.rm=TRUE)==pm.value) { Upper <- FALSE
  } else if (max(pheno, na.rm=TRUE)==pm.value) { Upper <- TRUE
  } else {
     pheno[pheno==pm.value] <- floor(min(pheno, na.rm=TRUE))
     pm.value <- floor(min(pheno, na.rm=TRUE))
     cross$pheno[,pheno.col] <- pheno 
     Upper <- FALSE
     pheno <-cross$pheno[,pheno.col]
  }
  
  if (any(is.na(pheno))) {
        cross <- subset(cross, ind = (!is.na(pheno)))
        pheno <- pheno[!is.na(pheno)]
    }
                                                 
   num.ind <- length(pheno)
 
   if (!missing(n.perm) && n.perm > 0) {
        results <- matrix(ncol = 1, nrow = n.perm)
        for (i in 1:n.perm) {
            o <- sample(length(pheno))
            pheno <- pheno[o]
            cross$pheno[, pheno.col] <- pheno
            temp <- twopartCIM(cross, pheno.col = pheno.col, n.marcovar = n.marcovar,
                window = window, pm.value=pm.value, threshold=threshold, imp.method = imp.method,
                error.prob = error.prob, map.function = map.function, use.log=FALSE)
            results[i, 1] <- max(temp[, 3], na.rm = TRUE)
        }
        class(results) <- c("scanoneperm", "matrix")
        return(results)
    }
  
  window <- window/2
  geno <- pull.geno(cross)
  if (any(is.na(geno))){
        geno <- pull.geno(fill.geno(cross, method = imp.method,
            error.prob = error.prob, map.function = map.function))
  }
  
  n.geno <- 2
  num.cont <- sum(pheno!=pm.value)
  max.cov  <- floor(num.cont/2)
  if (max.cov < (n.marcovar+n.geno)){
      warning("Number of covariates reduced due to too few continuous observations")
      n.marcovar <- max(max.cov-n.geno, 0)
  }
  
  if (num.cont==1){
      warning("Only one continuous observation. Cannot do two-part test.")
      firstlod <- 0
  } else {
     if (n.marcovar==0){
       firstlod <- scanone(cross, pheno.col=pheno.col, model="2part", method="em",
                    maxit=maxit, tol=tol, upper=Upper)[,1:3] 
       colnames(firstlod)[-(1:2)] <- c("lod")
    } else {  
  
  sig.marks <- PMforwardSelect(geno, pheno, threshold=threshold, n.marcovar, pm.value)
  MarkerNames <- colnames(geno)

  if (sig.marks[1]==0){
    n.marcovar <- 0
  } else n.marcovar <- length(sig.marks)
    if (sig.marks[1]==0){
       firstlod <- scanone(cross, pheno.col=pheno.col, model="2part", method="em",
                    maxit=maxit, tol=tol, upper=Upper)[,1:3]
       colnames(firstlod)[-(1:2)] <- c("lod")
   } else {         
      sig.names <- MarkerNames[sig.marks]
      chrpos <- find.markerpos(cross, sig.names)    
      ac <- as.matrix(geno[,sig.marks])

      firstlod<- pmscan(cross, pheno.col=pheno.col, pm.value=pm.value, maxit=maxit, 
                   tol=tol, addcovar=ac, verbose=verbose, imp.method = imp.method,
                   error.prob = error.prob, map.function = map.function)
 
  for (i in seq(along = sig.names)) {

          if (length(sig.names)==1){
               templod <- scanone(cross, chr = chrpos[i, 1],pheno.col=pheno.col, model="2part", method="em",
                        maxit=maxit, tol=tol, upper=Upper)[,1:3]
               colnames(templod)[-(1:2)] <- c("lod")
          } else {
          useac <- as.matrix(ac[, -i])
          templod <- pmscan(cross, pheno.col=pheno.col, pm.value=pm.value, addcovar = useac, 
                  chr = chrpos[i, 1],  maxit=maxit, tol=tol, verbose=verbose, imp.method = imp.method,
                  error.prob = error.prob, map.function = map.function)
          }
          wh1 <- (firstlod[, 1] == chrpos[i, 1] & firstlod[,
              2] >= chrpos[i, 2] - window & firstlod[, 2] <= chrpos[i,
              2] + window)
          wh2 <- (templod[, 2] >= chrpos[i, 2] - window & templod[, 2] <=
              chrpos[i, 2] + window)
          firstlod[wh1, 3] <- templod[wh2, 3]
      }

      attr(firstlod, "marker.covar") <- sig.names
      attr(firstlod, "marker.covar.pos") <- chrpos
      u <- table(chrpos[, 1])

        if (any(u > 1)) {
          u <- names(u)[u > 1]
          for (j in u) {
              wh <- which(chrpos[, 1] == j)
              pos <- chrpos[wh,2]
              scanpos <- firstlod[firstlod[, 1] == j, 2]
              need2drop <- t(sapply(scanpos, function(a, b, d) as.numeric(abs(a -
                  b) <= d), pos, window))

              n2drop <- apply(need2drop, 1, sum)

              if (any(n2drop > 1)) {
                  pat2drop <- apply(need2drop, 1, paste, collapse = "")
                  thepat <- unique(pat2drop[n2drop > 1])
                  for (k in thepat) {
                    whpos <- which(pat2drop == k)
                    whpos2 <- (firstlod[, 1] == j & !is.na(match(firstlod[,
                      2], scanpos[whpos])))
                    windowmarks <-  wh[as.logical(need2drop[whpos[1],])]
                    if (length(windowmarks)==n.marcovar){
                       templod <- scanone(cross, chr = chrpos[i, 1],pheno.col=pheno.col, 
                                  model="2part", method="em", upper=Upper)[,1:3]
                       colnames(templod)[-(1:2)] <- c("lod")
                    } else {
                      tempac <- ac[, -windowmarks]                    
                      useac <- matrix(tempac, nrow=num.ind)
                      templod <- pmscan(cross, pheno.col = pheno.col,pm.value=pm.value,
                               addcovar = useac, chr = j, maxit=maxit, tol=tol, verbose=verbose,
                               imp.method = imp.method,error.prob = error.prob, map.function = map.function)
                    }
                    firstlod[whpos2, 3] <- templod[whpos, 3]
                  }
              }
          }
       }
      }
     }  
    }  
    firstscan <- firstlod
    return(firstscan)
}

################################################################################
# Forward selection procedure up to the user specified number of markers
################################################################################

PMforwardSelect <- function(geno, pheno, threshold, n.marcovar, pm.value=0){
    num.ind <- length(pheno)
    num.marks <- dim(geno)[2]
    contIndex <- which(pheno!=pm.value)
    contPheno <- pheno[contIndex]
    binary <- rep(0,num.ind)
    binary[-contIndex] <- 1
    contGeno <- geno[contIndex,]

    binNull <- logitreg(rep(1,num.ind), binary, intercept=FALSE)$value
    linNull <- deviance(lm(contPheno~1))
    NullDev <- binNull+linNull

    binDev.1 <- sapply(1:num.marks, function(x) logitreg(matrix(geno[,x]), binary)$value)
    linDev.1 <- sapply(1:num.marks, function(x) deviance(lm(contPheno~contGeno[,x])))
    TotalDev.1 <- binDev.1+linDev.1
    pvalues <- 1-pchisq((NullDev-TotalDev.1), df=2)
    ms <- which.min(TotalDev.1)
    
    if (pvalues[ms] > threshold) {
        ms <- 0
        return(ms)
    } else if (n.marcovar==1) {
      return(ms) 
    } else { 
 
    newDev <- TotalDev.1[ms]
    marks <- c(1:num.marks)

      for (j in 1:c(n.marcovar-1)){
        binDev <- sapply(marks[-ms], function(x) logitreg(cbind(geno[,ms],geno[,x]), binary)$value)
        linDev <- sapply(marks[-ms], function(x) deviance(lm(contPheno~contGeno[,ms]+contGeno[,x])))
        TotalDev <- binDev+linDev
        pvalues <- 1-pchisq((newDev-TotalDev), df=2)
        if (all(pvalues > threshold)) break
        newMS <- which.min(TotalDev)
        newDev <- TotalDev[newMS]
        ms <- c(ms, marks[-ms][newMS])
      }
    return(ms)
  } 
}


################################################################################
# Logit regression function from in Venables and Ripley 2002, see pages 197-199
################################################################################

logitreg <- function(x, y, wt = rep(1, length(y)),
               intercept = TRUE, start = rep(0, p), ...)
{
  if (!is.matrix(x)) x <- matrix(x)
  if (any(is.na(x))){
      NAs <- apply(matrix(apply(x, 1, is.na), nrow=dim(x)[2]),2,any)
      y <- y[!NAs]
      x <- x[!NAs,]
  }
  fmin <- function(beta, X, y, w) {
      p <- plogis(X %*% beta)
      -sum(2 * w * ifelse(y, log(p), log(1-p)))
  }
  gmin <- function(beta, X, y, w) {
      eta <- X %*% beta; p <- plogis(eta)
      -2 * matrix(w *dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% X
  }
  if(is.null(dim(x))) dim(x) <- c(length(x), 1)
  dn <- dimnames(x)[[2]]
  if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
  p <- ncol(x) + intercept
  if(intercept) {x <- cbind(1, x); dn <- c("(Intercept)", dn)}
  if(is.factor(y)) y <- (unclass(y) != 1)
  fit <- optim(start, fmin, gmin, X = x, y = y, w = wt,
               method = "BFGS", ...)
  invisible(fit)
}

