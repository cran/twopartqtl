################################################################################
#
# pmscan.R
#
# copyright (c) 2009, Sandra L. Taylor
# 
# first written Aug, 2009
#
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the twopartqtl package
# Contains: pmscan
#
################################################################################

pmscan <- function(cross, chr, pheno.col=1,pm.value=0, maxit=4000, tol=1e-4,
          addcovar, verbose=FALSE,imp.method = c("imp", "argmax"), error.prob=1e-4,
          map.function = c("haldane", "kosambi", "c-v", "morgan"), n.perm, use.log=FALSE){

 if (!missing(chr))   cross <- subset(cross, chr)
 
 if (any(pheno.col < 1 | pheno.col > nphe(cross)))
        stop("pheno.col values should be between 1 and the no. phenotypes")

 num.ind <- length(cross$pheno[,pheno.col])

 if (use.log){
        for (i in 1:num.ind){
          if(!is.na(cross$pheno[i,pheno.col])){
            if (cross$pheno[i,pheno.col]!=pm.value) cross$pheno[i,pheno.col] <- log(cross$pheno[i,pheno.col])
          }
        }
   }
  pheno <-cross$pheno[,pheno.col]
  
  if (any(is.na(pheno))) {
        cross <- subset(cross, ind = (!is.na(pheno)))    
        addcovar <- addcovar[!is.na(pheno),]
        pheno <- pheno[!is.na(pheno)]      
  }
                                                 
   num.ind <- length(pheno)

   if (!missing(n.perm) && n.perm > 0) {
        results <- matrix(ncol = 1, nrow = n.perm)
        for (i in 1:n.perm) {
            o <- sample(length(pheno))
            pheno <- pheno[o]
            cross$pheno[, pheno.col] <- pheno
            temp <- pmscan(cross, chr, pheno.col = pheno.col, pm.value=pm.value,
                maxit=maxit, tol=tol,  addcovar=addcovar, imp.method = imp.method,
                error.prob = error.prob, map.function = map.function, use.log=FALSE)
            results[i, 1] <- max(temp[, 3], na.rm = TRUE)
        }
        class(results) <- c("scanoneperm", "matrix")
        return(results)
    }
  
  n.ind <- nind(cross)
  n.chr <- nchr(cross)
  type <- class(cross)[1]
  n.addcovar <- ncol(addcovar)

  geno <- pull.geno(cross)
    if (any(is.na(geno))){
          geno <- pull.geno(fill.geno(cross, method = imp.method,
              error.prob = error.prob, map.function = map.function))
    }

  pointmass <- rep(0, length(pheno))
  pointmass[pheno==pm.value] <- 1
  pheno[pointmass==1] <- 99999

  results <- coeff <- NULL
  gamma0 <- gamma1 <- beta0 <- beta1 <- NULL

  llik0 <- NA
 for (i in 1:n.chr){
    chrtype <- class(cross$geno[[i]]) 
    sexpgm <- NULL

    if (!("prob" %in% names(cross$geno[[i]]))) {
          warning("First running calc.genoprob.")
          cross <- calc.genoprob(cross, step=1)
      }
      genoprob <- cross$geno[[i]]$prob
      n.pos <- ncol(genoprob)
    if ("map" %in% names(attributes(cross$geno[[i]]$prob)))
        map <- attr(cross$geno[[i]]$prob, "map")
     else {
          stp <- attr(cross$geno[[i]]$prob, "step")
          oe <- attr(cross$geno[[i]]$prob, "off.end")
          if ("stepwidth" %in% names(attributes(cross$geno[[i]]$prob)))
            stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
          else stpw <- "fixed"
          map <- create.map(cross$geno[[i]]$map, stp, oe,
            stpw)
      }
    if (is.matrix(map))
      map <- map[1, ]

    gen.names <- getgenonames(type, chrtype, "full", sexpgm, attributes(cross))
    n.gen <- length(gen.names)
    
    if (is.na(llik0)) {
      null.fit.binary <- logitreg(addcovar, pointmass)      
      nullgamma <- null.fit.binary$par   
      fit <- exp(addcovar%*%nullgamma[-1]+nullgamma[1])
      binary.fitted <- fit/(1+fit)
      llik0.binary  <- sum(pointmass*log10(binary.fitted)+(1-pointmass)*log10(1-binary.fitted))
      bin_start <- rep(nullgamma[1], n.gen)
      bin_start <- c(bin_start, nullgamma[-1])
  
      cont <- pheno[pointmass==0]
      n.ind.cont <- length(cont)
      cont.ac <- as.matrix(addcovar[pointmass==0,])
      null.fit.cont <- lm(cont~cont.ac)
      resid0 <- null.fit.cont$resid
      sig0 <- sqrt(sum((resid0^2)/n.ind.cont))
      llik0.cont <- sum(log10(dnorm(resid0, 0 ,sig0)))
      llik0 <- llik0.binary + llik0.cont
     }
   
    out <- .C("R_pmscan_covar", as.integer(n.ind), as.integer(n.ind.cont), as.integer(n.pos),
      as.integer(n.gen), as.double(genoprob), as.double(addcovar), as.double(cont.ac),
      as.integer(n.addcovar), as.double(pheno), as.double(cont), as.integer(pointmass),
      as.double(bin_start), result = as.double(rep(0, n.pos)),
      gamma0 = as.double(rep(0, n.pos)), gamma1 = as.double(rep(0, n.pos)),
      beta0 = as.double(rep(0, n.pos)), beta1 = as.double(rep(0, n.pos)),
      as.integer(maxit), as.double(tol), as.integer(verbose),PACKAGE="twopartqtl")
 
    z <- matrix(out$result, nrow = n.pos)
    z[z[,1]==0,] <- NA
    z[, 1] <- z[, 1] - llik0
    z[is.na(z[, 1]), 1] <- 0
    z <- z[, 1, drop = FALSE]    

    colnames(z)[1] <- "lod"
    w <- names(map)
    o <- grep("^loc-*[0-9]+", w)
    if (length(o) > 0)
        w[o] <- paste("c", names(cross$geno)[i], ".", w[o],
            sep = "")
    z <- data.frame(chr = rep(names(cross$geno)[i], length(map)),
        pos = as.numeric(map), z)
    rownames(z) <- w
    results <- rbind(results, z)
    
    s <- cbind(out$gamma0, out$gamma1, out$beta0, out$beta1)
    colnames(s) <- c("gamma0", "gamma1", "beta0", "beta1")
    s <- data.frame(chr = rep(names(cross$geno)[i], length(map)),
        pos = as.numeric(map), s)
    rownames(s) <- w
    coeff <- rbind(coeff, s)
        
  }     
    class(results) <- c("scanone", "data.frame")
    attr(results, "method") <- "em"
    attr(results, "type") <- type
    attr(results, "model") <- "two.part"
    return(results) 
} 

