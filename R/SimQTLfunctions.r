######################################################################
#
# SimQTLfunctions.R
#
# sim.cross.exp, sim.cross.exp.bc are expansions to sim.cross and sim.cross.bc
# written by Karl W Broam for R/q tl and modified by Sandra L. Taylor for the 
# twopartqtl package
# copyright (c) 2001-6, Karl W Broman
# last modified Sep, 2006
# first written Apr, 2001
#
# modifications copyright(c) 2009, Sandra L Taylor
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# rtruncnorm, rpointmass written by Sandra L Taylor for the twopartqtl package
# copyright (c) 2009 Sandra L Taylor
# first written Aug, 2009
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the twopartqtl package
# Contains: sim.cross.exp, sim.cross.exp.bc, rtruncnorm, rpointmass
#
######################################################################

sim.cross.exp <- function(map, model = NULL, n.ind = 100, type = c("bc"), 
    error.prob = 0, missing.prob = 0, keep.qtlgeno = TRUE, keep.errorind = TRUE, 
    m = 0, p = 0,  map.function = c("haldane", "kosambi", "c-f", "morgan"), dist, params, logNorm=FALSE){
    type <- match.arg(type)
    map.function <- match.arg(map.function)
    if (error.prob < 1e-50)
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    if (!is.null(model) && is.matrix(model))
        model <- model[order(model[, 1], model[, 2]), ]
    if (type == "bc")
        cross <- sim.cross.bc.exp(map, model, n.ind, error.prob,
            missing.prob, keep.errorind, m, p, map.function, dist, params, logNorm)
    else if (type == "f2")
        cross <- sim.cross.f2(map, model, n.ind, error.prob,
            missing.prob, keep.errorind, m, p, map.function)
    else cross <- sim.cross.4way(map, model, n.ind, error.prob,
        missing.prob, keep.errorind, m, p, map.function)
    qtlgeno <- NULL
    for (i in 1:nchr(cross)) {
        o <- grep("^QTL[0-9]+", colnames(cross$geno[[i]]$data))
        if (length(o) != 0) {
            qtlgeno <- cbind(qtlgeno, cross$geno[[i]]$data[,
                o, drop = FALSE])
            cross$geno[[i]]$data <- cross$geno[[i]]$data[, -o,
                drop = FALSE]
            if (is.matrix(cross$geno[[i]]$map))
                cross$geno[[i]]$map <- cross$geno[[i]]$map[,
                  -o, drop = FALSE]
            else cross$geno[[i]]$map <- cross$geno[[i]]$map[-o]
        }
    }
    if (keep.qtlgeno)
        cross$qtlgeno <- qtlgeno
    for (i in 1:nchr(cross)) storage.mode(cross$geno[[i]]$data) <- "integer"
    cross
}

sim.cross.bc.exp <- function(map, model, n.ind, error.prob, missing.prob,
                    keep.errorind, m, p, map.function, dist, params, logNorm){
if (map.function == "kosambi")
    mf <- mf.k
else if (map.function == "c-f")
    mf <- mf.cf
else if (map.function == "morgan")
    mf <- mf.m
else mf <- mf.h
if (any(sapply(map, is.matrix)))
    stop("Map must not be sex-specific.")
n.chr <- length(map)
if (is.null(model))
    n.qtl <- 0
else {
    if (!((!is.matrix(model) && length(model) >= 3) || (is.matrix(model) &&
        ncol(model) >= 3)))
        stop("Model must be a matrix with 3 columns (chr, pos and effect).")
    if (!is.matrix(model))
        model <- rbind(model)
    n.qtl <- nrow(model)
    if (any(model[, 1] < 0 | model[, 1] > n.chr))
        stop("Chromosome indicators in model matrix out of range.")
    model[, 2] <- model[, 2] + 1e-14
}
if (n.qtl > 0) {
    for (i in 1:n.qtl) {
        temp <- map[[model[i, 1]]]
        if (model[i, 2] < min(temp)) {
            temp <- c(model[i, 2], temp)
            names(temp)[1] <- paste("QTL", i, sep = "")
        }
        else if (model[i, 2] > max(temp)) {
            temp <- c(temp, model[i, 2])
            names(temp)[length(temp)] <- paste("QTL", i,
              sep = "")
        }
        else {
            j <- max((seq(along = temp))[temp < model[i,
              2]])
            temp <- c(temp[1:j], model[i, 2], temp[(j + 1):length(temp)])
            names(temp)[j + 1] <- paste("QTL", i, sep = "")
        }
        map[[model[i, 1]]] <- temp
    }
}
geno <- vector("list", n.chr)
names(geno) <- names(map)
n.mar <- sapply(map, length)
mar.names <- lapply(map, names)
chr.type <- sapply(map, function(a) ifelse(class(a) == "X",
    "X", "A"))
for (i in 1:n.chr) {
    thedata <- sim.bcg(n.ind, map[[i]], m, p, map.function)
    dimnames(thedata) <- list(NULL, mar.names[[i]])
    geno[[i]] <- list(data = thedata, map = map[[i]])
    class(geno[[i]]) <- chr.type[i]
    class(geno[[i]]$map) <- NULL
}
if (dist!="pointmass"){
    pheno <- switch(dist, norm = rnorm(n.ind, mean=params[1], sd=params[2]),
            lognorm = exp(rnorm(n.ind, mean=params[1], sd=params[2])),
            gamma = rgamma(n.ind, shape=params[1], scale=params[2]),
            truncNorm = rtruncnorm(n.ind, mean=params[1], sd=params[2], floor=params[3], logNorm),
            t = rt(n.ind, df=params),
            cauchy = rcauchy(n.ind, location=params[1], scale=params[2]),
            logistic = rlogis(n.ind, location=params[1], scale=params[2]),
            expon = rexp(n.ind, rate=params), uniform = runif(n.ind, min=params[1], max=params[2]),
            pois = rpois(n.ind, lambda=params))
}        
if (n.qtl > 0) {
    QTL.chr <- QTL.loc <- NULL
    for (i in 1:n.chr) {
        o <- grep("^QTL[0-9]+", mar.names[[i]])
        if (length(o) > 0) {
            QTL.chr <- c(QTL.chr, rep(i, length(o)))
            QTL.loc <- c(QTL.loc, o)
        }
    }
  if (dist != "pointmass"){
  # Note geno==1 is homozygote and 2 is heterozygote
    for (i in 1:n.qtl) {
        QTL.geno <- geno[[QTL.chr[i]]]$data[, QTL.loc[i]]
        pheno[QTL.geno == 2] <- pheno[QTL.geno == 2] + model[i,3]
    }
  }
  else {
     combParams <- matrix(nrow=n.ind,ncol=2)
     QTL.geno <- geno[[QTL.chr[1]]]$data[, QTL.loc[1]]
     combParams[QTL.geno==1,1] <- params[1]
     combParams[QTL.geno==1,2] <- params[2]
     combParams[QTL.geno==2,1] <- params[1]+model[1,4]
     combParams[QTL.geno==2,2] <- params[2]+model[1,3]
   
     if (n.qtl > 1){
      for (i in 2:n.qtl){
        QTL.geno <- geno[[QTL.chr[i]]]$data[, QTL.loc[i]]
        combParams[QTL.geno==2,1] <- combParams[QTL.geno==2,1]+model[i,4]
        combParams[QTL.geno==2,2] <- combParams[QTL.geno==2,2]+model[i,3]
      }
     }
    pheno <- sapply(1:n.ind, function(x) rpointmass(1,mix=combParams[x,1],mean=combParams[x,2],sd=params[3], logNorm))
  }    
}  

n.mar <- sapply(geno, function(a) length(a$map))
if (error.prob > 0) {
    for (i in 1:n.chr) {
        a <- sample(0:1, n.mar[i] * n.ind, repl = TRUE, prob = c(1 -
            error.prob, error.prob))
        geno[[i]]$data[a == 1] <- 3 - geno[[i]]$data[a ==
            1]
        if (keep.errorind) {
            errors <- matrix(0, n.ind, n.mar[i])
            errors[a == 1] <- 1
            colnames(errors) <- colnames(geno[[i]]$data)
            geno[[i]]$errors <- errors
        }
    }
}
if (missing.prob > 0) {
    for (i in 1:n.chr) {
        o <- grep("^QTL[0-9]+", mar.names[[i]])
        if (length(o) > 0)
            x <- geno[[i]]$data[, o]
        geno[[i]]$data[sample(c(TRUE, FALSE), n.mar[i] *
            n.ind, repl = TRUE, prob = c(missing.prob, 1 -
            missing.prob))] <- NA
        if (length(o) > 0)
            geno[[i]]$data[, o] <- x
    }
}
pheno <- data.frame(phenotype = pheno)
cross <- list(geno = geno, pheno = pheno)
class(cross) <- c("bc", "cross")
cross
}


rtruncnorm <- function(num.ind, mean, sd, floor,logNorm){
   data <- NULL
   for (i in 1:num.ind){
         RV <- floor-1
         while (RV <= floor){
             RV <- rnorm(1, mean=mean, sd)
         }
         if (logNorm) data[i] <- exp(RV)
         else data[i] <- RV
   }
   return(data)
}

rpointmass <- function(num.ind, mix, mean, sd,logNorm){
    data <- rep(NA,num.ind)
    for (i in 1:num.ind){
      if (runif(1) <= mix) data[i] = 0 
      else data[i] = rtruncnorm(1,mean,sd,floor=0,logNorm) 
     }
     return(data)
}
