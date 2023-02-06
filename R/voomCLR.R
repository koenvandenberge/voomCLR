## example
# set.seed(495212344)
# n <- 40 # sample size
# P <- 10 # number of cell types
# mu0 <- rnbinom(n=P, size=1/2, mu=400)
# mu0 # absolute counts in group 0
# beta <- rlnorm(n=P, meanlog = 0, sdlog=2) * # these are log-fold-changes
#   rbinom(n=P, size=1, prob=.15) *
#   sample(c(-1,1), size=P, replace=TRUE) # fold change on log scale
# mu1 <- exp(beta) * mu0 # because we want log(mu2/mu1) = beta
# relAbundances <- data.frame(g0=mu0/sum(mu0),
#                             g1=mu1/sum(mu1)) # relative abundance of absolute count
# # relative abundance information (observed data in typical experiment)
# Y0 <- rmultinom(n=10, size=1e4, prob=relAbundances$g0)
# Y1 <- rmultinom(n=10, size=1e4, prob=relAbundances$g1)
# Y <- cbind(Y0, Y1)
# group <- factor(rep(0:1, each=10))
# design <- model.matrix(~group)
# v <- voomCLR(counts = Y,
#              design = design,
#              lib.size = NULL)
# fit <- lmFit(v, design)
# fit <- applyBiasCorrection(fit)
# fit <- eBayes(fit)

## for interactively running the code:
# counts=Y
# design=design
# lib.size=NULL
# normalize.method="none"
# block=NULL
# correlation=NULL
# weights=NULL
# span=0.5
# plot=FALSE
# save.plot=FALSE
# varCalc="empirical"




voomCLR <- function(counts,
                    design=NULL,
                    lib.size=NULL,
                    normalize.method="none",
                    block=NULL,
                    correlation=NULL,
                    weights=NULL,
                    span=0.5,
                    plot=FALSE,
                    save.plot=FALSE,
                    varCalc="empirical")
  #	Linear modelling of count data with mean-variance modelling at the observation level.
  #	Creates an EList object for entry to lmFit() etc in the limma pipeline.
  #	Gordon Smyth and Charity Law
  #	Created 22 June 2011.  Last modified 1 May 2021.
  # Modified by Koen Van den Berge:
  # - CLR transformation instead of CPM
  # - Allow for analytical calculation of standard deviation
  # TODO:
  # - check difference between lib.size=1 and lib.size=geoMeans
  # - should we be using CLR (instead of CPM) on x-axis when fitting trend?
  # - check CLR backtransform
  # - limma-trend? Allows for intensity-dependent prior variance per gene.
  # - analytical weights are much different in magnitude as compared to using the trend.
{
  out <- list()
  counts <- as.matrix(counts)
  
  #	Check counts
  n <- nrow(counts)
  if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")
  m <- min(counts)
  if(is.na(m)) stop("NA counts not allowed")
  if(m < 0) stop("Negative counts not allowed")
  
  #	Check design
  if(is.null(design)) {
    design <- matrix(1,ncol(counts),1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  
  
  #	Fit linear model to CLR
  geoMeans <- exp(colMeans(log(counts+0.5)))
  #	Check lib.size
  if(is.null(lib.size)){
    lib.size <- geoMeans
  }
  y <- t(log(t(counts+0.5)/geoMeans))
  y <- normalizeBetweenArrays(y,method=normalize.method)
  fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
  if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)
  
  #	If no replication found, set all weight to 1
  NWithReps <- sum(fit$df.residual > 0L)
  if(NWithReps < 2L) {
    if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
    if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    out$design <- design
    if(is.null(out$targets))
      out$targets <- data.frame(lib.size=lib.size)
    else
      out$targets$lib.size <- lib.size
    return(new("EList",out))
  }
  
  #	Fit lowess trend to sqrt-standard-deviations by log-count-size
  # TODO: check CLR vs CPM in next line. The mean is being back-transformed.
  # sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)
  sx <- fit$Amean+mean(log(lib.size)) # note that lib.size=geoMeans
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts)==0
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx,sy,f=span)
  if(plot) {
    plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
    title("voom: Mean-variance trend")
    lines(l,col="red")
  }
  
  #	Make interpolating rule
  #	Special treatment of zero counts is now removed;
  #	instead zero counts get same variance as smallest gene average.
  #	l$x <- c(0.5^0.25, l$x)
  #	l$x <- c(log2(0.5), l$x)
  #	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
  #	var0 <- max(var0,1e-6)
  #	l$y <- c(var0, l$y)

  
  #	Find individual quarter-root fitted counts
  if(fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
  } else {
    fitted.values <- fit$coefficients %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
  fitted.logcount <- log2(fitted.count)
  
  #	Apply trend to individual observations
  if(varCalc=="empirical"){
    f <- approxfun(l, rule=2, ties=list("ordered",mean))
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
  } else if(varCalc=="analytical"){
    # note that lib.size=geoMeans
    geoMeanMat <- matrix(lib.size, nrow=nrow(counts), ncol=ncol(counts),
                         byrow=TRUE)
    yBarFit <- exp(fitted.values)*geoMeanMat
    # calculate group-wise mean in simple designs
    # yBarDiffs <- counts %*% design %*% solve(crossprod(design))
    # yBarDiffs[,2:ncol(yBarDiffs)] <- yBarDiffs[,-1] + yBarDiffs[,1]
    # yBar <- yBarDiffs
    # yBar[yBar==0] <- 1e-10
    # calculate variances based on Delta method approximation
    estVar <- ((nrow(counts)-1)/nrow(counts))^2 * (1/(yBarFit))
    # calculate weights and expand dimension to count dimension
    w <- 1/estVar
  }
  
  
  #	Output
  out$E <- y
  out$weights <- w
  out$design <- design
  if(is.null(out$targets))
    out$targets <- data.frame(lib.size=lib.size)
  else
    out$targets$lib.size <- lib.size
  if(save.plot) {
    out$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )")
    out$voom.line <- l
  }
  
  new("EList",out)
}


getMode <- function(beta, n){
  suppressMessages(mode <- modeest::mlv(sqrt(n) * beta, 
                                        method = "meanshift", kernel = "gaussian")/sqrt(n))
}

applyBiasCorrection <- function(fit){
  n <- nrow(fit$design)
  biasCorrCoef <- apply(fit$coefficients, 2, function(x){
    x - getMode(x, n)
  })
  fit$coefficients <- biasCorrCoef
  return(fit)
}


#### workflow
## requires two objects:
##  - design matrix as constructed using model.matrix
##  - count matrix where rows are populations and columns are samples
# library(limma)
# v <- voomCLR(counts = countMatrix,
#              design = design,
#              lib.size = 1)
# fit <- lmFit(v, design)
# fit <- applyBiasCorrection(fit)
# fit <- eBayes(fit)
# testing is custom: fill in based on hypothesis of interest.
# tt <- topTable(fit, ....)



## in the case of random effects
# you need one extra input: the covariate you'd like to add as random effect.
# here I assume this is in the 'patient' object.
# library(limma)
# v <- voomCLR(counts = countMatrix,
#              design = design,
#              lib.size = 1)
# cf <- duplicateCorrelation(v, design, block=patient)
# v <- voomCLR(counts = countMatrix,
#              design = design,
#              lib.size = 1,
#              block = patient, 
#              correlation = cf$consensus)
# cf <- duplicateCorrelation(v, design, block=patient)
# fit <- lmFit(v, design,
#              block = patient, 
#              correlation = cf$consensus)
# fit <- applyBiasCorrection(fit)
# fit <- eBayes(fit)