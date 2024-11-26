.getMode <- function(beta, n){
  suppressMessages(mode <- modeest::mlv(sqrt(n) * beta, 
                                        method = "meanshift", kernel = "gaussian")/sqrt(n))
  return(mode)
}


.estimateNBDispersion <- function(counts,
                                 design,
                                 normalization="TMM"){
  require(edgeR)
  d <- edgeR::DGEList(counts)
  if(normalization == "TMM") d <- edgeR::calcNormFactors(d)
  d <- edgeR::estimateDisp(d, design)
  return(d$tagwise.dispersion)
}


### NON-PARAMETRIC BOOTSTRAP
### for vector of coefs
.calcBias <- function(beta, ids, n){
  return(.getMode(beta[ids], n=n))
}
.nonParametricBootBeta <- function(beta, n, R=4e3){
  library(boot)
  bootRes <- boot(beta, # we are resampling rows in 'boot'
                  .calcBias, 
                  R=4e3, 
                  n=n)
  varMode <- var(bootRes$t[,1])
  return(varMode)
}
# ### for matrix of coefs
# .calcBias <- function(betaMat, ids, n){
#   apply(betaMat[,ids], 1, function(x){ # on the rows since transposed in input
#     .getMode(x, n)
#   })
# }
# .nonParametricBootBeta <- function(betaMat, n, R=4e3){
#   library(boot)
#   hlp <- boot(t(betaMat), # we are resampling rows in 'boot'
#               .calcBias, 
#               R=4e3, 
#               n=n)
#   varMode <- apply(hlp$t,2,var)
#   return(varMode)
# }


### PARAMETRIC BOOTSTRAP
.parametricBootstrap <- function(beta, design, sigma2, weights, n, R=4e3, L=NULL){
  
  ### for dev:
  # beta <- fit$coefficients
  # sigma2 <- sigma2post
  # weights <- v$weights
  # n <- nrow(design)
  # L <- contrastMatrix
  
  X <- design

  ## loop over populations to get variance-covariance
  covBetaList <- list()
  for(pp in 1:nrow(beta)){
    curW <- diag(weights[pp,])
    curCovBeta <- sigma2[pp] * solve(t(X) %*% curW %*% X) # agrees with topTable
    covBetaList[[pp]] <- curCovBeta
  }
  
  ## simulate for each population
  if(is.null(L)){
    simBetaList <- list()
    for(pp in 1:nrow(beta)){
      curSimBeta <- mixtools::rmvnorm(n=R, mu=beta[pp,], sigma=covBetaList[[pp]])
      simBetaList[[pp]] <- curSimBeta
    }
  } else {
    simBetaList <- list()
    for(pp in 1:nrow(beta)){
      curSimBeta <- mixtools::rmvnorm(n=R, mu=beta[pp,], sigma=t(L) %*% covBetaList[[pp]] %*% L)
      simBetaList[[pp]] <- curSimBeta
    }
  }
  
  ## calculate mode for each bootstrap iteration
  modeMat <- matrix(NA, nrow=R, ncol=ncol(beta))
  for(bb in 1:R){
    curBeta <- do.call(rbind, lapply(simBetaList, function(x) x[bb,]))
    curModes <- apply(curBeta, 2, function(x) .getMode(x, n))
    modeMat[bb,] <- curModes
  }
  
  ## calculate covariance
  covMat <- matrix(NA, nrow=nrow(beta), ncol=ncol(beta),
                   dimnames=dimnames(beta))
  for(pp in 1:length(simBetaList)){
    for(cc in 1:ncol(modeMat)){
      covMat[pp,cc] <- cov(simBetaList[[pp]][,cc], modeMat[,cc])
    }
  }
  
  ## variance of mode
  varMode <- apply(modeMat, 2, var)
  names(varMode) <- colnames(beta)
  
  return(list("varMode"=varMode,
              "covMode"=covMat))
}

