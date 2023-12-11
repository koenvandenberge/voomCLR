.getMode <- function(beta, n){
  suppressMessages(mode <- modeest::mlv(sqrt(n) * beta, 
                                        method = "meanshift", kernel = "gaussian")/sqrt(n))
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
  return(.getMode(betaMat[ids], n=n))
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


