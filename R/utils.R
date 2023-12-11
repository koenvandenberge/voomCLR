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
