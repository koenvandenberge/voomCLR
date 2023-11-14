.getMode <- function(beta, n){
  suppressMessages(mode <- modeest::mlv(sqrt(n) * beta, 
                                        method = "meanshift", kernel = "gaussian")/sqrt(n))
}


.estimateNBDispersion <- function(counts,
                                 design,
                                 normalization="TMM"){
  require(edgeR)
  d <- DGEList(counts)
  if(normalization == "TMM") d <- calcNormFactors(d)
  d <- estimateDisp(d, design)
  return(d$tagwise.dispersion)
}
