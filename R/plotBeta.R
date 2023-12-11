#' @include utils.R


#' @title Plot distribution of coefficients across features.
#' @name plotBeta
#' @export
plotBeta <- function(fit){
  for(cc in 1:ncol(fit$coefficients)){
    curBeta <- fit$coefficients[,cc]
    plot(density(curBeta), main=colnames(fit$coefficients)[cc])
    abline(v=.getMode(beta=curBeta, n=nrow(fit$design)),
           col="orange", lwd=2)
  }
}