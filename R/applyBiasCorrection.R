#' @include utils.R

#' @title Apply bias correction to estimated regression coefficients.
#' @name applyBiasCorrection
#' @export
applyBiasCorrection <- function(fit){
  par(mfrow=c(3,3),
      cex.lab = 1, cex.main = 1.2, cex.axis = 1, mar = c(2.5, 2.5, 
                                                         1.6, 1.1), 
      mgp = c(1.5, 0.5, 0))
  n <- nrow(fit$design)
  biasCorrCoef <- apply(fit$coefficients, 2, function(x){
    x - .getMode(x, n)
  })
  fit$coefficients <- biasCorrCoef
  par(mfrow=c(1,1))
  return(fit)
}