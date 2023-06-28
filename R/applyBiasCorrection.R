#' @include utils.R

#' @name applyBiasCorrection
#' @export
applyBiasCorrection <- function(fit){
  n <- nrow(fit$design)
  biasCorrCoef <- apply(fit$coefficients, 2, function(x){
    x - .getMode(x, n)
  })
  fit$coefficients <- biasCorrCoef
  return(fit)
}