#' @include utils.R

### for dev:
# eb=fit[ebcols]
# A=fit$Amean
# sigma2post=fit$s2.post
# design = fit$design
# fit=fit[c("coefficients","stdev.unscaled")]
# coef=coef
# number=number
# genelist=genelist
# adjust.method=adjust.method
# sort.by=sort.by
# resort.by=resort.by
# p.value=p.value
# lfc=lfc
# confint=confint
# bootstrap <- "parametric"
# voomWeights <- voomWeights


.toptableTBC <- function(fit,
                         n,
                         design=NULL,
                         coef=1,
                         number=10,
                         genelist=NULL,
                         A=NULL,
                         eb=NULL,
                         adjust.method="BH",
                         sort.by="p",
                         resort.by=NULL,
                         p.value=1,
                         lfc=0,
                         confint=FALSE,
                         bootstrap=FALSE,
                         voomWeights=NULL,
                         sigma2post=NULL,
                         ...){
  #	Summary table of top genes for a single coefficient
  #	Gordon Smyth
  #	Created 21 Nov 2002. Was called toptable() until 1 Feb 2018. Last revised 12 Apr 2020.

  #	Check fit
  fit$coefficients <- as.matrix(fit$coefficients)
  rn <- rownames(fit$coefficients)
  
  #	Check coef is length 1
  if(length(coef)>1) {
    coef <- coef[1]
    warning("Treat is for single coefficients: only first value of coef being used")
  }
  
  #	Ensure genelist is a data.frame
  if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ID=genelist,stringsAsFactors=FALSE)
  
  #	Check rownames
  if(is.null(rn)){
    rn <- 1:nrow(fit$coefficients)
    } else {
      if(anyDuplicated(rn)) {
        if(is.null(genelist))
          genelist <- data.frame(ID=rn,stringsAsFactors=FALSE)
        else
          if("ID" %in% names(genelist))
            genelist$ID0 <- rn
          else
            genelist$ID <- rn
          rn <- 1:nrow(fit$coefficients)
    }
  }
  
  #	Check sort.by
  sort.by <- match.arg(sort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","none"))
  if(sort.by=="M") sort.by="logFC"
  if(sort.by=="A" || sort.by=="Amean") sort.by <- "AveExpr"
  if(sort.by=="T") sort.by <- "t"
  if(sort.by=="p") sort.by <- "P"
  
  #	Check resort.by
  if(!is.null(resort.by)) {
    resort.by <- match.arg(resort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t"))
    if(resort.by=="M") resort.by <- "logFC"
    if(resort.by=="A" || resort.by=="Amean") resort.by <- "AveExpr"
    if(resort.by=="p") resort.by <- "P"
    if(resort.by=="T") resort.by <- "t"
  }
  
  #	Check A
  if(is.null(A)) {
    if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
  } else {
    if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
  }
  
  #	Check for lods component
  if(is.null(eb$lods)) {
    if(sort.by=="B") stop("Trying to sort.by B, but B-statistic (lods) not found in MArrayLM object",call.=FALSE)
    if(!is.null(resort.by)) if(resort.by=="B") stop("Trying to resort.by B, but B-statistic (lods) not found in MArrayLM object",call.=FALSE)
    include.B <- FALSE
  } else {
    include.B <- TRUE
  }
  
  # Bias calculation
  bias <- apply(fit$coefficients, 2, function(x){
    .getMode(x, n=n)
  })
  
  #	Extract statistics for table
  M <- fit$coefficients[,coef] - bias[coef]
  se_coef <- as.matrix(fit$coefficients / eb$t)[, coef]
  if(!is.null(bootstrap)){
    if(bootstrap == "nonparametric"){
      var_mode <- suppressWarnings(.nonParametricBootBeta(fit$coefficients[,coef], n))
      varCombined <- se_coef^2 + var_mode
      tstat <- as.matrix(M / sqrt(varCombined))
    } else if(bootstrap == "parametric"){
      var_mode <- .parametricBootstrap(beta = fit$coefficients, 
                                       design = design,
                                       sigma2 = sigma2post,
                                       weights = voomWeights,
                                       n=n)
      varCombined <- se_coef^2 + var_mode[coef]
      tstat <- as.matrix(M / sqrt(varCombined))
    }
  } else { # bootstrap is NULL
    tstat <- as.matrix(M / se_coef)
  }
  df_coef <- matrix(eb$df.total, dimnames = list(names(tstat)))
  P.Value <- as.matrix(2*pt(abs(tstat), df = df_coef, lower.tail = F))

  #	Apply multiple testing adjustment
  adj.P.Value <- p.adjust(P.Value,method=adjust.method)
  
  #	Thin out fit by p.value and lfc thresholds	
  if(p.value < 1 | lfc > 0) {
    sig <- (adj.P.Value <= p.value) & (abs(M) >= lfc)
    if(anyNA(sig)) sig[is.na(sig)] <- FALSE
    if(!any(sig)) return(data.frame())
    genelist <- genelist[sig,,drop=FALSE]
    M <- M[sig]
    A <- A[sig]
    tstat <- tstat[sig]
    P.Value <- P.Value[sig]
    adj.P.Value <- adj.P.Value[sig]
    # if(include.B) B <- B[sig]
    rn <- rn[sig]
  }
  
  #	Are enough rows left?
  if(length(M) < number) number <- length(M)
  if(number < 1) return(data.frame())
  
  #	Select top rows
  ord <- switch(sort.by,
                logFC=order(abs(M),decreasing=TRUE),
                AveExpr=order(A,decreasing=TRUE),
                P=order(P.Value,decreasing=FALSE),
                t=order(abs(tstat),decreasing=TRUE),
                # B=order(B,decreasing=TRUE),
                none=1:length(M)
  )
  top <- ord[1:number]
  
  #	Assemble output data.frame
  if(is.null(genelist))
    tab <- data.frame(logFC=M[top])
  else {
    tab <- data.frame(genelist[top,,drop=FALSE],logFC=M[top],stringsAsFactors=FALSE)
  }
  if(confint) {
    if(is.numeric(confint)) alpha <- (1+confint[1])/2 else alpha <- 0.975
    margin.error <- sqrt(eb$s2.post[top])*fit$stdev.unscaled[top,coef]*qt(alpha,df=eb$df.total[top])
    tab$CI.L <- M[top]-margin.error
    tab$CI.R <- M[top]+margin.error
  }
  if(!is.null(A)) tab$AveExpr <- A[top]
  tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],adj.P.Val=adj.P.Value[top])
  # if(include.B) tab$B <- B[top]
  rownames(tab) <- rn[top]
  
  #	Resort table
  if(!is.null(resort.by)) {
    ord <- switch(resort.by,
                  logFC=order(tab$logFC,decreasing=TRUE),
                  AveExpr=order(tab$AveExpr,decreasing=TRUE),
                  P=order(tab$P.Value,decreasing=FALSE),
                  t=order(tab$t,decreasing=TRUE)
    )
    tab <- tab[ord,]
  }
  
  tab
}






#' @name topTableBC
#' @title Table of Top Genes from Linear Model Fit
#' 
#' @description Extract a table of the top-ranked genes from a linear model fit.
#' @usage 
#' topTable(fit, coef = NULL, number = 10, genelist = fit$genes,
#'            adjust.method = "BH", sort.by = "B", resort.by = NULL,
#'            p.value = 1, fc = NULL, lfc = NULL, confint = FALSE)
#' topTableF(fit, number = 10, genelist = fit$genes,
#'             adjust.method="BH", sort.by="F",
#'             p.value = 1, fc = NULL, lfc = NULL)
#' topTreat(fit, coef = 1, sort.by = "p", resort.by = NULL, \dots)
#' @param fit list containing a linear model fit produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{mrlm}.
#'     For \code{topTable}, \code{fit} should be an object of class \code{MArrayLM} as produced by \code{lmFit} and \code{eBayes}.
#' @param coef column number or column name specifying which coefficient or contrast of the linear model is of interest. For \code{topTable}, can also be a vector of column subscripts, in which case the gene ranking is by F-statistic for that set of contrasts.
#' @param number maximum number of genes to list.
#' @param genelist data frame or character vector containing gene information.
#'     For \code{topTable} only, this defaults to \code{fit$genes}.
#' @param adjust.method method used to adjust the p-values for multiple testing.  Options, in increasing conservatism, include \code{"none"}, \code{"BH"}, \code{"BY"} and \code{"holm"}.
#'     See \code{\link{p.adjust}} for the complete list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.
#' @param sort.by
#'     character string specifying which statistic to rank the genes by.
#'     Possible values for \code{topTable} are \code{"logFC"}, \code{"AveExpr"}, \code{"t"}, \code{"P"}, \code{"p"}, \code{"B"} or \code{"none"}.
#'     (Permitted synonyms are \code{"M"} for \code{"logFC"}, \code{"A"} or \code{"Amean"} for \code{"AveExpr"}, \code{"T"} for \code{"t"} and \code{"p"} for \code{"P"}.)
#'     Possible values for \code{topTableF} are \code{"F"} or \code{"none"}.
#'     \code{topTreat} accepts the same values as \code{topTable} except for \code{"B"}.
#' @param resort.by
#'     character string specifying statistic to sort the selected genes by in the output data.frame.  Possibilities are the same as for \code{sort.by}.
#' @param p.valuecutoff value for adjusted p-values. Only genes with lower p-values are listed.
#' @param fc optional minimum fold-change required.
#' @param lfc
#'     optional minimum log2-fold-change required, equal to \code{log2(fc)}.
#'     \code{fc} and \code{lfc} are alternative ways to specify a fold-change cutoff and, if both are specified, then \code{fc} take precedence.
#'     If specified, then the results from \code{topTable}, \code{topTableF} or \code{topTreat} will include only genes with (at least one) absolute log-fold-change greater than \code{lfc}.
#'     This argument is not normally used with \code{topTreat}, which handles fold-change thresholds differently via the \code{treat} function.
#' @param confint logical, should confidence 95\% intervals be output for \code{logFC}?  Alternatively, can be a numeric value between zero and one specifying the confidence level required.
#' @param bootstrap logical, should bootstrapping be performed to take into account uncertainty of the estimation of the bias correction term?
#' @param dots other \code{topTreat} arguments are passed to \code{topTable}.
#' @return 
#'   A dataframe with a row for the \code{number} top genes and the following columns:
#'   \item{genelist}{one or more columns of probe annotation, if genelist was included as input}
#'   \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast (for \code{topTableF} there may be several columns of log-fold-changes)}
#'   \item{CI.L}{left limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{CI.R}{right limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
#'   \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as \code{Amean} in the \code{MarrayLM} object}
#'   \item{t}{moderated t-statistic (omitted for \code{topTableF})}
#'   \item{F}{moderated F-statistic (omitted for \code{topTable} unless more than one coef is specified)}
#'   \item{P.Value}{raw p-value}
#'   \item{adj.P.Value}{adjusted p-value or q-value}
#'   \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}
#'   
#'   If \code{fit} had unique rownames, then the row.names of the above data.frame are the same in sorted order.
#'   Otherwise, the row.names of the data.frame indicate the row number in \code{fit}.
#'   If \code{fit} had duplicated row names, then these are preserved in the \code{ID} column of the data.frame, or in \code{ID0} if \code{genelist} already contained an \code{ID} column.
#' @details
#'   These functions summarize the linear model fit object produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{mrlm} by selecting the top-ranked genes for any given contrast, or for a set of contrasts.
#'   \code{topTable} assumes that the linear model fit has already been processed by \code{\link{eBayes}}.
#'   \code{topTreat} assumes that the fit has been processed by \code{\link{treat}}.
#'   
#'   If \code{coef} has a single value, then the moderated t-statistics and p-values for that coefficient or contrast are used.
#'   If \code{coef} takes two or more values, the moderated F-statistics for that set of coefficients or contrasts are used.
#'   If \code{coef} is left \code{NULL}, then all the coefficients or contrasts in the fitted model are used, except that any coefficient named \code{(Intercept)} will be removed.
#'   
#'   The p-values for the coefficient/contrast of interest are adjusted for multiple testing by a call to \code{\link[stats]{p.adjust}}.
#'   The \code{"BH"} method, which controls the expected false discovery rate (FDR) below the specified value, is the default adjustment method because it is the most likely to be appropriate for microarray studies.
#'   Note that the adjusted p-values from this method are bounds on the FDR rather than p-values in the usual sense.
#'   Because they relate to FDRs rather than rejection probabilities, they are sometimes called q-values.
#'   See \code{help("p.adjust")} for more information.
#'   
#'   Note, if there is no good evidence for differential expression in the experiment, that it is quite possible for all the adjusted p-values to be large, even for all of them to be equal to one.
#'   It is quite possible for all the adjusted p-values to be equal to one if the smallest p-value is no smaller than \code{1/ngenes} where \code{ngenes} is the number of genes with non-missing p-values.
#'   
#'   The \code{sort.by} argument specifies the criterion used to select the top genes.
#'   The choices are: \code{"logFC"} to sort by the (absolute) coefficient representing the log-fold-change; \code{"A"} to sort by average expression level (over all arrays) in descending order; \code{"T"} or \code{"t"} for absolute t-statistic; \code{"P"} or \code{"p"} for p-values; or \code{"B"} for the \code{lods} or B-statistic.
#'   
#'   Normally the genes appear in order of selection in the output table.
#'   If a different order is wanted, then the \code{resort.by} argument may be useful.
#'   For example, \code{topTable(fit, sort.by="B", resort.by="logFC")} selects the top genes according to log-odds of differential expression and then orders the selected genes by log-ratio in decreasing order.
#'   Or \code{topTable(fit, sort.by="logFC", resort.by="logFC")} would select the genes by absolute log-fold-change and then sort them from most positive to most negative.
#'   
#'   Toptable output for all probes in original (unsorted) order can be obtained by \code{topTable(fit,sort="none",n=Inf)}.
#'   However \code{\link{write.fit}} or \code{\link{write}} may be preferable if the intention is to write the results to a file.
#'   A related method is \code{as.data.frame(fit)} which coerces an \code{MArrayLM} object to a data.frame.
#'   
#'   By default \code{number} probes are listed.
#'   Alternatively, by specifying \code{p.value} and \code{number=Inf}, all genes with adjusted p-values below a specified value can be listed.
#'   
#'   The arguments \code{fc} and \code{lfc} give the ability to filter genes by log-fold change, but see the Note below.
#'   This argument is not available for \code{topTreat} because \code{treat} already handles fold-change thresholding in a more sophisticated way.
#'   
#'   The function \code{topTableF} is scheduled for removal in a future version of limma.
#'   It is equivalent to \code{topTable} with \code{coef=NULL}.
#' @note
#'   Although \code{topTable} enables users to set both p-value and fold-change cutoffs, the use of fold-change cutoffs is not generally recommended.
#'   If the fold changes and p-values are not highly correlated, then the use of a fold change cutoff can increase the false discovery rate above the nominal level.
#'   Users wanting to use fold change thresholding are usually recommended to use \code{treat} and \code{topTreat} instead.
#'   
#'   In general, the adjusted p-values returned by \code{adjust.method="BH"} remain valid as FDR bounds only when the genes remain sorted by p-value.
#'   Resorting the table by log-fold-change can increase the false discovery rate above the nominal level for genes at the top of resorted table.
#' @seealso
#'   An overview of linear model and testing functions is given in \link{06.LinearModels}.
#'   See also \code{\link[stats]{p.adjust}} in the \code{stats} package.
#' 
#' @author Gordon Smyth. Adapted by Lucas Beerland, Koen Van den Berge, Lieven Clement.
#' 
#' @examples
#' library(limma)
#' set.seed(495212344)
#' n <- 40 # sample size
#' P <- 10 # number of cell types
#' mu0 <- rnbinom(n=P, size=1/2, mu=400)
#' mu0 # absolute counts in group 0
#' beta <- rlnorm(n=P, meanlog = 0, sdlog=2) * # these are log-fold-changes
#'   rbinom(n=P, size=1, prob=.15) *
#'   sample(c(-1,1), size=P, replace=TRUE) # fold change on log scale
#' mu1 <- exp(beta) * mu0 # because we want log(mu2/mu1) = beta
#' relAbundances <- data.frame(g0=mu0/sum(mu0),
#'                             g1=mu1/sum(mu1)) # relative abundance of absolute count
#' # relative abundance information (observed data in typical experiment)
#' Y0 <- rmultinom(n=10, size=1e4, prob=relAbundances$g0)
#' Y1 <- rmultinom(n=10, size=1e4, prob=relAbundances$g1)
#' Y <- cbind(Y0, Y1)
#' rownames(Y) <- paste0("celltype",1:10)
#' colnames(Y) <- paste0("sample",1:20)
#' group <- factor(rep(0:1, each=10))
#' design <- model.matrix(~group)
#' v <- voomCLR(counts = Y,
#'              design = design,
#'              lib.size = NULL,
#'              plot = TRUE)
#' fit <- lmFit(v, design)
#' fit <- eBayes(fit)
#' ttNoBoot <- topTableBC(fit, coef=2, sort.by="none", n=Inf)
#' ttNonParBoot <- topTableBC(fit, coef=2, bootstrap="nonparametric", sort.by="none", n=Inf)
#' ttParBoot <- topTableBC(fit, coef=2, bootstrap="parametric", voomWeights=v$weights, sort.by="none", n=Inf)
#' @export
topTableBC <- function(fit,
                       coef=NULL,
                       number=10,
                       genelist=fit$genes,
                       adjust.method="BH",
                       sort.by="p",
                       resort.by=NULL,
                       p.value=1,
                       fc=NULL,
                       lfc=NULL,
                       confint=FALSE,
                       bootstrap=FALSE,
                       voomWeights=NULL){
    #	Summary table of top genes, object-orientated version
    #	Gordon Smyth
    #	4 August 2003.  Last modified 20 Aug 2022.

  
  ### for dev:
  # library(limma)
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
  # rownames(Y) <- paste0("celltype",1:10)
  # colnames(Y) <- paste0("sample",1:20)
  # group <- factor(rep(0:1, each=10))
  # design <- model.matrix(~group)
  # v <- voomCLR(counts = Y,
  #              design = design,
  #              lib.size = NULL,
  #              plot = TRUE)
  # fit <- lmFit(v, design)
  # fit <- eBayes(fit)
  # fit <- fit
  # coef <- 2
  # number <- 10
  # genelist=fit$genes
  # eb <- NULL
  # adjust.method <- "BH"
  # sort.by <- "M"
  # resort.by <- NULL
  # p.value <- 1
  # fc <- NULL
  # lfc <- NULL
  # confint <- FALSE
  # bootstrap <- NULL
  # voomWeights <- v$weights
  
    #	Check fit
    if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM object")
    if(is.null(fit$t) && is.null(fit$F)) stop("Need to run eBayes or treat first")
    if(is.null(fit$coefficients)) stop("coefficients not found in fit object")
    if(confint && is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
  
    # Check weights when bootstrapping
    if(isTRUE(bootstrap=="parametric")){
      if(is.null(voomWeights)){
        stop("Please provide the voom weights when using the parametric bootstrap approach.")
      }
    }
    
    #	Check coef
    if(is.null(coef)) {
      if(is.null(fit$treat.lfc)) {
        coef <- 1:ncol(fit)
        cn <- colnames(fit)
        if(!is.null(cn)) {
          i <- which(cn == "(Intercept)")
          if(length(i)) {
            coef <- coef[-i]
            message("Removing intercept from test coefficients")
          }
        }
      } else
        coef <- ncol(fit)
    }
  
  
  #	Check adjust.method
  if(is.null(adjust.method)) adjust.method <- "BH"
    
    #	Set log2-fold-change cutoff
    if(is.null(fc)) {
      if(is.null(lfc)) lfc <- 0
    } else {
      if(fc < 1) stop("fc must be greater than or equal to 1")
      lfc <- log2(fc)
    }
    
    #	If testing for multiple coefficients, call low-level topTable function for F-statistics
    if(length(coef)>1) {
      stop("Work in progress. Still need to implement bias-corrected F-tests.")
      # if(!is.null(fit$treat.lfc)) stop("Treat p-values can only be displayed for single coefficients")
      # coef <- unique(coef)
      # if(length(fit$coefficients[1,coef]) < ncol(fit)) fit <- fit[,coef]
      # if(sort.by=="B") sort.by <- "F"
      # return(.topTableF(fit,number=number,genelist=genelist,adjust.method=adjust.method,sort.by=sort.by,p.value=p.value,lfc=lfc))
    }
    
    #	Call low-level topTable function for t-statistics
    n <- nrow(fit$design)
    fit <- unclass(fit)
    ebcols <- c("t","p.value","lods", "df.total")
    if(confint) ebcols <- c("s2.post","df.total",ebcols)
    tt <- .toptableTBC(fit=fit[c("coefficients","stdev.unscaled")],
               coef=coef,
               number=number,
               genelist=genelist,
               A=fit$Amean,
               eb=fit[ebcols],
               adjust.method=adjust.method,
               sort.by=sort.by,
               resort.by=resort.by,
               p.value=p.value,
               lfc=lfc,
               confint=confint,
               n=n,
               bootstrap=bootstrap,
               design=fit$design,
               voomWeights=voomWeights,
               sigma2post=fit$s2.post)
    return(tt)
    }
