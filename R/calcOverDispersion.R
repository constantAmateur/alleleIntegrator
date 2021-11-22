#' Calculates burst kinetic over-dispersion
#'
#' The naive expectation is that without any copy number changes, expression should be equally distributed between alleles, with counts following a binomial distribution.  This basic expectation is incorrect for three reasons: errors, allele specific expression, and the burst like nature of transcription.  The model accounts for transcriptional bursts, where the tendency of each allele to have a burst of transcription not perfect correlated with the other allele, by using a beta-binomial model.  The extra dispersion introduced by this beta-binomial model recognises that even in CN neutral cells with no errors or allele specific bias in expression, we will still observe deviations from 50/50 reasonably often.
#'
#' This function calibrates exactly how much extra dispersion to expect, by considering cells which are very likely to be normal and fitting a beta-binomial model.
#'
#' @note The expectation is that \code{\link{calcASE}} has been run on \code{gCnts} before this function is called.  If it has not, \code{gCnts$matASE} must be defined (e.g. set to 0.5 for all) and \code{marginaliseASE} set to \code{FALSE}.
#'
#' @inheritParams findImprinting 
#' @param marginaliseASE Should we marginalise over the posterior distirbution for each gene's allele specific expression?
#' @param ... Passed to \code{\link{likelihoodModel}}.
#' @return The best fit for the over-dispersion.
#' @export
#' @importFrom VGAM dbetabinom
#' @importFrom stats optimise
calcOverDispersion = function(gCnts,normalCells,marginaliseASE=FALSE,...){
  if(missing(normalCells)){
    if(!is.null(gCnts$isNorm)){
      normalCells = unique(gCnts$cellID[which(gCnts$isNorm)])
    }else{
      message("No normal cells specified, assuming all cells are normal")
      normalCells = unique(gCnts$cellID)
    }
  }
  #Keep only those that are normal
  gCnts = gCnts[gCnts$cellID %in% normalCells]
  #Define the model
  message("Calculating optimal over-dispersion (be patient):") 
  optFun = function(overDisp) {
    cat('.')
    sum(likelihoodModel(gCnts$matCount,gCnts$patCount,0.5,overDisp,gCnts$errRate,gCnts$matASE,gCnts$matASE_postAlpha,gCnts$matASE_postBeta,marginaliseASE=marginaliseASE,verbose=FALSE))
  }
  fit = optimise(optFun,c(0,1))
  cat('\n')
  #Return the best value
  return(fit$minimum)
}
