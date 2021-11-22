#' Calculates the likelihood of each mutation in each cell
#'
#' Given a set of somatic mutations, calculates the likelihood of each cell containing a mutation under the null model (mutation not present + error rate) and the alternative model (mutation present at frequency specified).  The default is to assume that the mutation is on one of two copies (i.e., frequency 50%), but if the copy number configuration is known a more precise frequency can be given.  In general this tends not to matter as the posterior probability is based on the difference in log likelihood between the null and alt models.  If more mutant reads than expected by chance are present, that will be enough to favour the alt genotype even if the frequency is not correct.
#'
#' @inheritParams calcLikelihoods
#' @param mCnts Somatic mutations GRanges.  Usually the output of \code{\link{loadCanPipeMuts}}.
#' @param defaultErr Error rate to use when not specified by the 'errRate' column in \code{mCnts}.
#' @param mutFreq Frequency at which we expect to see a mutation in tumour cells.  Must be of length 1 or the same length as \code{mCnts}.
#' @importFrom VGAM dbetabinom
#' @export
calcMutationLikelihoods = function(mCnts,overDisp,defaultErr=0.01,mutFreq=0.5){
  errRate = getErrorRate(mCnts,defaultErr)
  #The alternate model
  tot = mCnts$altCount + mCnts$refCount
  mCnts$llAlt = -dbetabinom(mCnts$altCount,tot,mutFreq,overDisp,log=TRUE)
  mCnts$llNeut = -dbetabinom(mCnts$altCount,tot,errRate,log=TRUE)
  return(mCnts)
}

