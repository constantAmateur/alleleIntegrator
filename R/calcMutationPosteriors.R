#' Calculates posterior probability of each mutation
#'
#' Given a set of mutations, calculates the posterior probability of each cell containing ALL mutations.  Likelihoods on which this is based are calculated using \code{\link{calcMutationLikelihoods}}, which will be run by this function if it has not been run previously.
#'
#' @inheritParams calcMutationLikelihoods
#' @return A table giving the posterior probability of the mutant genotype for each cell.
#' @export
calcMutationPosteriors = function(mCnts,overDisp,defaultErr=0.01,mutFreq=0.5){
  #Do we need to calculate likelihoods?
  if(!all(c('llAlt','llNeut') %in% colnames(mcols(mCnts))))
    mCnts = calcMutationLikelihoods(mCnts,overDisp,defaultErr,mutFreq)
  #Now group across all cells
  ans = aggregate(cbind(llAlt,llNeut,altCount,refCount) ~ cellID,data=mcols(mCnts),FUN=sum)
  ans$llDiff = ans$llAlt-ans$llNeut
  ans$postTum = (1+exp(ans$llDiff))**-1
  ans$postNorm = 1-ans$postTum
  return(ans)
}
