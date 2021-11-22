#' Calculates probability of prespecified CN states being present
#'
#' Given a set of pre-specified CN states, such as copy number aberations in a tumour, calculates the posterior probability of each state, in each cell, for each segment.  This is a convenience function that calls \code{\link{calcStateProbs}} with parameters suited to identifying copy number changes in tumours.
#'
#' The main thing to be specified is the GRanges object defining the copy number segments.  By default, the saved version specified by the user when running \code{\link{phaseSNPsFromCN}} is attempted to be used.  The model needs to know the expected major allele frequency for abberant and normal copy number states for each segment, which are assumed to be stored in the GRanges object with column names 'tumFrac' and 'normFrac'.  But they can also be specified by either giving a different column name, or explicitly specifying a value.  When values specified, they must be either a single value (if the same for all segments) or a vector of the same length as \code{segs}.
#'
#' It is often the case that you are not interested in the individual segment probabilities for each cell, and just care about the joint probability of each cell having all the abbarent CN states.  The obvious example for this is the identification of tumour cells.  In this case, you can set \code{globalOnly=TRUE} and this function will return a table of global abbarent state probabilities for each cell.
#'
#' @inheritParams calcStateProbs
#' @param segs GRanges object specifying segments with different copy number states.
#' @param abbFrac The major allele frequency expected for the abberant state in each segment.  See details.
#' @param normFrac The major allele frequency expected for the normal state in each segment.  See details.
#' @param globalOnly Don't return segment specific calls, just the global model.
#' @param ... Passed to \code{\link{calcStateProbs}}
#' @return A GRanges object with each entry being a cell/CN segment combination with columns reporting the probability of the abbarent state (\code{probAbberant}).  If \code{globalOnly=TRUE} a data.frame is returned instead with probabilities for the abbarent state for each cell at a combined genome wide level.
#' @export
abbSegProb = function(gCnts,overDisp,segs=gCnts@metadata$segs,abbFrac=gCnts@metadata$segs$tumFrac,normFrac = 0.5,calcJointProb=TRUE,globalOnly=FALSE,...){
  #Can't have global model without global probabilities
  if(globalOnly)
    calcJointProb=TRUE
  #Check segs is the right sort of thing
  if(!inherits(segs,'GRanges'))
    stop("segs is not a valid GRanges object.")
  #Validate the abbFrac
  if(is.character(abbFrac)){
    if(!abbFrac %in% colnames(mcols(segs)))
      stop(sprintf("abbFrac points to column %s of segs which does not exist",abbFrac))
    abbFrac = mcols(segs)[,abbFrac]
  }
  if(!length(abbFrac) %in% c(1,length(segs)))
    stop(sprintf("abbFrac is of length %d.  It should be either length 1 or the length of segs (%d)",length(abbFrac),length(segs)))
  segs$abbFrac = abbFrac
  #Validate the normFrac
  if(is.character(normFrac)){
    if(!normFrac %in% colnames(mcols(segs)))
      stop(sprintf("normFrac points to column %s of segs which does not exist",normFrac))
    normFrac = mcols(segs)[,normFrac]
  }
  if(!length(normFrac) %in% c(1,length(segs)))
    stop(sprintf("normFrac is of length %d.  It should be either length 1 or the length of segs (%d)",length(normFrac),length(segs)))
  segs$normFrac = normFrac
  #Now run the thing
  out = calcStateProbs(gCnts,overDisp,regions=segs,CNStates=c('abbFrac','normFrac'),calcJointProb=calcJointProb,...)
  #Add a couple of columns that are most likely of interest
  out$probAbberant = out$postProb_abbFrac
  out$removeN_probAbberant = out$removeN_postProb_abbFrac
  if(globalOnly){
    out = data.frame(mcols(out[seqnames(out)=='genomeWide',]))
    out$regionID=NULL
    rownames(out)=out$cellID
  }
  return(out)
}
