#' Determine error rate for groups of counts
#' 
#' Assuming the genotype is correct, then only REF/ALT alleles should be observed and anything other than that is an error.  The number that are neither, sets a limit of the error rate.  This only measures sequencing errors, not getting the wrong allele due to ambient contamination, or mis-mapping.  There is also seldom sufficient power to detect error rates on a SNP by SNP basis.  So instead, SNPs are binned in some way (typically into intronic, exonic, and inter-genic) and the error rate for each group is calculated.
#'
#' @inheritParams findImprinting
#' @param groupBy How to group together things in \code{gCnts}.
#' @return A vector giving the error rate by each group.
#' @export
calcErrorRate = function(gCnts,groupBy=gCnts$regionType){
  if(is.null(groupBy) || length(groupBy)!=length(gCnts))
    stop("groupBy must be the same length as gCnts")
  errRates = aggregate(cbind(Tot,errors=Tot-refCount-altCount) ~ groupBy,FUN=sum,data=cbind(mcols(gCnts),groupBy))
  errRates$errRate = errRates$errors/errRates$Tot
  #This is the number that go from REF/ALT to one of the other bases.  Some errors will just jump between REF and ALT.  To get the true error rate need to account for this.
  errRates$errRate = errRates$errRate / (2/3)
  gCnts$errRate = errRates$errRate[match(groupBy,errRates$groupBy)]
  return(gCnts)
}
