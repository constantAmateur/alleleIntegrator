#' Calculates the allele specific expression prior
#'
#' The allele specific expression of each gene (or SNP) is estimated using a beta prior on a binomial likelihood.  As we don't expect any particular preference in allelic ratio, the mean of this beta prior should be 0.5.  This function sets the dispersion of this prior by fitting a beta distribution with mean 0.5 to the allelic ratios of high expression genes across cells highly likely to be normal.  That is, it looks at genes with extremely high expression where we can very confidently determine the allelic ratio and uses this to determine how much allele specific expression to expect for those genes where we don't have as much evidence. 
#'
#' We use the parameterisation of the beta distribution by mean mu and kappa = alpha + beta.  It is kappa that is optimised for and returned by this function.  If \code{\link{findImprinting}} has been run, genes marked as imprinted are not used. 
#'
#' @note Failure to properly specify a set of sensible normal cells (via \code{normalCells}) will lead to over-estimation of prior and reduced power to call copy number changes.
#'
#' @inheritParams findImprinting
#' @param imprintCut Anything within this value of 0 or 1 are ignored.  Beta likelihood is not defined at 0 or 1.
#' @param expCut Consider genes with more than this many reads across all normal cells to be "highly expressed".
#' @return The optimum value of kappa for the beta distribution prior.
#' @importFrom stats dbeta optimise
#' @importFrom Matrix rowSums
calcPriorASE = function(gCnts,normalCells,gns=NULL,imprintCut=0.01,expCut=400){
  #Validate genes
  if(is.null(gns))
    gns = gnsFromAnnotation(gCnts)
  #Get imprinted genes if it hasn't been run already
  if(!'imprinted' %in% colnames(mcols(gCnts))){
    gCnts = findImprinting(gCnts,normalCells,gns=gns)
  }
  imprinted = gCnts$imprinted
  #Add in extra regions for the SNPs that aren't in gns.  They're probably not expressed at a high enough level, but none the less
  extras = subsetByOverlaps(unique(gCnts),gns,invert=TRUE)
  gns = c(gns,extras)
  mcols(gns) = NULL
  gnCnts = aggregateByRegions(gCnts[!imprinted & gCnts$cellID %in% normalCells],gns,assays=c('matCount','patCount'))
  gnCnts = buildCountMatricies(gnCnts)
  #Get total counts
  mCounts = rowSums(gnCnts$matCount)
  pCounts = rowSums(gnCnts$patCount)
  mFrac = mCounts/(mCounts+pCounts)
  #Get the allelic ratios for our high coverage ones
  w = which(mCounts+pCounts>expCut)
  if(length(w)<10)
    stop(sprintf("Only %d genes have more than %d counts.  Consider reducing expCut or manually specifying kappa.",length(w),expCut))
  x = mFrac[w]
  #Ones that look genuinely imprinted are dropped
  x = x[x>imprintCut & x<1-imprintCut]
  #Do the fitting
  optFun = function(kappa,x) -sum(dbeta(x,0.5*kappa,0.5*kappa,log=TRUE))
  #0 shunts all the probability to 0 or 1.  Anything above 1,000 is pretty much a delta function at 0.5.
  fit = optimise(optFun,c(0,1e4),x=x) 
  #Anything below 2 (flat PDF) is probably a failure.
  kappa = fit$minimum
  if(kappa < 2)
    stop(sprintf("Calculated Kappa = %g.  Kappa below 2 indicates allele specific expression is more likely than balance.  Something has probably gone wrong.",kappa))
  #Return the value of kappa
  return(kappa)
}
