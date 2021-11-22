#' Calculates the allele specific expression of genes specified
#'
#' Determines what the default allelic balance of expression is for each gene, under the assumption of normality.  That is, given a set of cells that are very likely to be normal, the expression for each allele is aggregated across all cells and a posterior probability for the allelic ratio is calculated. Counts are modelled using a binomial distribution with beta prior, where the beta prior is pramaterised by it's mean (set to 0.5) and kappa = alpha + beta.  The shape of this posterior can be determined from the data, using \code{\link{calcPriorASE}}, although the default value of \code{priorKappa} should do a reasonable job.
#'
#' For SNPs that don't overlap any of the genes provided in \code{gns}, a SNP level allele specific expression value is calculated using the same method.
#'
#' Where a SNP is overlapped by more than one gene, the SNP is allocated the most extreme ASE value (i.e., biggest deviation from 0.5).  This behaviour can be changed by setting \code{overlapMax}.
#'
#' If \code{gns} is NULL, gene information will be loaded from the annotation stored in \code{gCnts} from \code{\link{annotateSNPs}}.
#'
#' @inheritParams findImprinting
#' @param normalCells Cells that you are very confident are normal.  If not specified, stored values used if present or all genes are set to have an expected allele fraction of 0.5
#' @param priorKappa Value of kappa used for the beta prior.  If NULL, calculated using \code{\link{calcPriorASE}}.
#' @param overlapMax When two genes overlap a SNP, take the most extreme value.  If FALSE, takes the least extreme value.
#' @param ... Passed to \code{\link{calcPriorASE}}.
#' @return The input \code{gCnts} object, with an extra column giving an estimate of the maternal allele expression fraction, using posterior likelihood.  The alpha and beta values of the posterior ASE distribution are also added.
#' @importFrom Matrix rowSums
#' @export
calcASE = function(gCnts,normalCells,gns=NULL,priorKappa=NULL,overlapMax=TRUE,...){
  #Turn off ASE
  if(missing(normalCells)){
    if(!is.null(gCnts$isNorm)){
      normalCells = unique(gCnts$cellID[which(gCnts$isNorm)])
    }else{
      message("No normal cells specified, all genes assumed to have no ASE.")
      normalCells = unique(gCnts$cellID)
      priorKappa = Inf
    }
  }
  #Validate genes
  if(is.null(gns))
    gns = gnsFromAnnotation(gCnts)
  #Get prior
  if(is.null(priorKappa))
    priorKappa = calcPriorASE(gCnts,normalCells,gns,...)
  message(sprintf("Using beta prior with mean 0.5 and kappa %g",priorKappa))
  #Build the matrix by gene
  gnsPlus = c(gns,subsetByOverlaps(unique(gCnts),gns,invert=TRUE))
  mcols(gnsPlus) = NULL
  gnCnts = aggregateByRegions(gCnts[gCnts$cellID %in% normalCells],gnsPlus,assays=c('matCount','patCount'))
  gnCnts = buildCountMatricies(gnCnts)
  #Get total counts
  mCounts = rowSums(gnCnts$matCount)
  pCounts = rowSums(gnCnts$patCount)
  #Get the posterior alpha/beta
  dd = data.frame(row.names=names(mCounts),
                  geneName = names(mCounts),
                  postAlpha = 0.5*priorKappa + mCounts,
                  postBeta = 0.5*priorKappa + pCounts)
  #Calculate the maximum posterior probability value
  dd$matFrac = (dd$postAlpha-1)/(dd$postAlpha+dd$postBeta-2)
  #Do a check for bias
  tmp = sum(dd$matFrac>0.5)/sum(dd$matFrac!=0.5)
  if(abs(tmp-0.5)>0.2)
    warning(sprintf("%d of %d have ASE greater than 0.5.  Your normal cells are probably not normal.  Consider running without normalCells.",sum(dd$matFrac>0.5),nrow(dd)))
  #Set them back into the gCnts object
  o = findOverlaps(gCnts,gnsPlus)
  dd = dd[names(gnsPlus)[subjectHits(o)],]
  dd$gCntsIdx = queryHits(o)
  dd = dd[order(dd$gCntsIdx,abs(dd$matFrac-0.5)*ifelse(overlapMax,-1,1)),]
  #Keep one per gCnts entry
  dd = dd[!duplicated(dd$gCntsIdx),]
  #Save the results
  gCnts$matASE = dd$matFrac
  gCnts$matASE_postAlpha = dd$postAlpha
  gCnts$matASE_postBeta = dd$postBeta
  #Set the ones for which there was no information to the prior
  w = is.na(gCnts$matASE)
  gCnts$matASE[w] = 0.5
  gCnts$matASE_postAlpha[w] = priorKappa*0.5
  gCnts$matASE_postBeta[w] = priorKappa*0.5
  return(gCnts)
}
