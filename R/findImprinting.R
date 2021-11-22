#' Find imprinted regions
#'
#' Given allele specific expression at het SNPs and gene definitions, identifies genes or individual SNPs that have only one allele expressed across all normal cells.  These regions could either be due to imprinting or mis-called heterozygous SNPs (see note).
#'
#' Imprinted genes are determined using a binomial test, against the null of strong allele specific expression (but not imprinting).  The strength of the null hypothesis is determined by \code{nullDiff}.  This sets how far from 0.5 the expected allele frequency should be for the null distribution.  For a gene that appears to be imprinted with no maternal expression, the null tests that the maternal fraction  is less than \code{0.5 - nullDif}.  For a gene with no paternal expression, it would test \code{0.5+nullDiff}.
#'
#' If genes are specified by \code{gns} then an extra column, geneImprinted, will be added to the output indicating if there is evidence of imprinting at the gene level.
#'
#' @note A gene appearing in this list is either evidence of imprinting **or** evidence that the heterozygous SNP in this region is not actually heterozygous.  For this reason, genes with higher values of \code{nSupportingSNPs} are more reliable.  If all you want to do is call CN segments, it doesn't matter, just remove everything that passes.  If you are interested in the biology, look closely at the SNPs underlying the call and treat those with only 1 or 2 SNPs with suspicion.  Even when many SNPs support a gene, make sure that the call is not being driven by high coverage of one SNP with a trivial contribution from the others.
#'
#' @note This differs from the results of \code{\link{calcASE}} in that the allele specific expression prior is calculated excluding imprinted genes.  The ideal workflow is to first run \code{findImprinting} and pass the resulting imprinted genes to \code{\link{calcASE}}.
#' 
#' @param gCnts Allele specific counts at individual loci in single cell data.  Produced by \code{\link{getAllelicExpression}}.
#' @param normalCells IDs of cells that can be considered normal.  If missing, use stored labelling of normal cells or consider all cells normal.
#' @param gns A GRanges object giving the coordinates of genes to use.  If NULL, load from annotation if \code{\link{annotateSNPs}} has been run.
#' @param useGenes Should we use genes?  If FALSE will consider each SNP independently and ignore gene information.
#' @param nullDiff The deviation from 0.5 to use for the maternal allele expression fraction under the null hypothesis.  E.g. 0.3 implies a null of 0.2 or 0.8.  See details.
#' @param FDR False discover rate to accept.
#' @return A modified version of \code{gCnts} with columns indicating imprinting.  A table with details underlying the calls is also added to \code{gCnts@metadata$imprinted}.
#' @importFrom Matrix rowSums
#' @export
findImprinting = function(gCnts,normalCells,gns=NULL,useGenes=TRUE,nullDiff=0.3,FDR=0.05){
  if(missing(normalCells)){
    #Check if already exist
    if(!is.null(gCnts$isNorm)){
      normalCells = unique(gCnts$cellID[which(gCnts$isNorm)])
    }else{
      message("No normal cells specified, assuming all cells are normal")
      normalCells = unique(gCnts$cellID)
    }
  }
  if(useGenes & is.null(gns))
    gns = gnsFromAnnotation(gCnts)
  #Get the evidence at the SNP level
  mat = buildCountMatricies(gCnts[gCnts$cellID %in% normalCells],assays=c('matCount','patCount'))
  mCnts = rowSums(mat$matCount)
  pCnts = rowSums(mat$patCount)
  #If we have genes, do that as well
  if(useGenes){
    gnCnts = aggregateByRegions(gCnts[gCnts$cellID %in% normalCells],gns,assays=c('matCount','patCount'))
    gnCnts = buildCountMatricies(gnCnts)
    mCnts = c(mCnts,rowSums(gnCnts$matCount))
    pCnts = c(pCnts,rowSums(gnCnts$patCount))
  }
  tot = mCnts+pCnts
  #Construct the test
  w = mCnts>pCnts
  pVals = rep(NA,length(tot))
  pVals[w] = pbinom(pCnts[w],tot[w],0.5-nullDiff)
  pVals[!w] = pbinom(mCnts[!w],tot[!w],0.5-nullDiff)
  qVals = p.adjust(pVals,method='BH')
  out = data.frame(row.names = names(tot),
                   regionID = names(tot),
                   matCnts = mCnts,
                   patCnts = pCnts,
                   totCnts = tot,
                   matFrac = mCnts/tot,
                   pVal = pVals,
                   FDR = qVals,
                   imprinted = qVals<FDR)
  #Get the number of supporting SNPs.
  out$nSupportingSNPs = 1
  if(!is.null(gns)){
    snpCnts = unique(gCnts[gCnts$cellID %in% normalCells])
    o = findOverlaps(snpCnts,gns)
    snpCnts = table(names(gns)[subjectHits(o)])
    out$nSupportingSNPs = snpCnts[out$regionID]
  }
  out = out[order(out$FDR),]
  gCnts@metadata$imprinted=out
  #Add the highlights to the GRanges thing
  if(is.null(gns)){
    gnsExtra = unique(gCnts)
  }else{
    gnsExtra = c(gns,unique(gCnts))
  }
  m = match(names(gnsExtra),out$regionID)
  gnsExtra = gnsExtra[!is.na(m)]
  m = m[!is.na(m)]
  out = out[m[!is.na(m)],]
  o = findOverlaps(gnsExtra,gCnts)
  out = out[queryHits(o),]
  out$gCntsIdx = subjectHits(o)
  #Keep the lowest FDR for each entry (gene or SNP level)
  out = out[order(out$gCntsIdx,out$pVal),]
  out = out[!duplicated(out$gCntsIdx),]
  gCnts$imprintingFDR = NA
  gCnts$imprinted = FALSE
  gCnts$imprinted[out$gCntsIdx] = out$imprinted
  gCnts$imprintingFDR[out$gCntsIdx] = out$FDR
  return(gCnts)
}
