#' Filter cells
#'
#' Counting of alleles using \code{\link{getAllelicExpression}} treats all barcodes as equal.  Usually we are not interested in most barcodes and only those that pass QC.  This function filters out useless barcodes, incorporates cluster information if available, and applies a few other sensible filters to produce a set of high quality counts.
#'
#' If \code{clusterIDs} is provided, as either a list or named vector mapping cluster IDs to cell IDs, the permitted cell IDs will be derived from this object if not explicitly specified by \code{passCellIDs}.  The specification format is pretty flexible, so any kind of named vector or list that maps cellIDs to clusterIDs will probably work.
#'
#' If \code{normIDs} are provided this information will be stored in the returned object and also used to calculated imprinted genes if needed and \code{dropImprinted=TRUE}.
#'
#' @param phCnts The counts to be filtered.  Usually produced by \code{\link{getAllelicExpression}}
#' @param clusterIDs Map of cluster IDs to cell IDs.  Either a vector or list, named by cluster ID, with entries being cell IDs.
#' @param passCellIDs The IDs of the cells to keep.  Derived from \code{clusterIDs} if given and this is set to \code{NULL}.
#' @param normIDs The IDs of the cells or clusters to mark as likely normal. If NULL, nothing stored.
#' @param regionsToKeep What types of regions should be kept?
#' @param genesToDrop Drop SNPs in these genes.  The defaults are genes that almost always have some hard to predict allele specific expression (either imprinting, or a bias) and are best ignored.
#' @param dropImprinted Drop imprinted regions.  Attempts to identify them using \code{\link{findImprinting}} if not already calculated.
#' @param dropInsane Drop SNPs that have failed sanity checks.
#' @param dropUninformative Drop SNPs marked as uninformative (i.e., couldn't be phased).
#' @param dropZeroCoverage Drop entries with no coverage of maternal or paternal allele.
#' @param segs CN segments to aggregate information for.  Only used if \code{verbose=2}.
#' @param verbose Report summary stats.
#' @return A GRanges object with filtering done.
#' @importFrom stats quantile
#' @export
filterCells = function(phCnts,clusterIDs=NULL,passCellIDs=NULL,normIDs=NULL,regionsToKeep=c('Intronic','Exonic'),genesToDrop=biasedGenes,dropImprinted=TRUE,dropInsane=TRUE,dropUninformative=TRUE,dropZeroCoverage=TRUE,segs=phCnts@metadata$segs,verbose=TRUE){
  #Check and fixup filtering options as needed
  if(is.null(phCnts$geneName) && is.null(phCnts$geneID) && length(genesToDrop)>0)
    stop("Gene name and ID not available in phCnts, cannot filter genes.  Gene names/IDs must be stored in columns named 'geneName' or geneID'") 
  #Toggle off options that we can't apply
  if(dropInsane & is.null(phCnts$passSanity))
    dropInsane=FALSE
  if(dropUninformative & is.null(phCnts$informative))
    dropUninformative=FALSE
  if(!is.null(clusterIDs)){
    #If it's not a list, convert it to a list
    if(!is.list(clusterIDs))
      clusterIDs = split(as.character(clusterIDs),names(clusterIDs))
    #Now make sure it's the right way around.  That is, the list should be named by clusters
    nomMatch = sum(rep(names(clusterIDs),lengths(clusterIDs)) %in% phCnts$cellID)
    entMatch = sum(unlist(clusterIDs) %in% phCnts$cellID)
    #Standardise on list named by clusterID where entries are cellIDs 
    if(nomMatch > entMatch)
      clusterIDs = split(rep(names(clusterIDs),lengths(clusterIDs)),unlist(clusterIDs))
    #Check that we don't have conflicting annotations
    if(any(duplicated(unlist(clusterIDs))))
      stop("The same cellID is allocated to mulitple clusters.  Cluster allocations must be unique.")
    #Construct passCellIDs if not explicitly provided
    if(is.null(passCellIDs))
      passCellIDs = unique(unlist(clusterIDs))
  }
  #Ensure that some of our IDs actually match
  if(!is.null(passCellIDs)){
    if(!any(passCellIDs %in% phCnts$cellID))
      stop("None of passCellIDs found in count object phCnts.")
  }
  #Incorporate cluster information
  phCnts$clusterID = rep(names(clusterIDs),lengths(clusterIDs))[match(phCnts$cellID,unlist(clusterIDs))]
  #Store the normal cell designation
  if(!is.null(normIDs)){
    if(all(normIDs %in% phCnts$clusterID)){
      phCnts$isNorm = phCnts$clusterID %in% normIDs
    }else{
      phCnts$isNorm = phCnts$cellID %in% normIDs
    }
  }
  #Check if imprinting needs calculating
  if(dropImprinted && is.null(phCnts$imprinted)){
    message("Calculating imprinted genes.")
    phCnts = findImprinting(phCnts)
  }
  #Do the filtering
  #Regions
  w = phCnts$regionType %in% regionsToKeep
  if(verbose)
    message(sprintf('Keeping %s (%.01f%%) entries after region type filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  #Genes
  if(!is.null(phCnts$geneID)){
    badIdxs = rep(seq_along(phCnts),lengths(phCnts$geneID))[(unlist(phCnts$geneID) %in% genesToDrop)]
    w = w & !(seq_along(phCnts) %in% badIdxs)
  }
  if(!is.null(phCnts$geneName)){
    badIdxs = rep(seq_along(phCnts),lengths(phCnts$geneName))[(unlist(phCnts$geneName) %in% genesToDrop)]
    w = w & !(seq_along(phCnts) %in% badIdxs)
  }
  if(verbose)
    message(sprintf('Keeping %s (%.01f%%) entries after gene filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  #Imprinted
  if(dropImprinted){
    w = w & !phCnts$imprinted
    if(verbose)
      message(sprintf('Keeping %s (%.01f%%) entries after imprinting filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  }
  #Cells
  if(!is.null(passCellIDs)){
    w = w & (phCnts$cellID %in% passCellIDs)
    if(verbose)
      message(sprintf('Keeping %s (%.01f%%) entries after cell filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  }
  #Basic filters
  if(dropInsane){
    w = w & phCnts$passSanity
    if(verbose)
      message(sprintf('Keeping %s (%.01f%%) entries after sanity filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  }
  if(dropUninformative){
    w = w & phCnts$informative
    if(verbose)
      message(sprintf('Keeping %s (%.01f%%) entries after informative filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  }
  if(dropZeroCoverage){
    w = w & (phCnts$matCount + phCnts$patCount)>0
    if(verbose)
      message(sprintf('Keeping %s (%.01f%%) entries after coverage filter',prettyNum(sum(w),big.mark=','),sum(w)/length(w)*100))
  }
  gCnts = phCnts[w,]
  #Summary stats
  #To work out how many informative
  gCnts$hasCov=ifelse(gCnts$matCount+gCnts$patCount>0,1,0)
  globCnts = aggregateByLists(gCnts,assays=c('matCount','patCount','hasCov'),regionList=1)
  globCnts$totCounts = globCnts$matCount + globCnts$patCount
  gCnts$hasCov=NULL
  quants=c(0,.25,.5,.75,1)
  if(verbose){
    message(sprintf('After filtering there are %s cells with:',prettyNum(nrow(globCnts),big.mark=',')))
    x = quantile(globCnts$hasCov,quants)
    message(sprintf("  Quantiles of SNPs covered: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
    x = quantile(globCnts$totCounts,quants)
    message(sprintf("  Quantiles of total counts: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
    x = quantile(globCnts$matCount,quants)
    message(sprintf("  Quantiles of maternal counts: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
    x = quantile(globCnts$patCount,quants)
    message(sprintf("  Quantiles of paternal counts: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
  }
  if(verbose>1 & !is.null(segs)){
    message('##########################\n# Segment specific stats #\n##########################')
    cCnts = aggregateByRegions(gCnts,segs,assays=c('matCount','patCount','hasCov'))
    cCnts$totCounts = cCnts$matCount + cCnts$patCount
    for(i in seq_along(segs)){
      lab = names(segs)[i]
      seg = segs[i]
      tmp = subsetByOverlaps(cCnts,seg)
      message(sprintf('For segment %s the stats are:',lab))
      x = quantile(tmp$hasCov,quants)
      message(sprintf("  Quantiles of SNPs covered: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
      x = quantile(tmp$totCounts,quants)
      message(sprintf("  Quantiles of total counts: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
      x = quantile(tmp$matCount,quants)
      message(sprintf("  Quantiles of maternal counts: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
      x = quantile(tmp$patCount,quants)
      message(sprintf("  Quantiles of paternal counts: %s",paste(prettyNum(x,big.mark=','),' (',names(x),')',sep='',collapse = ' ')))
    }
  }
  return(gCnts)
}
