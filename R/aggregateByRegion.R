#' Collapse counts by regions
#'
#' Given a number of assays (e.g. counts) at a series of regions (e.g. SNPs), summarises assays by groupings of regions (e.g. genes).
#'
#' Regions can be given as either a GRanges or GRangesList object.  If it is a GRangesList object, as created by \code{\link[GenomicFeatures]{exonsBy}} for example, the output will coerce this into a GRanges object using \code{range} (unless \code{returnMatricies} is TRUE).
#'
#' When a SNP overlaps with multiple groupings there are multiple approaches to aggregation.  Current options (passed to \code{overlapBehaviour}) are 'drop' which ignores any SNP that maps to more than one grouping, or 'all', which allocates the SNPs counts to all groupings that overlap it.
#'
#' @inheritParams aggregateByLists
#' @param regions GRanges or GRangesList object defining the regions in which to aggregate counts.
#' @param overlapBehaviour When the same SNP overlaps multiple groupings, what should we do.  Options are 'all' or 'drop'.  See details.
#' @return A GRanges object based on \code{regions}, with aggregated counts for the supplied assays as columns.
#' @seealso aggregateByClusters aggregateByLists
#' @export
#' @importFrom S4Vectors queryHits subjectHits
aggregateByRegions = function(gCnts,regions,assays,overlapBehaviour=c('all','drop')){
  #Validate params
  overlapBehaviour = match.arg(overlapBehaviour)
  if(!all(assays %in% colnames(mcols(gCnts))))
    stop("Not all assays found in gCnts object.")
  if(is.null(names(regions)))
    names(regions) = as.character(regions)
  regionNames = names(regions)
  if(any(duplicated(regionNames)))
    stop("Names of regions must be unique")
  if(is.null(names(gCnts)))
    names(gCnts) = as.character(gCnts)
  snpIDs = names(gCnts)
  cellIDs = gCnts$cellID
  o = findOverlaps(gCnts,regions,ignore.strand=TRUE)
  #Should we drop anything
  regSplit = names(regions)[subjectHits(o)]
  if(overlapBehaviour=='drop'){
    multi = queryHits(o)[duplicated(queryHits(o))]
    regSplit[queryHits(o) %in% multi]=NA
  }
  #Make into a list
  tmp = vector(mode='list',length(gCnts))
  regSplit = split(regSplit,queryHits(o))
  tmp[as.numeric(names(regSplit))]=regSplit
  out = aggregateByLists(gCnts,assays,gCnts$cellID,tmp)
  #Build the output base
  if(inherits(regions,'GRangesList')){
    base = range(regions,ignore.strand=TRUE)
    if(any(lengths(base)>1))
      stop("Cannot compress input regions to single range.  Do ranges span chromosomes?")
    base = unlist(base)
  }else{
    base = regions
  }
  ans = base[match(out$regionID,names(regions)),]
  mcols(ans) = cbind(mcols(ans),out)
  #Retain clusterID if it's there
  if(!is.null(gCnts$clusterID))
    ans$clusterID = gCnts$clusterID[match(ans$cellID,gCnts$cellID)]
  return(ans)
}
