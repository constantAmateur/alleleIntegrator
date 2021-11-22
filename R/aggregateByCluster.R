#' Collapse counts by clusters of cells
#'
#' Given a number of assays (e.g. counts) at a series of regions (e.g. SNPs), summarises assays by groupings of cells (e.g. clusters).
#'
#' Clusters can be given as either a named list, or a named vector.  If it is a list, each list entry should contain cellIDs to be grouped.  When it is a vector, it should be a vector of cellIDs with cluster IDs as the names.
#'
#' When a cell occurs in multiple multiple groupings there are multiple approaches to aggregation.  Current options (passed to \code{overlapBehaviour}) are 'drop' which ignores any cell that maps to more than one grouping, or 'all', which allocates the cell counts to all groupings that contain it.
#'
#' @inheritParams aggregateByLists
#' @param clusters List or named vector defining how cells are to be grouped.  If length matches \code{gCnts} assumed to directly indicate which clusters each row belongs to.  If a vector, names should be cellIDs and values cluster names.  If a list, list names should be cluster names.
#' @param preserve Columns in \code{gCnts} that should be copied across unchanged to the aggregated object.  Anything not present in gCnts in ignored.  E.g. geneID
#' @param overlapBehaviour When the same SNP overlaps multiple groupings, what should we do.  Options are 'all' or 'drop'.  See details.
#' @return A GRanges object based on \code{gCnts}, with aggregated assays
#' @seealso aggregateByRegions aggregateByLists
#' @export
aggregateByClusters = function(gCnts,clusters=gCnts$clusterID,assays=c('A','C','G','T','Tot','altCount','refCount','matCount','patCount'),preserve=c('REF','ALT','geneID','imprintingFDR','imprinted','errRate','matASE','matASE_postAlpha','matASE_postBeta','regionType','altIsMum','geneID','geneName'),overlapBehaviour=c('all','drop')){
  #Validate params
  overlapBehaviour = match.arg(overlapBehaviour)
  if(is.null(clusters))
    stop("Invalid cluster specification")
  if(!all(assays %in% colnames(mcols(gCnts))))
    stop("Not all assays found in gCnts object.")
  if(is.null(names(gCnts)))
    names(gCnts) = as.character(gCnts)
  #Format the clusters
  if(length(clusters)==length(gCnts)){
    groupBy = clusters
  }else{
    if(!is.list(clusters)){
      if(is.null(names(clusters)))
        stop("If clusters is not a list, it must be a vector of cellIDs with cluster IDs as names.")
      clusters = split(names(clusters),clusters)
    }
    clusterNames = names(clusters)
    #Convert to gCnts format list
    clusters = split(rep(clusterNames,lengths(clusters)),unlist(clusters))
    groupBy = clusters[gCnts$cellID]
  }
  #Drop duplicates if that's what we're doing
  if(overlapBehaviour=='drop'){
    groupBy[lengths(groupBy)>1]=NA
  }
  #Do the collapsing
  out = aggregateByLists(gCnts,assays,groupBy,names(gCnts))
  #Reformat as GRanges again
  base = match(out$regionID,names(gCnts))
  base = gCnts[base]
  #Strip back to essentials
  mcols(base) = NULL
  mcols(base) = cbind(mcols(base),out)
  #Add in any of the extras
  m = match(base,gCnts)
  preserve = intersect(preserve,colnames(mcols(gCnts)))
  if(length(preserve)>0)
    mcols(base) = cbind(mcols(base),mcols(gCnts)[m,preserve])
  return(base)
}
