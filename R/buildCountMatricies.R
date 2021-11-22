#' Constructs allele level count matricies
#'
#' Given allele level counts in single cell data, summarised at some genomic level (SNP, exon, gene), summarises counts into matricies.  These matricies have rows as genomic regions, columns as cell barcodes, and entries as the different allele level summaries specified in \code{assays}.
#'
#' Any region that has no information in the single cell data is dropped.
#'
#' @param counts GRanges object summarising allele level expression.  Typically the output of \code{\link{getAllelicExpression}} or \code{\link{aggregateByRegions}}.
#' @param cellIDs For each row in \code{counts}, what cell does it belong to?  Must be of the same length as \code{counts}.
#' @param regionIDs For each row in \code{counts}, how should we label this region?  \code{as.character(counts)} is slow and often less useful, but should always work.
#' @param assays Which allele level expresion counts to convert to matricies.  Must exist as columns in \code{counts}.
#' @return Region by cell sparse matricies containing counts for assays provided.  Returned as a list, unless only one assay given, in which case the matrix is returned.
#' @export
#' @importFrom Matrix sparseMatrix
buildCountMatricies = function(counts,cellIDs = counts$cellID, regionIDs = counts$regionID, assays=c('matCount','patCount')){
  #Check input
  if(!all(assays %in% colnames(mcols(counts))))
    stop("Not all assays found in counts.")
  if(length(cellIDs)!=length(counts) || length(regionIDs) != length(counts))
    stop("Length of cellIDs and regionIDs does not match length of counts.")
  #Get row/column names
  cellIDsUnique = unique(cellIDs)
  regionIDsUnique = unique(regionIDs)
  outs=list()
  for(assay in assays){
    o = which(!is.na(mcols(counts)[,assay]))
    #Construct the giant matrix
    i = match(regionIDs[o],regionIDsUnique)
    j = match(cellIDs[o],cellIDsUnique)
    tmp = sparseMatrix(i=i,
                       j=j,
                       x=mcols(counts[o])[,assay],
                       dims=c(length(regionIDsUnique),length(cellIDsUnique)),
                       dimnames=list(regionIDsUnique,cellIDsUnique)
                       )
    outs[[assay]] = tmp
  }
  #Don't return list, if only one entry
  if(length(assays)==1)
    outs = outs[[1]]
  return(outs)
}
