#' Convert positions between reference genomes
#'
#' Given a set of locations (SNPs, SNVs, etc) generated using one reference genome, uses the chain file located at \code{chainPath} to convert them to another reference system.
#'
#' @param loci Locations to be lifted over.  Should be a GRanges object where each entry has width 1.
#' @param chainPath The path to the lift over chain to convert between coordinates.
#' @return \code{loci} converted to new genome annotation.
#' @export
#' @importFrom rtracklayer import.chain liftOver
changeGenomeVersion = function(loci,chainPath){
  ch = rtracklayer::import.chain(chainPath)
  out = unlist(rtracklayer::liftOver(loci,ch))
  #Ensure the meta-data is lifted too.
  out@metadata = loci@metadata
  return(out)
}

