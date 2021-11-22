#' Get filtered cellIDs from 10X data
#'
#' Given BAM files and labels, as passed to \code{\link{getAllelicExpression}}, assumes the BAM files sit in a cellranger output folder and loads those barcodes that are in the filtered list.
#'
#' @param bams The 10X BAM files.
#' @param labels The label for each 10X bam file.
#' @return The barcodes that pass 10X filters.
#' @export
getCellBarcodes10X = function(bams,labels=names(bams)){
  #First try and get v3 IDs.
  out = list()
  for(i in seq_along(bams)){
    bam=bams[i]
    lab=labels[i]
    dir = list.files(dirname(bam))
    isV2 = 'filtered_gene_bc_matrices' %in% dir
    isV3 = 'filtered_feature_bc_matrix' %in% dir
    if(!isV2 && !isV3)
      stop("BAM file %s does not sit in a valid cellranger output directory.",bam)
    #Now load the filenames
    if(isV2){
      bcodes = file.path(dirname(bam),'filtered_gene_bc_matrices')
      bcodes = file.path(bcodes,list.files(bcodes)[1],'barcodes.tsv')
      bcodes = read.table(bcodes,sep='\t',header=FALSE)
    }
    if(isV3){
      bcodes = file.path(dirname(bam),'filtered_feature_bc_matrix','barcodes.tsv.gz')
      bcodes = read.table(bcodes,sep='\t',header=FALSE)
    }
    #Add label
    out[[i]] = paste0(lab,'_',bcodes[,1])
  }
  out = unlist(out,use.names=FALSE)
  return(out)
}


#' Get N distinct colours
#'
#' Gets N distinct colours.  Taken from here https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#'
#' @param n Number of colours
#' @importFrom grDevices hcl
getColoursN = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Get colours for a factor
#'
#' Given a factor, or something that can be made a factor, gets a vector of colours.  Named by unique values.
#'
#' @param x factor or thing to be turned into factor.
#' @return Vector of colours.
#' @export
getColours = function(x){
  #Ensure it's a factor
  x = factor(x)
  #Get the colours to use
  cols = getColoursN(nlevels(x))
  names(cols) = levels(x)
  return(cols[as.numeric(x)])
}

#' Get shape for a factor
#'
#' Given a factor, or something that can be made a factor, gets a vector of shapes.  Named by unique values.  Recycled if more than 25.
#'
#' @param x factor or thing to be turned into factor.
#' @return Vector of shapes.
#' @export
getShapes = function(x){
  #Ensure it's a factor
  x = factor(x)
  #Get the colours to use
  pchs = (seq(0,nlevels(x)-1) %% 25)+1
  names(pchs) = levels(x)
  return(pchs[as.numeric(x)])
}

#' Get number of cols/roms for plot grid
#'
#' Given a desired aspect ratio and a number of plots, decides on how many rows/columns are needed.
#'
#' @param nPlots Number of plots
#' @param r Aspect ratio.
#' @param nCol Number of columns.  Worked out if not fixed.
#' @param nRow Number of rows  Worked out if not fixed.
getRowColNums = function(nPlots,r=4/3,nCol=NULL,nRow=NULL){
  if(is.null(nCol) & is.null(nRow)){
    nCol = sqrt(nPlots*r)
    nRow = sqrt(nPlots/r)
    #Now we need to make an integer...
    nCol = ceiling(nCol)
    nRow = ceiling(nPlots/nCol)
  }
  if(is.null(nCol))
    nRow = ceiling(nPlots/nCol)
  if(is.null(nRow))
    nCol = ceiling(nPlots/nRow)
  if(nRow==1)
    nCol=nPlots
  if(nCol==1)
    nRow=nPlots
  return(c(nRow,nCol))
}

#' Build genes from annotation
#'
#' Given a \code{gCnts} object with gene information (usually from having run \code{\link{annotateSNPs}}), build a GRanges object with all the genes.
#'
#' @param gCnts GRanges object to construct genes object from.  Usually something on which \code{\link{annotateSNPs}} has been run.
#' @param geneID Vector or list of same length as \code{gCnts} which will be used to build output.
#' @return GRanges with gene ranges.
gnsFromAnnotation = function(gCnts,geneID = gCnts$geneID){
  if(is.null(geneID))
    stop("geneID not found in gCnts")
  gCnts$geneID = geneID
  gns = gCnts[!is.na(gCnts$geneID)]
  gns = split(gns[rep(seq_along(gns),lengths(gns$geneID))],unlist(gns$geneID))
  gns = unlist(range(GRangesList(gns)))
  gns$geneID = names(gns)
  return(gns)
}

