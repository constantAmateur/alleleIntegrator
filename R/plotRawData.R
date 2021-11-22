#' Plots a high level summary of raw data
#'
#' Summarises data at the cell (and cluster if available) level by either chromosome or segment specified.  The default plot shows how each cluster of cells behave in each segment.  Each point represents a cell, the boxplot shows the distribution, and the big red points show the value treating each cluster as one cell.  The colour scheme shows low/medium/high coverage for each cell/segment combination.
#'
#' If \code{returnData=TRUE}, this function will return two tables that give the summary statistics for each cell/segment combination (\code{cellCnts}) and cluster/segment combination (\code{clustCnts}).
#'
#' @param gCnts GRanges object with single cell level count data.  Usually output of \code{\link{getAllelicExpression}}.
#' @param clusterID  Vector of the same length as \code{gCnts} indicating how to group cells into clusters.  If NULL, all cells treated as one cluster.
#' @param segs Genomic segments within which data will be summarised.  If NULL, summarised by chromosome.
#' @param tgtMAF The target maternal allele fraction for each segment.  Must be same length as \code{segs} or NULL, in which case it will default to 0.5 everywhere.
#' @param returnData Return data as well as plot object.
#' @param showPlot Don't just return the plot object, explicitly plot it.
#' @return Either a ggplot object or if \code{returnData=TRUE} a list containing the cell and cluster summary tables and the ggplot object for the plot.
#' @importFrom ggplot2 ggplot aes geom_hline geom_jitter geom_boxplot facet_wrap ylim ylab theme element_blank element_line element_text guides guide_legend geom_point scale_colour_manual xlab
#' @export
plotRawData = function(gCnts,clusterID=gCnts$clusterID,segs=gCnts@metadata$segs,tgtMAF = segs$tumFrac,returnData=FALSE,showPlot=returnData){
  #Store cluster and segment IDs in main object
  if(is.null(clusterID))
    clusterID = rep('All',length(gCnts))
  gCnts$clusterID = clusterID
  if(is.null(segs)){
    segs = unique(as.character(seqnames(gCnts)))
    segs = GRanges(segs,IRanges(1,1e9))
    names(segs) = as.character(seqnames(segs))
  }
  gCnts$segID = names(segs)[findOverlaps(gCnts,segs,select='first')]
  #Generate counts by cell and segments
  cCnts = aggregateByRegions(gCnts,segs,assays=c('matCount','patCount'))
  #Put everything together if no clusters given
  cCnts$clusterID = gCnts$clusterID[match(cCnts$cellID,gCnts$cellID)]
  cCnts$totCount = cCnts$patCount+cCnts$matCount
  cCnts$MAF = cCnts$matCount/cCnts$totCount
  #And by cluster and segment
  clCnts = aggregateByLists(gCnts,assays=c('matCount','patCount'),gCnts$clusterID,gCnts$segID)
  colnames(clCnts) = gsub('^cellID$','clusterID',colnames(clCnts))
  clCnts$totCount = clCnts$patCount + clCnts$matCount
  clCnts$MAF = clCnts$matCount/clCnts$totCount
  #Segment ratios.  This is just for plotting
  tmp = data.frame(mcols(segs))
  tmp$regionID = names(segs)
  #Draw target line at 0.5 if tumour fraction not given
  if(is.null(tgtMAF))
    tgtMAF = rep(0.5,nrow(tmp))
  tmp$tumFrac = tgtMAF
  #Here's a version that is aware of the coverage
  cCnts$covBins = cut(cCnts$totCount,breaks=c(0,10,20,Inf))
  gg = ggplot(data.frame(mcols(cCnts)),aes(clusterID,MAF)) +
    geom_hline(data=tmp,aes(yintercept=tumFrac),colour='red',linetype='dashed') +
    geom_jitter(aes(colour=covBins),size=0.2,alpha=1/1,height=0) +
    geom_boxplot(outlier.shape=NA,alpha=1/100,lwd=1/4) +
    facet_wrap(~regionID) +
    ylim(0,1) + 
    ylab('Maternal allele frequency') +
    xlab('Cluster') + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
    guides(colour = guide_legend(override.aes = list(size=2),title='Coverage')) +
    geom_point(data=clCnts,colour='red') + 
    scale_colour_manual(values=c('#1b9e77','#d95f02','#7570b3'))
  if(showPlot)
    plot(gg)
  #Drop plotting column
  cCnts$covBins=NULL
  #Format output
  if(returnData){
    return(list(cellCnts = cCnts,clustCnts=clCnts,plot=gg))
  }else{
    return(gg)
  }
}
