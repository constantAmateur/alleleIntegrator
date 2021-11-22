#' Plot normalised counts split by some variable
#'
#' This function will plot one bar plot per grouping, which is normalised to sum to 1.  Each bar is split into groups of different colours and further sub-divisions indicated by black boxes.  The default is to show some genomic grouping (typically het SNPs or genes), colour the bar by maternal/paternal counts, and indicate contributions from individual cells with boxes.
#'
#' This seemed like a good idea at the time.  I'm sure there was some point to this, but I can't remember what it was.
#' 
#' @inheritParams aggregateByLists
#' @param assayCols Colours for assays.  Must be the same length as \code{assayCols}.
#' @param transposeMatrix Should count matricies built from assays be transposed?  This swaps the meaning of rows and columns. 
#' @param orderBy How to order columns.  Should be either the name of an assay or 'counts' for total counts.
#' @param countMax Function to apply to vector of transformed counts to get the upper limit for plotting.
#' @param countTrans What transformation should be applied to the counts?
#' @param countTransLab What label to give the count legend?
#' @param xlab Label for x-axis.
#' @param plotNames Should we plot the names of the columns?  If NULL, plots if fewer than 20 entries.
#' @param mar Margins for plot.
#' @param lwdInside The width for lines separating columns.
#' @param lwdBox The width for lines separating parts of columns.
#' @param showLegends Plot legends to right of plot?
#' @return A plot, obs.
#' @export
#' @importFrom circlize colorRamp2
#' @importFrom graphics legend
plotSplit = function(gCnts,assays=c('matCount','patCount'),assayCols=c('darkred','purple'),transposeMatrix=FALSE,orderBy='counts',countMax=max,countTrans=log10,countTransLab='(log10)',xlab=ifelse(transposeMatrix,'cell','region'),plotNames=NULL,mar=c(5,4,4,5),lwdInside=0.1,lwdBox=0.5,showLegends=TRUE){
  mtx = buildCountMatricies(gCnts,assays=assays)
  if(transposeMatrix)
    mtx = lapply(mtx,t)
  nRows = nrow(mtx[[1]])
  nCols = ncol(mtx[[1]])
  nCounts = colSums(do.call(rbind,mtx))
  if(is.null(plotNames))
    plotNames = nCols<20
  #Order each block
  mtx = lapply(mtx,apply,2,sort)
  #Collapse them into groups and normalise
  mm = as.matrix(do.call(rbind,lapply(mtx,colSums)))
  mm = t(t(mm)/colSums(mm))
  #Expand them out and normalise
  mtx = as.matrix(do.call(rbind,mtx))
  mtx = t(t(mtx)/colSums(mtx))
  #Order by the first one
  if(orderBy == 'counts'){
    o = order(nCounts)
  }else{
    o = order(mm[orderBy,])
  }
  mtx = mtx[,o]
  mm = mm[,o]
  nCounts = nCounts[o]
  #Do the plotting
  #Make space at top for track and right for legend
  par(mar=mar)
  barplot(mm,
          space=0,
          axisnames=plotNames,
          las=2,
          border = NA,
          xlab=xlab,
          ylab = 'Count fraction',
          col = assayCols)
  for(i in seq(ncol(mtx))){
    lines(c(i,i),c(0,1),lwd=lwdInside)
    for(y in unique(cumsum(mtx[,i]))){
      lines(c(i-1,i),c(y,y),lwd=lwdBox)
    }
  }
  #Add bar at top showing counts
  cMap = colorRamp2(c(0,countMax(countTrans(nCounts))),c('white','black'))
  for(i in seq(ncol(mtx)))
    rect(i-1,1.01,i,1.05,
         col=cMap(countTrans(nCounts[i])),
         border=NA,
         xpd=NA)
  #Add legends
  if(showLegends){
    xLeft = nCols*1.01
    #Main legend
    legend(xLeft,0.3,
           legend=assays,
           fill=assayCols,
           bty='n',
           title='Count\nType',
           xpd=NA
           )
    #Legned for count bar
    cPts = c(0,countMax(countTrans(nCounts))/2,countMax(countTrans(nCounts)))
    legend(xLeft,0.8,
           legend=round(cPts,1),
           fill = cMap(cPts),
           bty='n',
           xpd=NA,
           title=paste0('# Counts\n',countTransLab)
           )
  }
}
