#' Plot coverage and BAF
#'
#' Takes a \code{GRanges} object for which coverage and BAF have been previously calculated using \code{\link{generateCoverageAndBAF}} and plots the coverage, BAF, and ASCAT output (if calculated).
#'
#' @param hSNPs Output of \code{\link{generateCoverageAndBAF}}
#' @param useLogR Should logR be plotted instead of raw coverage?
#' @param linesBAF Draw horizontal lines at these BAFs.
#' @param alpha Transparency for points used in plotting.
#' @param minCoverage Don't plot points with coverage less than this.
#' @param plotASCAT Add a third panel showing the ASCAT results.
#' @param useCorrectedLogR If GC corrected logR values are available (from ASCAT) and \code{useLogR=TRUE} plot corrected values. 
#' @param chrsToPlot Which chromosomes to plot.
#' @return Nothing, but makes a pretty plot.
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom stats lowess
#' @export
plotCoverageAndBAF = function(hSNPs,useLogR=TRUE,linesBAF=c(),alpha=max(0.002,min(1,1e5/length(hSNPs))),minCoverage=10,plotASCAT = !is.null(hSNPs@metadata$ascat.output),useCorrectedLogR=TRUE,chrsToPlot=c(1:22,'X')) {
  #Filter to just the ones that we trust
  filt = hSNPs[hSNPs$coverage>=minCoverage,]
  #Work out if ASCAT has been run
  if(plotASCAT & is.null(filt@metadata$ascat.output)){
    warning("ASCAT results not found, plotting segments disabled")
    plotASCAT=FALSE
  }
  #Plot coverage and BAF in one plot
  if(plotASCAT){
    layout(matrix(1:3,ncol=1),heights=c(2,2,1))
  }else{
    par(mfrow=c(2,1))
  }
  par(mar=c(2.1,4.1,1.1,1.1))
  #Work out the chromosome boundaries
  chrs = chrsToPlot
  chrLens = seqlengths(filt)
  tmp = sapply(split(start(filt),as.character(seqnames(filt))),max)
  chrLens[is.na(chrLens)] = tmp[names(chrLens)[is.na(chrLens)]]
  chrLens = as.numeric(chrLens[chrs])
  if(useLogR){
    covMin=-1
    covMax=1
    if(useCorrectedLogR && !is.null(filt$correctedLogR)){
      dat = filt$correctedLogR
    }else{
      dat = filt$logR
    }
  }else{
    covMin = 0
    covMax = quantile(filt$coverage,0.99)*1.25
    #Add some vertical jitter to help visualisation
    dat = filt$coverage + runif(length(filt),min=-0.5,max=0.5)
  }
  x = start(filt) +cumsum(c(0,chrLens))[match(as.character(seqnames(filt)),chrs)]
  plot(x,dat,
       col=rgb(0,0,0,alpha=alpha),
       cex=0.01,
       las=2,
       xaxt='n',
       ylab=ifelse(useLogR,'logR','Coverage'),
       xlab='',
       xaxs='i',
       yaxs='i',
       ylim=c(covMin,covMax),
       xlim=c(1,sum(chrLens))
       )
  axis(side=1, at=cumsum(c(0,chrLens[-length(chrLens)]))+chrLens/2, labels = chrs)
  #axis(side=2, at=c(0,covMax,covMax*2),labels=c(covMax,0.5,1),las=1)
  abline(v=cumsum(chrLens),col='lightgrey')
  #Add loess smoothing for each chromosome.
  for(chr in chrs){
    w = as.character(seqnames(filt))==chr
    xx = x[w]
    y = dat[w]
    #Lowess fit of straight lines (degree=1) on 10% of data (span = 0.1)
    lo = lowess(y~xx,f=0.1)
    lines(lo,col='red')
  }
  #BAF plot
  if(!plotASCAT)
    par(mar=c(5.1,4.1,1.1,1.1))
  plot(x,filt$BAF,
       col=ifelse(!filt$isHet,rgb(0,255/255,0,alpha=alpha/1),rgb(0,0,0,alpha=alpha)), #Mark homozygous in darkgreen
       cex=0.01,
       las=2,
       xaxt='n',
       yaxt='n',
       xlab=ifelse(plotASCAT,'','Chromosomes'),
       ylab='BAF',
       xaxs='i',
       yaxs='i',
       ylim=c(0,1),
       xlim=c(1,sum(chrLens))
       )
  axis(side=1, at=cumsum(c(0,chrLens[-length(chrLens)]))+chrLens/2, labels = chrs)
  axis(side=2, at=c(0,0.5,1),labels=c(0,0.5,1),las=1)
  abline(v=cumsum(chrLens),col='lightgrey')
  if(length(linesBAF)>0)
    abline(h=linesBAF,col='black',lwd=1)
  #Add lowess smooths of het SNPs
  for(chr in chrs){
    #Only want to fit to het SNPs
    w = which(as.character(seqnames(filt))==chr & filt$isHet)
    #Split in two based on BAF
    isUp = filt$BAF[w]>0.5
    isDown = filt$BAF[w]<0.5
    y = filt$BAF[w[isUp]]
    xx = x[w[isUp]]
    lines(lowess(y~xx,f=0.1),col='darkblue')
    y = filt$BAF[w[isDown]]
    xx = x[w[isDown]]
    lines(lowess(y~xx,f=0.1),col='darkred')
  }
  #Add ASCAT results if you have them
  if(plotASCAT){
    par(mar=c(5.1,4.1,1.1,1.1))
    plot(0,
         type='n',
         las=2,
         ylim=c(0,3),
         xlim=c(1,sum(chrLens)),
         xaxs='i',
         ylab='nCopies',
         xlab='Chromosomes',
         yaxt='n',
         xaxt='n')
    segs = filt@metadata$ascat.output$segments
    segs$xStart = segs$startpos +cumsum(c(0,chrLens))[match(segs$chr,chrs)]
    segs$xEnd = segs$endpos +cumsum(c(0,chrLens))[match(segs$chr,chrs)]
    segs$nMajor = pmin(segs$nMajor,3)
    segs$nMinor = pmin(segs$nMinor,3)
    for(i in seq(nrow(segs))){
      #Offsets so they don't overlap if they're the same
      lines(c(segs$xStart[i],segs$xEnd[i]),c(segs$nMajor[i],segs$nMajor[i])+.1,col='red',lwd=2)
      lines(c(segs$xStart[i],segs$xEnd[i]),c(segs$nMinor[i],segs$nMinor[i])+.0,col='green',lwd=2)
    }
    axis(side=1, at=cumsum(c(0,chrLens[-length(chrLens)]))+chrLens/2, labels = chrs)
    axis(side=2, at=c(0,1,2,3),labels=c(0,1,2,3),las=1)
    abline(v=cumsum(chrLens),col='lightgrey')
  }
}
