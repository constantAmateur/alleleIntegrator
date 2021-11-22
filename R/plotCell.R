#' Plots the evidence for one cell in detail.
#'
#' Plots summary of evidence for an individual cell.  SNPs are ordered such that those that provide the most weight are on the right.  Shape shows if the SNP is in an intron,exon, or intergenic.  The red/black line shows the tumour/normal BAF for each chrommosome.  The dots and error bars show the estimate and confidence interval for the BAF of this SNP aggregating across cells marked as normal.  The dotted line shows the cumulative BAF aggregating across SNPs from the left.
#'
#' @inheritParams calcStateProbs
#' @param cellID The ID of the cell to plot evidence for.
#' @param regions GRanges regions to plot independently.  If NULL, plot everything split by chromosome.  Can also be a vecetor of chromosomes to plot.
#' @param geneIDs List or vector of same length as \code{gCnts}.
#' @param orderBy Used to order SNPs in plot.  Anything in \code{gCnts} is valid and start with '-' to reverse order.  'pos' is genomic position.
#' @param spaceBy How to space things.  Anything in \code{gCnts} is valid.  'pos' is genomic position.  If NULL, space evenly.  Will space on the value of \code{spaceBy} unless \code{spaceBy} starts with a '-' in which case space on the difference.
#' @param pointShape Vector of same length as \code{gCnts} giving the pch value to use for each point.
#' @param pointColour Vector of same length as \code{gCnts} giving the colour to use for each point.
#' @param plotState The copy number state to plot.  Requires that \code{\link{calcLikelihoods}} have been run on \code{gCnts}.  If NULL, state not plotted.
#' @param plotGeneNames Include gene names below.
#' @param plotGeneLines Put separating lines at the boundary of each gene?
#' @param plotLegends Show legends?
#' @param ncol Number of columns to use in array of chromosomes.
#' @param nrow Number of columns to use in array of chromosomes.
#' @return A lovely lovely plot.
#' @export
#' @importFrom stats qbeta
#' @importFrom IRanges IRanges
#' @importFrom grDevices rgb
#' @importFrom graphics abline axis barplot lines par plot.new points rect segments text layout title
plotCell = function(gCnts,cellID,regions=gCnts@metadata$segs,geneIDs=gCnts$geneName,orderBy='pos',spaceBy=NULL,pointShape = getShapes(ifelse(gCnts$altIsMum,'matIsALT','matIsREF')),pointColour = getColours(gCnts$regionType),plotState=NULL,plotGeneNames=FALSE,plotGeneLines=FALSE,plotLegends=TRUE,ncol=NULL,nrow=NULL){
  #Param validation
  if(is.null(plotState)){
    stateFlag=FALSE
  }else{
    stateFlag=TRUE
    if(stateFlag && !any(grepl(paste0('^nLL_',plotState),colnames(mcols(gCnts))))){
      warning(sprintf("Likelihood for model %s not found in gCnts.  Has calcLikelihoods been run?",plotState))
      stateFlag=FALSE
    }
  }
  aseFlag=TRUE
  if(!all(c('matASE','matASE_postAlpha','matASE_postBeta') %in% colnames(mcols(gCnts)))){
    warning("No allele specific expression data available.")
    aseFlag=FALSE
  }
  if(!is.null(geneIDs) && length(geneIDs)!=length(gCnts))
    stop("Length of geneIDs must match gCnts")
  if(!cellID %in% gCnts$cellID)
    stop("Invalid cell ID")
  if(length(pointShape)==1)
    pointShape = rep(pointShape,length(gCnts))
  if(length(pointColour)==1)
    pointColour = rep(pointColour,length(gCnts))
  if(length(pointShape)!=length(gCnts))
    stop("Length of pointShape and gCnts must match")
  if(length(pointColour)!=length(gCnts))
    stop("Length of pointColour and gCnts must match")
  gCnts$pointShape = pointShape
  gCnts$pointColour = pointColour
  #Turn into character as needed
  gCnts$gns = geneIDs
  gCnts$gns[lengths(gCnts$gns)>1] = lapply(gCnts$gns[lengths(gCnts$gns)>1],paste,collapse=',')
  gCnts$gns = as.character(unlist(gCnts$gns))
  #Get the cell specific cnts
  cCnts = gCnts[gCnts$cellID == cellID]
  if(is.null(regions)){
    regions=unique(as.character(seqnames(cCnts)))
    regions = regions[order(suppressWarnings(as.numeric(regions)))]
  }
  #If it's a character, make it into GRanges
  if(!inherits(regions,'GRanges'))
    regions = GRanges(regions,IRanges(1,1e9))
  if(is.null(names(regions)))
    names(regions) = as.character(regions)
  #Keep only the ones with data to plot.
  regions = subsetByOverlaps(regions,cCnts)
  cCnts = subsetByOverlaps(cCnts,regions)
  #Decide how to plot things
  nPlots = length(regions)
  #Add an extra plot for legends if we need it
  nPlots = nPlots + plotLegends
  if(nPlots==0)
    stop("No data found overlapping segments specified.")
  pDims = getRowColNums(nPlots,nRow=nrow,nCol=ncol)
  #Construct the layout
  pMat = matrix(seq(2,(pDims[1]*pDims[2])+1),nrow=pDims[1],ncol=pDims[2],byrow=TRUE)
  pMat = rbind(rep(1,pDims[2]),pMat)
  #Fraction of plot taken up by title
  tmp = .05
  tmp = tmp/(1-tmp)*pDims[1]
  #Alter widths slightly to make room for labels
  widths = rep(1,pDims[2])
  #Make left column a bit wider for wider plot
  widths[1] = 1.15
  #Make right column narrower if there's only one and it's the legend column
  if(pDims[1]==1 & plotLegends)
    widths[length(widths)]=0.5
  layout(pMat,heights=c(tmp,rep(1,pDims[1])),widths=widths)
  par_old = par()
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.5,0.5,cellID,cex=2,font=2)
  for(i in seq_along(regions)){
    hh = ifelse(plotGeneNames,5,1)
    #Put more space for the one on the very left
    if((i-1)%%pDims[2]==0){
      par(mar=c(hh,3.5,2,2))
    }else{
      par(mar=c(hh,1,2,2))
    }
    region = regions[i]
    xx = subsetByOverlaps(cCnts,region)
    dd = mcols(xx)
    dd$tot = dd$matCount + dd$patCount
    dd$pos = start(xx)
    dd$ratio = dd$matCount/dd$tot
    #Decide how to order them.
    oo = orderBy[gsub('^-','',orderBy) %in% colnames(dd)]
    oo = do.call(order,lapply(oo,function(e) ifelse(grepl('^-',e),-1,1)*dd[,gsub('^-','',e)]))
    dd = dd[oo,]
    #Decide what to space by
    if(is.null(spaceBy)){
      dd$x = seq(0,nrow(dd)-1)
    }else{
      #Ensure that ordering can't be changed by the spacing
      if(grepl('^-',spaceBy)){
        dd$x = cumsum(c(0,abs(diff(dd[,gsub('^-','',spaceBy)]))))
      }else{
        dd$x = cumsum(abs(dd[,spaceBy]))
      }
    }
    dd$x = dd$x/max(dd$x)
    #Make the plot
    plot(0,0,
         type='n',
         xlim=c(0,1),
         ylim=c(0,1),
         xaxt='n',
         xlab='',
         ylab='',
         main=names(regions)[i],
         las=2,
         frame.plot=FALSE)
    title(ylab="maternal Frac", line=2.2)
    #Normal fraction
    abline(h=0.5)
    #Put box and whisker range for ASE
    if(aseFlag){
      lwd = min(.01,1/length(xx))
      dd$high = qbeta(0.95,dd$matASE_postAlpha,dd$matASE_postBeta)
      dd$low = qbeta(0.05,dd$matASE_postAlpha,dd$matASE_postBeta)
      segments(dd$x-lwd/2,dd$high,dd$x+lwd/2,dd$high,col='darkblue')
      segments(dd$x,dd$high,dd$x,dd$low,col='darkblue',lwd=0.3)
      segments(dd$x-lwd/2,dd$low,dd$x+lwd/2,dd$low,col='darkblue')
      points(dd$x,dd$matASE,pch=18,cex=1.0,col='darkblue')
    }
    #The progressive estimate
    lines(dd$x,cumsum(dd$matCount)/cumsum(dd$tot),lty=2,col='darkgreen')
    #Add count labels
    text(dd$x,1.05,
         labels=dd$tot,
         adj=c(0,0.5),
         cex=0.5,
         srt=90,
         xpd=NA)
    #And shaded line
    for(i in seq(nrow(dd))){
      if(i==1){
        lft=-0.05
      }else{
        lft = mean(dd$x[c(i-1,i)])
      }
      if(i==nrow(dd)){
        rht=1.05
      }else{
        rht = mean(dd$x[c(i,i+1)])
      }
      frac = dd$tot[i]/max(dd$tot)
      #Below 0.2 the line is white and hard to see.  So scale so this is the floor
      #c(0,1) -> c(0.2,1) implies (frac+x)/(1+x) where x = 0.2/(1-0.2)
      minShade = 0.1
      minShade = minShade/(1-minShade)
      frac = (frac+minShade)/(1+minShade)
      lines(c(lft,rht),c(1.025,1.025),lwd=2,col=rgb(1-frac,1-frac,1-frac))
    }
    #The observations for this cell
    points(dd$x,dd$ratio,pch=dd$pointShape,cex=1,col=dd$pointColour)
    #Add gene boxes
    if(!is.null(geneIDs)){
      aa=TRUE
      tmp = rle(dd$gns)
      for(ii in seq_along(tmp$values)){
        tgt = tmp$values[ii]
        idx = which(dd$gns==tgt)
        #Work out limits.
        if(min(idx)==1){
          lft = -0.05
        }else{
          lft = mean(dd$x[c(min(idx),min(idx)-1)])
        }
        if(max(idx)==nrow(dd)){
          rht = 1.05
        }else{
          rht = mean(dd$x[c(max(idx),max(idx)+1)])
        }
        rect(lft,-.1,rht,-0.05,xpd=NA,col=ifelse(aa,'lightgrey','darkgrey'))
        #Add name
        if(plotGeneNames)
          text(mean(c(lft,rht)),-.1,labels=tgt,srt=270,xpd=NA,adj=c(0,0.5))
        #Add lines
        if(plotGeneLines)
          abline(v=c(lft,rht),lwd=0.1)
        aa=!aa
      }
    }
    #Add model track
    if(stateFlag){
      tmp = as.matrix(dd[,grep('^nLL_',colnames(dd))])
      tmp = tmp -tmp[,paste0('nLL_',plotState)]
      tmp = do.call(cbind,lapply(seq(ncol(tmp)),function(e) cumsum(tmp[,e])))
      dd$stateProb = 1/rowSums(exp(-tmp))      
      lines(dd$x,dd$stateProb,lty=2,col='darkorange')
    }
  }
  #Make legends.  I can't make legend function do what I want so....
  if(plotLegends){
    par(mar=c(0,0,0,0))
    plot.new()
    titleBuff=0.03
    text(0.5,1,'Shape',adj=0,xpd=NA,font=2)
    #Place things evenly in gap
    tmp = seq(1-titleBuff,0.66+titleBuff,length.out=length(unique(pointShape))+1)
    tmp = tmp[-length(tmp)]+diff(tmp)/2
    points(rep(0,length(tmp)),tmp,pch=unique(pointShape))
    text(rep(0.2,length(tmp)),tmp,names(pointShape)[match(unique(pointShape),pointShape)],adj=0)
    text(0.5,0.66,'Colour',adj=0,font=2)
    tmp = seq(0.66-titleBuff,0.33+titleBuff,length.out=length(unique(pointColour))+1)
    tmp = tmp[-length(tmp)]+diff(tmp)/2
    points(rep(0,length(tmp)),tmp,pch=19,col=unique(pointColour))
    text(rep(0.2,length(tmp)),tmp,names(pointColour)[match(unique(pointColour),pointColour)],adj=0)
    text(0.5,0.33,'Lines',adj=0,font=2)
    tmp = seq(0.33-titleBuff,0,length.out=2+1)
    tmp = tmp[-length(tmp)]+diff(tmp)/2
    lines(c(0,0.15),c(tmp[1],tmp[1]),lty=2,col='darkgreen')
    lines(c(0,0.15),c(tmp[2],tmp[2]),lty=2,col='darkorange')
    text(rep(0.2,length(tmp)),tmp,c('matFrac',paste0('prob_',plotState)),adj=0)
  }
}

