#' Calculates posterior probability for a set of copy number states in genomic regions
#'
#' Given a set of regions, and a set of candidate copy number states to test, calculates the posterior probability of each copy number state, in each cell, in each region.  
#'
#' There are really two common use cases for this function.  In the first, the user has no prior expectation of what CN changes might be present in which cells.  The defaults are tuned to this scenario.  In this case, all plausibly detectable CN states are specified and posteriors calculated in each chromosome indepedently.
#'
#' The other common use case is when a series of pre-defined copy number segments are supplied.  In this case each region has its own copy number state and \code{CNStates} should specify which columns of \code{regions} contain the expected BAF for the tumour and normal state.
#'
#' By default, this will simply calculate the posterior probabilities with all available data.  However, if \code{removeN} is set to something greater than 0 a second calculation is performed with \code{removeN} bits of data removed per cell.  The logic here is that we want to know if the posterior probability is highly dependent on just a few regions.  If so, then it may be an artefact of the allele specific expression in those regions.  The default will try and remove the \code{removeN} genes that have the largest contribution to the posterior probability for each cell.  However, you can remove the top \code{removeN} SNPs or some other region by specifying a different grouping to \code{removeGroup}.
#'
#' @inheritParams calcLikelihoods
#' @param regions GRanges object giving regions to calculate probabilities over.  Can also be a character vector specifying chromosomes.  If NULL, set to all chromosomes in \code{gCnts}.
#' @param CNStates What copy number states to test in each region.  Can be specified as a named vector (the default), in which case all regions are tested for the same CN states and the meaning is as in \code{\link{calcLikelihoods}}.  Alternatively, a pair of names can be given, specifying the names of columns in \code{regions} containing the allele fraction for tumour and normal respectively.
#' @param removeGroup Group things by this vector/list when removing.  Typically genes, but can be anything (position, whatever).  If NULL, no grouping, just remove SNPs.
#' @param removeN How many things to remove to check that we still have a significant result.
#' @param calcJointProb Include a genome-wide estimate of state probabilities as well?
#' @param verbose Be verbose?
#' @param ... Passed to \code{\link{calcLikelihoods}}.
#' @return A data.frame giving the posterior probabilities of each region being tumour and normal.  A genome wide posterior probability is also calculated and returned for each cell.
#' @importFrom IRanges subsetByOverlaps IRanges
#' @importFrom utils txtProgressBar setTxtProgressBar head
#' @importFrom GenomeInfoDb seqlevels
#' @export
calcStateProbs = function(gCnts,overDisp,regions=NULL,CNStates = c(pLoss=1,mGain=2/3,diploid=1/2,pGain=1/3,mLoss=0),removeGroup=gCnts$geneID,removeN=0,calcJointProb=FALSE,verbose=TRUE,...){
  #Convert regions to GRanges
  if(is.null(regions))
    regions = seqlevels(gCnts)
  if(is.character(regions))
    regions = GRanges(regions,IRanges(1,1e9))
  if(is.null(names(regions)))
    names(regions) = as.character(regions)
  if(any(duplicated(names(regions))))
    stop("Regions must be uniquely named.")
  #Save removeGroup and default to SNPs if none given.
  gCnts$removeGroup = removeGroup
  if(is.null(removeGroup))
    gCnts$removeGroup = as.character(gCnts)
  #Keep only the counts we need
  gCnts = subsetByOverlaps(gCnts,regions)
  #Decide if we're in tumour mode or not
  if(length(CNStates)==2 && all(is.character(CNStates)) && all(CNStates %in% colnames(mcols(regions)))){
    modelNames = paste0('nLL_',CNStates)
    tumourMode=TRUE
  }else{
    #Validate input
    if(!all(is.numeric(CNStates)) || is.null(names(CNStates)) || any(duplicated(names(CNStates))))
      stop("Invalid CNStates (see help).  Expected either two column names for regions, or a named numeric vector, where names are unique.")
    modelNames = paste0('nLL_',names(CNStates))
    tumourMode=FALSE
  }
  #Do each region independently.  This allows for regions to potentially overlap and for CNstates to vary.
  #Calculate each likelihoods one region at a time.
  ans = list()
  for(iReg in seq_along(regions)){
    region = regions[iReg]
    if(verbose)
      message(sprintf("Calculating state probabilities in region %s",names(region)))
    if(tumourMode){
      statesLocal = unlist(mcols(regions)[iReg,CNStates])
    }else{
      statesLocal = CNStates
    }
    lCnts = calcLikelihoods(subsetByOverlaps(gCnts,region),overDisp,statesLocal,...)
    #Collapse into regions
    out = aggregateByRegions(lCnts,region,assays=c('matCount','patCount',modelNames))
    out$totCount = out$matCount + out$patCount
    #Force storage as the correct type
    out$nInfSNPs = as.integer(table(lCnts$cellID)[out$cellID])
    out$nInfGroups = sapply(split(unlist(lCnts$removeGroup),rep(lCnts$cellID,lengths(lCnts$removeGroup))),function(e) length(unique(e)))[out$cellID]
    #And the posterior probabilities.  Do this in log space to enhance accuracy
    nLLs = as.matrix(mcols(out)[,modelNames])
    pp = nLLs
    for(i in seq(ncol(pp))){
      pp[,i] = rowSums(exp(nLLs[,i]-nLLs))**-1
    }
    colnames(pp) = gsub('nLL_','postProb_',colnames(pp))
    mcols(out) = cbind(mcols(out),pp)
    #Some other useful values
    out$maxPostProb = apply(pp,1,max)
    out$mostLikelyState = gsub('^postProb_','',colnames(pp)[apply(pp,1,which.max)])
    out$matFrac = out$matCount/out$totCount
    #Keep the interesting ones in interesting order
    toKeep = c(colnames(region),'regionID','cellID','clusterID','matFrac','matCount','patCount','totCount','nInfSNPs','nInfGroups',modelNames,colnames(pp),'maxPostProb','mostLikelyState')
    toKeep = intersect(toKeep,colnames(mcols(out)))
    out = out[,toKeep]
    #Do the removal if we're going to
    if(removeN>0){
      #Decide which things to drop
      #if(verbose)
      #  message(sprintf("Calcuating which %d entries to remove from each call...",removeN))
      #Aggregate it by the groups we're going to remove
      mtx = aggregateByLists(lCnts,modelNames,cellList = lCnts$cellID,regionList = lCnts$removeGroup)
      #Get the baseline.  Easier to work with as a data.frame
      full = as.data.frame(mcols(out))
      #full = setNames(full$pp,full$cellID)
      marks = unique(mtx$cellID)
      #if(verbose)
      #  pb = txtProgressBar(0,length(marks),style=1)
      #Decide which things to remove.  Select them to try and kill the strongest signal.
      #This isn't perfect.  In the case of overlapping regions it will remove any count that is in the top N for any region.  But that's unlikely to come up so eh.  I guess the ideal solution here would be to run each region independently.
      toDrop = list()
      nDropped = rep(0,length(marks))
      names(nDropped) = marks
      for(mark in marks){
        x = mtx[mtx$cellID==mark,,drop=FALSE]
        rids = x$regionID
        #Work out what the likelihoods would be without each gene
        x = t(unlist(full[full$cellID==mark,modelNames])-t(x[,modelNames,drop=FALSE]))
        #Convert them to probabilities
        x = exp(-x)/rowSums(exp(-x))
        #What is the change for the one that was the winner before we started fucking around and dropping data
        pDiff = x[,paste0('nLL_',unlist(full[full$cellID==mark,'mostLikelyState']))] - unlist(full[full$cellID==mark,'maxPostProb'])
        #Decide which ones to drop.  Only drop things that are going to make the call worse
        o = order(pDiff)
        oo = head(o,removeN)
        oo = oo[pDiff[oo]<0]
        nDropped[mark] = length(oo)
        if(length(oo)>0){
          rids = rids[oo]
          #Construct the naughty combination of cellID and removeGroup that we'll want to drop
          toDrop[[mark]] = paste0(mark,'___',rids)
        }
        #Progress
        #if(verbose)
        #  setTxtProgressBar(pb,i)
        i=i+1
      }
      #if(verbose)
      #  close(pb)
      #Now drop anything that is on the naughty list of cell/removeGroup combinations
      toDrop = unlist(toDrop,use.names=FALSE)
      #Need to do some expansion and contraction here too.
      toDropIdx = paste0(rep(lCnts$cellID,lengths(lCnts$removeGroup)),'___',unlist(lCnts$removeGroup))
      toDropIdx = which(sapply(relist(toDropIdx %in% toDrop,lCnts$removeGroup),any))
      if(length(toDropIdx)>0)
        lCnts = lCnts[-toDropIdx]
      #Do the sexy recursion
      #if(verbose)
      #  message(sprintf("Recalculating after removing %d top genes (or similar)",removeN))
      tmp = calcStateProbs(lCnts,overDisp,regions=region,CNStates = statesLocal,removeGroup=lCnts$removeGroup,removeN=0,calcJointProb=FALSE,verbose=FALSE,...)
      tmp = as.data.frame(mcols(tmp))
      #Do this cautiously, as the merge step is prone to breaking
      m = match(tmp$cellID,out$cellID)
      for(nom in colnames(tmp)){
        if(nom %in% c('regionID','cellID'))
          next
        #Create target 
        tgt = paste0('removeN_',nom)
        mcols(out)[,tgt]=NA
        mcols(out)[m,tgt] = tmp[,nom]
      }
      #Add some summary of the difference metrics
      out$nDropped = nDropped[out$cellID]
      out$nDropped[is.na(out$nDropped)]=0
      out$matFracDiff = out$matFrac - out$removeN_matFrac
      out$consistentState = out$mostLikelyState == out$removeN_mostLikelyState
    }
    ans[[iReg]] = out
  }
  #Merge back into one thing
  ans = unlist(GRangesList(ans))
  #Create genome wide assay if we want that.  I really hate this part, it would probably be better to exclude it but it is useful
  if(calcJointProb){
    #Aggregate things sensibly 
    toAgg = c('matCount','patCount','totCount','nInfSNPs','nInfGroups',grep('^nLL_',colnames(mcols(ans)),value=TRUE))
    if(removeN>0)
      toAgg = c(toAgg,paste0('removeN_',toAgg))
    toAgg = intersect(toAgg,colnames(mcols(ans)))
    tmp = aggregate(as.data.frame(mcols(ans))[,toAgg,drop=FALSE],list(cellID=ans$cellID),FUN=sum,na.rm=TRUE)
    #Format it in the same way as the main bits
    base = GRanges(rep('genomeWide',nrow(tmp)),IRanges(1,1e9))
    for(nom in colnames(mcols(ans))){
      mcols(base)[,nom]=NA
      if(nom %in% toAgg)
        mcols(base)[,nom] = tmp[,nom]
    }
    base$cellID = tmp$cellID
    #Fill in the derived values.  This very much violates the don't repeat yourself rule, but it's not worth refactoring the code for
    for(pre in c('','removeN_')){
      if(pre=='removeN_' & removeN<=0)
        next
      nLLs = as.matrix(mcols(base)[,paste0(pre,modelNames)])
      for(i in seq(ncol(pp)))
        mcols(base)[,gsub('nLL_','postProb_',colnames(nLLs)[i])] = rowSums(exp(nLLs[,i]-nLLs))**-1
      #Some other useful values
      nLLs = as.data.frame(mcols(base)[,gsub('nLL_','postProb_',colnames(nLLs))])
      mcols(base)[,paste0(pre,'maxPostProb')] = apply(nLLs,1,max)
      mcols(base)[,paste0(pre,'mostLikelyState')] = gsub(paste0('^',pre,'postProb_'),'',colnames(nLLs)[apply(nLLs,1,which.max)])
      mcols(base)[,paste0(pre,'matFrac')] = mcols(base)[,paste0(pre,'matCount')]/mcols(base)[,paste0(pre,'totCount')]
    }
    #Do some final calculations if we've removed things
    if(removeN>0){
      base$matFracDiff = base$matFrac - base$removeN_matFrac
      base$consistentState = base$mostLikelyState == base$removeN_mostLikelyState
      base$nDropped = base$nInfGroups - base$removeN_nInfGroups
    }
    #Update regionID if it's there
    if(!is.null(base$regionID))
      base$regionID='genomeWide'
    #Update clusterID if it's there
    if(!is.null(ans$clusterID))
      base$clusterID = ans$clusterID[match(base$cellID,ans$cellID)]
    #Finally merge it back in
    ans = suppressWarnings(c(ans,base))
  }
  return(ans)
}
