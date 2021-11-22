#' Phases heterozygous SNPs in a patient using CN segments
#'
#' Given heterozygous SNPs in an individual, as called using \code{\link{findHetSNPs}}, and a series of segments where a copy number change has occurred, checks genotype in Tumour BAM files at these locations.  This information is used to identify which allele is major and minor at each SNP.  The major allele is arbitrarily labelled as "maternal" for consistency with downstream functions.
#'
#' The phasing of SNPs relies on the major allele having a detectable and statistically significant deviation from the other allele.  As such, the larger the allelic imbalance, the more SNPs will be phased.  Note that for tumour samples, the allelic imbalance is determined by the nature of the copy number change and the amount of normal contamination in the tumour biopsy.  For subtle shifts (e.g. 3/2 CN ratio, 80% contamination), it may not be possible to phase many SNPs and those SNPs that are phased may not be particularly reliable.
#'
#' If \code{useEM=FALSE}, a binomial test to exclude a null hypothesis of allele ratio 0.5 is done and any SNP with false discover rate less than \code{FDR} are phased.  This is almost, but not exactly, equivalent to the default approach, which is to use an EM algorithm to model the SNPs as a mixture of two populations with binomial likelihood.  In this mode, (the default, \code{useEM=TRUE}), the allelic ratio is calculated within each segment and then used to estimate the probability of each SNP belonging to each allele.  This allows a check that the segment is consistent with a CN change being present in the tumour BAM file.  In this mode \code{FDR} is interpreted as how far the allele probability needs to be from 1 to phase a SNP.  For example, \code{FDR=0.01} implies that SNPs with major allele probabilities \code{0.99} and greater or \code{0.01} and lower will be used.  In EM mode the \code{stateProbs} column in the output represents the major allele probability, while when \code{useEM=FALSE} it represents the binomial test p-value.
#'
#' Note that both BAMs (the one used to generate \code{hSNPs} and \code{tBAM}) must have been mapped to the same reference genome, \code{refGenome}.
#'
#' Segments are sanity checked for the presence of an allelic imbalance.  If fewer than \code{minPhasable} SNPs can phased, the entire segment is flagged is likely to have failed. 
#'
#' @inheritParams findHetSNPs
#' @param hSNPs The heterozygous SNP locations.  See \code{\link{findHetSNPs}}.
#' @param segs A GRanges object giving the genomic coordinates where copy number changes that produce allelic imbalance are present.
#' @param tBAM BAM file containing tumour genotype.
#' @param outPath File name used to storing genotyping calculations.  If left at NULL, results are not saved and cannot be reloaded in subsequent runs.
#' @param minPhasable If the fraction of phasable SNPs in a segment is less than this, flag the segment as likely improperly called.
#' @param useEM Use the expectation-maximisation method to phase SNPs.  See details.
#' @param alleleCounterParams Parameters to be used when running allele counter.
#' @param plotMixtures If TRUE and using EM mode produces one plot for each CN segment showing the mixture solutions.
#' @param verbose Be verbose?  FALSE/0 for not at all, TRUE or 1 for a bit, 2 for very.
#' @param ... Does nothing.
#' @return \code{hSNPs} but with extra columns indicating results of phasing.
#' @importFrom stats runif
#' @importFrom graphics hist
#' @export
phaseSNPsFromCN = function(hSNPs,segs,refGenome,tBAM,outPath=NULL,FDR=0.05,minPhasable = 2*FDR,errRate=0.01,useEM=TRUE,alleleCounterParams=list(f=3,F=3852,m=20,q=35),nParallel=1,plotMixtures=TRUE,verbose=TRUE,...){
  #Ensure segments are named uniquely
  if(is.null(names(segs)))
    names(segs) = as.character(segs)
  if(any(duplicated(names(segs))))
    names(segs) = make.unique(names(segs))
  ##############################
  # Run alleleCounter as needed
  #Set the non-optional allelecounter parameters
  alleleCounterParams$x=FALSE
  alleleCounterParams$bams = tBAM
  #alleleCounterParams$tgtLoci = subsetByOverlaps(hSNPs,segs) #Only bother phasing SNPs that overlap CN changes
  #Just do it everywhere, the time cost is pretty small and it's easier to re-run
  alleleCounterParams$tgtLoci = hSNPs
  alleleCounterParams$refGenome = refGenome
  alleleCounterParams$outputs = outPath
  alleleCounterParams$nParallel = nParallel
  tCnts = do.call(alleleCounter,alleleCounterParams)[[1]]
  ####################
  # Determine phasing
  bases=c('A','C','G','T')
  #Get counts
  cnts = mcols(tCnts)
  tCnts$altCount = (as.matrix(cnts[,bases])[cbind(seq(nrow(cnts)),match(cnts$ALT,bases))])
  tCnts$refCount = (as.matrix(cnts[,bases])[cbind(seq(nrow(cnts)),match(cnts$REF,bases))])
  #Match to hSNPs
  m = match(tCnts,hSNPs)
  hSNPs$altCountTum = hSNPs$refCountTum = hSNPs$totCountTum = hSNPs$altIsMum = hSNPs$matProb = hSNPs$passSanity = NA
  hSNPs$informative = FALSE
  hSNPs$altCountTum[m] = tCnts$altCount
  hSNPs$refCountTum[m] = tCnts$refCount
  hSNPs$totCountTum[m] = tCnts$Tot
  hSNPs$passSanity[m] = TRUE
  #Use the EM algorithm to validate the CN change and phase
  if(useEM){
    #For each segment, work out the ideal mixture using EM.  
    for(ii in seq_along(segs)){
      if(verbose)
        message(sprintf("Fitting mixture model for segment %s",names(segs[ii])))
      tmp = subsetByOverlaps(hSNPs,segs[ii])
      aC = tmp$altCountTum
      tC = tmp$totCountTum
      ####The constrained version.  Not very useful for detecting bad configs.
      ###rho = runif(k)
      ###tau = 0.5
      ###i = 0
      ###Q = -Inf
      ###while(TRUE){
      ###  #EStep.  In log space for stability
      ###  states = cbind(dbinom(aC,tC,rho,log=TRUE),dbinom(aC,tC,1-rho,log=TRUE))
      ###  states = t(t(states)+log(c(tau,1-tau)))
      ###  #Another numerical stability fix.  Subtraction in log space equivalent to multiplication by a constant which changes nothing cebause normalisation
      ###  stateProbs = exp(states-apply(states,1,max))
      ###  stateProbs = stateProbs/rowSums(stateProbs)
      ###  #Mstep
      ###  Qnew = sum(stateProbs[,1]*(log(tau)+dbinom(aC,tC,rho,log=TRUE)) +
      ###             stateProbs[,2]*(log(1-tau) + dbinom(aC,tC,1-rho,log=TRUE)))

      ###  p=stateProbs[,1]
      ###  #This is the analytic solution for the minimum of Q with respect to rho and tau
      ###  rho = sum(p*aC+(1-p)*(tC-aC))/sum(tC)
      ###  tau = mean(stateProbs[,1])
      ###  #Termination condition
      ###  dQ = (Qnew-Q)
      ###  Q = Qnew
      ###  if(verbose>1)
      ###    message(sprintf("  %03d: rho=%g, tau=%g, dQ = %g",i,rho,tau,dQ))
      ###  #Check termination condition
      ###  i = i+1
      ###  if(abs(dQ)<1e-6 | i==5000)
      ###    break
      ###}
      #Full version for arbitrary number of mixtures without constraint
      k=2 #Number of mixtures.
      n=length(tC)
      rhos = runif(k)
      taus = rep(1/k,k)
      i=0
      Q=-Inf
      while(TRUE){
        #EStep
        ####Linear space
        ###states = do.call(cbind,lapply(rhos,function(e) dbinom(aC,tC,e)))
        ###states = t(t(states)*taus)
        ###stateProbs = states/rowSums(states)
        #Do in log space for stability
        ll = matrix(dbinom(rep(aC,k),rep(tC,k),rep(rhos,each=n),log=TRUE),nrow=n,ncol=k)
        states = t(t(ll)+log(taus))
        #states = do.call(cbind,lapply(rhos,function(e) dbinom(aC,tC,e,log=TRUE)))
        #states = t(t(states)+log(taus))
        #Another numerical stability fix.  Subtraction in log space equivalent to multiplication by a constant which changes nothing cebause normalisation
        stateProbs = exp(states-apply(states,1,max))
        stateProbs = stateProbs/rowSums(stateProbs)
        #Mstep
        Qnew = sum(stateProbs*(rep(log(taus),each=n)+ll))
        #Qnew = sum(unlist(lapply(seq(k),function(e) sum(stateProbs[,e]*(log(taus[e]) + dbinom(aC,tC,rhos[e],log=TRUE))))))
        #This is the analytic solution for the minimum of Q with respect to rho
        rhos = colSums(stateProbs*aC)/colSums(stateProbs*tC)
        #Likewise tau
        taus = colMeans(stateProbs)
        dQ = (Qnew-Q)
        Q = Qnew
        if(verbose>1)
          message(sprintf("  %03d: rho=%g, tau=%g, dQ = %g",i,rhos[1],taus[1],dQ))
        #Check termination condition
        i = i+1
        if(abs(dQ)<1e-6 | i==5000)
          break
      }
      #Organise them so that the first one is the lower one
      o = order(rhos)
      rhos = rhos[o]
      taus = taus[o]
      stateProbs = stateProbs[,o]
      #Plot if we need to
      if(plotMixtures){
        mixCols = c('#FF000080','#0000FF80')
        breaks = seq(0,1,length.out=100)
        a = hist((aC/tC)[stateProbs[,1]>0.5],breaks=breaks,plot=FALSE)
        b = hist((aC/tC)[stateProbs[,1]<0.5],breaks=breaks,plot=FALSE)
        yMax = max(c(a$counts,b$counts))
        plot(a,
             border=FALSE,
             xlim=c(0,1),
             col=mixCols[1],
             xlab='BAF',
             main = sprintf('seg %s, mixes %.02f (%.01f%%) & %.02f (%.01f%%)',names(segs[ii]),rhos[1],100*taus[1],rhos[2],100*taus[2]),
             ylim=c(0,yMax)
             )
        plot(b,
             border=FALSE,
             xlim=c(0,1),
             col=mixCols[2],
             add=TRUE
             )
        abline(v=rhos,col=mixCols,lty=2)
      }
      #Do sanity checks that CN change is detected
      if(abs(sum(rhos)-1)>0.2 || max(taus)>0.65){
        stop(sprintf("Mixture solutions for CN segment %s inconsistent with CN change.  Found mixture with means %.02f (%.01f%%) and %.02f (%.01f%%)",names(segs[ii]),rhos[1],100*taus[1],rhos[2],100*taus[2]))
      }else if(verbose){
        message("  Tumour DNA consistent with CN change at this segment.")
      }
      mm = match(tmp,hSNPs)
      #Maternal arbitrarily defined as major allele
      hSNPs$matProb[mm] = stateProbs[,2]
    }
    hSNPs$informative = !is.na(hSNPs$matProb) & (0.5-abs(hSNPs$matProb-0.5)) < FDR
  }else{
    #Simple p-value based model
    pVals = pbinom(pmax(tCnts$altCount,tCnts$refCount)-1,tCnts$Tot,0.5,lower.tail=FALSE)
    qVals = p.adjust(pVals,method='BH')
    hSNPs$matProb[m] = pVals
    hSNPs$informative[m] = qVals<FDR
  }
  #Do a sanity check that we can detect both alleles.  A location is heterozygous if there are enough counts for both the ref and alt to reject the null hypothesis of random errors.
  if(errRate>0){
    pVals = pbinom(pmin(hSNPs$altCountTum,hSNPs$refCountTum)-1,hSNPs$totCountTum,errRate,lower.tail=FALSE)
    qVals = p.adjust(pVals,method='BH')
    hSNPs$passSanity[which(qVals>=FDR)] = FALSE
  }
  hSNPs$altIsMum = ifelse(hSNPs$informative,
                          hSNPs$altCountTum > hSNPs$refCountTum, #Designate the major allele as "Mum"
                          NA)
  #Sanity check each segment. 
  if(verbose){
    message('###########\n# Summary #\n###########')
    message(sprintf('Of %s heterozygous SNPs:',prettyNum(length(hSNPs),big.mark=',')))
  }
  o = findOverlaps(hSNPs,segs)
  for(i in seq(0,length(segs))){
    #How many are there
    qH = queryHits(o)[subjectHits(o)==i]
    lab = names(segs[i])
    if(i==0){
      qH = unique(queryHits(o))
      lab = 'global'
    }
    nSNPs = length(qH)
    nPhased = sum(hSNPs$informative[qH],na.rm=TRUE)
    if(verbose){
      message(sprintf('  %s (%.01f%%) fall in CN segment %s, of which:',prettyNum(nSNPs,big.mark=','),nSNPs/length(hSNPs)*100,lab))
      message(sprintf('    %s (%.01f%%) could be phased, of which:',prettyNum(nPhased,big.mark=','),nPhased/nSNPs*100))
      nSane =sum(hSNPs$passSanity[qH] & hSNPs$informative[qH],na.rm=TRUE)
      message(sprintf('      %s (%.01f%%) passed sanity checks, of which:',prettyNum(nSane,big.mark=','),nSane/nPhased*100))
      nAltMum = sum(hSNPs$altIsMum[qH] & hSNPs$passSanity[qH] & hSNPs$informative[qH],na.rm=TRUE)
      message(sprintf('        %s (%.01f%%) have the maternal allele as ALT',prettyNum(nAltMum,big.mark=','),nAltMum/nSane*100))
    }
    if(i!=0 && nPhased/nSNPs < minPhasable){
      hSNPs$passSanity[qH]=FALSE
      hSNPs$altIsMum[qH]=NA
      warning(sprintf("Only %s of %s SNPs (%.02f%%) could be phased in copy number segment %s.  This segment may not contain a real CN change.",prettyNum(nPhased,big.mark=','),prettyNum(nSNPs,big.mark=','),nPhased/nSNPs*100,names(segs[i])))
    }
  }
  #Explicitly record which allele is maternal (major allele)/paternal (minor allele)
  hSNPs$matAllele = ifelse(hSNPs$altIsMum,hSNPs$ALT,hSNPs$REF)
  hSNPs$patAllele = ifelse(hSNPs$altIsMum,hSNPs$REF,hSNPs$ALT)
  #Add extra alias for clarity
  hSNPs$majorAllele = hSNPs$matAllele
  hSNPs$minorAllele = hSNPs$patAllele
  #Store segments
  hSNPs@metadata$segs = segs
  return(hSNPs)
}
