#' Phases heterozygous SNPs in a patient using parents genotype
#'
#' Given heterozygous SNPs in an individual, as called using \code{\link{findHetSNPs}}, checks genotype in Maternal and/or Paternal BAM files at these locations.  This information is used to identify which allele is maternal and paternal at each SNP.
#'
#' To be phased, one of the alleles must be completely abscent in one of the parents.  This is determined by a binomial test against heterozygosity with FDR of \code{FDR}, plus a requirement that the allele frequency not exceed \code{nullAF}.
#'
#' Where both parents are available, additionaly sanity checks are performed to make sure the trio is consistent.  This tests for presence of one allele in each parent, as measured by the allele frequency \code{notNullAF}.
#'
#' Note that all three BAMs (the one used to generate \code{hSNPs}, \code{mBAM}, and \code{pBAM}) must have been mapped to the same reference genome, \code{refGenome}.
#'
#' @param hSNPs The heterozygous SNP locations.  See \code{\link{findHetSNPs}}.
#' @param refGenome Reference genome used to map \code{mBAM} and \code{pBAM}.
#' @param mBAM BAM file containing maternal genotype.
#' @param pBAM BAM file containing paternal genotype.
#' @param outBase Base name used to construct files storing genotyping calculations.  If left at NULL, results are not saved and cannot be reloaded in subsequent runs.
#' @param FDR False discovery rate for declaring an allele not present in the parental data.
#' @param nullAF As well as passing the FDR criteria, an allele must have frequency less than this to be declared abscent in the parent.
#' @param notNullAF As well as passing the FDR criteria, an allele must have frequency greater than this to be declared present in the parent.
#' @param alleleCounterParams Parameters to be used when running allele counter.
#' @param nParallel How many threads.
#' @param ... Does nothing.
#' @return \code{hSNPs} but with extra columns indicating results of phasing.
#' @export
phaseSNPsFromParents = function(hSNPs,refGenome,mBAM=NULL,pBAM=NULL,outBase=NULL,FDR=0.05,nullAF=0.05,notNullAF=0.3,alleleCounterParams=list(f=3,F=3852,m=20,q=35),nParallel=1,...){
  mum = !is.null(mBAM)
  dad = !is.null(pBAM)
  if(!mum & !dad)
    stop("At least one parent must be given")
  if(is.null(outBase))
    warning("outBase not specified.  Results of genotyping calculations not saved to disk and cannot be automatically reloaded on subsequent runs.")
  ##############################
  # Run alleleCounter as needed
  if(!is.null(outBase)){
    outputs = paste(outBase,c('maternal','paternal')[c(mum,dad)],'countAtHetSNPs.tsv',sep='_')
  }else{
    outputs = NULL
  }
  #Set the non-optional allelecounter parameters
  alleleCounterParams$x=FALSE
  alleleCounterParams$bams = c(mBAM,pBAM)
  alleleCounterParams$tgtLoci = hSNPs
  alleleCounterParams$refGenome = refGenome
  alleleCounterParams$outputs = outputs
  alleleCounterParams$nParallel = nParallel
  pCnts = do.call(alleleCounter,alleleCounterParams)
  ####################
  # Determine phasing
  bases=c('A','C','G','T')
  #Can we phase this SNP?
  informative = rep(FALSE,length(hSNPs))
  #Is the alternate mum?
  altIsMum = rep(NA,length(hSNPs))
  #Are the data consistent?
  sanityFail = rep(FALSE,length(hSNPs))
  for(i in seq(2)){
    if(c(mum,dad)[i]){
      #Load count object
      cnts =  mcols(pCnts[[min(i,length(pCnts))]])
      cnts$BAF = (as.matrix(cnts[,bases])[cbind(seq(nrow(cnts)),match(cnts$ALT,bases))])/cnts$Tot
      cnts$AAF = (as.matrix(cnts[,bases])[cbind(seq(nrow(cnts)),match(cnts$REF,bases))])/cnts$Tot
      #Save.
      tmp = cnts[,c(bases,'Tot','BAF','AAF')]
      colnames(tmp) = paste0(c('maternal','paternal')[i],'_',colnames(tmp))
      mcols(hSNPs) = cbind(mcols(hSNPs),tmp)
      #Work out if there's evidence of heterozygous absence
      pValNoRef = ifelse(cnts$Tot==0,1,pbinom(cnts$AAF*cnts$Tot,cnts$Tot,0.5))
      pValNoAlt = ifelse(cnts$Tot==0,1,pbinom(cnts$BAF*cnts$Tot,cnts$Tot,0.5))
      noRef = p.adjust(pValNoRef,method='BH') < FDR & cnts$AAF< nullAF
      noAlt = p.adjust(pValNoAlt,method='BH') < FDR & cnts$BAF< nullAF
      #If it's missing both we're in trouble
      sanityFail[noRef & noAlt] = TRUE
      #It's informative if either REF or ALT is completely abscent in Mum or Dad
      informative = informative | noRef | noAlt
      #If there's no REF allele in Mum, the alt must be from dad.
      w = is.na(altIsMum)
      #Update where we have both bits of data first
      if(any(!w)){
        #What would we set it to if it wasn't already set?
        ww = which(!w)
        tmp = rep(NA,length(ww))
        tmp[noRef[ww]] = c(TRUE,FALSE)[i]
        tmp[noAlt[ww]] = c(FALSE,TRUE)[i]
        #Now check for consistency and NA those that are inconsistent (others already set)
        www = ww[which(tmp!=altIsMum[ww])]
        sanityFail[www]=TRUE
        altIsMum[www] = NA
      }
      #Now set those not already set
      altIsMum[which(w & noRef)] = c(TRUE,FALSE)[i]
      altIsMum[which(w & noAlt)] = c(FALSE,TRUE)[i]
    }
  }
  #Do extra checks if we have both
  if(mum & dad){
    #Make sure that both REF and ALT are present
    hasBoth = (hSNPs$maternal_BAF > notNullAF & hSNPs$paternal_AAF > notNullAF) | (hSNPs$maternal_AAF > notNullAF & hSNPs$paternal_BAF > notNullAF)
    #Where data is missing from one, it's equivalent to not being able to check
    hasBoth[is.na(hasBoth)]=TRUE
    #Where both can't be found.  Consider this a bad SNP.
    sanityFail[which(!hasBoth)] = TRUE
  }
  #Now save the information and the phasing
  hSNPs$informative = informative
  hSNPs$passSanity = !sanityFail
  hSNPs$altIsMum = NA
  #Exclude those that fail the sanity checks
  hSNPs$altIsMum[which(informative & !sanityFail)] = altIsMum[which(informative & !sanityFail)]
  #Explicitly record which allele is maternal/paternal
  hSNPs$matAllele = ifelse(hSNPs$altIsMum,hSNPs$ALT,hSNPs$REF)
  hSNPs$patAllele = ifelse(hSNPs$altIsMum,hSNPs$REF,hSNPs$ALT)
  return(hSNPs)
}
 
