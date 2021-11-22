#' Generate coverage (logR) and BAF from BAM file
#'
#' The most "raw" way to look at copy number status of a BAM file is to calculate the coverage and alternate allele frequency at heterozygous SNPs.  This function uses alleleCounter to generate this from a BAM file, for a set of heterozygous SNPs.  You can specify the SNPs however you like, but usually this would either be using a pre-specified panel of SNPs, or using the \code{\link{findHetSNPs}} function to identify potential heterozygous SNPs directly.
#'
#' Will optionally run ASCAT to infer CN segments (in unmatched mode) and produces a basic plot to visualise the raw data.  If ASCAT is run, there is some extra information that is needed.  These values almost all relate to guessing if a given SNP is heterozygous in the germline.  The allele frequency \code{AF} is used as a prior and posterior probabilities of heterozygosity and homozygosity (in the germline) is determined using the provided error rate (for homozygosity) and maximal abbarent cell fraction (\code{abbCellFrac} for heterozygosity).  Any SNP with posterior probability of homozygosity greater than \code{minHomProb} or observed BAF less than \code{minHetAF} will be coded as homozygous.  This errs on the side of marking things as homozygous if the data is not clear, which is the correct thing to do.
#'
#' @inheritParams phaseSNPsFromCN
#' @param BAM The BAM file to calculate coverage and BAF on.
#' @param hSNPs The heterozygous SNPs to investigate.  Defaults to 1000 genomes SNPs.  Must be a GRanges object with REF and ALT columns.
#' @param AF The population allele frequency of the alternate allele.  If NULL, defaults to 0.05
#' @param gcFile ASCAT GC correction file for the SNPs provided.  If given, used to GC correct logR.
#' @param chrs Which chromosomes to plot and/or run ASCAT on.
#' @param errRate Assumed sequencing error rate.
#' @param abbCellFrac Assumed fraction of abberant cells for determining heterozygous SNPs.  Higher values are more stringent.
#' @param minHomProb Any SNP with posterior probability of homozygosity in the germline greater than this will be code as homozygous.
#' @param minHetAF Any SNP with BAF below this will be coded as homozygous.
#' @param runASCAT Should we run ASCAT?
#' @return The hSNPs object with extra columns, particularly total coverage and BAF.  If ASCAT has been run, columns associated with the fit are also added and the full ASCAT output is stored at \code{hSNPs@metadata$ascat.output}.
#' @importFrom stats median
#' @export
generateCoverageAndBAF = function(BAM,refGenome,hSNPs=SNPs1k,AF=hSNPs$AF,outPath=NULL,gcFile=NULL,alleleCounterParams=list(f=3,F=3852,m=20,q=35),chrs=c(1:22,'X'),errRate=0.01,abbCellFrac=0.95,minHomProb=0.05,minHetAF=0.1,runASCAT=FALSE,nParallel=1){
  ###############################
  # Preliminary checks and setup
  if(is.null(AF))
    AF=0.05
  hSNPs$AF = AF
  ##############################
  # Run alleleCounter as needed
  #Set the non-optional allelecounter parameters
  alleleCounterParams$x=FALSE
  alleleCounterParams$bams = BAM
  alleleCounterParams$tgtLoci = hSNPs
  alleleCounterParams$refGenome = refGenome
  alleleCounterParams$outputs = outPath
  alleleCounterParams$nParallel = nParallel
  rawCnts = do.call(alleleCounter,alleleCounterParams)[[1]]
  ###################################
  # Extract the values we care about
  rawCnts$refCnts = as.matrix(mcols(rawCnts)[,c('A','C','G','T')])[cbind(seq_along(rawCnts$REF),match(rawCnts$REF,c('A','C','G','T')))]
  rawCnts$altCnts = as.matrix(mcols(rawCnts)[,c('A','C','G','T')])[cbind(seq_along(rawCnts$ALT),match(rawCnts$ALT,c('A','C','G','T')))]
  rawCnts$sourceBAM = BAM
  rawCnts$coverage = rawCnts$refCnts + rawCnts$altCnts
  rawCnts$BAF = rawCnts$altCnts/rawCnts$coverage
  chrs = intersect(chrs,unique(as.character(seqnames(rawCnts[rawCnts$coverage>0]))))
  #Determine which ones are likely heterozygous in the germline
  #The model is that we assume an error rate, that by default is conservatively chosen to be a high value of 0.01 (see this paper https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1659-6).  Then we need to specify what the highest purity we expect is.  Obviously if it's pure tumour, then LOH will result in complete loss of one allele and we cannot tell which are heterozygous and which are homozygous.  The default value is 95% purity and the model assumes a loss of one allele in the aberant cell (i.e. haploid).
  #We also incoporate the allele frequency, if known, as a prior.  Otherwise a uniform prior of 5% allele frequency is used.
  #Convert everything to the less than 0.5 case
  pp = ifelse(rawCnts$AF>0.5,1-rawCnts$AF,rawCnts$AF)
  x = ifelse(rawCnts$altCnts>rawCnts$altCnts/rawCnts$coverage,rawCnts$coverage-rawCnts$altCnts,rawCnts$altCnts)
  pHom = (1-pp*(1-pp))*dbinom(x,rawCnts$coverage,errRate)
  pHet = pp*(1-pp)*dbinom(x,rawCnts$coverage,(1-abbCellFrac)/(2-abbCellFrac))
  #Workaround for numerical issues
  pHom[is.na(pHom)] = 0
  pHet[is.na(pHet)] = 0
  normFac = pHet + pHom
  pHet = pHet/normFac
  pHom = pHom/normFac
  rawCnts$isHet = !(pHom>minHomProb | abs(rawCnts$BAF-0.5)>(0.5-minHetAF))
  #Now calculate logR.  Because the definition is normalised to the median, dividing by any constant here will change nothing.
  rawCnts$logR = log(rawCnts$coverage) - median(log(rawCnts$coverage))
  #Run ascat if required
  if(runASCAT){
    message("Fitting CN states using unmatched ASCAT")
    #Make temporary object to manipulate, dropping zeros
    tmp = rawCnts[rawCnts$coverage>0,]
    nom = paste0('BAM',1)
    #Drop anything that is not a target chromosome
    chrsLocal = intersect(chrs,unique(as.character(seqnames(tmp))))
    tmp = tmp[as.character(seqnames(tmp)) %in% chrsLocal,]
    #Reorder as required
    tmp = tmp[order(factor(as.character(seqnames(tmp)),levels=chrsLocal),start(tmp)),]
    #Make position object
    SNPpos = data.frame(row.names=as.character(tmp),chrs=as.character(seqnames(tmp)),pos=start(tmp))
    ch = split(seq_along(tmp),as.character(seqnames(tmp)))[chrs]
    names(ch) = NULL
    #Construct in ASCAT format
    inDat = list(Tumor_LogR = matrix(tmp$logR,ncol=1,dimnames=list(as.character(tmp),nom)),
                 Tumor_BAF = matrix(tmp$BAF,ncol=1,dimnames=list(as.character(tmp),nom)), 
                 Tumor_LogR_segmented = NULL, 
                 Tumor_BAF_segmented = NULL, 
                 Germline_LogR = NULL, 
                 Germline_BAF = NULL, 
                 SNPpos = SNPpos,
                 ch = ch,
                 chr = ASCAT:::split_genome(SNPpos),
                 chrs = chrsLocal,
                 samples = nom,
                 gender = 'XX',
                 sexchromosomes = 'XX',
                 failedarrays = NULL
                 )
    #Do we have a GC content file.  If so, use it to correct
    if(!is.null(gcFile))
      inDat = ASCAT::ascat.GCcorrect(inDat,gcFile)
    #ascat.plotRawData(inDat,'~')
    isHom = list(germlinegenotypes = matrix(!tmp$isHet,ncol=1,dimnames=dimnames(inDat$Tumor_BAF)),failedarrays=NULL)
    inDat = ASCAT::ascat.aspcf(inDat, ascat.gg=isHom)
    #ascat.plotSegmentedData(inDat,'~')
    ascat.output = ASCAT::ascat.runAscat(inDat,gamma=1)
    #ascat.plotSunrise(ascat.output$distance_matrix[[1]],ascat.output$ploidy,ascat.output$aberrantcellfraction)
    #Include 
    rawCnts@metadata$ascat.output = ascat.output
    rawCnts$nA = NA
    rawCnts$nB = NA
    m = match(rownames(ascat.output$nA),as.character(rawCnts))
    rawCnts$nA[m] = ascat.output$nA[,1]
    rawCnts$nB[m] = ascat.output$nB[,1]
    if(!is.null(gcFile)){
      rawCnts$correctedLogR = NA
      rawCnts$correctedLogR[m] = inDat$Tumor_LogR[,1]
    }
  }
  return(rawCnts)
}
