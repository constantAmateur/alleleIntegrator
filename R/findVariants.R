#' Use BCFtools to find variants in a BAM file
#'
#' Find any base with decent coverage and BAF within a range and save them a VCF. If an output filename is given via \code{outVCF}, an earlier result will be loaded if it exists.
#'
#' This is mostly just a wrapper around external calls to bcftools, which obviously need to be installed for this to work.
#'
#' @param bam Path to BAM file to call variants from.
#' @param refGenome Path to reference genome used to map BAM.
#' @param minCoverage Any location with coverage below this will not be reported.
#' @param BAF_lim Any SNP with BAF outside this range will be discarded.  Vector of length tow, e.g. \code{BAF_lim=c(0.2,0.8)}.
#' @param outVCF The path of the VCF file containing SNPs to be created.
#' @param chrsToProcess Which chromosomes to process?  If NULL, all chromosomes in the BAM header are used.  Any chromosomes not present in the BAM header will be dropped.  If none found, tries prepending 'chr' to this before failing.
#' @param minMapQual Minimum map quality to pass to -q flag of bcftools mpileup.
#' @param minBaseQual Minimum base quality to pass to -Q flag of bcftools mpileup.
#' @param minVarQual Filter variants with quality scores below this value.
#' @param bin Path to BCF tools binary.
#' @param skipIfExists Skip running the SNP caller if the VCF already exists.
#' @param nParallel How many threads to use.  Note that the final variants called depends (very mildly) on the genomic chunks processed together.  As such, it is not advisable to set this much larger than the number of chromosomes.
#' @param nChunks When running in parallel, split into this many chunks per parallel thread.
#' @param ... Absorb unused parameters.
#' @return A GRanges object containing all variants matching specified parameters.
#' @import GenomicRanges
#' @importFrom Rsamtools BamFile
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom VariantAnnotation readVcf info
#' @importFrom DelayedArray rowRanges
#' @importFrom utils write.table
#' @export
findVariants = function(bam,refGenome,minCoverage,BAF_lim,outVCF=NULL,chrsToProcess=c('X',1:22),minMapQual=35,minBaseQual=20,minVarQual=225,bin='bcftools',skipIfExists=TRUE,nParallel=1,nChunks=1,...){
  if(is.null(outVCF))
    warning("output file not specified, so results cannot be reused without recalculation.  These calculations are time consuming, so it is advisable to save them somewhere.")
  #Check BAF_lim set sensibly
  if(any(BAF_lim<0 | BAF_lim>1) || length(BAF_lim)!=2 || BAF_lim[1]>BAF_lim[2])
    stop("Invalid BAF_lim.")
  #Get list of chromosomes to process
  chrLens = seqinfo(BamFile(bam))
  chrsBAM = seqnames(chrLens)
  if(!is.null(chrsToProcess))
    chrs = chrsToProcess[chrsToProcess %in% chrsBAM]
  if(length(chrs)==0){
    #Check if we can fix by stripping 'chr'
    chrsToProcess = gsub('^chr','',chrsToProcess)
    if(!is.null(chrsToProcess))
      chrs = chrsToProcess[chrsToProcess %in% chrsBAM]
    #Or adding it in
    chrsToProcess = paste0('chr',chrsToProcess)
    if(!is.null(chrsToProcess))
      chrs = chrsToProcess[chrsToProcess %in% chrsBAM]
    #Still no good :(
    if(length(chrs)==0)
      stop("No chromosomes found to process!")
  }
  #Break it up into pieces
  chrLens = as.data.frame(chrLens[chrs])
  chrLens$cumLen = cumsum(as.numeric(chrLens$seqlengths))
  nChunksTot = ifelse(nParallel==1,1,nParallel*nChunks)
  breaks = seq(0,max(chrLens$cumLen),length.out=nChunksTot+1)
  #Construct GRanges for each
  chunks = lapply(seq(nChunksTot),function(e) {
                    start = breaks[e]+1
                    stop = breaks[e+1]
                    startChr = min(which(start<=chrLens$cumLen))
                    stopChr = min(which(stop<=chrLens$cumLen))
                    start = start - c(0,chrLens$cumLen)[startChr]
                    stop = stop - c(0,chrLens$cumLen)[stopChr]
                    if(startChr==stopChr){
                      GRanges(rownames(chrLens)[startChr],IRanges(start,stop))
                    }else{
                      nBits = stopChr-startChr+1
                      GRanges(rownames(chrLens)[seq(startChr,stopChr)],IRanges(start=c(start,rep(1,nBits-1)),end=c(chrLens$seqlengths[seq(startChr,stopChr-1)],stop)))
                    }})
  chunkStr = sapply(chunks,function(e) paste(as.character(e),collapse=','))
  #Define base command.  Need to split in two for this to work downstream.  Quality thresholds are set to alleleCounter defaults
  baseCommand = sprintf("%s mpileup -d 10000 -q %d -Q %d -Ou -f %s -r %%s %s | %s call -cv -V indels | %s filter -Oz -o %%s -e ",bin,minMapQual,minBaseQual,refGenome,bam,bin,bin)
  baseFilter = sprintf("'%%SUM(DP4)<%d || (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) < %g || (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > %g'",minCoverage,BAF_lim[1],BAF_lim[2])
  #Check if we need to actually run anything
  if(!skipIfExists || is.null(outVCF) || !file.exists(outVCF)){
    #Processs one chunk at a time
    tFiles = rep(NA,nChunksTot)
    cmds = rep(NA,nChunksTot)
    for(i in seq_along(chunkStr)){
      tFiles[i] = tempfile()
      cmds[i] = paste0(sprintf(baseCommand,chunkStr[i],tFiles[i]),baseFilter)
    }
    #Actually run things
    rFile = tempfile()
    write.table(cmds,rFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
    #Actually run things
    system(sprintf("parallel -j %d < %s",nParallel,rFile))
    #All done, concatenate to local home.
    output = tempfile()
    system(sprintf('%s concat -o %s %s',bin,output,paste(tFiles,collapse = ' ')))
    #Move to final home if this exists
    if(!is.null(outVCF))
      file.copy(output,outVCF,overwrite=TRUE)
    #Read it in
    tmp = readVcf(outVCF)
    #Cleanup
    unlink(c(tFiles,rFile,output))
  }else{
    #Just read it in
    tmp = readVcf(outVCF)
  }
  #Process into final format
  vcf = rowRanges(tmp)
  cnts = as.data.frame(as.matrix(info(tmp)$DP4))
  colnames(cnts) = c('fREF','rREF','fALT','rALT')
  mcols(vcf) = cbind(mcols(vcf)[,c('REF','ALT','QUAL')],cnts)
  vcf$refCnt = rowSums(cnts[,1:2])
  vcf$altCnt = rowSums(cnts[,3:4])
  vcf$total = rowSums(cnts)
  vcf = vcf[lengths(vcf$ALT)==1]
  #Convert REF/ALT to character
  vcf$REF = as.character(vcf$REF)
  vcf$ALT = as.character(unlist(vcf$ALT))
  #Drop anything that's not length 1
  vcf = vcf[nchar(vcf$REF)==1 & nchar(vcf$ALT)==1]
  #Or anything with NAs
  vcf = vcf[!is.na(vcf$total),]
  #Or anything below quality threshold
  vcf = vcf[vcf$QUAL >= minVarQual,]
  return(vcf)
}
