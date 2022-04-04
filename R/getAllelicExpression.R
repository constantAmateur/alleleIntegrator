#' Get allele specific expression from RNA BAM files at base pair resolution
#'
#' Given a set of individual nucleotide genomic coordinates, interrogates all locations in a set of RNA bam files and generates a cell specific allele count.  The output is \code{loci}, with extra columns added recording the number of counts for A,C,G,T, the refence/alternative alleles (as encoded by REF/ALT) and the maternal/paternal allele (if known).  The labels for each are determined by the \code{labels} param and derived from the BAM file name by default.  \code{loci} must contain columns named \code{REF} and \code{ALT} specifying the REF/ALT allele.
#'
#' Each entry is given a regionID, which is given in the format chr:pos_refAllele/altAllele with the maternal allele (where known) coded as a lower case letter.
#'
#' @param loci GRanges object containing heterozygous SNP locations to measure.  Must have REF and ALT.
#' @param refGenome Reference genome used to align BAM files.
#' @param bams Aligned BAM files
#' @param labels Unique label for each BAM.  Name of \code{bams} used if available.  If not, generated automatically.
#' @param outputs File names to store output for each BAM file in. See \code{\link{alleleCounter}}.
#' @param minCounts Any SNP that doesn't have at least this many counts from both alleles across all cells is marked as \code{passSanity=FALSE}.
#' @param assayType The type of expression assay run.  Used to determine the default parameters passed to \code{\link{alleleCounter}}.  These can always be over-riden by directly specifying allelelCounter params yourself.
#' @param matAlleleCol Column in \code{loci} indicating which base is maternally derived.
#' @param patAlleleCol Column in \code{loci} indicating which base is maternally derived.
#' @param segs CN segments.  Used for summary stats if present and \code{verbose=2}.
#' @param errRate Assumed sequencing error rate.  Only used for summary stats.
#' @param verbose Be verbose?  0/FALSE reports nothing, 1/TRUE reports a minimal summary, 2 for a detailed summary. 
#' @param nParallel Number of threads to use.
#' @param ... Passed to alleleCounter.
#' @return The input GRanges object \code{loci}, but expanded to have one line per cell in the expression data and with columns added to give the number of counts in that cell, at that location in the expression data.
#' @export
getAllelicExpression = function(loci,refGenome,bams,labels=names(bams),outputs=NULL,minCounts=0,assayType=c('10X','SS2','Other'),matAlleleCol='matAllele',patAlleleCol='patAllele',segs=loci@metadata$segs,errRate=0.02,verbose=TRUE,nParallel=1,...){
  assayType = match.arg(assayType)
  if(is.null(labels) || any(duplicated(labels))){
    warning("No valid BAM labels found.  Setting to generic Sample1, Sample2, Sample3, etc.")
    labels = paste0('Sample',seq_along(bams))
  }
  params = list(bams=bams,
                refGenome=refGenome,
                tgtLoci=loci,
                outputs=outputs,
                nParallel=nParallel)
  if(assayType=='SS2'){
    params$x=FALSE
    params$f=3
    params$F=3852
    params$m=20
    params$q=35
  }
  theDots = list(...)
  for(nom in names(theDots))
    params[[nom]] = theDots[[nom]]
  cnts = do.call(alleleCounter,params)
  #Construct the giant GRanges of everything
  out=list()
  for(i in seq_along(bams)){
    lab = labels[i]
    dat = cnts[[i]]
    #If it's been run in 10X mode, need to construct barcodes
    if(assayType=='10X'){
      dat$cellID = paste0(lab,'_',dat$barcode)
    }else{
      dat$cellID = lab
    }
    #Record sample
    dat$scSource = bams[i]
    dat$sample = lab
    out[[i]]=dat
  }
  out = do.call(c,out)
  #Calculate counts for different summaries of data
  bases=c('A','C','G','T')
  vars = c(bases,'Tot','altCount','refCount')
  tmp = as.matrix(mcols(out)[,bases])
  out$altCount = tmp[cbind(seq(length(out)),match(out$ALT,bases))]
  out$refCount = tmp[cbind(seq(length(out)),match(out$REF,bases))]
  if(matAlleleCol %in% colnames(mcols(out))){
    out$matCount = tmp[cbind(seq(length(out)),match(mcols(out)[,matAlleleCol],bases))]
    vars = c(vars,'matCount')
  }else{
    out$matCount = as.numeric(NA)
  }
  if(patAlleleCol %in% colnames(mcols(out))){
    out$patCount = tmp[cbind(seq(length(out)),match(mcols(out)[,patAlleleCol],bases))]
    vars = c(vars,'patCount')
  }else{
    out$patCount = as.numeric(NA)
  }
  #Sort it by position, then sample
  out = out[order(out,out$sample)]
  #Make labels for each SNP.  Format is chr:pos_REF/ALT with maternal allele coded as lower case.
  if(matAlleleCol %in% colnames(mcols(out))){
    tmp = ifelse(is.na(mcols(out)[,matAlleleCol]),
                 paste(out$REF,out$ALT,sep='/'),
                 ifelse(mcols(out)[,matAlleleCol] == out$REF,
                        paste(tolower(out$REF),out$ALT,sep='/'),
                        paste(out$REF,tolower(out$ALT),sep='/')
                        )
                 )
    out$regionID = paste0(as.character(out),'_',tmp)
  }else{
    out$regionID = paste0(as.character(out),'_',out$REF,'/',out$ALT)
  }
  names(out) = out$regionID
  #Do the sanity check filter
  if(minCounts>0){
    tmp = buildCountMatricies(out,assays=c('refCount','altCount'))
    failIDs = rownames(tmp$refCount)[rowSums(tmp$refCount)<minCounts | rowSums(tmp$altCount)<minCounts]
    if(length(failIDs)>0)
      out$passSanity[which(out$regionID %in% failIDs)]=FALSE
  }
  #Summary stats about fit
  #Do summary stats by segment?  Need mat/pat columns if so
  if(length(segs)>0 && !all(is.na(out$matCount)) && !all(is.na(out$patCount))){
    for(i in seq(0,length(segs))){
      if(i==0){
        src = loci
        dst = out
        lab = 'Global'
      }else{
        src = subsetByOverlaps(loci,segs[i])
        dst = subsetByOverlaps(out,segs[i])
        lab =names(segs)[i]
      }
      nSNPs = sum(src$informative & src$passSanity)
      if(verbose)
        message(sprintf('There are %s informative SNPs in segment %s, of which:',prettyNum(nSNPs,big.mark=','),lab))
      #Collapse and discard cell specific info
      gCnts = dst[(dst$matCount + dst$patCount)>0 & dst$informative & dst$passSanity,]
      cCnts = aggregateByClusters(gCnts,rep(1,length(gCnts)),assays=c('matCount','patCount'))
      cCnts$totCount = cCnts$matCount + cCnts$patCount
      nCov = length(cCnts)
      if(verbose)
        message(sprintf('  %s (%.01f%%) have some coverage, of which:',prettyNum(nCov,big.mark=','),nCov/nSNPs*100))
      #Coverage summary
      tmp = quantile(cCnts$totCount,seq(0,1,.1))
      if(verbose)
        message(sprintf('    Coverage quantiles are: %s',paste(tmp,' (',names(tmp),')',sep='',collapse = ' ')))
      #Breakdown by type
      types = unique(cCnts$regionType)
      for(type in types){
        nType = sum(cCnts$regionType==type)
        if(verbose)
          message(sprintf('    %s (%.01f%%) are in %s',prettyNum(nType,big.mark=','),nType/nCov*100,type))
      }
      #Can we detect both alleles
      pVals = pbinom(pmin(cCnts$matCount,cCnts$patCount)-1,cCnts$totCount,errRate,lower.tail=FALSE)
      if(verbose)
        message(sprintf('    %s (%.01f%%) look heterozygous at 0.05 p-value cut-off.',prettyNum(sum(pVals<0.05),big.mark=','),sum(pVals<0.05)/nCov*100))
      w = which(cCnts$totCount>10)
      nHet = sum(pVals[w]<0.05)
      if(verbose){
        message(sprintf('    %s (%.01f%%) have high coverage (>10 UMIs), of which:',prettyNum(length(w),big.mark=,','),length(w)/nCov*100))
        message(sprintf('      %s (%.01f%%) look heterozygous at 0.05 p-value cut-off.',prettyNum(sum(pVals[w]<0.05),big.mark=','),sum(pVals[w]<0.05)/length(w)*100))
      }
      #Do sanity check that thing match up
      if(i==0 && nHet/length(w)<0.5){
        warning(sprintf("A very low fraction of high coverage SNPs have both alleles detected %s (%.01f%%).  This likely indicates a sample or genome mismatch between single cell and DNA data.",prettyNum(nHet,big.mark=','),nHet/length(w)*100))
      }
      #Only progress to segment specific summary if very verbose
      if(verbose<=1)
        break
    }
  }
  return(out)
}
