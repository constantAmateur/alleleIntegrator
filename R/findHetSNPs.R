#' Find heterozygous single nucleotide polymorphisms
#'
#' Uses bcftools to find heterozygous SNPs.  See \code{\link{findVariants}}.
#'
#' All variants with reasonable coverage (\code{minCoverage}) and BAF consistent with a het SNP (\code{BAF_lim}) are identified.  These are then further filtered and any with a significant deviation from 0.5 BAF are removed.  Significant deviation is defined by binomial FDR less than \code{FDR} and 0.5 - \code{minDeviation} < BAF < 0.5 + \code{minDeviation}.
#'
#' The defaults assume that \code{bam} comes from normal tissue.  If instead this comes from tumour tissue (or some other type of mixed genotype tissue), you should change some parameters.  The difference with calling SNPs from tumour sequencing SNPs will not have BAF of 0.5, so the \code{minDeviation} filter is not appropriate and should be set to \code{minDeviation = 0.5}.  You may also want to adjust `BAF_lim` to be more permissive and/or adjust \code{errRate}.
#'
#' @param bam Path to BAM file to call variants from.
#' @param refGenome Path to reference genome used to map BAM.
#' @param outVCF The path of the VCF file containing SNPs to be created.
#' @param minCoverage Any location with coverage below this will not be reported.
#' @param BAF_lim Any SNP with BAF outside this range will be discarded.
#' @param FDR False discover rate to accept.
#' @param minDeviation Minimum amount that the BAF must deviate from 0.5 before we reject it.
#' @param errRate Assumed error rate for sanity check that location is heterozygous.  This check is disabled if \code{errRate=0}
#' @param nParallel How many threads to use.
#' @param ... Other things passed to \code{\link{findVariants}}
#' @return A GRanges object containing heterozygous SNPs.
#' @export
#' @importFrom stats pbinom p.adjust
findHetSNPs = function(bam,refGenome,outVCF=NULL,minCoverage=10,errRate=0.01,BAF_lim=c(0.2,0.8),FDR=0.05,minDeviation=0.1,nParallel=1,...){
  #Call SNPs
  snps = findVariants(bam,refGenome,outVCF=outVCF,BAF_lim=c(0.2,0.8),minCoverage=10,nParallel=nParallel,...)
  snps$BAF = snps$altCnt/snps$total
  #Exclude those that stray too far from the expected 50/50
  pVals = pmin(pbinom(snps$altCnt-1,snps$total,0.5,lower.tail=FALSE),
               pbinom(snps$altCnt,snps$total,0.5,lower.tail=TRUE)
               )
  qVals = p.adjust(pVals,method='BH')
  #Reject those with significant pValue, and effect size we consider reasonable.  That is, if it's within minDeviation of 50/50, we keep it no matter what.
  snps = snps[qVals>=FDR | abs(snps$BAF-0.5)<=minDeviation,]
  #Keep only those where we can exclude getting the observed number of alt/ref counts by chance at the specified error rate.
  if(errRate>0){
    pVals = pbinom(pmin(snps$altCnt,snps$refCnt)-1,snps$total,errRate,lower.tail=FALSE)
    qVals = p.adjust(pVals,method='BH')
    snps = snps[qVals<FDR,]
  }
  return(snps)
}
