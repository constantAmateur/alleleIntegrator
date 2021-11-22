#' Load Mutations from CanPipe
#'
#' Loads somatic point mutations from Caveman output file (flag.vcf.gz) and format them as needed.  Only point mutations with FILTER=='PASS' are kept.
#'
#' @param cavemanOutput Path to caveman mutation output file (the one ending flag.vcf.gz).
#' @return GRanges file with mutations.
#' @export
loadCanPipeMuts = function(cavemanOutput){
  muts = readVcf(cavemanOutput)
  out = rowRanges(muts)
  #Which ones to keep
  w = out$FILTER=='PASS' & lengths(out$REF)==1 & lengths(out$ALT)==1
  out = out[w]
  muts = muts[w]
  #Make ALT/REF characters
  out$ALT = as.character(unlist(out$ALT))
  out$REF = as.character(out$REF)
  #Get the extra information
  mat = cbind((geno(muts)$FAZ + geno(muts)$RAZ)[,'NORMAL'],
              (geno(muts)$FCZ + geno(muts)$RCZ)[,'NORMAL'],
              (geno(muts)$FGZ + geno(muts)$RGZ)[,'NORMAL'],
              (geno(muts)$FTZ + geno(muts)$RTZ)[,'NORMAL'])
  out$normRefCnt = mat[cbind(seq_along(out),match(out$REF,c('A','C','G','T')))]
  out$normAltCnt = mat[cbind(seq_along(out),match(out$ALT,c('A','C','G','T')))]
  mat = cbind((geno(muts)$FAZ + geno(muts)$RAZ)[,'TUMOUR'],
              (geno(muts)$FCZ + geno(muts)$RCZ)[,'TUMOUR'],
              (geno(muts)$FGZ + geno(muts)$RGZ)[,'TUMOUR'],
              (geno(muts)$FTZ + geno(muts)$RTZ)[,'TUMOUR'])
  out$tumRefCnt = mat[cbind(seq_along(out),match(out$REF,c('A','C','G','T')))]
  out$tumAltCnt = mat[cbind(seq_along(out),match(out$ALT,c('A','C','G','T')))]
  #The VAFs
  out$normVAF = geno(muts)$PM[,'NORMAL']
  out$tumVAF = geno(muts)$PM[,'TUMOUR']
  #Get rid of the pointless columns
  out$QUAL = NULL
  out$paramRangeID = NULL
  return(out)
}
