#' Load SNPs from CanPipe
#'
#' Loads Caveman SNPs from canpipe file and format them as needed. 
#'
#' @param cavemanOutput Path to caveman SNP output.
#' @return GRanges SNP file.
#' @importFrom VariantAnnotation geno
#' @export
loadCanPipeSNP = function(cavemanOutput){
  snps = readVcf(cavemanOutput)
  out = rowRanges(snps)
  out$normGT = geno(snps)$GT[,'NORMAL']
  out$tumGT = geno(snps)$GT[,'TUMOUR']
  #Which ones to keep
  w = out$normGT=='0|1' & lengths(out$REF)==1 & lengths(out$ALT)==1
  out = out[w]
  snps = snps[w]
  #Make ALT/REF characters
  out$ALT = as.character(unlist(out$ALT))
  out$REF = as.character(out$REF)
  #Get the matrix of counts
  mat = cbind((geno(snps)$FAZ + geno(snps)$RAZ)[,'NORMAL'],
              (geno(snps)$FCZ + geno(snps)$RCZ)[,'NORMAL'],
              (geno(snps)$FGZ + geno(snps)$RGZ)[,'NORMAL'],
              (geno(snps)$FTZ + geno(snps)$RTZ)[,'NORMAL'])
  out$refCnt = mat[cbind(seq_along(out),match(out$REF,c('A','C','G','T')))]
  out$altCnt = mat[cbind(seq_along(out),match(out$ALT,c('A','C','G','T')))]
  out$total = rowSums(mat)
  #Get rid of the pointless columns
  out$QUAL = NULL
  out$FILTER = NULL
  out$paramRangeID = NULL
  return(out)
}

