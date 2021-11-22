#' 1000 genomes SNPs hg19
#'
#' The 1000 genomes SNPs Phase3, mapped to hg19 is the starting point.  This is filtered to include only biallelic SNPs with a population allele frequency of at least 5%.  This is all squashed into a \code{GRanges} object with minimal metadata.
#'
#' The following gives a complete description of how \code{SNPs1k} was generated.
#' \itemize{
#' \item  #Download the starting data
#' \item  curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/|grep 'ALL.chr'  |sed 's/.a.href="\(.*\)">.*/\1/g' > tmp.txt
#' \item  for file in `cat tmp.txt`; do  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$file; done
#' \item  #Now do the processing in R
#' \item library(VariantAnnotation)
#' \item tgts = list.files('.')
#' \item tgts = tgts[grepl('^ALL.*vcf.gz$',tgts)]
#' \item #Process them one chromosome at a time
#' \item param = ScanVcfParam(info = c('AF','EAS_AF','EUR_AF','AFR_AF','AMR_AF','SAS_AF','VT'), geno=NA)
#' \item for(tgt in tgts){
#' \item   message(sprintf("Loading chromosome %s",tgt))
#' \item   vcf = readVcfAsVRanges(tgt,'hg19',param=param)
#' \item   #Do some sanity checks
#' \item   good = nchar(ref(vcf))==1 & nchar(alt(vcf))==1 & width(vcf)==1
#' \item   #More length sanity checks
#' \item   good = rep(TRUE,length(vcf))
#' \item   for(nom in colnames(mcols(vcf)))
#' \item     good = good & lengths(mcols(vcf)[,nom])==1
#' \item   vcf = vcf[good]
#' \item   vcf = vcf[unlist(vcf$VT)=='SNP']
#' \item   bad = duplicated(start(vcf))
#' \item   if(any(bad))
#' \item     vcf = vcf[!start(vcf) %in% start(vcf)[bad]]
#' \item   #Construct GRanges
#' \item   tmp = mcols(vcf)
#' \item   for(nom in colnames(tmp))
#' \item     tmp[,nom] = unlist(tmp[,nom])
#' \item   #Drop uninteresting ones
#' \item   tmp = tmp[,!colnames(tmp) %in% c('QUAL','VT')]
#' \item   vcf = GRanges(as.character(seqnames(vcf)),IRanges(start(vcf),width=1),REF=ref(vcf),ALT=alt(vcf),seqinfo=seqinfo(vcf))
#' \item   mcols(vcf) = cbind(mcols(vcf),tmp)
#' \item   #Save big version
#' \item   chr = gsub('^ALL.(chr.*?)\\..*','\\1',tgt)
#' \item   chrLen = seqlengths(vcf)[unique(as.character(seqnames(vcf)))]
#' \item   message(sprintf("  Full thing has %s SNPs, that's one per %d bases",prettyNum(length(vcf),big.mark=','),round(chrLen/length(vcf))))
#' \item   saveRDS(vcf,paste0('SNPs1kBiallelic.All.',chr,'.RDS'))
#' \item   #Apply AF filter.  Exclude the high end as well
#' \item   vcf = vcf[abs(vcf$AF-0.5)<(0.5-0.01),]
#' \item   message(sprintf("  Filtered to 0.01 has %s SNPs, that's one per %d bases",prettyNum(length(vcf),big.mark=','),round(chrLen/length(vcf))))
#' \item   saveRDS(vcf,paste0('SNPs1kBiallelic.AF.01.',chr,'.RDS'))
#' \item   vcf = vcf[abs(vcf$AF-0.5)<(0.5-0.05),]
#' \item   message(sprintf("  Filtered to 0.05 has %s SNPs, that's one per %d bases",prettyNum(length(vcf),big.mark=','),round(chrLen/length(vcf))))
#' \item   saveRDS(vcf,paste0('SNPs1kBiallelic.AF.05.',chr,'.RDS'))
#' \item   vcf = vcf[abs(vcf$AF-0.5)<(0.5-0.1),]
#' \item   message(sprintf("  Filtered to 0.10 has %s SNPs, that's one per %d bases",prettyNum(length(vcf),big.mark=','),round(chrLen/length(vcf))))
#' \item   saveRDS(vcf,paste0('SNPs1kBiallelic.AF.10.',chr,'.RDS'))
#' \item   vcf = vcf[abs(vcf$AF-0.5)<(0.5-0.2),]
#' \item   message(sprintf("  Filtered to 0.20 has %s SNPs, that's one per %d bases",prettyNum(length(vcf),big.mark=','),round(chrLen/length(vcf))))
#' \item   saveRDS(vcf,paste0('SNPs1kBiallelic.AF.20.',chr,'.RDS'))
#' \item }
#' \item #Load in the AF restricted data and merge
#' \item SNPs1k = unlist(GRangesList(lapply(grep('SNPs1kBiallelic.AF.10.chr.*.RDS',list.files('.'),value=TRUE),readRDS)))
#' \item save(SNPs1k,file='SNPs1k.RData',compress='xz')
#' }
#' 
#' @docType data
#'
#' @usage data(SNPs1k)
#'
#' @source \href{http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/}{1000GenomesFTP}
#'
#' @format An object of class \code{"GRanges"}; see \code{\link[GenomicRanges]{GRanges}}.
#'
#' @keywords datasets
"SNPs1k"
