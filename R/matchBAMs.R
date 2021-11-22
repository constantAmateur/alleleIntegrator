#' Checks if BAM files have matching genotype. 
#'
#' Calculates identity by similarity for a series of BAM files. Looks for locations that support the ALT allele in the target regions in a reference sample.  Then counts at those locations in all BAMs to see if the genotypes match.  Intended use is with a \code{regions} parameter containing exons.
#'
#' If not given, the SNPs are calculated from a reference sample, where heterozygous SNPs are identified.  This only really works well when you have exome or genome sequencing data, otherwise you are better off providing a pre-calculated set of SNPs.
#'
#' @inheritParams alleleCounter
#' @param BAMs The BAM files to genotype.  Names are used as labels.
#' @param refGenomes The FASTA file for the genome of each sample.
#' @param snps A GRanges object with SNPs to inspect.  Must contain REF and ALT columns.
#' @param outputs Store output SNP counts to prevent recalculation.  Vector of length to match \code{BAMs}.
#' @param liftOvers If one or more of the BAMs uses a different coordinate system to the reference, lift the coordinates over with these files.  Vector of length BAMs, with NAs for no-liftover.
#' @param is10X If a BAM is 10X, set to TRUE so it can be counted properly.
#' @param usesChr Do any samples have a 'chr' prefix on chromosome names.
#' @param doPlot Make a plot of the results?
#' @param colPal Colour scheme to use.  Passed to \code{\link[RColorBrewer]{brewer.pal}}.
#' @param ... Pass to \code{\link{alleleCounter}}.
#' @return A list with all the different bits of information.  Probably the main bit you care about is \code{ibs$ibs}
#' @importFrom SNPRelate snpgdsCreateGeno snpgdsOpen snpgdsLDpruning snpgdsIBS
#' @importFrom stats dbinom setNames
#' @importFrom ComplexHeatmap draw Heatmap
#' @importFrom GenomeInfoDb renameSeqlevels
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @export
matchBAMs = function(BAMs,refGenomes,snps=ExAcSNPs,outputs=NULL,liftOvers=NULL,is10X=FALSE,usesChr=FALSE,nParallel=1,doPlot=TRUE,colPal='Greens',...){
  if(is.null(names(BAMs)))
    names(BAMs) = paste0('BAM_number',seq_along(BAMs))
  w = which(names(BAMs)=='')
  if(length(w)>0)
    names(BAMs)[w] = paste0('BAM_number',seq(length(w)))
  if(any(duplicated(names(BAMs))))
    names(BAMs) = make.unique(names(BAMs))
  if(is.null(liftOvers))
    liftOvers=NA
  if(length(liftOvers)==1)
    liftOvers = rep(liftOvers,length(BAMs))
  if(length(is10X)==1)
    is10X = rep(is10X,length(BAMs))
  if(length(refGenomes)==1)
    refGenomes = rep(refGenomes,length(BAMs))
  if(length(usesChr)==1)
    usesChr = rep(usesChr,length(BAMs))
  if(!is.null(outputs) && length(outputs)!=length(BAMs))
    stop("outputs must have same length as BAMs")
  #Check that everything is what it should be
  if(all(c('REF','ALT') %in% colnames(snps)))
    stop("snps must have REF and ALT")
  if(!is.character(snps$REF))
    stop("REF and ALT must be a character vector")
  if(!is.character(snps$ALT))
    stop("REF and ALT must be a character vector")
  refSNPs = snps[nchar(snps$REF)==1 & nchar(snps$ALT)==1]
  #Obliterate chr if it's there
  refSNPs = renameSeqlevels(refSNPs,setNames(gsub('^chr','',seqlevels(refSNPs)),seqlevels(refSNPs)))
  #Useful for matching when using liftover
  refSNPs$snpID = paste0('SNP',seq_along(refSNPs))
  #Add default AF if not given
  if(!'AF' %in% colnames(refSNPs))
    refSNPs$AF=0.5
  #Convert NAs to string NA so we can work with them in a standard way
  liftOvers[is.na(liftOvers)]='NA'
  #Do each split differently
  snpCnts = list()
  for(tgtLiftOver in unique(liftOvers)){
    for(tgt10X in c(TRUE,FALSE)){
      for(tgtUsesChr in c(TRUE,FALSE)){
        for(refGenome in unique(refGenomes)){
          #Check if there are any matches for this combinatio
          w = liftOvers==tgtLiftOver & is10X==tgt10X & usesChr==tgtUsesChr & refGenomes==refGenome
          if(any(w)){
            message(sprintf("Processing %d samples",sum(w)))
            #What are the outputs
            if(is.null(outputs)){
              tgtOuts = NULL
            }else{
              tgtOuts = outputs[w]
            }
            #Do we need to do the liftover?
            if(tgtLiftOver!='NA'){
              ch = rtracklayer::import.chain(tgtLiftOver)
              tgtSNPs = unlist(rtracklayer::liftOver(refSNPs,ch))
            }else{
              tgtSNPs = refSNPs
            }
            if(tgtUsesChr){
              tgtSNPs = renameSeqlevels(tgtSNPs,setNames(paste0('chr',seqlevels(tgtSNPs)),gsub('chr','',seqlevels(tgtSNPs))))
            }
            if(tgt10X){
              out = alleleCounter(BAMs[w],refGenome,tgtSNPs,outputs=tgtOuts,nParallel=nParallel,...)
            }else{
              out = alleleCounter(BAMs[w],refGenome,tgtSNPs,outputs=tgtOuts,x=FALSE,f=3,F=3852,m=20,q=35,nParallel=nParallel,...)
            }
            if(tgt10X){
              #We don't actually care about the cellID column for 10X, so aggregate over it
              out = lapply(out,function(e) {
                             tmp = mcols(e)
                             tmp$mark = as.character(e)
                             tmp = aggregate(cbind(A,C,G,T,Tot) ~ mark,data=tmp,FUN=sum)
                             m = match(tmp$mark,as.character(e))
                             e = e[m]
                             e$barcode = NULL
                             e$A = tmp$A
                             e$C = tmp$C
                             e$G = tmp$G
                             e$T = tmp$T
                             e$Tot = tmp$Tot
                             e
                 })
            }
            #Put them where they should go
            for(i in seq_along(out)){
              snpCnts[[(names(BAMs)[w])[i]]] = out[[i]]
            }
          }
        }
      }
    }
  }
  #Make into a giant comparison matrix
  cnts = matrix(NA,
                nrow=length(refSNPs),
                ncol=length(BAMs),
                dimnames = list(refSNPs$snpID,names(BAMs))
                )
  altCnts = cnts
  #Fill in each
  for(i in seq_along(snpCnts)){
    jj = match(names(snpCnts)[i],colnames(cnts))
    ii = match(snpCnts[[i]]$snpID,rownames(cnts))
    #Get the Alt allele frequency for each
    bases=c('A','C','G','T')
    tmp = mcols(snpCnts[[i]])
    tmp = as.matrix(tmp[,bases])
    cnts[ii,jj] = rowSums(tmp)
    altCnts[ii,jj] = tmp[cbind(seq_along(ii),match(snpCnts[[i]]$ALT,bases))]
  }
  #What is the evidence for AA,Aa,aa?
  pAA = dbinom(altCnts,cnts,0.1,log=TRUE)
  pAa = dbinom(altCnts,cnts,0.5,log=TRUE)
  paa = dbinom(altCnts,cnts,0.9,log=TRUE)
  #Work out the normalised prob of aa, the most interesting one
  pp = data.frame(pAA=as.vector(pAA),
                  pAa=as.vector(pAa),
                  paa=as.vector(paa))
  pp = exp(pp)/rowSums(exp(pp),na.rm=TRUE)
  pp$sampleID = rep(colnames(cnts),each=nrow(cnts))
  pp$snpID = rep(rownames(cnts),ncol(cnts))
  pp$snpAF = rep(refSNPs$AF,ncol(cnts))
  #Make call of the number of Alt alleles
  pp$numAlt = rep(3,nrow(pp))
  pp$numAlt[pp$paa>0.9]=2
  pp$numAlt[pp$pAa>0.9]=1
  pp$numAlt[pp$pAA>0.9]=0
  #Make this into a genotype matrix
  gtMat = matrix(3,
                 nrow=length(refSNPs),
                 ncol=length(BAMs),
                 dimnames = list(refSNPs$snpID,names(BAMs))
                 )
  #Fill it in
  gtMat[cbind(match(pp$snpID,rownames(gtMat)),match(pp$sampleID,colnames(gtMat)))]=pp$numAlt
  #Make the necessary object for relatdness
  gdsFile = tempfile()
  snpgdsCreateGeno(gdsFile,
                   genmat = gtMat,
                   sample.id=colnames(gtMat),
                   snp.id=rownames(gtMat),
                   snp.chromosome = as.character(seqnames(refSNPs)),
                   snp.position = start(refSNPs),
                   snp.allele = paste0(refSNPs$ALT,'/',refSNPs$REF),
                   snpfirstdim=TRUE,
                   )
  gds = snpgdsOpen(gdsFile)
  gdsSub = snpgdsLDpruning(gds,ld.threshold=0.2)
  ibs = snpgdsIBS(gds,num.thread=nParallel)
  #Add row/col names to ibs thing
  rownames(ibs$ibs) = colnames(ibs$ibs) = names(BAMs)
  if(doPlot){
    #Parameter checking
    colFun = suppressWarnings(brewer.pal(100,colPal))
    colFun = colorRamp2(seq(0.5,1,length.out=length(colFun)),colFun)
    hm = Heatmap(ibs$ibs,
                 col=colFun,
                 name = 'IBS',
                 show_row_names=TRUE,
                 show_column_names=TRUE,
                 show_row_dend=FALSE,
                 show_column_dend=FALSE
                 )
    draw(hm)
  }
  #Return everything
  return(list(ibs=ibs,pp=pp,snpCnts=snpCnts,refSNPs))
}

