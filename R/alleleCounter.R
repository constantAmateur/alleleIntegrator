#' Runs allele Counter in parallel over a set of input and output files
#'
#' Calculates the pileup using alleleCounter binary for a set of BAM files and target locations in parallel.  All parameters to alleleCounter can be either length 1 or the length of \code{bams}, in which case each parameter is matched to the corresponding BAM file.
#'
#' If \code{outputs} is specified, the results will be saved to a file and reused instead of re-running alleleCounter if called again.
#'
#' The default alleleCounter parameters are tuned for 10X output.
#'
#' @param bams The BAM files to get allele counts from. 
#' @param refGenome The reference genome each BAM file was mapped using.
#' @param tgtLoci Locations to interogate as a GRanges object.  Either one per sample or if length 1 the same loci are used for all samples.
#' @param outputs The name of the file to write allele counts for each sample.  Lengths of outputs and bams must match.  If NULL, a temp file is used and deleted.
#' @param f Allelecounter -f param
#' @param F Allelecounter -F param
#' @param x Include the -x option?
#' @param d Include the -d option?
#' @param m Allelecounter -m param
#' @param q Allelecounter -q param
#' @param bin The location of the alleleCounter binary.
#' @param nParallel Number of processors to use.  Should be a multiple of number of samples.
#' @param nChunks When running in parallel, split into this many chunks per parallel thread.
#' @param skipIfExists If \code{outputs} already exists.  Load results from there instead of calculating from scratch.
#' @return A list of GRanges objects, with each entry giving counts of A,C,G,T (and total) at the target loci provided.
#' @importFrom utils read.table read.delim write.table
#' @export
alleleCounter = function(bams,refGenome,tgtLoci,outputs=NULL,f=0,F=0,x=TRUE,d=TRUE,m=20,q=200,bin='alleleCounter',nParallel=1,nChunks=4,skipIfExists=TRUE){
  if(is.null(outputs))
    warning("output files not specified, so results cannot be reused without recalculation.  These calculations are time consuming, so it is advisable to save them somewhere.")
  #make tgtLoci a list
  if(!is.list(tgtLoci))
    tgtLoci = list(tgtLoci)
  #Check things are the correct length
  if(!(length(tgtLoci) %in% c(1,length(bams))) ||
     !(length(refGenome) %in% c(1,length(bams))) || 
     (!is.null(outputs) & length(outputs)!=length(bams)) || 
     !(length(f) %in% c(1,length(bams))) ||
     !(length(F) %in% c(1,length(bams))) ||
     !(length(x) %in% c(1,length(bams))) ||
     !(length(d) %in% c(1,length(bams))) ||
     !(length(m) %in% c(1,length(bams))) ||
     !(length(q) %in% c(1,length(bams))))
    stop("Length of all parameters must be 1, or match bams length")
  #Build a data.frame of inputs
  params = data.frame(callIdx = seq_along(bams),
                      bam = bams,
                      loci = if(length(tgtLoci)==1) rep(1,length(bams)) else seq_along(bams), # Point all bams at the same loci file, or their unique one.
                      refGenome = refGenome,
                      f=f,
                      F=F,
                      x=x,
                      d=d,
                      m=m,
                      q=q)
  #Construct the base string
  baseCommand = paste0(bin, ' -l %s -b %s -o %s -r %s -f %d -F %d -m %d -q %d') 
  #Make this more parallelisable by chopping into nParallel * x pieces, such that each piece isn't too small.
  nChunksTot = ifelse(nParallel==1,1,nParallel*nChunks)
  #Split loci into chromosomes and save
  lociPaths = list()
  for(i in seq_along(tgtLoci)){
    #Get the loci for this BAM
    tmp = data.frame(chr=as.character(seqnames(tgtLoci[[i]])),pos=start(tgtLoci[[i]]))
    #Try and split into nChunksTot parts.  Ensure that each chunk has at least ~1000 loci
    nEff = min(nChunksTot,ceiling(nrow(tmp)/1000))
    splitIdx = seq(nrow(tmp)) %/% (nrow(tmp)/nEff)
    #Fix the allocation of last entry
    splitIdx[splitIdx>=nEff] = nEff-1
    dd = data.frame(chunkIdx = seq(nEff),
                    tFile = sapply(seq(nEff),function(e){
                                     tFile = tempfile()
                                     write.table(tmp[splitIdx == e-1,,drop=FALSE],tFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
                                     tFile}),
                    lociIdx = i)
    lociPaths[[i]] = dd
  }
  #Loop over samples and construct all commands to run
  toRun=list()
  for(i in seq_along(bams)){
    #Check if we actually need to run this one.
    if(skipIfExists && !is.null(outputs) && file.exists(outputs[i]))
      next
    #Base the construction on the loci data.frame
    dd = lociPaths[[params$loci[i]]]
    dd$outputs = sapply(seq(nrow(dd)),function(e) tempfile()) # Create per-chunk temporary output files
    #Build the final commands to run
    dd$callIdx = params$callIdx[i]
    dd$bam = params$bam[i]
    dd$refGenome = params$refGenome[i]
    dd$f = params$f[i]
    dd$F = params$F[i]
    dd$x = params$x[i]
    dd$d = params$d[i]
    dd$m = params$m[i]
    dd$q = params$q[i]
    dd$cmd = sprintf(baseCommand,dd$tFile,dd$bam,dd$outputs,dd$refGenome,dd$f,dd$F,dd$m,dd$q)
    dd$cmd = paste0(dd$cmd,ifelse(dd$x,' -x',''),ifelse(dd$d,' -d',''))
    toRun[[length(toRun)+1]] = dd
  }
  toRun = do.call(rbind,toRun)
  if(length(toRun)>0){
    #Write and run in parallel
    rFile = tempfile()
    write.table(toRun$cmd,rFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
    #Actually run things
    system(sprintf("parallel -j %d < %s",nParallel,rFile))
  }
  #Load everything
  AChead = c('chr','pos','A','C','G','T','Tot')
  ACheadx = c('chr','pos','barcode','A','C','G','T','Tot')
  weFailed=FALSE
  outs=list()
  for(i in seq_along(bams)){
    if(skipIfExists && !is.null(outputs) && file.exists(outputs[i])){
      tmp = read.table(outputs[i],sep='\t',header=TRUE)
      if(ncol(tmp) != ifelse(params$x[i],length(ACheadx),length(AChead)) || !all(colnames(tmp) == if(params$x[i]) ACheadx else AChead)){
        warning(sprintf("Output for BAM file %s at %s, does not match expected format.",bams[i],outputs[i]))
        weFailed=TRUE
      }
    }else{
      #Have to load one at a time.
      tmp = toRun[toRun$callIdx == i,,drop=FALSE]
      tmp = lapply(tmp$outputs,read.delim,sep='\t')
      tmp = do.call(rbind,tmp)
      colnames(tmp) = if(params$x[i]) ACheadx else AChead
      #Should we save it?
      if(!is.null(outputs))
        write.table(tmp,outputs[i],sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
    #Keep loci as the basis
    out = tgtLoci[[params$loci[i]]]
    #Match results back to source.  This looks different for 10X mode.
    if(params$x[i]){
      #In 10X mode, we may have 0, 1, or many lines of output per loci in tgtLoci (depending on how many cell barcodes exist at this location)
      if(nrow(tmp)==0){
       warning(sprintf("Output for BAM file %s is empty.  alleleCounter may not have run properly for this sample.",bams[i]))
       out = out[NULL,]
      }else{
        m = match(paste0(tmp$chr,':',tmp$pos),as.character(out))
        #Sanity check
        if(any(is.na(m))){
          warning(sprintf("Output for BAM file %s does not match input loci.",bams[i]))
          weFailed=TRUE
        }
        #Construct the output
        out = out[m,]
      }
      mcols(out) = cbind(mcols(out),tmp[,-c(1,2)]) #Don't need to keep chr/pos, we used that to match
    }else{
      #In non-10X mode we must have one line of output per line of input.
      m = match(as.character(out),paste0(tmp$chr,':',tmp$pos))
      #Sanity checks
      if(length(out)!=nrow(tmp) || any(is.na(m))){
        warning(sprintf("Output for BAM file %s does not match input loci.",bams[i]))
        weFailed=TRUE
      }
      #Construct the output
      mcols(out) = cbind(mcols(out),tmp[m,-c(1,2)]) #Don't need to keep chr/pos, we used that to match
    }
    outs[[i]] = out
  }
  #Cleanup
  if(length(toRun)>0){
    unlink(toRun$inputs)
    unlink(rFile)
  }
  for(tmp in lociPaths){
    unlink(tmp$tFile)
  }
  if(weFailed)
    stop("Not all output loaded successfully.")
  return(outs)
}
