#' Runs allele Counter in parallel over a set of input and output files
#'
#' Calculates the pileup using alleleCounter binary for a set of BAM files and target locations in parallel.  All parameters to alleleCounter can be either length 1 or the length of \code{bams}, in which case each parameter is matched to the corresponding BAM file.
#'
#' If \code{outputs} is specified, the results will be saved to a file and reused instead of re-running alleleCounter if called again.
#'
#' The default alleleCounter parameters are tuned for 10X output.
#'
#' Note that there are two levels of possible parallel execution.  Firstly, when multiple BAM files are specified, they can be processed simultaneously.  Secondly, the loci can be split into groups (e.g. by chromosome) and each group can be processed simultaneously.  Rather than split by chromosome, which will result in very uneven numbers of loci in each group, this function instead splits the input loci up into \code{nChunks} roughly equally sized groups.  Note increasing \code{nChunks} may actually compromise perfomance by creating too many threads trying to read from the same file at once.
#'
#' The parallel execution is further controlled by the \code{nMaxSim} parameter.  When processing a large number of BAM files simultaneously, this parameter will split them up into groups of size \code{nMaxSim} or smaller.  The purpose of this is so that the output is progressively written across the course of a job rather than in one big lump at the end.  As such, this parameter is ignored with \code{outputs} is set to NULL.
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
#' @param autoChr Try to automatically strip or prepend 'chr' to make loci and BAM files match.
#' @param bin The location of the alleleCounter binary.
#' @param nParallel Number of processors to use.  Should be a multiple of number of samples.
#' @param nChunks When running in parallel, split into this many chunks.  Set this to roughly the number of threads that can simultaneously read a file.
#' @param nMaxSim Don't do more than too many samples simultaneously.  This prevents computation being wasted if a big job fails at the end.
#' @param skipIfExists If \code{outputs} already exists.  Load results from there instead of calculating from scratch.
#' @param ignoreOutputs Suppress printing of outputs from alleleCounter binary.
#' @return A list of GRanges objects, with each entry giving counts of A,C,G,T (and total) at the target loci provided.
#' @importFrom utils read.table read.delim write.table
#' @export
alleleCounter = function(bams,refGenome,tgtLoci,outputs=NULL,f=0,F=0,x=TRUE,d=TRUE,m=20,q=200,autoChr=TRUE,bin='alleleCounter',nParallel=1,nChunks=8,nMaxSim=6,skipIfExists=TRUE,ignoreOutputs=TRUE){
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
  #Don't do a huge number of files at the same time.  This prevents having to redo large calculations on failure of the final files.
  #There's no point in this if we're not saving the outputs
  if(!is.null(outputs) && length(bams)>nMaxSim){
    out = list()
    #Loop over groups
    for(w in split(seq_along(bams),(seq_along(bams)-1) %/% nMaxSim)){
      out = c(out,alleleCounter(bams=bams[w],
                                refGenome = if(length(refGenome)==1){refGenome}else{refGenome[w]},
                                tgtLoci = if(length(tgtLoci)==1){tgtLoci}else{tgtLoci[w]},
                                outputs = outputs[w],
                                f = if(length(f)==1){f}else{f[w]},
                                F = if(length(F)==1){F}else{F[w]},
                                x = if(length(x)==1){x}else{x[w]},
                                d = if(length(d)==1){d}else{d[w]},
                                m = if(length(m)==1){m}else{m[w]},
                                q = if(length(q)==1){q}else{q[w]},
                                autoChr=autoChr,
                                bin=bin,
                                nParallel=nParallel,
                                nChunks=nChunks,
                                nMaxSim=nMaxSim,
                                skipIfExists=skipIfExists
                                ))
    }
    return(out)
  }
  if(autoChr){
    #tgtLoci may now become non-unique due to chr matching
    if(length(tgtLoci)==1)
      tgtLoci = tgtLoci[rep(1,length(bams))]
    #Do the bams have chr?
    bamsHaveChr = sapply(bams,function(e) any(grepl('^chr',seqlevels(BamFile(e)))))
    lociHaveChr = sapply(tgtLoci,function(e) any(grepl('^chr',seqlevels(e))))
    needChanging = which(bamsHaveChr!=lociHaveChr)
  }else{
    needChanging=integer()
  }
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
  #Make this more parallelisable by chopping into nParallel * x pieces, such that each piece isn't too large
  nChunksTot = ifelse(nParallel==1,1,nChunks)
  #Split loci into chromosomes and save
  lociPaths = list()
  for(i in seq_along(tgtLoci)){
    #Get the loci for this BAM
    tmp = data.frame(chr=as.character(seqnames(tgtLoci[[i]])),pos=start(tgtLoci[[i]]))
    #Add in (or remove) chr as needed
    if(i %in% needChanging){
      if(lociHaveChr[i]){
        tmp$chr = gsub('^chr','',tmp$chr)
      }else{
        tmp$chr = paste0('chr',tmp$chr)
      }
    }
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
    system(sprintf("parallel -j %d < %s",nParallel,rFile),ignore.stderr=ignoreOutputs,ignore.stdout=ignoreOutputs)
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
      #Convert it back if needed
      if(i %in% needChanging){
        if(lociHaveChr[i]){
          tmp$chr = paste0('chr',tmp$chr)
        }else{
          tmp$chr = gsub('^chr','',tmp$chr)
        }
      }
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
