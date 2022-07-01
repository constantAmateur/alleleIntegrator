#' Runs at startup
#'
#' Checks if singularity is installed on system and sets the binary paths appropriately.
#'
#' @param libname Library name.
#' @param pkgname Package name.
#' @return Nothing.
.onLoad = function(libname,pkgname){
  #Check for singularity
  hasSing = system('which singularity',ignore.stdout=TRUE,ignore.stderr=TRUE)==0
  hasAC = system('which alleleCounter',ignore.stdout=TRUE,ignore.stderr=TRUE)==0
  hasBCF = system('which bcftools',ignore.stdout=TRUE,ignore.stderr=TRUE)==0
  hasParallel = system('which parallel',ignore.stdout=TRUE,ignore.stderr=TRUE)==0
  #message(sprintf("At load, singularity=%d, alleleCounter=%d, bcftools=%d, and parallel=%d",hasSing,hasAC,hasBCF,hasParallel))
  if(hasParallel){
    .parBin = 'parallel'
  }else{
    warning('Binary "parallel", not found on system.  Do you need to run "apt install parallel"?\nalleleIntegrator\'s parallel execution relies on this binary, will run on just one core unless installed.') 
    .parBin = NULL
  }
  if(hasSing){
    .acBin = system.file('singularity','alleleCounter',package='alleleIntegrator')
    .bcfBin = system.file('singularity','bcftools',package='alleleIntegrator')
  }else{
    warning("Singularity not found on system.  External calls will rely on alleleCounter and bcftools being installed.")
    #Assume they're installed on the system
    .acBin = 'alleleCounter'
    .bcfBin = 'bcftools'
    #Print more frantic warnings if they're not found
    if(!hasAC)
      warning("alleleCounter binary not found.  Most of alleleIntegrator package will not work!")
    if(!hasBCF)
      warning("bcftools binary not found.  Variant and SNP finding functions will not work!")
  }
  #Now store things
  assign('.acBin',.acBin,envir = parent.env(environment()))
  assign('.bcfBin',.bcfBin,envir = parent.env(environment()))
  assign('.parBin',.parBin,envir = parent.env(environment()))
}
