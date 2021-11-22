#' Calculate effective allelic imbalance
#'
#' Given an error rate, copy number ratio, and allele specific expression, calculates the effective allelic imbalance that you would expect to observe in the data.
#'
#' @param CNState The copy number state, expressed as a ratio of number of copies of the maternal allele out of the total number of copies.
#' @param matASE The allele specific expression bias, expressed as the number of counts expected from the maternal allele in a diploid cell.
#' @param errRate The error rate.
#' @return The expected allelic imbalance for this set of parameters.
calcEffectiveImbalance = function(CNState,matASE,errRate){
  #First get the expected allelic ratio, taking into account both ASE and CN state
  effRat = (CNState*matASE)/(CNState*matASE+(1-CNState)*(1-matASE))
  #Account for errors
  effRat = effRat*(1-2*errRate)+errRate
  return(effRat)
}

#' Likelihood model
#'
#' Calcualtes the likelihood model, with or without marginalisation, ASE, etc.  For user facing wrapper see \code{\link{calcLikelihoods}}.
#' 
#' @inheritParams calcEffectiveImbalance
#' @param obsMat Vector of observed maternal counts.
#' @param obsPat Vector of observed paternal counts.
#' @param overDisp Over-dispersion parameter for beta-binomial distribution.
#' @param alpha Alpha parameter of the prior on the allele specific expression ratio (\code{matASE}).
#' @param beta Beta parameter of the prior on the allele specific expression ratio (\code{matASE}).
#' @param defaultErr Default error rate to use if none given.
#' @param marginaliseASE Marginalise or not?
#' @param minVarASE If marginalising and the variance in the ASE is less than this, just use the central estimate.
#' @param verbose Show progress?
#' @param ... Extra things passed to \code{integrate} if marginalising ASE.
#' @return The negative log-likelihood for each observation.
#' @importFrom VGAM dbetabinom
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats dbeta integrate
likelihoodModel = function(obsMat,obsPat,CNState,overDisp,errRate,matASE,alpha,beta,defaultErr=0.01,marginaliseASE=FALSE,minVarASE=1e-4,verbose=marginaliseASE,...){
  tot = obsMat+obsPat
  if(is.null(errRate)){
    warning(sprintf("Error rate not specified.  Using default value of %.02f",errRate))
    errRate = defaultErr
  }
  if(marginaliseASE){
    #Most sensible way to expand out length 1 vectors that need expanding.
    dd = data.frame(k = obsMat,
                    n = tot,
                    CNState = CNState,
                    errRate = errRate,
                    rho = overDisp,
                    alpha = alpha,
                    beta = beta,
                    betaVar = (alpha*beta)/((alpha+beta)**2*(1+alpha*beta))
                    )
    if(verbose)
      pb = txtProgressBar(0,length(obsMat),style=1)
    #Define the integrand
    ff = function(matASE,k,n,CNState,errRate,rho,alpha,beta) dbetabinom(k,n,calcEffectiveImbalance(CNState,matASE,errRate),rho)*dbeta(matASE,alpha,beta)
    dd$ans = NA
    #Marginalise one at a time.
    for(i in seq(nrow(dd))){
      #Don't integrate if beta is too spiky.  There's both no point and integration will likely fail
      if(dd$betaVar[i]<minVarASE){
        dd$ans[i] = dbetabinom(dd$k[i],dd$n[i],calcEffectiveImbalance(dd$CNState[i],matASE[i],dd$errRate[i]),dd$rho[i])
      }else{
        #Do the integration
        tmp = integrate(ff,lower=0,upper=1,k=dd$k[i],n=dd$n[i],CNState=dd$CNState[i],errRate=dd$errRate[i],rho=dd$rho[i],alpha=dd$alpha[i],beta=dd$beta[i],...)
        dd$ans[i] = tmp$value
      }
      if(i%%1000==0 && verbose)
        setTxtProgressBar(pb,i)
    }
    if(verbose)
      close(pb)
    return(-log(dd$ans))
  }else{
    #Get the imbalance we'd expect to observe for this model
    matEff = calcEffectiveImbalance(CNState,matASE,errRate)
    #Calculate the actual likelihood
    return(-dbetabinom(obsMat,tot,matEff,overDisp,log=TRUE))
  }
}

#' Calcualtes the log-likelihood of various copy number states
#'
#' Given allelic counts in individual cells (calculated by \code{getAllelicExpression}), determines the log-likelihood of each cell having a range of copy number configurations at each SNP.  By default, a range of the most plausible copy number states are calculated.  That is, maternal loss, paternal gain, diploid, maternal gain, and paternal loss.  However, any combination can be specified using the \code{CNStates} parameter.
#'
#' The copy number states are given by the ratio of the number of maternal copies to total copies.  So maternal loss would be coded as \code{1/(0+1) = 1}, maternal gain \code{2/(1+2) = 0.33}, diploid \code{1/(1+1) = 0.5}, etc.  The names for \code{CNStates} are used to name output columns, so it is highly advisable to name them sensibly.  For instance, \code{CNStates=c(diploid=0.5,pLoH = 1)} will produce output columns \code{nLL_diploid} and \code{nLL_pLoH}.
#'
#' The naiive model would be to compute a likelihood for each copy number state based on the ratio specified in \code{CNStates} using a binomial model.  This is insufficient for three reasons.  Firstly, sequencing contains errors.  This is most important in regions of loss of heterozygosity (\code{CNStates=0 or 1}), where without errors even a single read mapping to the lost allele would give probability zero.  The more accurately the error rate for each SNP can be specified, the more the likelihoods will represent reality.  Error rates are loaded from \code{gCnts} from the \code{errRate} column and set to \code{defaultErr} when not found.
#'
#' Secondly, some genes have a persistent bias for one allele over another.  These genes will deviate from the 50/50 allele ratio of expression even in diploid regions.  When calculating the likelihood for the diploid state (\code{CNStates=0.5}), the allele specific expression of genes is taken into account.  This information is taken from the \code{matASE} column of \code{gCnts}, which gives the fraction of maternal counts expected in diploid cells for each SNP.  These values are usually calculated by running \code{\link{calcASE}} on \code{gCnts}.  This procedure will actually produce a probability distribution of likely allele specific expression values.  If \code{marginaliseASE} is set to \code{TRUE} the likelihood is obtained by marginalising over this distribution rather than using the most likely estimate.  Note that marginalising is slower and usually not worth the extra effort.  In the extreme case of imprinted genes, only one allele is expressed in all circumstances and so no information is left to determine the copy number state.  Therefore, imprinted genes are usually filtered out and allele specific expression is calculated by \code{\link{calcASE}}.  To prevent allele specific expression being considered at all, set \code{useASE=FALSE}.
#'
#' The final reason the naiive model is insufficient is that transcription occurs in bursts on each allele and burst times are not perfectly correlated across alleles.  This means that more deviation from the ratio given in \code{CNStates} is expected that a binomial model would expect.  To account for this, we use a beta-binomial model, where the over-dispersion is parameterised by \code{overDispersion} and usually calcualted from the data using \code{\link{calcOverDispersion}}.
#'
#' @inheritParams calcOverDispersion
#' @param overDisp The beta-binomial over-dispersion parameter.  Usually calculated by \code{\link{calcOverDispersion}}.
#' @param CNStates Specify which copy number states to test.  States are specified as ratios of number of maternal copies to total copies.  Names are used to name output columns.
#' @param useASE Should allele specific expression be incorporated into the model?
#' @param verbose Be verbose?
#' @param ... Passed to \code{\link{likelihoodModel}}.
#' @return \code{gCnts} with extra columns giving the negative log-likelihood for the model under each copy number state.  Extra columns named as \code{paste0('nLL_',names(CNStates))}.
#' @importFrom VGAM dbetabinom
#' @export
calcLikelihoods = function(gCnts,overDisp,CNStates=c(pLoss=1,mGain=2/3,diploid=1/2,pGain=1/3,mLoss=0),useASE=TRUE,marginaliseASE=FALSE,verbose=marginaliseASE,...){
  if(is.null(names(CNStates)) || any(duplicated(names(CNStates))))
    stop("CN States must be named uniquely.")
  if(useASE & !'matASE' %in% colnames(mcols(gCnts)))
    stop("Allele specific expression estimates not found so not used.  Run calcASE first or set useASE=FALSE.")
  if(useASE & marginaliseASE & !all(c('matASE_postAlpha','matASE_postBeta') %in% colnames(mcols(gCnts))))
    stop("Columns matASE_postAlpha and/or matASE_postBeta not found.  Marginalisation of ASE requires knowledge of distribution.")
  #Loop over the states to test
  for(nom in names(CNStates)){
    tgt = sprintf('nLL_%s',nom)
    if(verbose)
      message(sprintf("Calculating likelihood for state %s",nom))
    #Get the states to use.  Either a vector or a number
    tmp = unlist(CNStates[nom],use.names=FALSE)
    #Should we ignore ASE?
    if(!useASE){
      mcols(gCnts)[,tgt] = likelihoodModel(gCnts$matCount,gCnts$patCount,tmp,overDisp,gCnts$errRate,0.5,marginaliseASE=FALSE,verbose=verbose)
    }else{
      mcols(gCnts)[,tgt] = likelihoodModel(gCnts$matCount,gCnts$patCount,tmp,overDisp,gCnts$errRate,gCnts$matASE,gCnts$matASE_postAlpha,gCnts$matASE_postBeta,marginaliseASE=marginaliseASE,verbose=verbose)
    }
  }
  #Store map of CN state names to values
  gCnts@metadata$CNStates = CNStates
  return(gCnts)
}
