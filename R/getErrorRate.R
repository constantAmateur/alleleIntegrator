#' Get error rate
#'
#' Retrieves error rate previously calculated using \code{calcErrorRate} and falls back on a default if not present.
#'
#' @inheritParams findImprinting
#' @param defaultErr The error rate to use if none can be found.
#' @return A vector of error rates.
getErrorRate = function(gCnts,defaultErr=0.01){
  #Do we have it stored
  if(!'errRate' %in% colnames(mcols(gCnts))){
    warning(sprintf("Error rate not found.  Defaulting to uniform error of %g",defaultErr))
    return(rep(defaultErr,length(gCnts)))
  }else{
    #Error rate should never be 0
    ans = gCnts$errRate
    ans[ans==0] = defaultErr
    return(ans)
  }
}
