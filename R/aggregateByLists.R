#' Collapse counts by region and cluster
#'
#' Given a count summarisation GRanges object with a cell like and region like element, aggregate values by vectors/lists of grouping factors.  This is typically something like cluster ID and/or genes, but this function works generally.
#'
#' To aggregate on just cell-like or region-like variables, set the other entry to \code{1}.  That is, \code{cellList=1} and \code{regionList=gCnts$regionID} will aggregate by \code{regionID} ignoring any \code{cellIDs}.
#'
#' @inheritParams findImprinting
#' @param assays Columns in \code{gCnts} to aggregate over.
#' @param cellList Cell-like thing to aggregate over.  List or vector of same length as \code{gCnts}.
#' @param regionList Cell-like thing to aggregate over.  List or vector of same length as \code{gCnts}.
#' @param aggFun Function used for aggregation.
#' @return A data.frame giving the combination of cellID and regionID from \code{cellList} and \code{regionList} and aggregated information in \code{assays}.
#' @export
aggregateByLists = function(gCnts,assays,cellList=gCnts$cellID,regionList=gCnts$regionID,aggFun=sum){
  gCnts$cellList = cellList
  gCnts$regionList = regionList
  #Expand by cell like
  tmp = gCnts[rep(seq_along(gCnts),lengths(gCnts$cellList))]
  tmp$cellList = unlist(gCnts$cellList)
  #Expand by region like
  t2 = tmp[rep(seq_along(tmp),lengths(tmp$regionList))]
  t2$regionList = unlist(tmp$regionList)
  #Make the mark
  ans = as.data.frame(mcols(t2)[,assays])
  ans = aggregate(ans,by=list(cellID = t2$cellList,regionID = t2$regionList),FUN=aggFun)
  return(ans)
}
