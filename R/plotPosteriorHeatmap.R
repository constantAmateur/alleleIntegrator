#' A region by cell heatmap of posterior probabilities
#'
#' Takes the output of \code{\link{calcStateProbs}} and produces a heatmap showing the posterior probabilities for each cell in each region.  The rows and columns of the heatmap are constructed from \code{factor(pp$cellID)} and \code{factor(pp$regionID)}.  The order can be changed, and extras such as row/column splits provided, by pre-setting these columns of \code{pp} to factors and supplying entries based on their levels.
#'
#' This is essentialy just a smart wrapper around the \code{\link[ComplexHeatmap]{Heatmap}} function, so you can pass any of the same arguments to this function.
#'
#' @param pp Posterior probabilities produces by \code{\link{calcStateProbs}}.
#' @param plotType Either 'nLL' to plot the negative log-likelihoods or 'postProb' for posterior probabilities.
#' @param refState If type 'nLL', plot differences in log-likelihood relative to this state.  Ignored otherwise or if NULL.
#' @param plotStates Plot posterior probability for these states.  NULL uses all states if >2 states or first if 2. 
#' @param useRemoveGroup Should stats be extracted from the group where N have been left out?
#' @param zLims Truncate plot values to this range.  If \code{plotType} is 'postProb', transformed to on the interval (0,1) with an inverse logit transformation.
#' @param colPal Colour scheme to use.  Passed to \code{\link[RColorBrewer]{brewer.pal}}.
#' @param groupCellsBy A vector of the same length as \code{pp} indicating how to group cells.
#' @param ... Extra parameters passed to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @return A plot, obvs.
#' @importFrom ComplexHeatmap Heatmap 
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @export
plotPosteriorHeatmap = function(pp,plotType=c('postProb','nLL'),refState=NULL,plotStates=NULL,useRemoveGroup=FALSE,zLims=c(-5,5),colPal=ifelse(plotType=='postProb','Greens','BrBG'),groupCellsBy=pp$clusterID,...){
  #Parameter checking
  plotType = match.arg(plotType)
  #Construct the base string to search for stuff
  searchBase=paste0(ifelse(useRemoveGroup,'removeN_',''),plotType,'_')
  modelNames = grep(paste0('^',searchBase),colnames(mcols(pp)),value=TRUE)
  modelNames = gsub(paste0('^',searchBase),'',modelNames)
  #Check refState makes sense
  if(!is.null(refState) && !refState %in% plotStates)
    stop(sprintf("refState %s not found in states: %s",refState,paste(plotStates,collapse=',')))
  if(!is.null(plotStates)){
    if(!all(plotStates %in% modelNames))
      stop(sprintf("Not all plotStates found.  These are the available states: %s",paste(modelNames,collapse=', ')))
  }else{
    if(length(modelNames)>2){
      plotStates = modelNames
    }else{
      plotStates = modelNames[1]
      if(plotType=='nLL')
        refState = modelNames[2]
    }
  }
  #Build colour scale
  colFun = suppressWarnings(brewer.pal(100,colPal))
  if(plotType=='postProb')
    zLims = (1+exp(-zLims))**-1
  colFun = colorRamp2(seq(zLims[1],zLims[2],length.out=length(colFun)),colFun)
  #Construct matrix
  cellIDs = factor(pp$cellID)
  regionIDs = factor(pp$regionID)
  mtx = matrix(nrow=nlevels(cellIDs),
               ncol = nlevels(regionIDs)*length(plotStates),
               dimnames = list(levels(cellIDs),rep(levels(regionIDs),each=length(plotStates))))
  #Load the data
  for(i in seq_along(plotStates)){
    plotState = plotStates[i]
    idxs = cbind(cellIDs,regionIDs)
    idxs[,2] = (idxs[,2]-1)*length(plotStates)+i
    #Make changes needed
    tmp = mcols(pp)[,paste0(searchBase,plotState)]
    if(!is.null(refState) & plotType=='nLL')
      tmp = mcols(pp)[,paste0(searchBase,refState)] - tmp
    mtx[idxs] = tmp
  }
  stateLabs = rep(plotStates,nlevels(regionIDs))
  rowSmall = nlevels(cellIDs)<30
  colSmall = nlevels(regionIDs)<30
  #Construct a thing to group cells by if that's given
  if(!is.null(groupCellsBy))
    groupCellsBy = groupCellsBy[match(levels(cellIDs),pp$cellID)]
  params = list(matrix=mtx,
                col=colFun,
                name=ifelse(plotType=='nLL','nLL','prob'),
                column_split = factor(colnames(mtx),levels=levels(regionIDs)), #Split be regions
                row_title_rot = 0, #If we do supply splits, we almost certainly want this
                row_split = groupCellsBy,
                show_row_dend = FALSE, #Dendograms are pretty useless mostly
                show_column_dend = FALSE,
                cluster_column_slices = !colSmall, #Allow any pre-set ordering specified in factor levels to stand for small numbers
                cluster_columns = FALSE,
                column_labels = stateLabs, #Set labels for inside slices 
                cluster_rows = !rowSmall,
                show_row_names=rowSmall, #Only show names if they're aren't a bajillion of them
                show_column_names=colSmall
                )
  #Allow the dots to overwrite default parameters.
  theDots = list(...)
  for(nom in names(theDots))
    params[[nom]] = theDots[[nom]]
  hm = do.call(Heatmap,params)
  return(hm)
}
