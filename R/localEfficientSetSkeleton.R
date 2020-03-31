#' @export
localEfficientSetSkeleton = function(grid, integration="fast") {
  less = list()
  
  cat('1/3 Finding critical points ...\n')
  critical = getCriticalPointsCellCPP(grid$mo.grad, grid$so.grad, grid$div, grid$dims)
  
  sinks = critical$sinks

  ccs = connectedComponentsGrid(sinks, grid$dims)
  valid.ccs.ids = (ccs != 0 & !(ccs %in% as.numeric(names(table(ccs)))[table(ccs) < 4]))
  valid.ccs = ccs[valid.ccs.ids]
  valid.sinks = sinks[valid.ccs.ids]
  
  cat('2/3 Integrating vector field ...\n')
  if (integration %in% c("fast")) {
    integrated = computeCumulatedPathLengths(grid$dec.space, grid$mo.grad, valid.sinks, fix.diagonals = T, prec.vector.length = 0, prec.norm = 0)
  } else {
    integrated = integrateVectorField(grid$mo.grad, grid$dims, valid.sinks)
  }
  
  less$height = integrated$height
  
  integration.sinks = sort(unique(integrated$last.visited))

  ccs = connectedComponentsGrid(integration.sinks, grid$dims)
  valid.ccs.ids = (ccs != 0 & !(ccs %in% as.numeric(names(table(ccs)))[table(ccs) < 4]))
  valid.ccs = ccs[valid.ccs.ids]
  valid.sinks = integration.sinks[valid.ccs.ids]
  
  less$sinks = valid.sinks

  cat('3/3 Calculating basins of attraction ...\n')
  sink.to.basin = rep(-1, prod(grid$dims))
  sink.to.basin[valid.sinks] = valid.ccs
  
  basins = sapply(1:length(integrated$last.visited), function(i) {
    ifelse(integrated$last.visited[i] != -1, sink.to.basin[integrated$last.visited[i]], -1)
  })
  
  less$basins = basins
  
  ridges = changeOfBasin(basins, grid$dims)
  ridges = union(ridges, which(basins == -1))
  
  less$ridges = ridges
  
  return(less)
}