#' @export
generateDesign = function(fn, points.per.dimension=NULL, step.size=NULL) {
  upper = smoof::getUpperBoxConstraints(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  n = smoof::getNumberOfParameters(fn)
  
  if (!is.null(step.size) & !is.null(points.per.dimension)) {
    stop('Not both step.size and points.per.dimension can be set!')
  }
  
  l = list()
  
  if (!is.null(step.size)) {
    # Use step.size
    for (i in 1:n) {
      x.i = c(paste0("x", i))
      l[[x.i]] = seq(lower[i], upper[i], step.size)
    }
  } else {
    # Use points.per.dimension
    # Or calculate a default
    if (length(points.per.dimension) == 0L) {
      # default: approximately 1 Million points in total
      # n = 2: 1001 per dimension
      # n = 3: 101 per dimension
      points.per.dimension = round((1e6 ** (1/n))) + 1
      
      points.per.dimension = rep(points.per.dimension, n)
    } else if (length(points.per.dimension) == 1L) {
      points.per.dimension = rep(points.per.dimension, n)
    }
    
    for (i in 1:n) {
      x.i = c(paste0("x", i))
      l[[x.i]] = seq(lower[i], upper[i], length.out = points.per.dimension[i])
    }
  }
  
  expand.grid(l)
}
