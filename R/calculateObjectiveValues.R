#' @export
calculateObjectiveValues = function(points, fn) {
  n = smoof::getNumberOfObjectives(fn)
  names = sapply(1:n, function(i) paste0("y", i))
  
  obj.space = apply(points, 1, fn)
  obj.space = as.data.frame(t(obj.space))
  names(obj.space) = names
  
  return(obj.space)
}
