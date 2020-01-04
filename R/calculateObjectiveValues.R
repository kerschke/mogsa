#' @export
calculateObjectiveValues = function(points, fn, parallelize = FALSE) {
  n = smoof::getNumberOfObjectives(fn)
  names = paste0("y", 1:n)
  
  if (parallelize) {
    # setup should occur somewhere else? see documentation for future::plan
    future::plan(future::multisession, workers = (future::availableCores() - 1))
    obj.space = future.apply::future_apply(points, 1, fn)
    future::plan(future::sequential)
  } else {
    obj.space = apply(points, 1, fn)
  }
  
  obj.space = t(obj.space)
  colnames(obj.space) = names
  
  return(obj.space)
}
