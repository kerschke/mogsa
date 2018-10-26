#' Display MOGSA's Optimization Path as a HTML or LaTeX table.
#'
#' @description
#'   Convenience function for displaying the optimization path of MOGSA as
#'   a prettified (HTML or LaTeX) table.
#'
#' @param opt.path [\code{\link{data.frame}}]\cr
#'   Optimization path returned by \code{\link{runMOGSA}}.
#' @param format [\code{\link{character}(1L)}]\cr
#'   Type of table being generated; either \code{"html"} (default) or \code{"latex"}.
#' @param ... [any]\cr
#'   Further arguments to be passed to \code{kable}.
#' @return Returns MOGSA's optimization path either as HTML (default) or LaTeX table.
#' @examples
#' ## Example 1:
#' # (i) Define two single-objective test problems (fn1, fn2) and a starting point (ind):
#' fn1 = function(x) sum((x - c(2, 0))^2)
#' fn2 = function(x) sum((x - c(0, 1))^2)
#' ind = c(1, 1)
#' 
#' # (ii) Run MOGSA:
#' mogsa.result = runMOGSA(ind = ind, fn1 = fn1, fn2 = fn2)
#'
#' # (iii) Display the optimization path of MOGSA as HTML (or LaTeX) table:
#' displayOptPath(mogsa.result)
#' @export
displayOptPath = function(opt.path, format = "html", ...) {
  assertDataFrame(opt.path, any.missing = FALSE, min.rows = 1L, min.cols = 5L)
  assertSubset(c("x1", "x2", "iterations", "basin", "type"), choices = colnames(opt.path))
  assertChoice(format, choices = c("html", "latex"), null.ok = FALSE)
  rownames(opt.path) = NULL
  op = opt.path[, c("iterations", "x1", "x2")]
  colnames(op)[1L] = "Iterations"
  kable.result = kable(x = op, format = format, ...)
  # the following three formats are not able to handle further highlighting
  kable.result = kable_styling(kable_input = kable.result, bootstrap_options = "striped", full_width = TRUE)
  for (current.basin in sort(unique(opt.path[, "basin"]))) {
    # highlight all basins
    basin.index = which(opt.path$basin == current.basin)
    kable.result = group_rows(kable_input = kable.result,
      group_label = sprintf("Basin %i", current.basin),
      basin.index[1L], basin.index[length(basin.index)])
    if (current.basin == 1) {
      # an initial point should only exist in the first phase
      kable.result = group_rows(kable_input = kable.result,
        group_label = "Initialization", 1L, 1L)
    }
    # find the positions of the local search results (within the current basin)
    ls.index = basin.index[opt.path$type[basin.index] == "local search"]
    if (length(ls.index) > 0) {
      kable.result = group_rows(kable_input = kable.result,
        group_label = "Local Search", ls.index[1L], ls.index[length(ls.index)])
    }
    # find the positions of the exploration steps (within the current basin)
    exp.index = basin.index[opt.path$type[basin.index] == "exploration"]
    if (length(exp.index) > 0) {
      kable.result = group_rows(kable_input = kable.result,
        group_label = "Exploration", exp.index[1L], exp.index[length(exp.index)])
    }
    # list the external points within the current basin
    ext.index = basin.index[opt.path$type[basin.index] == "external"]
    if (length(ext.index) > 0) {
      kable.result = group_rows(kable_input = kable.result,
        group_label = "External", ext.index[1L], ext.index[length(ext.index)])
    }
  }
  return(kable.result)
}
