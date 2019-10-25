#' Performs a weighted bi-section optimization.
#'
#' @description
#'   Weighted version of the bisection optimization method. Given two points \code{x1} and \code{x2} on opposite
#'   sides of the optimum, this optimizer iteratively splits the interval [\code{x1}, \code{x2}] into two parts
#'   [\code{x1}, \code{x.new}] and [\code{x.new}, \code{x2}] and proceeds with the interval, whose boundaries are
#'   still located on opposite sides of the optimum. Instead to the \href{https://en.wikipedia.org/wiki/Bisection_method}{classical bisection method},
#'   where \code{x.new} is the arithmetic mean of \code{x1} and \code{x2}, this version uses the lengths of the
#'   bi-objective gradients in \code{x1} and \code{x2} to compute a more promising cut-point \code{x.new}.
#'
#' @param x1 [\code{\link{numeric}(d)}]\cr
#'   d-dimensional individual located on one side of the (bi-objective) optimum.
#' @param x2 [\code{\link{numeric}(d)}]\cr
#'   d-dimensional individual located on the opposite side (w.r.t. x1) of the (bi-objective) optimum.
#' @template arg_fni
#' @template arg_gri
#' @template arg_precgrad
#' @template arg_precnorm
#' @param max.steps [\code{\link{integer}(1L)}]
#'   Maximum number of allowed bi-section steps to reach an optimum. The default is \code{1000L}.
#' @template arg_lower
#' @template arg_upper
#' @return [\code{\link{list}(4L)}]\cr
#'   List containing a matrix (\code{opt.path}) with the individuals along the optimization path,
#'   the corresponding number of function evaluations (\code{fn.evals}), the single-objective
#'   gradients of the last individual (\code{gradient.list}) and a flag, indicating whether the
#'   optimizer found a local optimum.
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(2, 0))^2)
#' fn2 = function(x) sum((x - c(0, 1))^2)
#' fn = function(x) return(c(fn1(x), fn2(x)))
#' 
#' # Visualize locally efficient set, i.e., the "area" where we ideally want to find a point:
#' plot(c(2, 0), c(0, 1), type = "o", pch = 19,
#'   xlab = expression(x[1]), ylab = expression(x[2]), las = 1, asp = 1)
#' text(2, 0, "Optimum of fn1", pos = 2, offset = 1.5)
#' text(0, 1, "Optimum of fn2", pos = 4, offset = 1.5)
#' 
#' # Place two points x1 and x2 on opposite sides of the bi-objective optimum:
#' x1 = c(1, 1)
#' x2 = c(0.5, 0)
#' points(rbind(x1, x2), pch = 19, type = "o", lty = "dotted")
#' text(rbind(x1, x2), labels = c("x1", "x2"), pos = 4)
#' 
#' # Optimize using weighted bisection optimization:
#' opt.path = performWeightedBisectionOptimization(x1 = x1, x2 = x2, fn = fn)$opt.path
#' 
#' # Visualize the optimization path:
#' points(opt.path)
#' 
#' # Highlight the found local efficient point (= local optimum w.r.t. both objectives):
#' n = nrow(opt.path)
#' points(opt.path[n, 1], opt.path[n, 2], pch = 4, col = "red", cex = 2)
#' text(opt.path[n, 1], opt.path[n, 2], "Found Local Efficient Point", pos = 4, offset = 1.5)
#' @export
performWeightedBisectionOptimization = function(x1, x2, fn, g1 = NULL, g2 = NULL,
  prec.grad = 1e-6, prec.norm = 1e-6, max.steps = 1000L, lower, upper) {

  if (missing(lower)) {
    lower = ifelse(x1 < x2, x1, x2) - prec.grad
  }
  if (missing(upper)) {
    upper = ifelse(x1 > x2, x1, x2) + prec.grad
  }
  ## initialize optimization path
  opt.path = rbind(x1, x2)

  ## FIXME: currently only fn1 and fn2 supported (hard-coded)
  p = 2L
  sp = seq_len(p)
  fn.evals = matrix(0L, nrow = 2L, ncol = p)
  d = length(x1)
  ## approximate MO gradients (if not provided) per objective
  if (is.null(g1)) {
    # TODO fix fn1/fn2 split here
    g = estimateGradientBothDirections(
      fn = fn, ind = x1, prec.grad = prec.grad, lower = lower, upper = upper)
    
    v1 = normalizeVectorCPP(
      vec = g[1,],
      prec = prec.norm
    )
    v2 = normalizeVectorCPP(
      vec = g[2,],
      prec = prec.norm
    )
    
    fn.evals[1L, seq_len(p)] = p * d
    g1 = -(v1 + v2)
  }
  if (is.null(g2)) {
    g = estimateGradientBothDirections(
      fn = fn, ind = x2, prec.grad = prec.grad, lower = lower, upper = upper)
    
    v1 = normalizeVectorCPP(
      vec = g[1,],
      prec = prec.norm
    )
    v2 = normalizeVectorCPP(
      vec = g[2,],
      prec = prec.norm
    )
    
    fn.evals[2L, seq_len(p)] = p * d
    g2 = -(v1 + v2)
  }

  vl1 = computeVectorLengthCPP(g1)
  vl2 = computeVectorLengthCPP(g2)

  found.optimum = FALSE
  ## for a maximum of max.steps iterations cut the interval [x1, x2] into two parts
  ## using the lengths of the bi-objective gradients as weights
  for (i in seq_len(max.steps)) {
    ## compute weighted split-point
    stepTox2 = (x2 - x1) * (vl1 / (vl1 + vl2))
    if (computeVectorLengthCPP(stepTox2) < prec.norm / 100) {
      ## if the step is "too short", we can stop here
      break
    }
    x.new = x1 + stepTox2

    ## add split-point to optimization path
    opt.path = rbind(opt.path, x.new)
    fn.evals = rbind(fn.evals, 0L)

    ## estimate single-objective gradients in x.new
    i = nrow(fn.evals)
    
    g = estimateGradientBothDirections(
      fn = fn, ind = x.new, prec.grad = prec.grad, lower = lower, upper = upper)
    
    v1 = normalizeVectorCPP(
      vec = g[1,],
      prec = prec.norm
    )
    v2 = normalizeVectorCPP(
      vec = g[2,],
      prec = prec.norm
    )
    
    fn.evals[i, 1L] = fn.evals[i, 1L] + p * d
    fn.evals[i, 2L] = fn.evals[i, 2L] + p * d

    if (computeVectorLengthCPP(v1) == 0 || computeVectorLengthCPP(v2) == 0) {
      ## local optimum for one of the objectives
      found.optimum = TRUE
      break
    }

    ## compute corresponding bi-objective gradient (and its length)
    g.new = -(v1 + v2)
    vl.new = computeVectorLengthCPP(g.new)

    if (vl.new < prec.norm) {
      ## if the length of the computed gradient is "close enough" to the
      ## optimum, a local efficient point was found
      found.optimum = TRUE
      break
    } else {
      ## otherwise, decide on the interval to proceed with
      angle = computeAngleCPP(g1, g.new, prec.norm)
      if (angle < 90) {
        ## g1 and g.new show in a "similar" direction
        ## continue with [x.new, x2]
        x1 = x.new
        g1 = g.new
        vl1 = vl.new
      } else {
        ## g1 and g.new show in an "opposite" direction
        ## continue with [x1, x.new]
        x2 = x.new
        g2 = g.new
        vl2 = vl.new
      }
    }
  }

  rownames(opt.path) = NULL
  rownames(fn.evals) = NULL
  return(list(
    opt.path = opt.path,
    fn.evals = fn.evals,
    gradient.list = list(g1 = -v1, g2 = -v2),
    found.optimum = found.optimum))
}
