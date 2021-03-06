#' Approximate gradient of a single-objective function.
#'
#' @description
#'   \code{estimateGradient} approximates the gradient of a single-objective function
#'   \code{fn} in position \code{x} based on small changes in one direction per
#'   dimension of \code{x}.
#'   \deqn{\nabla f(\mathbf{x}) = \frac{f(\mathbf{x} + \boldsymbol{\varepsilon}) - f(\mathbf{x})}{||\boldsymbol{\varepsilon}||}}
#'   The gradient approximation is performed separately for each dimension of the search space.
#'   That is, the i-th element of the gradient results from a slight step (step size: \code{prec.grad})
#'   of the i-th element of \code{x}, while all remaining elements of \code{x} are not altered.\cr
#'   
#'   In contrast to the aforementioned function, \code{estimateGradientBothDirections}
#'   performs a small step in both directions (of the i-th element) of \code{x}.
#'   \deqn{\nabla f(\mathbf{x}) = \frac{f(\mathbf{x} + \boldsymbol{\varepsilon}) - f(\mathbf{x} - \boldsymbol{\varepsilon})}{2 \cdot ||\boldsymbol{\varepsilon}||}}
#'
#' @note 
#'   This function basically is a slightly modified version of
#'   \code{numDeriv::grad(..., method = "simple")} and was mainly developed for
#'   internal usage to have a speed up over the aforementioned version. However,
#'   this speed up will only have an effect, if you call this function very
#'   frequently and if you turn off the sanity checks.\cr
#'
#'   ATTENTION: Only turn off the sanity checks (\code{check.data = FALSE}),
#'   if you can ensure that all input parameters are provided in the correct format.
#'
#' @param fn [\code{\link{function}}]\cr
#'   Single-objective function, whose gradient will be approximated.
#' @template arg_ind
#' @param side [\code{\link{logical}(d)}]\cr
#'   Logical vector of the same length as \code{x}, stating per element whether to
#'   approximate the gradient into the positive (\code{TRUE}) or negative (\code{FALSE})
#'   direction of \code{x}. The default is \code{rep(TRUE, length(x))}.
#' @template arg_precgrad
#' @template arg_checkdata
#' @template arg_lower
#' @template arg_upper
#' @param ... Further arguments to be passed to \code{fn}.
#' @return [\code{\link{numeric}(d)}]
#' @name estimateGradient
#' @rdname estimateGradient
#' @examples
#' fn = function(x) sum(x^2)
#' estimateGradientSingleDirection(fn, c(0.3, 0.5))
#' estimateGradientSingleDirection(fn, c(0.3, 0.5), side = rep(FALSE, 2L))
#' estimateGradientBothDirections(fn, c(0.3, 0.5))
NULL

#' @rdname estimateGradient
#' @export
estimateGradientSingleDirection = function(fn, ind, side = NULL, prec.grad = 1e-6, check.data = TRUE, lower, upper, ...) {
  
  if (missing(lower)) {
    lower = ind - prec.grad
  }
  if (missing(upper)) {
    upper = ind + prec.grad
  }
  if (check.data) {
    assertNumeric(ind, any.missing = FALSE, null.ok = FALSE)
    assertFunction(fn, null.ok = FALSE)
  }
  f = fn(ind, ...)
  d = length(ind)
  if (is.null(side)) {
    side = rep(TRUE, d)
  }
  ## automatically adjust sides ensuring to stay in bounds
  side[ind < lower + prec.grad] = TRUE
  side[ind > upper - prec.grad] = FALSE
  if (check.data) {
    assertNumber(f, na.ok = FALSE, null.ok = FALSE)
    assertLogical(side, any.missing = FALSE, len = d)
    assertNumber(prec.grad, na.ok = FALSE, lower = 0, null.ok = FALSE)
  }
  precision.vector = ifelse(side, prec.grad, -prec.grad)
  grad = vapply(seq_len(d), function(i) {
    x = ind
    x[i] = x[i] + precision.vector[i]
    # scale by precision.vector[i] rather than prec.grad to account for the direction of the difference in the numerator
    (fn(x) - f) / precision.vector[i]
  }, double(1L))
  return(grad)
}

#' @rdname estimateGradient
#' @export
estimateGradientBothDirections = function(fn, ind, prec.grad = 1e-4, check.data = TRUE, lower, upper, ...) {
  
  if (missing(lower)) {
    lower = ind - prec.grad
  }
  if (missing(upper)) {
    upper = ind + prec.grad
  }
  if (check.data) {
    assertNumeric(ind, any.missing = FALSE, null.ok = FALSE)
    assertFunction(fn, null.ok = FALSE)
    assertNumber(prec.grad, na.ok = FALSE, lower = 0, null.ok = FALSE)
  }
  d = length(ind)
  if (any(ind < lower + prec.grad) || any(ind > upper - prec.grad)) {
    ## for those values, where the individual is to close to the lower bound,
    ## compute the gradient in positive direction (otherwise in negative)
    grad1 = estimateGradientSingleDirection(fn = fn, ind = ind,
      side = rep(FALSE, d), prec.grad = prec.grad, check.data = FALSE,
      lower = lower, upper = upper)
    ## for those values, where the individual is to close to the upper bound,
    ## compute the gradient in negative direction (otherwise in positive)
    grad2 = estimateGradientSingleDirection(fn = fn, ind = ind,
      side = rep(TRUE, d), prec.grad = prec.grad, check.data = FALSE,
      lower = lower, upper = upper)
    grad = (grad1 + grad2) / 2L
  } else {
    grad = vapply(seq_len(d), function(i) {
      x.pos = ind
      x.neg = ind
      x.pos[i] = x.pos[i] + prec.grad
      x.neg[i] = x.neg[i] - prec.grad
      (fn(x.pos) - fn(x.neg)) / (2 * prec.grad)
    }, double(1L))
  }
  return(grad)
}
