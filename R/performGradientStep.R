#' Perform multi-objective gradient descent step.
#'
#' @description
#'   Computes the bi-objective gradient for the two objectives \code{fn1} and \code{fn2}
#'   in position \code{ind}. The bi-objective gradient is the combined vector of the
#'   two single gradients. Note that the step size depends on the length of the combined
#'   gradient vector and thus automatically decreases when approaching an efficient set.
#'
#' @note 
#'   ATTENTION: Only turn off the sanity checks (\code{check.data = FALSE}),
#'   if you can ensure that all input parameters are provided in the correct format.
#'
#' @template arg_ind
#' @template arg_fni
#' @template arg_grad
#' @template arg_scalestep
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_lower
#' @template arg_upper
#' @template arg_checkdata
#' @return [\code{\link{numeric}}(d) | \code{\link{NULL}}]\cr
#'   Returns \code{NULL} if \code{ind} is a local efficient point. Otherwise a numeric
#'   vector of length \code{d} will be returned, which shows the result of a downhill
#'   step in the direction of the multi-objective gradient.
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(0.2, 1))^2)
#' fn2 = function(x) sum(x)
#'
#' # Perform a gradient step:
#' performGradientStep(c(0.3, 0.5), fn1, fn2)
#'
#' # Here, we have found the optimum of fn1:
#' performGradientStep(c(0.2, 1), fn1, fn2)
#' @export
performGradientStep = function(ind, fn1, fn2,
  gradient.list = list(g1 = NULL, g2 = NULL), scale.step = 0.5,
  prec.grad = 1e-6, prec.norm = 1e-6, prec.angle = 1e-4,
  lower, upper, check.data = TRUE) {

  assertNumeric(ind, any.missing = FALSE, null.ok = FALSE)
  d = length(ind)
  if (check.data) {
    assertFunction(fn1)
    assertFunction(fn2)
    assertNumber(scale.step, lower = 0, finite = TRUE, null.ok = FALSE)
    assertNumber(prec.grad, lower = 0, finite = TRUE, null.ok = FALSE)
    assertNumber(prec.norm, lower = 0, finite = TRUE, null.ok = FALSE)
    assertNumber(prec.angle, lower = 0, upper = 180, null.ok = FALSE)
    if (missing(lower)) {
      lower = rep(-Inf, d)
    } else if ((length(lower) == 1L) & (d != 1L)) {
      lower = rep(lower, d)
    }
    if (missing(upper)) {
      upper = rep(Inf, d)
    } else if ((length(upper) == 1L) & (d != 1L)) {
      upper = rep(upper, d)
    }
    assertNumeric(lower, len = d, any.missing = FALSE, null.ok = FALSE)
    assertNumeric(upper, len = d, any.missing = FALSE, null.ok = FALSE)
    assertTRUE(all(lower <= upper))
    assertList(gradient.list, types = c("numeric", "null", "integerish", "double"), any.missing = TRUE, min.len = 1L, null.ok = FALSE)
    if (is.null(gradient.list)) {
      stop("The gradient list needs to consist of p list elements.")
    } else {
      for (k in seq_along(gradient.list)) {
        gradient.list.element = gradient.list[[k]]
        if (!is.null(gradient.list.element)) {
          assertNumeric(gradient.list.element, len = d, any.missing = FALSE, null.ok = FALSE)
        }
      }
    }
  }

  p = length(gradient.list)
  fn.evals = matrix(0L, nrow = 1L, ncol = p)
  names(fn.evals) = sprintf("fn%i.evals", seq_len(p))

  for (k in seq_along(gradient.list)) {
    ## iterate over all gradients from the gradient list
    if (is.null(gradient.list[[k]])) {
      ## if a gradient is not existent, estimate it
      gi = -estimateGradientBothDirections(fn = list(fn1, fn2)[[k]],
        ind = ind, prec.grad = prec.grad, check.data = FALSE, lower = lower, upper = upper)
      gi = normalizeVectorCPP(gi, prec = prec.norm)
      gradient.list[[k]] = gi
      fn.evals[1L, k] = p * d
      if (all(gi == 0)) {
        ## if the gradient is all zero, we can already stop here as
        ## we have found a single-objective local optimum
        return(list(
          offspring = NULL,
          fn.evals = fn.evals,
          gradient.list = gradient.list
        ))
      }
    }
  }

  ## FIXME: this part currently only works for p = 2
  ## if one is between the peaks of fn1 and fn2, move towards a point along
  ## their efficient set
  g1 = gradient.list[[1L]]
  g2 = gradient.list[[2L]]
  angle = computeAngleCPP(vec1 = g1, vec2 = g2, prec = prec.norm)
  if (abs(180 - angle) < prec.angle) {
    # if the angle between both gradients is (approximately) 180 degree,
    # this has to be a local efficient point
    return(list(
      offspring = NULL,
      fn.evals = fn.evals,
      gradient.list = gradient.list
    ))
  } else {
    # FIXME: analyze, which of the following two approaches works better
    # g = scale.step * normalizeVectorCPP(vec = g1 + g2, prec = prec.norm) * (180 - angle) / 180
    g = scale.step * (g1 + g2)
    offspring = assureBoundsCPP(ind = ind, g = g, lower = lower, upper = upper)
    ## note that the gradients of the offspring itself have not (yet) been computed, so
    ## fn.evals does not need to be updated
    return(list(
      offspring = offspring,
      fn.evals = fn.evals,
      gradient.list = gradient.list
    ))
  }
}
