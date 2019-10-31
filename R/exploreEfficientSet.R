#' Explore locally efficient set.
#'
#' @description
#'   Explores a locally efficient set separately along the gradients of the two objectives
#'   \code{fn1} and \code{fn2}. Both explorations initially start in point \code{ind}. Such
#'   explorations obviously are only reasonable, if \code{ind} is a locally efficient point.
#'   See \code{\link{findLocallyEfficientPoint}} for further details.
#'
#' @note 
#'   ATTENTION: Only turn off the sanity checks (\code{check.data = FALSE}),
#'   if you can ensure that all input parameters are provided in the correct format.
#'
#' @template arg_ind
#' @template arg_fni
#' @template arg_grad
#' @template arg_max_no_steps_exploration
#' @template arg_explorationstep
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_lower
#' @template arg_upper
#' @template arg_checkdata
#' @template arg_show_info
#' @return [\code{\link{list}(4)}]\cr
#'   Returns a list with four elements. The first one provides the points that were found
#'   to be part of the efficient set.\cr
#'   The second element lists the corresponding number of performed function evaluations.\cr
#'   The third and fourth elements are lists, containing information on the edges of the
#'   efficient set. This information consists of a \code{\link{logical}} indicating whether
#'   it could be confirmed that the found efficient set definitely is only a local efficient
#'   set. Important note: a value of \code{FALSE} does not imply that the found efficient set
#'   has to be globally efficient, as it could also be a multi-objective trap. In addition
#'   those two sublists also contain the actual positions of the evaluated points outside
#'   the efficient set, their single-objective gradients as well as the number of function
#'   evaluations (per objective) performed to compute these gradients.\cr
#'   Note that these (up to two) points (\code{end1$external} and \code{end2$external})
#'   are reasonable starting points for a gradient descent towards a better efficient set.
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(0.2, 1))^2)
#' fn2 = function(x) 2 * x[1]^2 - x[1] * x[2] + 0.3 * x[2]^2
#'
#' # c(0.2, 1) is obviously a local optimum of fn1, so let's explore the efficient set from there:
#' exploreEfficientSet(c(0.2, 1), fn1, fn2, max.no.steps.exploration = 50L, exploration.step = 0.05)
#' @export
exploreEfficientSet = function(ind, fn, gradient.mat,
  max.no.steps.exploration = 400L, exploration.step = 0.2, prec.grad = 1e-6,
  prec.norm = 1e-6, prec.angle = 1e-4, check.data = TRUE, show.info = TRUE) {

  assertNumeric(ind, any.missing = FALSE, null.ok = FALSE)
  if (check.data) {
    assertFunction(fn)
    assertNumber(exploration.step, lower = 0, finite = TRUE, null.ok = FALSE)
    assertNumber(prec.grad, lower = 0, finite = TRUE, null.ok = FALSE)
    assertNumber(prec.norm, lower = 0, finite = TRUE, null.ok = FALSE)
    # assertList(gradient.mat, types = c("numeric", "null", "integerish", "double"),
    #   any.missing = TRUE, min.len = 1L, null.ok = FALSE)
    # if (is.null(gradient.mat)) {
    #   stop("The gradient list needs to consist of p list elements.")
    # }
    assertLogical(show.info, any.missing = FALSE, len = 1L, null.ok = FALSE)
    # ensure that ind is locally efficient
    if (show.info && !is.null(performGradientStep(ind = ind, fn = fn, gradient.mat = gradient.mat,
      prec.grad = prec.grad, prec.norm = prec.norm, prec.angle = prec.angle,
      check.data = FALSE)$offspring)) {
      catf("ATTENTION: The exploration of the efficient set was started from a locally non-efficient point.")
    }
  }

  ## initialize
  d = getNumberOfParameters(fn)
  p = getNumberOfObjectives(fn)
  lower = getLowerBoxConstraints(fn)
  upper = getUpperBoxConstraints(fn)

  fn.evals = matrix(0L, nrow = 1L, ncol = 1L)
  opt.path = matrix(ind, nrow = 1L)

  ## FIXME: hard-coded for p = 2
  ## single-objective gradients of the very first individual
  if (nrow(gradient.mat) >= 2) {
    g1 = gradient.mat[1,]
    g2 = gradient.mat[2,]
  }

  # compute gradients if missing
  if (is.null(g1) || is.null(g2)) {
    g = -estimateGradientBothDirections(fn = fn, ind = ind,
      prec.grad = prec.grad, check.data = FALSE,
      lower = lower, upper = upper)
    g1 = g[1,]
    g2 = g[2,]
    
    fn.evals[1L,1L] = p * d
  }
  
  grad1 = g1 = normalizeVectorCPP(vec = g1, prec = prec.norm)

  grad2 = g2 = normalizeVectorCPP(vec = g2, prec = prec.norm)

  ## initialize external points per objective, their gradient lists,
  ## as well as their costs
  external1 = external2 = NULL
  ext1.gradient.mat = rbind(g1, g2)
  ext2.gradient.mat = ext1.gradient.mat
  fn.evals.ext1 = fn.evals.ext2 = NULL

  ## explore along fn1
  is.local1 = FALSE
  for (i in seq_len(max.no.steps.exploration)) {
    parent = opt.path[nrow(opt.path),]
    offspring = assureBoundsCPP(ind = parent,
      g = exploration.step * normalizeVectorCPP(vec = grad1, prec = prec.norm),
      lower = lower, upper = upper)
    if (identical(parent, offspring)) {
      ## if the offspring and parent are identical, we can't move further into
      ## the direction of fn1
      break
    }

    ## if they are not identical, we have to figure out, whether the offspring
    ## is (a) a point along the efficient set, (b) a point outside the
    ## efficient set (pointing back towards it), or (c) a point in an adjacent
    ## basin
    g.off = -estimateGradientBothDirections(fn = fn, ind = offspring,
      prec.grad = prec.grad, check.data = FALSE,
      lower = lower, upper = upper)
    
    g1.off = g.off[1,]
    g1.off = normalizeVectorCPP(vec = g1.off, prec = prec.norm)
    
    g2.off = g.off[2,]
    g2.off = normalizeVectorCPP(vec = g2.off, prec = prec.norm)
    
    if (computeVectorLengthCPP(g1.off) == 0 || computeVectorLengthCPP(g2.off) == 0) {
      ## if the offspring's gradient in the direction of fn1 or fn2 is 0, we must have
      ## reached a local optimum and can't proceed
      opt.path = rbind(opt.path, offspring)
      
      ## costs for computing g.off, attributed to this step in this case
      fn.evals = rbind(fn.evals, 0L)
      fn.evals[nrow(fn.evals),1L] = p * d
      
      break
    }

    angle.offspring = computeAngleCPP(vec1 = g1.off, vec2 = g2.off, prec = prec.norm)
    angle.grad1 = computeAngleCPP(g1.off, grad1, prec = prec.norm)
    if ((angle.offspring > 90) & (angle.grad1 < 90)) {
      ## as long as the offspring's gradients show in opposite directions
      ## (> 90 degree), as well as the offspring's current and its preceeding
      ## individual's gradients in the direction of fn1 (angle.grad1) show in
      ## similar directions (< 90 degree), we can move on
      opt.path = rbind(opt.path, offspring)
      
      ## costs for computing g.off, attributed to this step in this case
      fn.evals = rbind(fn.evals, 0L)
      fn.evals[nrow(fn.evals),1L] = p * d

      ## replace the "old" gradient in direction of fn1 by the current
      ## one and continue with the next iteration of the loop
      grad1 = g1.off
      next
    } else {
      ## now, the offspring must have been an external point and we can store
      ## it as such
      external1 = offspring
      
      # costs for gradient, attributed to external point
      fn.evals.ext1 = matrix(0L, nrow = 1L, ncol = 1L)
      fn.evals.ext1[1L,] = p * d
      names(fn.evals.ext1) = c("fn.evals")
      
      if (angle.grad1 > 90) {
        ## if the direction of the gradients (of fn1) have changed between the
        ## offspring and its predecessor, we must have left the efficient set
        ## (but are still in the same basin)
        break
      } else {
        ## otherwise, we must even have passed a ridge into a new basin
        is.local1 = TRUE
        ext1.gradient.mat = g.off
        break
      }
    }
  }

  ## perform the same steps to explore along fn2
  is.local2 = FALSE
  for (i in seq_len(max.no.steps.exploration)) {
    if (i == 1) {
      ## let's *start* the exploration with the original individual
      parent = opt.path[1L,]
    } else {
      parent = opt.path[nrow(opt.path),]
    }
    offspring = assureBoundsCPP(ind = parent,
      g = exploration.step * normalizeVectorCPP(vec = grad2, prec = prec.norm),
      lower = lower, upper = upper)

    if (identical(parent, offspring)) {
      ## if the offspring and parent are identical, we can't move further into
      ## the direction of fn2
      break
    }
    
    g.off = -estimateGradientBothDirections(fn = fn, ind = offspring,
                                            prec.grad = prec.grad, check.data = FALSE,
                                            lower = lower, upper = upper)

    g1.off = g.off[1,]
    g1.off = normalizeVectorCPP(vec = g1.off, prec = prec.norm)
    
    g2.off = g.off[2,]
    g2.off = normalizeVectorCPP(vec = g2.off, prec = prec.norm)

    if (computeVectorLengthCPP(g2.off) == 0 || computeVectorLengthCPP(g1.off) == 0) {
      ## if the offspring's gradient in the direction of fn2 is 0, we must have
      ## reached a local optimum of fn2 and can't proceed
      opt.path = rbind(opt.path, offspring)
      
      ## costs for computing g.off, attributed to this step in this case
      fn.evals = rbind(fn.evals, 0L)
      fn.evals[nrow(fn.evals),] = p * d
      
      break
    }

    angle.offspring = computeAngleCPP(vec1 = g1.off, vec2 = g2.off, prec = prec.norm)
    angle.grad2 = computeAngleCPP(g2.off, grad2, prec = prec.norm)
    if ((angle.offspring > 90) & (angle.grad2 < 90)) {
      ## as long as the offspring's gradients show in opposite directions
      ## (> 90 degree), as well as the offspring's current and its preceeding
      ## individual's gradients in the direction of fn2 (angle.grad2) show in
      ## similar directions (< 90 degree), we can move on
      opt.path = rbind(opt.path, offspring)
      
      ## costs for computing g.off, attributed to this step in this case
      fn.evals = rbind(fn.evals, 0L)
      fn.evals[nrow(fn.evals),] = p * d

      ## replace the "old" gradient in direction of fn2 by the current
      ## one and continue with the next iteration of the loop
      grad2 = g2.off
      next
    } else {
      ## now, the offspring must have been an external point and we can store
      ## it as such
      external2 = offspring
      
      # costs for gradient, attributed to external point
      fn.evals.ext2 = matrix(0L, nrow = 1L, ncol = 1L)
      fn.evals.ext2[1L,] = p * d
      names(fn.evals.ext2) = c("fn.evals")
      
      if (angle.grad2 > 90) {
        ## if the direction of the gradients (of fn2) have changed between the
        ## offspring and its predecessor, we must have left the efficient set
        ## (but are still in the same basin)
        break
      } else {
        ## otherwise, we must even have passed a ridge into a new basin
        is.local2 = TRUE
        ext2.gradient.mat = g.off
        break
      }
    }
  }

  rownames(opt.path) = NULL
  rownames(fn.evals) = NULL

  return(list(
    opt.path = opt.path,
    fn.evals = fn.evals,
    end1 = list(
      is.local = is.local1,
      external = external1,
      gradient.mat = ext1.gradient.mat,
      evals = fn.evals.ext1
    ),
    end2 = list(
      is.local = is.local2,
      external = external2,
      gradient.mat = ext2.gradient.mat,
      evals = fn.evals.ext2
    )
  ))
}
