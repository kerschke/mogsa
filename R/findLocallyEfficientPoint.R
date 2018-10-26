#' Find locally efficient point.
#'
#' @description
#'   Performs sequential downhill steps in the direction of the multi-objective gradient
#'   until a locally efficient point was found. The latter possess a multi-objective
#'   gradient of length zero (see \code{\link{performGradientStep}} for details) as the
#'   normalized gradients of its objectives cancel each other out.\cr
#'   Note that once the local search detects that it crossed the optimum (i.e., the
#'   multi-objective gradients of succeeding iterations show in opposite directions),
#'   it will perform a weighted bisection optimization to find a local efficient point.
#'   See \code{\link{performWeightedBisectionOptimization}} for further details on the
#'   latter.
#'
#' @note 
#'   ATTENTION: Only turn off the sanity checks (\code{check.data = FALSE}),
#'   if you can ensure that all input parameters are provided in the correct format.
#'
#' @template arg_ind
#' @template arg_fni
#' @template arg_grad
#' @template arg_max_no_steps_ls
#' @template arg_scalestep
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_ls_method
#' @template arg_lower
#' @template arg_upper
#' @template arg_checkdata
#' @template arg_show_info
#' @template arg_allow_restarts
#' @return [\code{\link{list}(3L)}]\cr
#'   Returns a list, whose first element (\code{opt.path}) provides the matrix
#'   of the optimization path, which has been evaluated on the way to a locally
#'   efficient point. As the 1st row represents the initial individual
#'   (\code{ind}), the i-th row corresponds to the (i-1)-th individual.\cr
#'   The second element of the list (\code{fn.evals}) contains a matrix with
#'   the corresponding number of function evaluations per objective. As a
#'   result, this matrix consists of \code{p} columns and as many rows as
#'   \code{opt.path}.\cr
#'   The final list element (\code{gradient.list}) provides another list, whose
#'   \code{p} elements provide the single-objective gradients of the last
#'   individual from \code{opt.path}. This information might be used to safe
#'   function evaluations when proceeding with that particular individual
#'   (e.g., when exploring a potential existing efficient set).
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(2, 0))^2)
#' fn2 = function(x) sum((x - c(0, 1))^2)
#' 
#' # Visualize locally efficient set, i.e., the "area" where we ideally want to find a point:
#' plot(c(2, 0), c(0, 1), type = "o", pch = 19,
#'   xlab = expression(x[1]), ylab = expression(x[2]), las = 1, asp = 1)
#' text(2, 0, "Optimum of fn1", pos = 2, offset = 1.5)
#' text(0, 1, "Optimum of fn2", pos = 4, offset = 1.5)
#' 
#' # Place two points x1 and x2 on opposite sides of the bi-objective optimum:
#' x.start = c(0.3, 0.5)
#' points(rbind(x.start), pch = 19, type = "o", lty = "dotted")
#' text(rbind(x.start), labels = c("start"), pos = 2)
#'
#' # Visualize path of bisection search in blue
#' res.bisection = findLocallyEfficientPoint(c(0.3, 0.5), fn1, fn2, ls.method = "bisection")
#' points(res.bisection$opt.path, pch = 20, lty = 2, type = "o", col = "blue")
#'
#' # Visualize path of multi-objective local search in red
#' res.mo.ls = findLocallyEfficientPoint(c(0.3, 0.5), fn1, fn2, ls.method = "mo-ls")
#' points(res.mo.ls$opt.path, pch = 20, lty = 2, type = "o", col = "red")
#' @export
findLocallyEfficientPoint = function(ind, fn1, fn2, gradient.list = list(g1 = NULL, g2 = NULL),
  max.no.steps.ls = 500L, scale.step = 0.5, prec.grad = 1e-6, prec.norm = 1e-6, prec.angle = 1e-4,
  ls.method = "both", lower, upper, check.data = TRUE, show.info = TRUE, allow.restarts = TRUE) {

  d = length(ind)
  # perform sanity checks
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
    assertChoice(ls.method, c("both", "bisection", "mo-ls"), null.ok = FALSE)
    assertLogical(show.info, any.missing = FALSE, len = 1L, null.ok = FALSE)
    assertList(gradient.list, types = c("numeric", "null", "integerish", "double"), any.missing = TRUE, min.len = 1L, null.ok = FALSE)
    if (is.null(gradient.list)) {
      stop("The gradient list needs to consist of p list elements.")
    }
    assertLogical(show.info, any.missing = FALSE, len = 1L, null.ok = FALSE)
  }

  ## Initialization:
  p = length(gradient.list)
  sp = seq_len(p)
  fn.evals = matrix(0L, nrow = 1L, ncol = p)
  colnames(fn.evals) = sprintf("fn%i.evals", seq_len(p))
  opt.path = matrix(ind, nrow = 1L)
  restart.counter = 0L
  j = 0L # counter for the iterations since last restart

  while (TRUE) {
    i = nrow(opt.path)
    j = j + 1L

    ## generate offspring
    individual = opt.path[nrow(opt.path),]
    gradient.step = performGradientStep(ind = individual, fn1 = fn1, fn2 = fn2, gradient.list = gradient.list,
      scale.step = scale.step, prec.grad = prec.grad, prec.norm = prec.norm, prec.angle = prec.angle,
      lower = lower, upper = upper, check.data = check.data)
    offspring = gradient.step$offspring

    ## update function evaluations
    fn.evals[i, sp] = fn.evals[i, sp] + gradient.step$fn.evals[1L, sp]

    if (is.null(offspring)) {
      ## if there is no offspring (= NULL), the current
      ## individual is already a locally efficient point
      gradient.list = gradient.step$gradient.list
      break
    }

    if (computeVectorLengthCPP(individual - offspring) < prec.norm) {
      if (allow.restarts) {
        ## there was basically no improvement in the step from
        ## the individual to its offspring --> resample
        offspring = runif(n = d, min = lower, max = upper)
        opt.path = rbind(opt.path, offspring)
        fn.evals = rbind(fn.evals, 0L)
        gradient.list = vector(mode = "list", length = p)
        names(gradient.list) = sprintf("g%i", seq_len(p))
        j = 0L
        i = nrow(opt.path)
        if (i == (max.no.steps.ls + 1L)) {
          break
        } else {
          next
        }
      } else {
        break
      }
    }

    ## now, we can safely add the freshly generated offspring to opt.path and
    ## update the function evaluations
    opt.path = rbind(opt.path, offspring)
    fn.evals = rbind(fn.evals, 0L)
    i = nrow(opt.path)

    if (i == (max.no.steps.ls + 1L)) {
      ## if we have reached the maximum number of allowed
      ## function evaluations, we can stop the loop here
      break
    }

    if (j == 1L) {
      ## if this is still the first loop after a restart, we can move on to the
      ## next iteration of the loop; but should first reset the gradient list
      if (i == 2L) {
        gradient.list = vector(mode = "list", length = p)
        names(gradient.list) = sprintf("g%i", seq_len(p))
      }
      next
    }

    ## now, let's check whether the algorithm remains on
    ## the same side of an optimum
    x1 = opt.path[i - 2L,]
    x2 = opt.path[i - 1L,]
    x3 = offspring
    v1 = x2 - x1
    v2 = x3 - x2
    angle = computeAngleCPP(vec1 = v1, vec2 = v2, prec = prec.norm / 1000)
    if (angle < 90) {
      ## if the angle between the two succeeding steps (x1 -> x2 and
      ## x2 -> x3) is smaller than 90 degree, the algorithm still
      ## moves towards the same efficient set --> use x3 as start
      ## for the next loop
      next
    }

    if (ls.method %in% c("both", "bisection")) {
      ## FIXME: double-check whether max.steps is set correctly

      bisect.opt.result = performWeightedBisectionOptimization(
        x1 = x1, x2 = x2, fn1 = fn1, fn2 = fn2,
        g1 = v1, g2 = v2, prec.grad = prec.grad,
        prec.norm = prec.norm, max.steps = max.no.steps.ls - i + 1L)
      gradient.list = bisect.opt.result$gradient.list
      offsprings = bisect.opt.result$opt.path

      ## the first two rows of bisect.opt.result contain information on x1 and
      ## x2 --> need to update them here properly
      i = nrow(fn.evals)
      fn.evals[i - 2L,] = fn.evals[i - 2L,] + bisect.opt.result$fn.evals[1L,]
      fn.evals[i - 1L,] = fn.evals[i - 1L,] + bisect.opt.result$fn.evals[2L,]
      fn.evals = rbind(fn.evals, bisect.opt.result$fn.evals[-(1:2),, drop = FALSE])
      opt.path = rbind(opt.path, offsprings[-(1:2),, drop = FALSE])
      i = nrow(opt.path)
      if (bisect.opt.result$found.optimum) {
        ## if the optimizer claimed to have
        ## found an optimum, stop here
        break
      }
      if (i == (max.no.steps.ls + 1L)) {
        ## if we have reached the maximum number of allowed
        ## function evaluations, we can stop the loop here
        break
      }
    }

    if (ls.method %in% c("both", "mo-ls")) {
      ## FIXME: double-check whether max.steps is set correctly
      mo.ls.opt.result = performMultiObjectiveLocalSearch(
        x1 = x1, x2 = x2, x3 = x3, fn1 = fn1, fn2 = fn2,
        prec.grad = prec.grad, prec.norm = prec.norm, prec.angle = prec.angle,
        scale.step = scale.step, lower = lower, upper = upper,
        max.steps = max.no.steps.ls - i + 1L)
      gradient.list = mo.ls.opt.result$gradient.list
      i = nrow(fn.evals)
      ## update last row of fn.evals as this one already contains
      ## x3, but without any function evaluations, whereas
      ## mo.ls.opt.result provides further information on x3
      ## within its first row
      fn.evals[i, sp] = fn.evals[i, sp] + mo.ls.opt.result$fn.evals[1L, sp]
      fn.evals = rbind(fn.evals, mo.ls.opt.result$fn.evals[-1L,, drop = FALSE])
      opt.path = rbind(opt.path, mo.ls.opt.result$opt.path[-1L,, drop = FALSE])
      i = nrow(opt.path)
      if (mo.ls.opt.result$found.optimum) {
        ## if the optimizer claimed to have
        ## found an optimum, stop here
        break
      }
      if (i == (max.no.steps.ls + 1L)) {
        ## if we have reached the maximum number of allowed
        ## function evaluations, we can stop the loop here
        break
      }
    }

    ## having reached this point implies that the local search solvers did not
    ## succeed in finding an optimum but there still is some budget left
    ## --> restart (if allowed)
    if (allow.restarts) {
      offspring = runif(n = d, min = lower, max = upper)
      opt.path = rbind(opt.path, offspring)
      fn.evals = rbind(fn.evals, 0L)
      gradient.list = vector(mode = "list", length = p)
      names(gradient.list) = sprintf("g%i", seq_len(p))
      j = 0L
      next
    } else {
      break
    }
  }

  if (show.info && (i == max.no.steps.ls + 1L)) {
    catf("The maximum number of allowed steps for the local search phase was reached. Enlarge the budget (max.no.steps.ls) to ensure local efficiency of the final observation.")
  }
  if (show.info && (restart.counter > 0L)) {
    BBmisc::catf("Performed %i restart%s while looking for a local efficient point",
      restart.counter, ifelse(restart.counter == 1, "", "s"))
  }
  rownames(opt.path) = NULL
  rownames(fn.evals) = NULL
  return(list(opt.path = opt.path, fn.evals = fn.evals, gradient.list = gradient.list))
}
