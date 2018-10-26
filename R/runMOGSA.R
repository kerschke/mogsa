#' Run MOGSA.
#'
#' @description
#'   Starting from a specific point \code{ind} in the search space, MOGSA will perform
#'   multiple iterations of (a) multi-objective gradient descent towards a local efficient
#'   point, and (b) exploration of the efficient set to which the recently found efficient
#'   point belongs. These two phases are alternatively executed until the maximum number
#'   of visited basins (\code{max.no.basins}) was reached or the algorithm has found
#'   a global efficient set (or at least a multi-objective trap).
#'
#' @note 
#'   ATTENTION: Only turn off the sanity checks (\code{check.data = FALSE}),
#'   if you can ensure that all input parameters are provided in the correct format.
#'
#' @template arg_ind
#' @template arg_biobj_fn
#' @template arg_fni
#' @template arg_max_no_basins
#' @template arg_max_no_steps_ls
#' @template arg_max_no_steps_exploration
#' @template arg_scalestep
#' @template arg_explorationstep
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_ls_method
#' @template arg_lower
#' @template arg_upper
#' @template arg_checkdata
#' @template arg_show_info
#' @template arg_allow_restarts
#' @return [\code{\link{matrix}}]\cr
#'   Returns MOGSA's optimization path, consisting of the visited individuals
#'   (\code{x1} and \code{x2}) and the visited basins (\code{basin}). Furthermore
#'   the sequential \code{iterations} performed per basin and phase of the algorithm
#'   (start, local search, exploration and external points), as well as the actual
#'   phases (\code{type}) are provided. The external points result from evaluations
#'   during the exploration phase, in which the algorithm actually evaluated a point
#'   outside of the efficient set.
#' @examples
#' ## Example 1:
#' # (i) Define two single-objective test problems (fn1, fn2) and a starting point (ind):
#' fn1 = function(x) sum((x - c(2, 0))^2)
#' fn2 = function(x) sum((x - c(0, 1))^2)
#' ind = c(1, 1)
#' 
#' # (ii) Visualize locally efficient set, i.e., the "area" (here: line), of globally
#' # non-dominated points:
#' plot(c(2, 0), c(0, 1), type = "o", pch = 19, lty = 2,
#'   xlim = c(-0.2, 2.2), ylim = c(-0.1, 1.1),
#'   xlab = expression(x[1]), ylab = expression(x[2]), las = 1, asp = 1)
#' text(2, 0, "Optimum of fn1", pos = 2, offset = 1.5)
#' text(0, 1, "Optimum of fn2", pos = 4, offset = 1.5)
#' 
#' # (iii) Run MOGSA:
#' mogsa.result = runMOGSA(ind = ind, fn1 = fn1, fn2 = fn2)
#' 
#' # (iv) Visualize the optimization path of MOGSA:
#' points(mogsa.result[, c("x1", "x2")],
#'   col = as.integer(mogsa.result$type) + 1L,
#'   pch = as.integer(mogsa.result$type) + 14L)
#' legend("topright", pch = 15:18, col = 2:5, legend = levels(mogsa.result$type))
#' @export
runMOGSA = function(ind, fn = NULL, fn1 = NULL, fn2 = NULL,
  max.no.basins = 15L, max.no.steps.ls = 600L, max.no.steps.exploration = 100L,
  scale.step = 0.5, exploration.step = 0.2,
  prec.grad = 1e-6, prec.norm = 1e-6, prec.angle = 1e-4, ls.method = "both",
  lower, upper, check.data = TRUE, show.info = TRUE, allow.restarts = TRUE) {

  # initialize single-objective test functions
  if (is.null(fn) & (is.null(fn1) || is.null(fn2))) {
    stop("Either provide (i) fn or (ii) fn1 and fn2.")
  }
  if (is.null(fn1) || is.null(fn2)) {
    assertFunction(fn)
    if (is.null(fn1)) {
      fn1 = function(...) fn(...)[1L]
    }
    if (is.null(fn2)) {
      fn2 = function(...) fn(...)[2L]
    }
  } else if (!is.null(fn)) {
    warning("As the single-objective functions fn1 and fn2 are provided,
      the bi-objective function fn will be ignored.")
  }

  # perform sanity checks
  d = length(ind)
  assertFunction(fn1)
  assertFunction(fn2)
  assertIntegerish(max.no.basins, lower = 1L, any.missing = FALSE, len = 1L)
  assertIntegerish(max.no.steps.ls, lower = 1L, any.missing = FALSE, len = 1L)
  assertIntegerish(max.no.steps.exploration, lower = 1L, any.missing = FALSE, len = 1L)
  mc.min = .Machine$double.eps
  assertNumber(scale.step, lower = mc.min, finite = TRUE, null.ok = FALSE)
  assertNumber(exploration.step, lower = mc.min, finite = TRUE, null.ok = FALSE)
  assertNumber(prec.grad, lower = mc.min, finite = TRUE, null.ok = FALSE)
  assertNumber(prec.norm, lower = mc.min, finite = TRUE, null.ok = FALSE)
  assertNumber(prec.angle, lower = mc.min, upper = 180, null.ok = FALSE)
  assertLogical(check.data, any.missing = FALSE, len = 1L, null.ok = FALSE)
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
  # each dimension should contain more than one point
  assertTRUE(all(lower < upper))
  assertChoice(ls.method, choices = c("both", "bisection", "mo-ls"), null.ok = FALSE)
  assertLogical(show.info, any.missing = FALSE, len = 1L, null.ok = FALSE)

  ## actual algorithm loop
  ## FIXME: so far, it is hard-coded for p = 2L
  p = 2L
  sp = seq_len(p)
  opt.path = matrix(ind, nrow = 1L)
  fn.evals = matrix(0L, nrow = 1L, ncol = p)
  gradient.list = vector(mode = "list", length = p)
  names(gradient.list) = sprintf("g%i", seq_len(p))
  ls.steps = exploration.steps = external.steps = NULL
  for (ctr in seq_len(max.no.basins)) {
    ind = opt.path[nrow(opt.path),]
    # perform (bi-objective and gradient-driven) local search
    if (show.info) {
      catf("-- Performing Local Search in Basin %i", ctr)
    }
    ls.opt.result = findLocallyEfficientPoint(ind = ind, fn1 = fn1, fn2 = fn2,
      gradient.list = gradient.list, max.no.steps.ls = max.no.steps.ls,
      scale.step = scale.step, prec.grad = prec.grad, prec.norm = prec.norm,
      prec.angle = prec.angle, ls.method = ls.method, lower = lower,
      upper = upper, check.data = check.data, show.info = show.info,
      allow.restarts = allow.restarts)

    ## update function evaluations
    i = nrow(fn.evals)
    fn.evals[i, sp] = fn.evals[i, sp] + ls.opt.result$fn.evals[1L, sp]
    fn.evals = rbind(fn.evals, ls.opt.result$fn.evals[-1L,, drop = FALSE])

    # append local search optimization path to the current one
    opt.path = rbind(opt.path, ls.opt.result$opt.path[-1L,, drop = FALSE])
    ls.steps = c(ls.steps, nrow(ls.opt.result$opt.path) - 1L)

    # explore the previously found local efficient set
    ind = opt.path[nrow(opt.path),]
    if (show.info) {
      catf("-- Exploring a Local Efficient Set in Basin %i", ctr)
    }
    exploration.result = exploreEfficientSet(ind = ind, fn1 = fn1, fn2 = fn2,
      gradient.list = ls.opt.result$gradient.list,
      max.no.steps.exploration = max.no.steps.exploration,
      exploration.step = exploration.step,
      prec.grad = prec.grad, prec.norm = prec.norm, prec.angle = prec.angle,
      lower = lower, upper = upper, check.data = check.data,
      show.info = show.info)

    # append opt.path from the exploration.result to the current opt.path;
    # in cases, in which the local efficiency of the set could be confirmed
    # (i.e., end1$is.local and end1$is.local are both TRUE), the eventually
    # found 1-2 observations (end1$external and end2$external) provide the
    # starting position in an adjacent (and better) basin for the next
    # iteration of the for-loop
    i = nrow(fn.evals)
    fn.evals[i, sp] = fn.evals[i, sp] + exploration.result$fn.evals[1L, sp]
    fn.evals = rbind(fn.evals, exploration.result$fn.evals[-1L,, drop = FALSE])
    opt.path = rbind(opt.path, exploration.result$opt.path[-1L,, drop = FALSE])
    if (!exploration.result$end1$is.local) {
      opt.path = rbind(opt.path, exploration.result$end1$external, exploration.result$end2$external)
      fn.evals = rbind(fn.evals, exploration.result$end1$evals, exploration.result$end2$evals)
      gradient.list = exploration.result$end2$gradient.list
    } else {
      opt.path = rbind(opt.path, exploration.result$end2$external, exploration.result$end1$external)
      fn.evals = rbind(fn.evals, exploration.result$end2$evals, exploration.result$end1$evals)
      gradient.list = exploration.result$end1$gradient.list
    }
    exploration.steps = c(exploration.steps, nrow(exploration.result$opt.path) - 1L)
    external.steps = c(external.steps,
      length(c(exploration.result$end1$external, exploration.result$end2$external)) / d)

    # if the exploration step did not cross a ridge, MOGSA might have found
    # a globally efficient set (or a multi-objective trap) and will terminate
    if (!exploration.result$end1$is.local & !exploration.result$end2$is.local) {
      break
    }
  }

  # store the entire opt.path in a data frame
  opt.path = as.data.frame(opt.path)
  colnames(opt.path) = sprintf("x%i", seq_len(d))
  fn.evals = apply(fn.evals, 2, cumsum)
  if (!is.matrix(fn.evals)) {
    fn.evals = matrix(fn.evals, ncol = p)
  }
  colnames(fn.evals) = sprintf("fn%i.evals", seq_len(p))
  rownames(opt.path) = rownames(fn.evals) = NULL
  opt.path = cbind(opt.path, fn.evals)

  # enhance data frame by information about the MOGSA-loops (closely related
  # to the number of visited basins), performed local search steps within each
  # basin (starting at 1 per basin) and the number of steps used for exploring
  # each basin's efficient set
  opt.path$basin = c(1, unlist(lapply(seq_len(ctr), function(i) {
    c(rep(i, ls.steps[i] + exploration.steps[i] + external.steps[i]))
  })))
  opt.path$iterations = c(1, unlist(lapply(seq_len(ctr), function(i) {
    c(seq_len(ls.steps[i]), seq_len(exploration.steps[i]), seq_len(external.steps[i]))
  })))
  opt.path$type = c("start", unlist(lapply(seq_len(ctr), function(i) {
    c(rep("local search", ls.steps[i]),
      rep("exploration", exploration.steps[i]),
      rep("external", external.steps[i]))
  })))
  opt.path$type = factor(opt.path$type, levels = c("start", "local search", "exploration", "external"))
  return(opt.path)
}
