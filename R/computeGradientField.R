#' Compute the multi-objective gradient vector for a set of points.
#'
#' @description
#'   Computes the multi-objective gradients for a matrix of \code{points}.
#'
#' @param points [\code{\link{matrix}}]\cr
#'   Matrix of points, for which the multi-objective gradient should be computed.
#'   Each row of the matrix will be considered as a separate point, thus the number
#'   of rows corresponds to the number of observations and the number of columns
#'   to the dimensionality of the search space.
#' @template arg_fni
#' @param fn3 [\code{\link{function}}]\cr
#'   The third objective used for computing the multi-objective gradient. If not provided
#'   (\code{fn3 = NULL}), the gradient field will only consider fn1 and fn2.
#' @template arg_scalestep
#' @template arg_precgrad
#' @template arg_precnorm
#' @template arg_precangle
#' @template arg_lower
#' @template arg_upper
#' @param parallelize [\code{\link{logical}(1L)}]\cr
#'   Should the computation of the gradient vectors be parallelized (with \code{parallel::mclapply})?
#'   The default is \code{FALSE}.
#' @return [\code{\link{matrix}}]\cr
#'   Returns \code{matrix} of multi-objective gradients. The i-th row of the matrix
#'   contains the multi-objective gradient vector of the i-th observation (= row)
#'   of the input matrix \code{points}.
#' @examples
#' # Define two single-objective test problems:
#' fn1 = function(x) sum((x - c(0.2, 1))^2)
#' fn2 = function(x) sum(x)
#'
#' # Create a grid of points, for which the gradients should be computed:
#' points = expand.grid(x1 = seq(0, 1, 0.01), x2 = seq(0, 1, 0.05))
#' gradient.field = computeGradientField(points, fn1, fn2)
#' @export
computeGradientField = function(points, fn, prec.grad = 1e-6,
                                prec.norm = 1e-6, prec.angle = 1e-4, parallelize = FALSE) {
  if (parallelize) {
    r = parallel::mclapply(seq_row(points), function(i) {
      ind = as.numeric(points[i,])
      return(calcMOGradient(ind, fn, prec.grad, prec.norm, prec.angle))
    })
  } else {
    r = lapply(seq_row(points), function(i) {
      ind = as.numeric(points[i,])
      return(calcMOGradient(ind, fn, prec.grad, prec.norm, prec.angle))
    })
  }
  
  return(as.matrix(do.call(rbind, r)))
}

#' @export
computeGradientFieldGrid = function(points, fn, obj.values = NULL, prec.norm = 1e-6, prec.angle = 1e-4) {
  obj = smoof::getNumberOfObjectives(fn)
  d = ncol(points) # number of input dimensions
  n = nrow(points) # total number of points
  # calculate dimensions of given field of points
  dims = c()
  for (j in 1:d) {
    dims = c(dims, length(unique(points[,j])))
  }
  
  if (is.null(obj.values)) {
    cat("Evaluating grid of objective values ...\n")
    fn.mat = calculateObjectiveValues(points, fn)
  } else {
    fn.mat = obj.values
  }

  step.sizes = getStepSizes(points)
  
  cat("Estimating single-objective gradients ...\n")
  
  grad.mat.1 = -gridBasedGradientCPP(fn.mat[,1], dims, step.sizes, prec.norm, prec.angle)
  cat("Finished objective 1\n")
  grad.mat.2 = -gridBasedGradientCPP(fn.mat[,2], dims, step.sizes, prec.norm, prec.angle)
  cat("Finished objective 2\n")
  
  if (obj == 2) {
    cat("Estimating multi-objective gradients ...\n")
    
    mo.grad.mat = getBiObjGradientGridCPP(grad.mat.1, grad.mat.2, prec.norm, prec.angle)
  }
  
  if (obj == 3) {
    grad.mat.3 = -gridBasedGradientCPP(fn.mat[,3], dims, step.sizes, prec.norm, prec.angle)
    cat("Finished objective 3\n")
    
    cat("Estimating multi-objective gradients ...\n")

    mo.grad.mat = getTriObjGradientGridCPP(grad.mat.1, grad.mat.2, grad.mat.3, prec.norm, prec.angle)
  }
  
  cat("Finished multiobjective gradient\n")
  return(mo.grad.mat)
}

calcMOGradient = function(ind, fn, prec.grad, prec.norm, prec.angle) {
  len = length(ind)
  
  g = -estimateGradientBothDirections(fn = fn, ind = ind, prec.grad = prec.grad, check.data = FALSE)
  
  if (nrow(g) < 3) {
    getBiObjGradientCPP(g[1,], g[2,], prec.norm, prec.angle)
  } else {
    getTriObjGradientCPP(g[1,], g[2,], g[3,], prec.norm, prec.angle)
  }
}

getStepSizes = function(df) {
  df.x = signif(df, 6)
  unique.x = lapply(df.x, unique)
  
  # we need to assume that the steps in each dimension are alway the same
  # e.g. step size in x1 is 0.01, x2 is 0.02 etc.
  step.sizes = c()
  
  for (i in 1:length(unique.x)) {
    diff.x = diff(sort(unique.x[[i]]))
    diff.x = signif(diff.x, 6)
    step.sizes = c(step.sizes, min(diff.x))
  }
  
  return(step.sizes)
}

# computeHighDimensionalGradientField = function(points, fn1, fn2,
#   scale.step = 0.5, prec.grad = 1e-6, prec.norm = 1e-6,
#   prec.angle = 1e-4, parallelize = FALSE, lower, upper) {
# 
#   len = ncol(points)
#   if (missing(lower)) {
#     lower = apply(points, 2, min) - prec.grad
#   }
#   if (missing(upper)) {
#     upper = apply(points, 2, max) + prec.grad
#   }
# 
#   if (parallelize) {
#     r1 = parallel::mclapply(seq_row(points), function(i) {
#       ind = as.numeric(points[i,])
#       g1 = -estimateGradientBothDirections(fn = fn1, ind = ind, prec.grad = prec.grad, check.data = FALSE)
#       g1 = normalizeVectorCPP(vec = g1, prec = prec.norm)
#       if (all(g1 == 0)) {
#         # if the gradient of fn1 is zero, this has to be a local efficient point
#         return(rep(0L, len))
#       }
# 
#       g2 = -estimateGradientBothDirections(fn = fn2, ind = ind, prec.grad = prec.grad, check.data = FALSE)
#       g2 = normalizeVectorCPP(vec = g2, prec = prec.norm)
#       if (all(g2 == 0)) {
#         # if the gradient of fn2 is zero, this has to be a local efficient point
#         return(rep(0L, len))
#       }
# 
#       angle = computeAngleCPP(vec1 = g1, vec2 = g2, prec = prec.norm)
#       if (abs(180 - angle) < prec.angle) {
#         # if the angle between both gradients is (approximately) 180 degree,
#         # this has to be a local efficient point
#         return(rep(0L, len))
#       } else {
#         return(g1 + g2)
#       }
#     })
#     return(as.matrix(do.call(rbind, r1)))
#   } else {
#     r2 = t(apply(points, 1, function(ind) {
#       g1 = -estimateGradientBothDirections(fn = fn1, ind = ind, prec.grad = prec.grad, check.data = FALSE)
#       g1 = normalizeVectorCPP(vec = g1, prec = prec.norm)
#       if (all(g1 == 0)) {
#         # if the gradient of fn1 is zero, this has to be a local efficient point
#         return(rep(0L, len))
#       }
#       
#       g2 = -estimateGradientBothDirections(fn = fn2, ind = ind, prec.grad = prec.grad, check.data = FALSE)
#       g2 = normalizeVectorCPP(vec = g2, prec = prec.norm)
#       if (all(g2 == 0)) {
#         # if the gradient of fn2 is zero, this has to be a local efficient point
#         return(rep(0L, len))
#       }
# 
#       angle = computeAngleCPP(vec1 = g1, vec2 = g2, prec = prec.norm)
#       if (abs(180 - angle) < prec.angle) {
#         # if the angle between both gradients is (approximately) 180 degree,
#         # this has to be a local efficient point
#         return(rep(0L, len))
#       } else {
#         return(g1 + g2)
#       }
#     }))
#     return(as.matrix(r2))
#   }
# }
