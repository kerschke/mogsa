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
computeGradientField = function(points, fn,
  scale.step = 0.5, prec.grad = 1e-6, prec.norm = 1e-6,
  prec.angle = 1e-4, parallelize = FALSE) {
  
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

calcMOGradient = function(ind, fn, prec.grad, prec.norm, prec.angle) {
  
  len = length(ind)
  
  g = -estimateGradientBothDirections(fn = fn, ind = ind, prec.grad = prec.grad, check.data = FALSE)
  g1 = normalizeVectorCPP(vec = g[1,], prec = prec.norm)
  if (all(g1 == 0)) {
    # if the gradient of fn1 is zero, this has to be a local efficient point
    return(rep(0L, len))
  }
  
  g2 = normalizeVectorCPP(vec = g[2,], prec = prec.norm)
  if (all(g2 == 0)) {
    # if the gradient of fn2 is zero, this has to be a local efficient point
    return(rep(0L, len))
  }
  
  angle1 = computeAngleCPP(vec1 = g1, vec2 = g2, prec = prec.norm)
  if (nrow(g) < 3) {
    if (abs(180 - angle1) < prec.angle) {
      # if the angle between both gradients is (approximately) 180 degree,
      # this has to be a local efficient point
      return(rep(0L, len))
    } else {
      return(g1 + g2)
    }
  } else {
    g3 = normalizeVectorCPP(vec = g[3,], prec = prec.norm)
    if (all(g3 == 0)) {
      # if the gradient of fn3 is zero, this has to be a local efficient point
      return(rep(0L, len))
    }
    
    angle2 = computeAngleCPP(vec1 = g1, vec2 = g3, prec = prec.norm)
    if (abs(180 - angle2) < prec.angle) {
      # if the angle between both gradients is (approximately) 180 degree,
      # this has to be a local efficient point
      return(rep(0L, len))
    }
    angle3 = computeAngleCPP(vec1 = g2, vec2 = g3, prec = prec.norm)
    if (abs(180 - angle3) < prec.angle) {
      # if the angle between both gradients is (approximately) 180 degree,
      # this has to be a local efficient point
      return(rep(0L, len))
    }
    if (abs(angle1 + angle2 + angle3 - 360) < prec.angle) {
      # if all gradients show in "opposite" directions, this has to be a local effient point
      return(rep(0L, len))
    } else {
      # otherwise go with the widest angle
      max.angle = max(c(angle1, angle2, angle3))
      if (angle1 == max.angle) {
        return(g1 + g2)
      } else if (angle2 == max.angle) {
        return(g1 + g3)
      } else {
        return(g2 + g3)
      }
    }
  }
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
