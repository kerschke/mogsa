## compute points along the curved connection between two local optima (from different objectives)
computePointsOnCurvedFunction = function(mu1, mu2, sig1, sig2, steps = 100L, sd.val = 3L) {
  n = max(steps, 100L)
  z = seq(1 / (2 * n), 1 - 1 / (2 * n), 1 / n)
  ## help function, which ensures that the points are spread rather evenly across the connection
  ks = 2^qnorm(z, sd = sd.val)
  res = lapply(ks, function(ki) as.numeric(mu1 - ki * solve(sig1 + ki * sig2) %*% sig2 %*% (mu1 - mu2)))
  res = unique(rbind(mu1, do.call(rbind, res), mu2))
  attr(res, "dimnames") = NULL
  return(res)
}

## compute linear function connecting mu1 and mu2
computeLinearFunction = function(mu1, mu2) {
  a = (mu2[2] - mu1[2]) / (mu2[1] - mu1[1])
  b = mu1[2] - a * mu1[1]
  fun = function(x) a * x + b
  return(list(a = a, b = b, fun = fun))
}

## add gradient vector arrows to a base R plot
addArrowsToPlot = function(points, gradients, fac = 0.01 / 2, length = 0.01, ...) {

  for (i in seq_row(points)) {
    st = points[i,]
    ziel = st + fac * gradients[i,]
    arrows(st[1], st[2], ziel[1], ziel[2], col = "grey", length = length, ...)
  }
}

# extractEdgePoints = function(fun1, fun2, fun3 = NULL, s1 = 1L, s2 = 1L, s3 = 1L,
#   shape1 = "sphere", shape2 = "sphere", shape3 = "sphere", eps.x = 1e-4) {
# 
#   elliptical = (shape1 != "sphere") | (shape2 != "sphere") | (shape3 != "sphere")
#   if (elliptical) {
#     eval(rPython::python.load(system.file("mpm2.py", package = "mogsa")), envir = .GlobalEnv)
#   }
# 
#   ## needed for the python function that extracts the covariances
#   n1 = nrow(smoof::getLocalOptimum(fun1)$params)
#   n2 = nrow(smoof::getLocalOptimum(fun2)$params)
#   if (!is.null(fun3)) {
#     n3 = nrow(smoof::getLocalOptimum(fun3)$params)
#   }
#   edgePoints = NULL
#   line.counter = 0L
#   for (i in seq_len(n1)) {
#     mu1 = as.numeric(smoof::getLocalOptimum(fun1)$params[i,])
#     if (elliptical) {
#       D1 = rPython::python.call("getCovarianceMatrices", n1, 2L, "random", s1, TRUE, shape1)
#       D1 = matrix(unlist(D1[[i]]), 2L)
#     }
#     for (j in seq_len(n2)) {
#       line.counter = line.counter + 1L
#       mu2 = as.numeric(smoof::getLocalOptimum(fun2)$params[j,])
#       len = sqrt(sum((mu1 - mu2)^2))
#       steps2 = ceiling(len / eps.x)
#       if (!elliptical) {
#         # if shapes are spherical
#         delta.x = abs(mu2[1L] - mu1[1L]) / steps2
#         ## compute lines between the centers of the spheres
#         z1 = seq(mu1[1L], mu2[1L], length.out = steps2 + 1L)
#         g1 = computeLinearFunction(mu1 = mu1, mu2 = mu2)$fun
#         z2 = g1(z1)
#         XX = cbind(x1 = z1, x2 = z2)
#       } else {
#         # if shapes are elliptical
#         D2 = rPython::python.call("getCovarianceMatrices", n2, 2L, "random", s2, TRUE, shape2)
#         D2 = matrix(unlist(D2[[j]]), 2L)
#         steps2 = as.integer(ceiling(steps2 * 1.2))
#         XX = computePointsOnCurvedFunction(mu1 = mu1, mu2 = mu2, sig1 = D1, sig2 = D2, steps = steps2)
#         colnames(XX) = c("x1", "x2")
#         rownames(XX) = NULL
#       }
#       edgePoints = rbind(edgePoints, cbind(XX, line.counter = line.counter))
#     }
#     if (!is.null(fun3)) {
#       for (j in seq_len(n3)) {
#         line.counter = line.counter + 1L
#         mu2 = as.numeric(smoof::getLocalOptimum(fun3)$params[j,])
#         len = sqrt(sum((mu1 - mu2)^2))
#         steps2 = ceiling(len / eps.x)
#         if (!elliptical) {
#           # if shapes are spherical
#           delta.x = abs(mu2[1L] - mu1[1L]) / steps2
#           ## compute lines between the centers of the spheres
#           z1 = seq(mu1[1], mu2[1], length.out = steps2 + 1L)
#           g1 = computeLinearFunction(mu1 = mu1, mu2 = mu2)$fun
#           z2 = g1(z1)
#           XX = cbind(x1 = z1, x2 = z2)
#         } else {
#           # if shapes are elliptical
#           D2 = rPython::python.call("getCovarianceMatrices", n3, 2L, "random", s3, TRUE, shape3)
#           D2 = matrix(unlist(D2[[j]]), 2L)
#           steps2 = as.integer(ceiling(steps2 * 1.2))
#           XX = computePointsOnCurvedFunction(mu1 = mu1, mu2 = mu2, sig1 = D1, sig2 = D2, steps = steps2)
#           colnames(XX) = c("x1", "x2")
#           rownames(XX) = NULL
#         }
#         edgePoints = rbind(edgePoints, cbind(XX, line.counter = line.counter))
#       }
#     }
#   }
#   if (!is.null(fun3)) {
#     for (i in seq_len(n2)) {
#       mu1 = as.numeric(smoof::getLocalOptimum(fun2)$params[i,])
#       if (elliptical) {
#         D1 = rPython::python.call("getCovarianceMatrices", n2, 2L, "random", s2, TRUE, shape2)
#         D1 = matrix(unlist(D1[[i]]), 2L)
#       }
#       for (j in seq_len(n3)) {
#         line.counter = line.counter + 1L
#         mu2 = as.numeric(smoof::getLocalOptimum(fun3)$params[j,])
#         len = sqrt(sum((mu1 - mu2)^2))
#         steps2 = ceiling(len / eps.x)
#         if (!elliptical) {
#           # if shapes are spherical
#           delta.x = abs(mu2[1L] - mu1[1L]) / steps2
#           ## compute lines between the centers of the spheres
#           z1 = seq(mu1[1L], mu2[1L], length.out = steps2 + 1L)
#           g1 = computeLinearFunction(mu1 = mu1, mu2 = mu2)$fun
#           z2 = g1(z1)
#           XX = cbind(x1 = z1, x2 = z2)
#         } else {
#           # if shapes are elliptical
#           D2 = rPython::python.call("getCovarianceMatrices", n3, 2L, "random", s3, TRUE, shape3)
#           D2 = matrix(unlist(D2[[j]]), 2L)
#           steps2 = as.integer(ceiling(steps2 * 1.2))
#           XX = computePointsOnCurvedFunction(mu1 = mu1, mu2 = mu2, sig1 = D1, sig2 = D2, steps = steps2)
#           colnames(XX) = c("x1", "x2")
#           rownames(XX) = NULL
#         }
#         edgePoints = rbind(edgePoints, cbind(XX, line.counter = line.counter))
#       }
#     }
#   }
#   return(edgePoints)
# }

# ## highlight edges (i.e., potential local efficient sets)
# addGGEdges = function(g, fun1, fun2, fun3 = NULL, s1, s2, s3 = NULL,
#   shape1, shape2, shape3 = NULL, eps.x = 1e-4, ...) {
# 
#   df.edges = extractEdgePoints(fun1 = fun1, fun2 = fun2, fun3 = fun3, s1 = s1, s2 = s2, s3 = s3,
#     shape1 = shape1, shape2 = shape2, shape3 = shape3, eps.x = eps.x)
#   g = g + geom_path(data = as.data.frame(df.edges), mapping = aes(x = x1, y = x2, group = line.counter), ...)
#   return(g)
# }
