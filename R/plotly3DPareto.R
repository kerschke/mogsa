#' @export
plotly3DPareto = function(x, fn, mode = "decision.space", impute.zero = T) {
  # x: columns x1,x2,x3,height
  # if include.objectives, also y1,y2(,y3) are required
  # fn: smoof function, 3 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  if (impute.zero) {
    ## impute heights of zero for log-scale visualizations
    z = x$height
    mz = min(z[z != 0])
    z[z == 0] = mz / 2
    x$height = z
  }
  
  obj.space = as.matrix(x[,grep("y.*", names(x))])
  dec.space = as.matrix(x[,grep("x.*", names(x))])
  
  dims = c()
  for (j in 1:ncol(dec.space)) {
    dims = c(dims, length(unique(dec.space[,j])))
  }
  
  nondomIndices = nondominated(obj.space, dims)
  
  x.nondom = x[nondomIndices,]
  
  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x₁'),
    yaxis = list(range = c(lower[2],upper[2]), title='x₂'),
    zaxis = list(range = c(lower[3],upper[3]), title='x₃')
  )
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x.nondom$y1),max(x.nondom$y1)), title='y₁'),
      yaxis = list(range = c(min(x.nondom$y2),max(x.nondom$y2)), title='y₂'),
      zaxis = list(range = c(min(x.nondom$y3),max(x.nondom$y3)), title='y₃')
    )
  }
  
  marker = list(
    color=~log(height),
    colorscale=plotlyColorscale(),
    cmin=log(min(x$height)),
    cmax=log(max(x$height))
  )
  
  if (mode == "both") {
    x.shared = highlight_key(x.nondom)
    p.decision = plotly3DParetoDecisionSpace(x.shared, fn, marker, scene="scene")
    p.objective = plotly3DParetoObjectiveSpace(x.shared, fn, marker, scene="scene2")
    
    domain.left = list(
      x=c(0,0.5),
      y=c(0,1)
    )
    decision.scene$domain = domain.left
    
    domain.right = list(
      x=c(0.5,1),
      y=c(0,1)
    )
    if (n == 3) {
      objective.scene$domain = domain.right
    } else {
      objective.scene = list(domain=domain.right)
    }
    
    subplot(p.decision, p.objective) %>% layout(
      title = "Decision and Objective Space",
      scene = decision.scene,
      scene2 = objective.scene
    ) %>% highlight(
      on="plotly_click",
      off="plotly_deselect",
      opacityDim = 0.5,
      color = "red"
    ) %>% hide_guides()
  } else if (mode == "decision.space") {
    plotly3DParetoDecisionSpace(x.nondom, fn, marker) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly3DParetoObjectiveSpace(x.nondom, fn, marker) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly3DScanObjectiveSpace(x.nondom, fn, marker)
    }
  }
  
}

#' @export
nondominated = function(fnMat, dims) {
  locallyNondom = locallyNondominatedCPP(fnMat, dims)
  nondom = ecr::nondominated(t(fnMat[locallyNondom,]))
  
  locallyNondom[nondom]
}

plotly3DParetoObjectiveSpace = function(x, fn, marker.style, scene="scene") {
  n = smoof::getNumberOfObjectives(fn)
  
  if (n == 2) {
    plot_ly(data = x,
            type="scatter",
            x=~y1,y=~y2,
            ids=ids,
            mode = "markers",
            marker = marker.style
    )
  } else if (n == 3) {
    plot_ly(data = x,
            type="scatter3d",
            x=~y1,y=~y2,z=~y3,
            scene = scene,
            mode = "markers",
            marker = marker.style
    )
  }
}

plotly3DParetoDecisionSpace = function(x, fn, marker.style, scene="scene") {
  plot_ly(data = x,
          type="scatter3d",
          x=~x1,y=~x2,z=~x3,
          scene = scene,
          mode = "markers",
          marker = marker.style
  )
}
