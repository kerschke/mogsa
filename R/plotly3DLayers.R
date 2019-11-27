#' @export
plotly3DLayers = function(x, fn, mode = "decision.space", max.quantile = 0.05, no.steps = 20) {
  # x: columns x1,x2,x3,height
  # if include.objectives, also y1,y2(,y3) are required
  # fn: smoof function, 3 dimensional decision space
  
  max.height = quantile(x$height, max.quantile)
  min.height = min(x$height)
  step.sizes = getStepSizes(x)
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  x.boundaries = c()
  
  df.max = calculateMaxDisplayHeight(x, max.height, include.diagonals = F)
  
  for (height in seq(min.height, max.height, (max.height - min.height) / (no.steps - 1))) {
    boundary = df.max[which(df.max$height <= height & df.max$max.height >= height),]
    
    boundary$frame = height
    x.boundaries = rbind(x.boundaries, boundary)
  }
  
  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x₁'),
    yaxis = list(range = c(lower[2],upper[2]), title='x₂'),
    zaxis = list(range = c(lower[3],upper[3]), title='x₃')
  )
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x.boundaries$y1),max(x.boundaries$y1)), title='y₁'),
      yaxis = list(range = c(min(x.boundaries$y2),max(x.boundaries$y2)), title='y₂'),
      zaxis = list(range = c(min(x.boundaries$y3),max(x.boundaries$y3)), title='y₃')
    )
  }
  
  if (mode == "both") {
    x.shared = highlight_key(x.boundaries)
    p.decision = plotly3DLayersDecisionSpace(x.shared, fn, scene="scene")
    p.objective = plotly3DLayersObjectiveSpace(x.shared, fn, scene="scene2")
    
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
    ) %>% hide_guides() %>% animation_opts(
      # does not work?
      frame = 5000
    )
  } else if (mode == "decision.space") {
    plotly3DLayersDecisionSpace(x.boundaries,fn) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly3DLayersObjectiveSpace(x.boundaries,fn) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly3DLayersObjectiveSpace(x.boundaries,fn)
    }
  }
  
}

plotly3DLayersObjectiveSpace = function(x, fn, scene="scene") {
  p = smoof::getNumberOfObjectives(fn)
  
  if (p == 2) {
    plot_ly(data = x,
            x=~y1,y=~y2,
            frame=~frame,
            ids=~paste(x1,x2,x3)
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
    )
  } else if (p == 3) {
    plot_ly(data = x,
            x=~y1,y=~y2,z=~y3,
            frame=~frame,
            scene=scene,
            ids=~paste(x1,x2,x3)
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
    )
  }
}

plotly3DLayersDecisionSpace = function(x, fn, scene="scene") {
  # TODO: adapt color scale
  plot_ly(data = x,
          x=~x1,y=~x2,z=~x3,
          scene=scene,
          frame = ~frame,
          ids=~paste(x1,x2,x3)
  ) %>% add_markers(
    color=~log(height+1)
  ) %>% animation_opts(
    frame = 1000
  )
}

calculateMaxDisplayHeight = function(df, max.height, include.diagonals = T) {
  # points are "dominated" by their neighbours from some point on
  # calculate this point here
  
  if (include.diagonals) {
    deltas = expand.grid(list(-1:1,-1:1,-1:1))
  } else {
    deltas = as.data.frame(matrix(data = c(0,0,0,1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1), ncol = 3, byrow=T))
  }
  
  sorted.df = df[with(df, order(x3,x2,x1)),]
  dims = apply(sorted.df[,c("x1","x2","x3")], 2, function(x) length(unique(x)))
  
  ids = which(sorted.df$height <= max.height)
  indices = lapply(ids, function(id) convertCellID2IndicesCPP(id, dims))
  
  sorted.df$max.height = Inf
  sorted.df[ids,]$max.height = lapply(indices, function(i) {
    neighbour.heights = apply(deltas, 1, function(d) {
      j = i+d
      if (any(j <= 0) | any(j > dims)) {
        0
      } else {
        id = convertIndices2CellIDCPP(j, dims)
        sorted.df[id,"height"]
      }
    })
    max(neighbour.heights)
  })
  
  sorted.df[with(sorted.df, order(height, decreasing = T)),]
}
