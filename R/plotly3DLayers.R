#' @export
plotly3DLayers = function(x, fn, mode = "decision.space", max.quantile = 0.05, no.steps = 20) {
  # x: columns x1,x2,x3,height
  # if include.objectives, also y1,y2(,y3) are required
  # fn: smoof function, 3 dimensional decision space
  
  max.height = quantile(x$height, max.quantile)
  min.height = min(x$height)
  x.boundaries = c()
  
  for (height in seq(min.height, max.height, max.height / no.steps)) {
    cat("Height", height, "\n")
    boundary = x[which(x$height <= height),]
    boundary = boundary[!enclosedPoints(boundary),]
    boundary$frame = height
    x.boundaries = rbind(x.boundaries, boundary)
  }
  
  if (mode == "both") {
    x.shared = highlight_key(x.boundaries)
    p.objective = plotly3DLayersObjectiveSpace(x.shared,fn)
    p.decision = plotly3DLayersDecisionSpace(x.shared,fn)
    
    subplot(p.decision, p.objective) %>% layout(
      title = "Decision and Objective Space",
      scene = list(domain=list(x=c(0,0.5),y=c(0,1)),
                   aspectmode='cube'),
      scene2 = list(domain=list(x=c(0.5,1),y=c(0,1)),
                    aspectmode='cube')
    ) %>% highlight(
      on="plotly_click",
      off="plotly_doubleclick",
      opacityDim = 0.5,
      color = "red"
    ) %>% hide_guides() %>% animation_opts(
      # does not work?
      frame = 5000,
      transition = 0
    )
  } else if (mode == "decision.space") {
    plotly3DLayersDecisionSpace(x.boundaries,fn)
  } else if (mode == "objective.space") {
    plotly3DLayersObjectiveSpace(x.boundaries,fn)
  }
  
}

plotly3DLayersObjectiveSpace = function(x, fn) {
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
    # TODO: Add Scene
    plot_ly(data = x,
            x=~y1,y=~y2,z=~y3,
            frame=~frame,
            scene="scene2",
            ids=~paste(x1,x2,x3)
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000,
      transition = 0
    )
  }
}

plotly3DLayersDecisionSpace = function(x, fn) {
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  scene = list(
    xaxis = list(range = c(lower[1],upper[1])),
    yaxis = list(range = c(lower[2],upper[2])),
    zaxis = list(range = c(lower[3],upper[3]))
  )
  
  # TODO: adapt color scale
  plot_ly(data = x,
          x=~x1,y=~x2,z=~x3,
          scene="scene",
          frame = ~frame,
          ids=~paste(x1,x2,x3)
  ) %>% add_markers(
    color=~log(height+1)
  ) %>% layout(
    scene = scene
  ) %>% animation_opts(
    frame = 1000,
    transition = 0
  )
}

enclosedPoints = function(df) {
  # TODO currently rather inefficient
  
  # df: has to include x1,x2,x3
  df.x = signif(df[,c("x1","x2","x3")], 6)
  unique.df = apply(x.fn[,c("x1","x2","x3")], 2, unique)
  
  # we need to assume that the steps in each dimension are alway the same
  # e.g. step size in x1 is 0.01, x2 is 0.02 etc.
  column.deltas = abs(apply(signif(apply(unique.df, 2, diff), 6), 2, min))
  
  column.deltas = column.deltas * 1.5 # 1.5 for stability when comparing later
  
  apply(df.x, 1, function(x) {
    sum((df.x$x1 <= x[1] + column.deltas[1]) &
        (df.x$x1 >= x[1] - column.deltas[1]) &
        (df.x$x2 <= x[2] + column.deltas[2]) &
        (df.x$x2 >= x[2] - column.deltas[2]) &
        (df.x$x3 <= x[3] + column.deltas[3]) &
        (df.x$x3 >= x[3] - column.deltas[3])
    ) == 27
  })
}

