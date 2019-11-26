#' @export
plotly3DScan = function(x, fn, mode = "decision.space") {
  # x: columns x1,x2,x3,height
  # if include.objectives, also y1,y2(,y3) are required
  # fn: smoof function, 3 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  decision.scene = list(
    aspectmode='cube',
    xaxis = list(range = c(lower[1],upper[1]), title='x₁'),
    yaxis = list(range = c(lower[2],upper[2]), title='x₂'),
    zaxis = list(range = c(lower[3],upper[3]), title='x₃')
  )
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x$y1),max(x$y1)), title='y₁'),
      yaxis = list(range = c(min(x$y2),max(x$y2)), title='y₂'),
      zaxis = list(range = c(min(x$y3),max(x$y3)), title='y₃')
    )
  }
  
  if (mode == "both") {
    x.shared = highlight_key(x)
    p.decision = plotly3DScanDecisionSpace(x.shared, fn, scene="scene")
    p.objective = plotly3DScanObjectiveSpace(x.shared, fn, scene="scene2")
    
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
    plotly3DScanDecisionSpace(x,fn) %>% layout(
      scene = decision.scene
    )
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly3DScanObjectiveSpace(x,fn) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly3DScanObjectiveSpace(x,fn)
    }
  }
  
}

plotly3DScanObjectiveSpace = function(x, fn, frame="x3", scene="scene") {
  n = smoof::getNumberOfObjectives(fn)
  
  if (frame == "x1") {
    frame = ~x1
    ids = ~paste(x2,x3)
  } else if (frame == "x2") {
    frame = ~x2
    ids = ~paste(x1,x3)
  } else if (frame == "x3") {
    frame = ~x3
    ids = ~paste(x1,x2)
  }
  
  if (n == 2) {
    plot_ly(data = x,
            x=~y1,y=~y2,
            ids=ids,
            frame=frame
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000
    )
  } else if (n == 3) {
    plot_ly(data = x,
            x=~y1,y=~y2,z=~y3,
            ids = ids,
            frame = frame,
            scene = scene
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000
    )
  }
}

plotly3DScanDecisionSpace = function(x, fn, frame="x3", scene="scene") {
  if (frame == "x1") {
    frame = ~x1
    ids = ~paste(x2,x3)
  } else if (frame == "x2") {
    frame = ~x2
    ids = ~paste(x1,x3)
  } else if (frame == "x3") {
    frame = ~x3
    ids = ~paste(x1,x2)
  }
  
  # TODO: adapt color scale
  plot_ly(data = x,
          x=~x1,y=~x2,z=~x3,
          frame = frame,
          ids = ids,
          scene = scene
  ) %>% add_markers(
    color=~log(height+1)
  ) %>% animation_opts(
    frame = 1000
  )
}

