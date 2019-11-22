#' @export
plotly3DScan = function(x, fn, mode = "decision.space") {
  # x: columns x1,x2,x3,height
  # if include.objectives, also y1,y2(,y3) are required
  # fn: smoof function, 3 dimensional decision space
  
  if (mode == "both") {
    x.shared = highlight_key(x)
    p.objective = plotly3DScanObjectiveSpace(x.shared,fn)
    p.decision = plotly3DScanDecisionSpace(x.shared,fn)
    
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
      frame = 5000
    )
  } else if (mode == "decision.space") {
    plotly3DScanDecisionSpace(x,fn)
  } else if (mode == "objective.space") {
    plotly3DScanObjectiveSpace(x,fn)
  }
  
}

plotly3DScanObjectiveSpace = function(x, fn, frame="x3") {
  p = smoof::getNumberOfObjectives(fn)
  
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
  
  if (p == 2) {
    plot_ly(data = x,
            x=~y1,y=~y2,
            ids=~paste(x1,x2),
            frame=frame
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000
    )
  } else if (p == 3) {
    # TODO: Add Scene
    plot_ly(data = x,
            x=~y1,y=~y2,z=~y3,
            ids=~paste(x1,x2),
            frame=frame,
            scene="scene2"
    ) %>% add_markers(
      color=~log(height+1)
    ) %>% animation_opts(
      frame = 1000
    )
  }
}

plotly3DScanDecisionSpace = function(x, fn, frame="x3") {
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
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
  
  scene = list(
    xaxis = list(range = c(lower[1],upper[1])),
    yaxis = list(range = c(lower[2],upper[2])),
    zaxis = list(range = c(lower[3],upper[3]))
  )
  
  # TODO: adapt color scale
  plot_ly(data = x,
          x=~x1,y=~x2,z=~x3,
          scene="scene",
          frame = frame,
          ids = ids
  ) %>% add_markers(
    color=~log(height+1)
  ) %>% layout(
    scene = scene
  ) %>% animation_opts(
    frame = 1000
  )
}
