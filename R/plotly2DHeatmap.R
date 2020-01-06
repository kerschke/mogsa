#' @export
plotly2DHeatmap = function(grid, fn, mode = "decision.space", impute.zero = T) {
  # grid: list of obj.space, dims, dec.space, step.sizes
  # fn: smoof function, 2 dimensional decision space
  
  n = smoof::getNumberOfObjectives(fn)
  lower = smoof::getLowerBoxConstraints(fn)
  upper = smoof::getUpperBoxConstraints(fn)
  
  if (impute.zero) {
    grid$height = imputeZero(grid$height)
  }
  
  if (n == 3) {
    objective.scene = list(
      aspectmode='cube',
      xaxis = list(range = c(min(x.nondom[,'y1']),max(x.nondom[,'y1'])), title='y₁'),
      yaxis = list(range = c(min(x.nondom[,'y2']),max(x.nondom[,'y2'])), title='y₂'),
      zaxis = list(range = c(min(x.nondom[,'y3']),max(x.nondom[,'y3'])), title='y₃')
    )
  }
  
  x = cbind.data.frame(grid$dec.space, grid$height, grid$obj.space)
  x = x[order(x$height, decreasing=F),] # relevant for obj.space
  
  marker = plotlyMarker(grid)
  
  if (mode == "both") {
    x.shared = highlight_key(x)
    p.decision = plotly2DHeatmapDecisionSpace(x.shared, fn, marker)
    p.objective = plotly2DHeatmapObjectiveSpace(x.shared, fn, marker)
    
    domain.left = list(
      x=c(0,0.5),
      y=c(0,1)
    )
    decision.scene = list(domain=domain.left)
    
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
      title = paste("Decision and Objective Space of", smoof::getName(fn)),
      scene = decision.scene,
      scene2 = objective.scene
    ) %>% highlight(
      on="plotly_click",
      off="plotly_deselect",
      opacityDim = 0.5,
      color = "red"
    ) %>% hide_guides()
  } else if (mode == "decision.space") {
    plotly2DHeatmapDecisionSpace(x, fn, marker) %>%
      toWebGL() %>%
      hide_guides()
  } else if (mode == "objective.space") {
    if (n == 3) {
      plotly2DHeatmapObjectiveSpace(x, fn, marker) %>% layout(
        scene = objective.scene
      )
    } else {
      plotly2DHeatmapObjectiveSpace(x, fn, marker)
    }
  }
  
}

plotly2DHeatmapObjectiveSpace = function(x, fn, marker.style) {
  n = smoof::getNumberOfObjectives(fn)
  
  if (n == 2) {
    plot_ly(data = x,
            type="scattergl",
            x=~y1,y=~y2,
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

plotly2DHeatmapDecisionSpace = function(x, fn, marker.style) {
  plot_ly(data = x,
          type="heatmap",
          x=~x1,y=~x2,z=~log(height),
          colorscale=plotlyColorscale()
  )
}
