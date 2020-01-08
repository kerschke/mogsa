source("functions.R")

ui <- fluidPage(
  h1("Mogsa Dashboard <Working Title>"),
  sidebarLayout(
    sidebarPanel(
      selectInput("fn.id", "Multi-objective function", c("Select a function"="",test.function.ids)),
      selectInput("plot.type", "Type of plot", c("Select a function first" = "")),
      numericInput("grid.size", "Grid size", 100, min=20, max=300, step=1),
      selectInput("space", "Select space to plot", c("Decision Space"="decision.space", "Objective Space"="objective.space", "Decision + Objective Space"="both")),
      actionButton("update.plot", "Update Plot")
    ),
    mainPanel(
      plotly::plotlyOutput(outputId = "plot", height = "100%")
    )
  )
)

server <- function(input, output, session) {
  observe({
    fn = test.functions[[as.numeric(input$fn.id)]]
    if (is.null(fn)) return()
    
    if (smoof::getNumberOfParameters(fn) == 2) {
      updateSliderInput("grid.size", session = session, value = 100, min=20, max=500, step=1)
      updateSelectInput(session = session, inputId = "plot.type", choices = list("Heatmap" = "heatmap"))
    } else {
      updateSliderInput("grid.size", session = session, value = 50, min=20, max=100, step=1)
      updateSelectInput(session = session, inputId = "plot.type", choices = list("Onion Layers" = "layers", "MRI Scan" = "scan", "Pareto Front" = "pareto"))
    }
  })
  
  get.plot = function() {
    fn = test.functions[[as.numeric(input$fn.id)]]
    if (is.null(fn)) {
      return(NULL)
    }
    
    d = smoof::getNumberOfParameters(fn)
    n = smoof::getNumberOfObjectives(fn)
    
    grid = mogsa::generateDesign(fn, input$grid.size)
    grid$obj.space = mogsa::calculateObjectiveValues(grid$dec.space, fn, parallelize = T)
    
    gradients.fn = mogsa::computeGradientFieldGrid(grid, fn)
    grid$height = mogsa::computeCumulatedPathLengths(grid$dec.space, gradients.fn, fix.diagonals = T)
    
    if (d == 2) {
      mogsa::plotly2DHeatmap(grid, fn, mode=input$space)
    } else {
      switch (input$plot.type,
        pareto = mogsa::plotly3DPareto(grid, fn, mode=input$space),
        layers = mogsa::plotly3DLayers(grid, fn, mode=input$space),
        scan = mogsa::plotly3DScan(grid, fn, mode=input$space),
        NULL # if plot.type is invalid
      )
    }
  }
  
  output$plot = plotly::renderPlotly(
    eventReactive(input$update.plot, {
      get.plot()
    })()
  )
}

shinyApp(ui, server)
