#select color
colors_ui = function(id, color, label = NULL) {
  ns = shiny::NS(id)
  shiny::fluidRow(
    column(width = 10,shiny::h4(label)),
    column(width = 2,
    shinyWidgets::colorPickr(
      inputId = ns('Color'),
      label = NULL,
      selected = color,
      swatches = c(
        'set1' = RColorBrewer::brewer.pal(name = 'Set1', n = 8),
        'set2' = RColorBrewer::brewer.pal(name = 'Set2', n = 8),
        'set3' = RColorBrewer::brewer.pal(name = 'Set3', n = 11),
        'dark2' = RColorBrewer::brewer.pal(name = 'Dark2', n = 8)
      ),
      useAsButton = T,
      width = '100%',
      interaction = list(
        hex = F,
        clear =
          F,
        rgba =
          F
      ),
      theme = 'nano',
      update = 'save',
      hideOnSave = T
    )
  ))
}

colors_server = function(id, input, output, session) {
  shiny::moduleServer(id =  id, function(input, output, session) {
    return(shiny::reactive(input$Color))
  })
}


# height and width of a plot

hw_plot_ui=function(id,height=10,width=15,up=T,right=F){
  ns = shiny::NS(id)

  shinyWidgets::dropdown(up = up,right = right,circle=T,
    shiny::h4('Plot saving options'),
    shiny::radioButtons(inputId = ns('unit'),label = 'Units',choices = c('cm','in'),inline = T,width = '100%',selected = 'cm'),
    shiny::numericInput(inputId = ns('height'),label = 'height',value = height),
    shiny::numericInput(inputId = ns('width'),label = 'width',value = width)
  )
}

hw_plot_server = function(id, input, output, session) {
  shiny::moduleServer(id = id, function(input, output, session) {
    #return parameters if anything change
      return(shiny::reactive({
        list(
          Unit = input$unit,
          height = input$height,
          width = input$width
        )
      }))
  })
}
