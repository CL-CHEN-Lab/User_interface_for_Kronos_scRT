BinRepProb_ui = function(id) {
  ns = shiny::NS(paste0('BinRepProb',id))
  shinydashboard::box(
    id = ns('box'),
    status = 'primary',
    title = id,
    width = 12,
    solidHeader = T,
    collapsible = T,
    align = 'center',
    shiny::fluidRow(shiny::plotOutput(ns('plot'),width = '100%'))
  )
}

BinRepProb_server = function(id, variabilityBR, out, save) {
  shiny::moduleServer(paste0('BinRepProb',id),
                      function(input,
                               output,
                               session,
                               VarBR = variabilityBR,
                               Save = save,
                               Out = out) {
                        shiny::observe({
                          #load required operators
                          `%>%` = tidyr::`%>%`
                          #plot function
                          q=Kronos.scRT::BinRepProbPlot(VarBR)
                          output$plot <- shiny::renderPlot({q})
                            if(Save){
                              if(!dir.exists(Out)){
                                dir.create(Out,recursive = T)
                              }
                              ggplot2::ggsave(q,filename = file.path(Out,
                                                                           paste0(id, '_BinRepProv.pdf')),
                                              device = grDevices::cairo_pdf,width = 14,height = 10)
                            }
                        })
                      })
}


