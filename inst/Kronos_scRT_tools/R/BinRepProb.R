BinRepProb_ui = function(id,title=NULL) {
  ns = shiny::NS(id)
  shinydashboard::box(
    id = ns('box'),
    status = 'primary',
    title = title,
    width = 12,
    solidHeader = T,
    collapsible = T,
    align = 'center',
    shiny::fluidRow(shiny::plotOutput(ns('plot'), width = '80%',height = 'auto'))
  )
}

BinRepProb_server = function(id, variabilityBR, out, colors, sizes, save,file_name) {
  shiny::moduleServer(id,
                      function(input,
                               output,
                               session,
                               VarBR = variabilityBR,
                               Save = save,
                               Colors = colors,
                               plot_sizes = sizes,
                               Out = out,
                               basename=file_name) {
                        shiny::observe({
                          #load required operators
                          `%>%` = tidyr::`%>%`
                          #plot function
                          q = Kronos.scRT::BinRepProbPlot(VarBR, colors = Colors)
                          output$plot <- shiny::renderPlot({
                            q
                          },
                          height = function() {
                            session$clientData[[paste0('output_', id, '-plot_width')]] * plot_sizes$height / plot_sizes$width
                          })

                          if (Save) {
                            if (!dir.exists(Out)) {
                              dir.create(Out, recursive = T)
                            }
                            ggplot2::ggsave(
                              q,
                              filename = file.path(Out,
                                                   paste0(basename, '_BinRepProb.pdf')),
                              device = grDevices::cairo_pdf,
                              width = plot_sizes$width,
                              height = plot_sizes$height,
                              units = plot_sizes$unit
                            )
                          }
                        })
                      })
}
