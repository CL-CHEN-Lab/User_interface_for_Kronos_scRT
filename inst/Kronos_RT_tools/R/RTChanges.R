RTChanges_ui = function(id) {
  ns <- shiny::NS(id)

  shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shinyjs::disabled(
        shiny::div(
          id = ns('dsave'),
      shiny::column(width = 3,
                        shinyFiles::shinySaveButton(
                          id =  ns('save'),
                          filetype = '.pdf',
                          label = 'Save',
                          title = 'Save',
                          style = 'width:100%;'
                        )
                      ),
      shiny::column(width = 1,
                      hw_plot_ui(ns('RTChanges'),up = F,height = 7,width = 8))
                    )),
      shiny::column(
        width = 3,
        shiny::radioButtons(
          inputId = ns("mode"),
          label = "Mode",
          choices = c("Normal", "Specific ΔRT"),
          width = '100%',
          selected = "Normal",
          inline = T
        )
      ),
      shiny::column(width = 3,
                    shinyjs::hidden(
                      shiny::sliderInput(
                        inputId = ns('dRT'),
                        label = 'ΔRT threshold',
                        min = 0,
                        max = 1,
                        value = 0.2,
                        step = 0.01
                      )
                    )),
      shiny::column(
        offset = 1,
        width = 1,
        shinyWidgets::dropdown(
          inputId = ns('color_dropdown_change'),
          colors_ui(id =ns('LTE'), "darkred", label = 'Late to Early'),
          colors_ui(id =ns('TE'), "red", label = 'to Earlier'),
          colors_ui(id =ns('TL'), "blue", label = 'to Later'),
          colors_ui(id =ns('ETL'), "darkblue", label = 'Early to Late'),
          shiny::fluidRow(
            shiny::actionButton(
              inputId = ns('apply_color_changes_change'),
              label = 'Apply',
              width = '100%'
            )
          ),
          status = 'primary',
          inline = T,
          icon = icon("palette", lib =
                        "font-awesome"),
          width = 300,
          right = T
        )
      ),
      shiny::column(
        width = 12,
        shiny::checkboxGroupInput(
          inputId = ns("Sample"),
          label = "Samples to keep",
          choices = '',
          width = '100%',
          inline = T
        )
      ),
      shiny::fluidRow(align = 'center',
                      shiny::column(
                        width = 12,
                        shiny::div(
                        shiny::plotOutput(ns('RT_changes_out'), width = '80%'),align='center'
                      )))
    )
  )

}


RTChanges_server = function(id, scRT,Colors) {
  shiny::moduleServer(id,
                      function(input,
                               output,
                               session,
                               RT = scRT,
                               ID =  id) {
                        #load required operators
                        `%>%` = tidyr::`%>%`

                        #set save button
                        shinyFiles::shinyFileSave(
                          input = input,
                          id = 'save',
                          roots = c(shinyFiles::getVolumes()(),
                                    Home = Sys.getenv("HOME")),
                          defaultRoot = 'Home',
                          allowDirCreate = T
                        )

                        #initialize
                        Plots = reactiveValues(p = NULL)
                        shiny::updateCheckboxGroupInput(
                          inputId = 'Sample',
                          choices = unique(RT$basename),
                          selected = unique(RT$basename),
                          inline = T
                        )
                        ### RT Changes
                        Colors_changes = shiny::reactiveValues(
                          colors = c(
                            'Late to Early' = 'darkred',
                            'to Earlier' = 'red',
                            'to Later' = 'blue',
                            'Early to Late' = 'darkblue'
                          )
                        )

                        shiny::observeEvent(input$apply_color_changes_change, {
                          if (input$apply_color_changes_change > 0) {
                            LTE = colors_server('LTE')
                            TE = colors_server('TE')
                            TL = colors_server('TL')
                            ETL = colors_server('ETL')

                            Colors_changes$colors = c(
                              LTE(),
                              TE(),
                              TL(),
                              ETL()
                            )
                          }
                        })

                        shiny::observe({
                          print(Colors_changes$colors)
                          print(input$Sample)
                          if (length(input$Sample) > 1) {
                            shinyjs::disable('dsave')
                            #select samples of interest
                            if (input$mode != 'Normal') {
                              shinyjs::show('dRT')
                              Plots$p = Kronos.scRT::RT_changes_plot(
                                RT %>%
                                  dplyr::filter(basename %in% input$Sample),
                                deltaRT = input$dRT,
                                 colors = Colors_changes$colors
                              )
                            } else{
                              shinyjs::hide('dRT')
                              Plots$p = Kronos.scRT::RT_changes_plot(
                                RT %>%
                                  dplyr::filter(basename %in% input$Sample) ,
                                 colors = Colors_changes$colors
                              )
                            }
                            shinyjs::enable('dsave')
                          } else{
                            shinyjs::disable('dsave')
                            Plots$p = NULL
                          }
                          output$RT_changes_out = shiny::renderPlot(
                            Plots$p,
                            height = function() {
                              session$clientData[[paste0('output_', ID, '-RT_changes_out_width')]]
                            }
                          )
                        })

                        shiny::observeEvent(input$save, {
                          shinyjs::disable('dsave')
                          if (!is.numeric(input$save)) {
                            path = shinyFiles::parseSavePath(
                              roots = c(shinyFiles::getVolumes()(),
                                        Home = Sys.getenv("HOME")),
                              selection = input$save
                            )
                            path = path$datapath

                            sizes=hw_plot_server('RTChanges')
                            sizes=sizes()

                            ggplot2::ggsave(
                              plot = Plots$p,
                              device = grDevices::cairo_pdf,
                              filename = path,
                              units = sizes$unit,
                              height =sizes$height,
                              width = sizes$width
                            )
                          }
                          shinyjs::enable('dsave')
                        })
                      })
}
