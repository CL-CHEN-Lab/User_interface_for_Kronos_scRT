Correlation_ui = function(id) {
  ns <- shiny::NS(paste0('Correlation', id))

  shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shiny::column(
        width = 3,
        shiny::actionButton(
          inputId = ns('save'),
          label = 'Save',
          width = '100%'
        )
      ),
      shiny::column(
        width = 3,
        shiny::radioButtons(
          inputId = ns("method"),
          label = "Correlation method",
          choices = c("Spearman", "Pearson", "Kendall"),
          width = '100%',
          selected = "Spearman",
          inline = T
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
                        shiny::plotOutput(ns('correlation_out'), width = '100%')
                      ))
    )
  )

}


Correlation_server = function(id, scRT, out) {
  shiny::moduleServer(paste0('Correlation', id),
                      function(input,
                               output,
                               session,
                               RT = scRT,
                               Out = out,
                               ID = paste0('Correlation', id)) {
                        #load required operators
                        `%>%` = tidyr::`%>%`

                        #initialise
                        Plots = reactiveValues(p = NULL)
                        shiny::updateCheckboxGroupInput(
                          inputId = 'Sample',
                          choices = unique(RT$basename),
                          selected = unique(RT$basename),
                          inline = T
                        )

                        shiny::observe({
                          if (length(input$Sample) > 1) {
                            shinyjs::disable('save')
                            #select samples of interest
                            Plots$p = Kronos.scRT::KCorr_plot(
                              RT %>%
                                dplyr::filter(basename %in% input$Sample),
                              method = stringr::str_to_lower(input$method)
                            )

                            shinyjs::enable('save')
                          } else{
                            shinyjs::disable('save')
                            Plots$p = NULL
                          }
                          output$correlation_out = shiny::renderPlot(
                            Plots$p,
                            height = function() {
                              session$clientData[[paste0('output_', ID, '-correlation_out_width')]]
                            }
                          )})

                          shiny::observeEvent(input$save, {
                            shinyjs::disable('save')
                            if (input$save > 0) {
                              if (!dir.exists(file.path(Out))) {
                                dir.create(file.path(Out), recursive = T)
                              }
                              ggplot2::ggsave(
                                plot = Plots$p,
                                device = grDevices::cairo_pdf,
                                filename = file.path(
                                  Out,
                                  paste0(
                                    input$method,
                                    '_correlation',
                                    '_Samples_',
                                    paste(input$Sample, collapse = '_'),
                                    '.pdf'
                                  )
                                )
                              )
                            }
                            shinyjs::enable('save')

                          })
                        })
                      }
