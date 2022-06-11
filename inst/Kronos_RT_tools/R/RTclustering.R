RTClustering_ui = function(id) {
  ns <- shiny::NS(paste0('RTClustering', id))

  shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shinyjs::disabled(shiny::column(
        width = 3,
        shiny::actionButton(
          inputId = ns('run'),
          label = 'Run',
          width = '100%'
        )
      )),
      shinyjs::disabled(shiny::div(
        id = ns('dsave'),
        shiny::column(
          width = 3,
          shinyFiles::shinySaveButton(
            id =  ns('psave'),
            filetype = '.pdf',
            label = 'Save Plot',
            title = 'Save Plot',
            style = 'width:100%;'
          )
        ),      shiny::column(
          width = 1,
          hw_plot_ui(ns('Hclusts'),up = F)
          ),
        shiny::column(
          width = 3,
          shinyFiles::shinySaveButton(
            id =  ns('tsave'),
            filetype = '.tsv',
            label = 'Save Clusters',
            title = 'Save Clusters',
            style = 'width:100%;'
          )
        )
      )),
      shiny::column(
        width = 1,
        offset = 1,
        shinyWidgets::dropdown(
          shiny::tags$h4("Heatmap Colors"),

          colors_ui(id = ns('min_hclust'), '#005095', label = 'Min value'),

          colors_ui(
            id = ns('midpoint_hclust'),
            '#000000',
            label = 'Mid value'
          ),

          colors_ui(id = ns('max_hclust'), "#a7001b", label = 'Max value'),
          status = 'primary',
          inline = T,
          icon = icon("palette", lib =
                        "font-awesome"),
          width = 300,
          right = T
        )
      )
    )
    ,
    shiny::fluidRow(
      shiny::column(
        width = 3,
        shiny::radioButtons(
          inputId = ns("mode"),
          label = "Crossing RT 0.5?",
          choices = c(TRUE, FALSE),
          width = '100%',
          selected = TRUE,
          inline = T
        )
      ),
      shiny::column(
        width = 3,
        shiny::sliderInput(
          inputId = ns('dRT'),
          label = 'ΔRT threshold',
          min = 0,
          max = 1,
          value = 0.1,
          step = 0.01
        )
      ),
      shiny::column(
        width = 3,
        shiny::radioButtons(
          inputId = ns("Cdef"),
          label = "number of clusters",
          choices = c('Auto', 'Manual'),
          width = '100%',
          selected = 'Auto',
          inline = T
        )
      ),
      shiny::column(width = 3,
                    shinyjs::hidden(
                      shiny::numericInput(
                        inputId = ns('nClust'),
                        label = 'number of clusters',
                        value = 2,
                        step = 1,
                        width = '100%'
                      )
                    ))
    ),
    shiny::fluidRow(shiny::column(
      width = 12,
      shiny::checkboxGroupInput(
        inputId = ns("Sample"),
        label = "Samples to keep",
        choices = '',
        width = '100%',
        inline = T
      )
    )),
    shiny::fluidRow(align = 'center',
                    shiny::column(
                      width = 12,
                      shiny::div(shiny::plotOutput(ns(
                        'RTClustering_out'
                      ), width = '80%'), align = 'center')
                    ))
  )

}


RTClustering_server = function(id, scRT, Colors) {
  shiny::moduleServer(paste0('RTClustering', id),
                      function(input,
                               output,
                               session,
                               RT = scRT,
                               ID = paste0('RTClustering', id)) {
                        #load required operators
                        `%>%` = tidyr::`%>%`

                        #set save button
                        shinyFiles::shinyFileSave(
                          input = input,
                          id = 'psave',
                          roots = c(shinyFiles::getVolumes()(),
                                    Home = Sys.getenv("HOME")),
                          defaultRoot = 'Home',
                          allowDirCreate = T
                        )

                        shinyFiles::shinyFileSave(
                          input = input,
                          id = 'tsave',
                          roots = c(shinyFiles::getVolumes()(),
                                    Home = Sys.getenv("HOME")),
                          defaultRoot = 'Home',
                          allowDirCreate = T
                        )


                        #initialize
                        Results = reactiveValues()
                        shiny::updateCheckboxGroupInput(
                          inputId = 'Sample',
                          choices = unique(RT$basename),
                          selected = unique(RT$basename),
                          inline = T
                        )
                        #set colors plot
                        Colors_HCluster = shiny::reactiveValues()

                        shiny::observeEvent(input$run, {
                          if (input$run > 0) {
                            min_hclust = colors_server(id = 'min_hclust')
                            midpoint_hclust = colors_server(id = 'midpoint_hclust')
                            max_hclust = colors_server(id = 'max_hclust')

                            Colors_HCluster$colors = c(min_hclust(),
                                                       midpoint_hclust(),
                                                       max_hclust())
                          }
                        })

                        #show manual option if needed
                        shiny::observeEvent(input$Cdef, {
                          if (input$Cdef != 'Auto') {
                            shinyjs::show('nClust')
                          } else{
                            shinyjs::hide('nClust')
                          }
                        })

                        shiny::observe({
                          if (length(input$Sample) > 1) {
                            shinyjs::enable('run')
                          } else{
                            shinyjs::disable('run')
                            shinyjs::disable('dsave')
                            Results$p = NULL
                          }
                        })



                        shiny::observeEvent(input$run, {
                          if (input$run > 0) {
                            shinyjs::disable('run')
                            shinyjs::disable('dsave')
                            #select samples of interest
                            if (input$Cdef == 'Auto') {
                              Results$p = Kronos.scRT::RT_clustering(
                                RT %>%
                                  dplyr::filter(basename %in% input$Sample),
                                deltaRT_th = input$dRT ,
                                CrossingRT = input$mode,
                                colors = Colors_HCluster$colors
                              )
                            } else{
                              Results$p = Kronos.scRT::RT_clustering(
                                RT %>%
                                  dplyr::filter(basename %in% input$Sample),
                                deltaRT_th = input$dRT ,
                                CrossingRT = input$mode,
                                n_clusters = input$nClust,
                                colors = Colors_HCluster$colors
                              )
                            }
                            shinyjs::enable('run')
                            shinyjs::enable('dsave')
                          }
                          output$RTClustering_out = shiny::renderPlot(
                            Results$p$plot,
                            height = function() {
                              session$clientData[[paste0('output_', ID, '-RTClustering_out_width')]]
                            }
                          )
                        })

                        shiny::observeEvent(input$psave, {
                          shinyjs::disable('dsave')
                          shinyjs::disable('run')
                          if (!is.numeric(input$psave)) {
                            path = shinyFiles::parseSavePath(
                              roots = c(shinyFiles::getVolumes()(),
                                        Home = Sys.getenv("HOME")),
                              selection = input$psave
                            )
                            path = path$datapath

                            sizes=hw_plot_server('Hclusts')
                            sizes=sizes()

                            ggplot2::ggsave(
                              plot = Results$p$plot,
                              device = grDevices::cairo_pdf,
                              filename = path,
                              units = sizes$unit,
                              height =sizes$height,
                              width = sizes$width
                            )
                          }
                          shinyjs::enable('dsave')
                          shinyjs::enable('run')
                        })

                        shiny::observeEvent(input$tsave, {
                          shinyjs::disable('dsave')
                          shinyjs::disable('run')
                          if (!is.numeric(input$tsave)) {
                            path = shinyFiles::parseSavePath(
                              roots = c(shinyFiles::getVolumes()(),
                                        Home = Sys.getenv("HOME")),
                              selection = input$tsave
                            )
                            path = path$datapath

                            readr::write_tsv(x = Results$p$clusters, file = path)
                          }
                          shinyjs::enable('dsave')
                          shinyjs::enable('run')

                        })
                      })
}
