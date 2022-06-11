Twidth_ui = function(id, title = NULL) {
  ns = shiny::NS(id)

  shinydashboard::box(
    height = '100%',
    shinyjs::useShinyjs(),
    id = ns('box'),
    status = 'primary',
    title = title,
    width = 12,
    solidHeader = T,
    collapsible = T,
    align = 'center',
    shiny::fluidPage(
      shiny::fluidRow(
        shiny::column(
          width = 3,
          shiny::selectInput(
            inputId = ns('Plot_type'),
            label = 'Type of plot',
            choices = c('Extended', 'Barplot'),
            selected = 'Extended',
            width = '100%'
          )
        ),
        shiny::column(
          width = 3,
          shiny::selectInput(
            inputId = ns('n_cat_tw'),
            label = 'Number of RT categories',
            choices = c('1', '2', '3', '5'),
            selected = '2',
            width = '100%'
          )
        ),
        shiny::column(
          width = 1,offset = 5,
          shinyWidgets::dropdown(
            inputId = 'color_dropdown_BinRep',
            colors_ui(ns('color'),'red','Plot color'),
            status = 'primary',
            inline = T,
            icon = icon("palette", lib =
                          "font-awesome"),
            width = 300,
            right = T
          )

        )
      ),
      shiny::fluidRow(shiny::plotOutput(ns('plot'), width = '100%')),
      shiny::fluidRow(
        shiny::column(
          width = 3,
          shiny::actionButton(
            inputId = ns('run'),
            label = 'Run',
            width = '100%'
          )
        ),
        shiny::column(
          width = 3,
          shiny::actionButton(
            inputId = ns('save_tw'),
            label = 'Save',
            width = '100%'
          )
        ),
        shiny::column(width = 1,
                      shiny::uiOutput(ns('HW_TW_UI')))
      )
    )
  )
}

Twidth_server = function(id,
                         variability,
                         GenomeAnnotation = dplyr::tibble(),
                         out,
                         file_name,
                         cores = 4) {
  shiny::moduleServer(id,
                      function(input,
                               output,
                               session,
                               Var = variability,
                               GA = GenomeAnnotation,
                               Out = out,
                               maxCores = cores,
                               basename = file_name,
                               ID = id) {
                        shiny::observeEvent(input$Plot_type, {
                          if (input$Plot_type == 'Extended') {
                            output$HW_TW_UI = shiny::renderUI({
                              hw_plot_ui(
                                session$ns('hw_TW'),
                                right = F,
                                up = T,
                                height = 10,
                                width = as.numeric(input$n_cat_tw) * 5
                              )
                            })
                          } else{
                            output$HW_TW_UI = shiny::renderUI({
                              hw_plot_ui(
                                session$ns('hw_TW'),
                                right = F,
                                up = T,
                                height = 5,
                                width = as.numeric(input$n_cat_tw) * 2.5
                              )
                            })
                          }
                        })



                        #load required operators
                        `%>%` = tidyr::`%>%`


                        #initialize
                        data = reactiveValues(
                          twidth = dplyr::tibble(),
                          twidth_fitted_data = dplyr::tibble(),
                          variability = dplyr::tibble(),
                          ncores = 1,
                          hw_TW_size = NULL,
                          rerun = T
                        )


                        shiny::observeEvent(input$run, {
                          shinyjs::disable('Plot_type')
                          shinyjs::disable('n_cat_tw')
                          shinyjs::disable('save_tw')
                          shinyjs::disable('run')

                          if (!is.null(input$n_cat_tw)) {
                            if (data$rerun) {
                              if (nrow(GA) > 0) {
                                data$variability = Kronos.scRT::TW_GenomeAnnotation(Variability = Var,
                                                                                    GenomeAnnotation = GA)
                              } else{
                                # prepare df to calculate variability over all bins as well
                                data$variability = Kronos.scRT::TW_RTAnnotation(Variability =
                                                                                  Var,
                                                                                RT_Groups = as.integer(input$n_cat_tw))
                              }

                              data$ncores = length(unique(data$variability$category))

                              #fit data
                              data$twidth_fitted_data = Kronos.scRT::Twidth_fit_data(
                                df = data$variability,
                                ncores = ifelse(data$ncores > maxCores,
                                                maxCores,
                                                data$ncores)
                              )
                              #calculate TW
                              data$twidth = Kronos.scRT::Twidth(data$twidth_fitted_data)



                              data$rerun = F
                            }
                          }
                        })


                        shiny::observeEvent(input$n_cat_tw, {
                          data$rerun = T
                        })

                        shiny::observeEvent(input$run, {

                          data$color=colors_server('color')
                          data$color=data$color()
                          data$hw_TW_size = hw_plot_server('hw_TW')
                          data$hw_TW_size = data$hw_TW_size()

                          #plot
                          if (!is.null(input$Plot_type) &
                              nrow(data$twidth) > 0) {
                            if (input$Plot_type == 'Extended') {
                              shinyjs::show('n_cat_tw')
                              data$plot = Kronos.scRT::Twidth_extended_plot(
                                Variability = data$variability,
                                Fitted_data = data$twidth_fitted_data,
                                Twidth = data$twidth,
                                Color = data$color
                              )

                              #render plot
                              output$plot = shiny::renderPlot({
                                data$plot + ggplot2::theme(aspect.ratio = data$hw_TW_size$height / data$hw_TW_size$width)
                              })
                            } else{
                              shinyjs::hide('n_cat_tw')

                              data$plot = Kronos.scRT::Twidth_barplot(Variability = data$variability,
                                                                        Twidth = data$twidth,
                                                                      Color = data$color)
                              output$plot = shiny::renderPlot({
                                data$plot + ggplot2::theme(aspect.ratio = data$hw_TW_size$height / data$hw_TW_size$width)
                              })

                            }
                            shinyjs::enable('Plot_type')
                            shinyjs::enable('n_cat_tw')
                            shinyjs::enable('run')
                            shinyjs::enable('save_tw')
                          }
                        })



                        shiny::observeEvent(input$save_tw, {
                          if (!dir.exists(Out)) {
                            dir.create(Out, showWarnings = FALSE)
                          }
                          data$hw_TW_size = hw_plot_server('hw_TW')
                          data$hw_TW_size = data$hw_TW_size()

                          if (input$save_tw > 0) {
                            ggplot2::ggsave(
                              plot = data$plot,
                              filename = file.path(
                                Out,
                                paste0(
                                  basename,
                                  '_twith_',
                                  data$ncores,
                                  '_',
                                  input$Plot_type,
                                  '_',
                                  '_categories.pdf'
                                )
                              ),
                              device = grDevices::cairo_pdf,
                              width = data$hw_TW_size$width,
                              height = data$hw_TW_size$height,
                              units = data$hw_TW_size$Unit
                            )
                            readr::write_tsv(x = data$twidth,
                                             file = file.path(
                                               Out,
                                               paste0(
                                                 basename,
                                                 '_twith_',
                                                 data$ncores,
                                                 '_',
                                                 input$Plot_type,
                                                 '_',
                                                 '_categories.txt'
                                               )
                                             ))
                          }
                        })
                      })
}
