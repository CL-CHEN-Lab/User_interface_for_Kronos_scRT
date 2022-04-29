Twidth_ui = function(id) {
  ns = shiny::NS(paste0('Twidth', id))

  shinydashboard::box(
    height = '100%',
    shinyjs::useShinyjs(),
    id = ns('box'),
    status = 'primary',
    title = id,
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
            selected = 'Extended'
          )
        ),
        shiny::column(
          width = 3,
          shiny::selectInput(
            inputId = ns('n_cat_tw'),
            label = 'Number of RT categories',
            choices = c('1', '2', '3', '5'),
            selected = '2'
          )
        )
      ),
      shiny::fluidRow(shiny::plotOutput(ns('plot'))),
      shiny::fluidRow(shiny::column(
        width = 3,
        shiny::actionButton(
          inputId = ns('save_tw'),
          label = 'Save',
          width = '100%'
        )
      ))

    )
  )
}

Twidth_server = function(id,
                         variability,
                         GenomeAnnotation = dplyr::tibble(),
                         out,
                         cores = 4) {
  shiny::moduleServer(paste0('Twidth', id),
                      function(input,
                               output,
                               session,
                               Var = variability,
                               GA = GenomeAnnotation,
                               Out = out,
                               maxCores = cores) {
                        #load required operators
                        `%>%` = tidyr::`%>%`

                        #initialise
                        data = reactiveValues(
                          twidth = dplyr::tibble(),
                          twidth_fitted_data = dplyr::tibble(),
                          variability = dplyr::tibble(),
                          ncores = 1
                        )
                        shiny::observe({
                          shinyjs::disable('Plot_type')
                          shinyjs::disable('n_cat_tw')
                          if (!is.null(input$n_cat_tw)) {
                            shinyjs::disable('save_tw')
                            if (nrow(GA) > 0) {


                              data$variability = Kronos.scRT::TW_GenomeAnnotation(Variability=Var,GenomeAnnotation=GA)
                            } else{
                              # prepare df to calculate variability over all bins as well
                              data$variability = Kronos.scRT::TW_RTAnnotation(Variability=Var,RT_Groups=as.integer(input$n_cat_tw))
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

                            shinyjs::enable('Plot_type')
                            shinyjs::enable('n_cat_tw')
                          }
                        })

                        shiny::observe({
                          #plot
                          if (!is.null(input$Plot_type) &
                              nrow(data$twidth) > 0) {
                            if (input$Plot_type == 'Extended') {
                              shinyjs::show('n_cat_tw')
                              data$plot = Kronos.scRT::Twidth_extended_plot(
                                Variability = data$variability,
                                Fitted_data = data$twidth_fitted_data,
                                Twidth = data$twidth
                              )
                            } else{
                              shinyjs::hide('n_cat_tw')

                              data$plot = Kronos.scRT::Twidth_barplot(Variability = data$variability,
                                                                 Twidth = data$twidth)
                            }

                            output$plot = shiny::renderPlot({
                              data$plot
                            })
                            shinyjs::enable('save_tw')
                          }
                        })
                        #render plot

                        shiny::observeEvent(input$save_tw, {
                          if (!dir.exists(Out)) {
                            dir.create(Out, showWarnings = FALSE)
                          }
                          if (input$save_tw > 0) {
                            ggplot2::ggsave(
                              plot = data$plot,
                              filename = file.path(
                                Out,
                                paste0(
                                  id,
                                  '_twith_',
                                  data$ncores,
                                  '_',
                                  input$Plot_type,
                                  '_',
                                  '_categories.pdf'
                                )
                              ),
                              device = grDevices::cairo_pdf
                            )
                            readr::write_tsv(x = data$twidth,
                                             file = file.path(
                                               Out,
                                               paste0(id, '_twith_',
                                                      data$ncores,
                                                      '_',
                                                      input$Plot_type,
                                                      '_', '_categories.txt')
                                             ))
                          }
                        })
                      })
}
