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
          offset = 3,
          shiny::selectInput(
            inputId = ns('n_cat_tw'),
            label = 'Number of RT categories',
            choices = c('1', '2', '3', '5'),
            selected = '1'
          )
        ),
        shiny::column(
          width = 3,
          shiny::actionButton(
            inputId = ns('save_tw'),
            label = 'Save',
            width = '100%'
          )
        )
      ),
      shiny::fluidRow(shiny::plotOutput(ns('plot'), width = '100%'))
    )
  )
}

Twidth_server = function(id, variability, out, cores = 4) {
  shiny::moduleServer(paste0('Twidth', id),
                      function(input,
                               output,
                               session,
                               Var = variability,
                               Out = out,
                               maxCores = cores) {
                        #load required operators
                        `%>%` = tidyr::`%>%`

                        #initialize
                        data = reactiveValues(
                          twidth = dplyr::tibble(),
                          twidth_fitted_data = dplyr::tibble(),
                          variability = dplyr::tibble()
                        )

                        shiny::observe({
                          shinyjs::disable('save_tw')
                          shinyjs::disable('n_cat_tw')

                          if (!is.null(input$n_cat_tw)) {
                            # prepare df to calculate variability over all bins as well
                            data$variability = Kronos.scRT::TW_RTAnnotation(Variability=Var,RT_Groups=as.integer(input$n_cat_tw))

                            #fit data
                            data$twidth_fitted_data = Kronos.scRT::Twidth_fit_data(
                              df = data$variability,
                              ncores = ifelse(
                                as.integer(input$n_cat_tw) > maxCores,
                                maxCores,
                                as.integer(input$n_cat_tw)
                              )
                            )

                            #calculate RT
                            data$twidth = Kronos.scRT::Twidth(data$twidth_fitted_data)

                            plot = Kronos.scRT::Twidth_extended_plot(Variability =data$variability,Fitted_data = data$twidth_fitted_data,data$twidth )

                            #render plot

                            output$plot = shiny::renderPlot({
                              plot
                            })
                            shinyjs::enable('save_tw')
                            shinyjs::enable('n_cat_tw')
                            shiny::observeEvent(input$save_tw, {
                              if (!dir.exists(Out)) {
                                dir.create(Out, recursive = T)
                              }

                              if (input$save_tw == 1) {
                                ggplot2::ggsave(
                                  plot = plot,
                                  filename = file.path(
                                    Out,
                                    paste0(id, '_twith_', input$n_cat_tw, '_categories.pdf')
                                  ),
                                  device = grDevices::cairo_pdf
                                )

                                readr::write_tsv(x = data$twidth,
                                                 file = file.path(
                                                   Out,
                                                   paste0(id, '_twith_', input$n_cat_tw, '_categories.txt')
                                                 ))

                              }
                              shinyjs::reset('save_tw')

                            })

                          }
                        })



                      })
}
