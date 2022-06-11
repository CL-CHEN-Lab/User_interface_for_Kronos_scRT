Correlation_ui = function(id) {
  ns <- shiny::NS(id)

  shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      shinyjs::disabled(shiny::div(
        id = ns('dsave'),
        shiny::column(
          width = 3,
          shinyFiles::shinySaveButton(
            id =  ns('Save_Corr'),
            filetype = '.pdf',
            label = 'Save',
            title = 'Save',
            style = 'width:100%;'
          )
        ),
        shiny::column(width = 1, hw_plot_ui(
          id = ns('save_size'),
          height = 5,
          width = 5,
          up = F
        ))
      )),
      shiny::column(
        width = 7,
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
        width = 1,
        shinyWidgets::dropdown(
          inputId = ns('color_dropdown_corr'),
          shiny::tags$h4("Diagonal"),
          colors_ui(
            id = ns('density_fill_corr'),
            "grey",
            label = 'Density'
          ),
          shiny::tags$h4("Upper Triangle"),
          colors_ui(id =ns('min_corr'), "#FDE725FF", label = 'Min value'),
          colors_ui(
            id =ns('midpoint_corr'),
            "#21908CFF",
            label = 'Mid value'
          ),
          colors_ui(id = ns('max_corr'), "#440154FF", label = 'Max Value'),
          colors_ui(id = ns('text_color_corr'), "red", label = 'Text Color'),
          shiny::numericInput(
            inputId = ns('text_size_corr'),
            label = shiny::h4('Text size'),
            value = 15,
            min = 1,
            step = 0.5
          ),
          shiny::tags$h4("Lower Triangle"),
          shiny::selectInput(
            inputId = ns('hex_color_palette_corr'),
            label = NULL,
            choices = c('viridis', 'inferno', 'magma', 'plasma', 'cividis'),
            selected = 'viridis',
            multiple = F,
            width = '100%'
          ),

          shiny::fluidRow(
            shiny::actionButton(
              inputId = ns('apply_color_changes_corr'),
              label = 'Apply',
              width = '100%'
            )
          ),
          circle = TRUE,
          status = 'primary',
          inline = T,
          icon = icon("palette", lib =
                        "font-awesome"),
          width = 300,
          right = T
        )
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
                          'correlation_out'
                        ), width = '80%'), align = 'center')
                      ))
    )
  )

}


Correlation_server = function(id, RT, aes) {
  shiny::moduleServer(id,
                      function(input,
                               output,
                               session,
                               rt = RT,
                               ID =  id) {
                        variables = shiny::reactiveValues()

                        #load required operators
                        `%>%` = tidyr::`%>%`
                        #set save button

                        shinyFiles::shinyFileSave(
                          input = input,
                          id = 'Save_Corr',
                          roots = c(shinyFiles::getVolumes()(),
                                    Home = Sys.getenv("HOME")),
                          defaultRoot = 'Home',
                          allowDirCreate = T
                        )

                        #initialize
                        #set aesthetics correlation
                        Plots = shiny::reactiveValues(
                          density_fill = 'grey',
                          correlation = c("#FDE725FF", "#21908CFF", "#440154FF"),
                          txt_color = "red",
                          font_size = 15,
                          hex_color_palette = 'viridis',
                          p = NULL
                        )

                        shiny::observeEvent(input$apply_color_changes_corr, {
                          if (input$apply_color_changes_corr > 0) {
                            density_fill_corr = colors_server('density_fill_corr')
                            Plots$density_fill = density_fill_corr()

                            text_color_corr = colors_server('text_color_corr')
                            Plots$txt_color = text_color_corr()

                            Plots$font_size = input$text_size_corr

                            min_corr = colors_server(id = 'min_corr')
                            midpoint_corr = colors_server(id = 'midpoint_corr')
                            max_corr = colors_server(id = 'max_corr')
                            Plots$correlation = c(min_corr(),
                                                  midpoint_corr(),
                                                  max_corr())

                            Plots$hex_color_palette = input$hex_color_palette_corr
                          }
                        })

                        shiny::updateCheckboxGroupInput(
                          inputId = 'Sample',
                          choices = unique(RT$basename),
                          selected = unique(RT$basename),
                          inline = T
                        )


                        shiny::observe({
                          if (length(input$Sample) > 1) {
                            shinyjs::disable('dsave')
                            #select samples of interest
                            Plots$p = Kronos.scRT::KCorr_plot(
                              rt %>%
                                dplyr::filter(basename %in% input$Sample),
                              method = stringr::str_to_lower(input$method),
                              density_fill = Plots$density_fill,
                              correlation_gradient = Plots$correlation,
                              correlation_text_color = Plots$txt_color,
                              correlation_text_size = Plots$font_size,
                              hex_color_palette = Plots$hex_color_palette

                            )

                            shinyjs::enable('dsave')
                          } else{
                            shinyjs::disable('dsave')
                            Plots$p = NULL
                          }
                          output$correlation_out = shiny::renderPlot(
                            Plots$p,
                            height = function() {
                              session$clientData[[paste0('output_', ID, '-correlation_out_width')]]
                            }
                          )
                        })

                        shiny::observeEvent(input$Save_Corr, {
                          shinyjs::disable('dsave')
                          if (!is.numeric(input$Save_Corr)) {
                            path = shinyFiles::parseSavePath(
                              roots = c(shinyFiles::getVolumes()(),
                                        Home = Sys.getenv("HOME")),
                              selection = input$Save_Corr
                            )
                            path = path$datapath

                            #recover sizes
                            variables$sizes = hw_plot_server('save_size')
                            print(variables$sizes())
                            variables$sizes = variables$sizes()

                            ggplot2::ggsave(
                              plot = Plots$p,
                              device = grDevices::cairo_pdf,
                              filename = path,
                              units = variables$sizes$unit,
                              height = variables$sizes$height,
                              width = variables$sizes$width
                            )
                          }
                          shinyjs::enable('dsave')
                        })
                      })
}
