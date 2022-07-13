#options
options(stringsAsFactors = FALSE,
        dplyr.summarise.inform=FALSE,
        warn = 1,
        scipen = 999)

#load operators
`%>%` = tidyr::`%>%`
`%dopar%` = foreach::`%dopar%`
`%do%` = foreach::`%do%`

#set theme for plots
ggplot2::theme_set(new = ggplot2::theme_bw())

#function to close app
jscode <- "shinyjs.closeWindow = function() { window.close(); }"


# Define UI
ui <- shinydashboard::dashboardPage(
  title = 'Kronos scRT',
  skin = 'blue',
  shinydashboard::dashboardHeader(title = shiny::span(
    shiny::img(src = 'KronosLogo.png', width = '100%')
  )),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      id = 'Sidebar',
      shinydashboard::menuItem(text = "Home",
                               tabName = "Home"),
      shinydashboard::menuItem(text = "Plots",
                               tabName = "Plots"),
      shinydashboard::menuItem(text = "Correlation",
                               tabName = "Correlation"),
      shinydashboard::menuItem(text = "RT Changes",
                               tabName = 'Changes'),
      shinydashboard::menuItem(text = "Hierarchical Clustering",
                               tabName = "HClust"),
      shinydashboard::menuItem(text = "Exit",
                               tabName = "Exit")
    )
  ),
  shinydashboard::dashboardBody(
    #use shinyjs
    shinyjs::useShinyjs(),

    #to close window
    shinyjs::extendShinyjs(text = jscode, functions = c("closeWindow")),

    #use CSS
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet",
                       type = "text/css",
                       href = "custom.css")
    ),

    #setting spinner
    shinybusy::add_busy_spinner(
      spin = "fading-circle",
      position = 'bottom-right',
      color = 'blue',
      height = '200px',
      width = '200px'
    ),

    shinydashboard::tabItems(
      shinydashboard::tabItem(
        tabName = 'Home',
        shiny::fluidRow(
          shiny::column(
            width = 3,
            shinyFiles::shinyFilesButton(
              id = 'LoadRT',
              label = 'RT',
              multiple = F,
              title = 'RT',
              style = 'width:100%;'
            )
          ),
          shiny::column(width = 7,
                        shiny::htmlOutput('RT_loaded')),
          shiny::column(
            width = 1,
            shinyWidgets::dropdown(
              inputId = 'color_dropdown',

              shiny::fluidPage(
                shiny::tags$h4("Color settings"),

                shiny::fluidRow(
                  shiny::uiOutput('Colors'),
                  shiny::actionButton(
                    inputId = 'apply_color_changes',
                    label = 'Apply',
                    width = '100%'
                  )
                )
              ),
              circle = TRUE,
              status = 'primary',
              inline = T,
              icon = icon("palette", lib = "font-awesome"),
              width = 300,
              right = T
            )
          )

        ),
        shiny::fluidRow(column(
          width = 6,
          offset = 3,
          shinyjs::hidden(shiny::plotOutput('distribution',width = 'auto',height = 'auto'))
        )),
        shiny::fluidRow(shinyjs::hidden(div(
          id = 'SdivB',

          shiny::column(
            width = 3,
              shinyFiles::shinySaveButton(
              id = 'Sdist',
              label = 'Save Plot',
              title =  'Save Plot',
              filetype = '.pdf',
              style = 'width:100%;'
            )
          ), shiny::column(width = 1,hw_plot_ui(id = 'Sdist',height = 5,width = 7)))
        ))
      ),
      shinydashboard::tabItem(
        tabName = "Plots",
        shiny::fluidPage(
          shinydashboard::box(
            title = 'Samples Options',
            solidHeader = T,
            status = 'primary',
            collapsible = T,
            collapsed = T,
            width = 12,
            shiny::fluidRow(
              shiny::column(width = 3,
                            shiny::h4('Sample')),
              shiny::column(width = 1,
                            shiny::h4('Color')),
              shiny::column(width = 7,
                            shiny::h4('Plotting Group',
                              bsplus::shiny_iconlink() %>%
                              bsplus::bs_embed_popover(title = 'If a sample has to be plotted in multiple groups, please divide each group with a semicolon.' , placement = 'right')))
            ),
            shiny::uiOutput('Samples_colors_group_plot'),
            shiny::fluidRow(shiny::column(
              width = 3,
              shiny::actionButton(
                inputId = 'Samples_colors_group_apply_change',
                label = 'Apply Changes',
                width = '100%'
              )
            ))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 2,
              shiny::selectInput(
                inputId = 'Chr__RTPlot',
                label = 'Chrmosome',
                choices =  list('Chrom'),
                selected = 'Chrom',
                multiple = F,
                width = '100%'
              )
            ),
            shiny::column(
              width = 10,
              shiny::sliderInput(
                width = '100%',
                inputId = 'range__RTPlot',
                label = 'Coodrinates',
                min = 0,
                max = 0,
                value = c(0, 0),
                dragRange = TRUE,
                post = 'Mb'
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(width = 10,offset = 1,
          shiny::plotOutput(
            'Grouped_RT_plot',
            brush = shiny::brushOpts(
              id = 'Grouped_RT_plot_brush',
              direction = 'x',
              resetOnNew = T
            ),
            dblclick = shiny::dblclickOpts(id = 'Grouped_RT_plot_dbclick'),width = 'auto',height = 'auto'
          ))),
          shiny::fluidRow( shinyjs::hidden(
            shiny::div( id = 'save_grouped_plot_div',
              shiny::column(
            width = 3,
                shinyFiles::shinySaveButton(
                  id = 'save_grouped_plot',
                  label = 'Save',
                  title = 'Save',
                  filetype = '.pdf',
                  style = 'width:100%;'
                )
              ),
            shiny::column(width = 1,
              shiny::uiOutput('save_grouped_plot')
            )
          ))
        )
      )),
      shinydashboard::tabItem(
        tabName = "Correlation",
        shiny::fluidPage(
          Correlation_ui('Correlation')
        )
      ),
      shinydashboard::tabItem(tabName = "Changes",
                              shiny::fluidPage(
                                RTChanges_ui('Changes')
                              )),
      shinydashboard::tabItem(tabName = "HClust",
                              shiny::fluidPage(
                                  RTClustering_ui('HClust')
                                )
                              )
    )
  )
)


server <- function(input, output, session) {
  varialbes = shiny::reactiveValues(
    path = NULL,
    data = dplyr::tibble(),
    dist_plot = NULL,
    basename = NULL,
    colors = NULL
  )

  color_module = shiny::reactiveValues(ui = NULL, server = NULL)
  Samples_colors_group_plot = shiny::reactiveValues(
    ui = NULL,
    server = NULL,
    Sample = NULL,
    Group = NULL,
    Color = NULL,
    plot = NULL,
    plot_backup = NULL,
    Use_color_all_selections =
      F
  )


  #stop app when the session ends
  session$onSessionEnded(function() {
    shiny::stopApp()
  })

  #define operator
  `%>%` = tidyr::`%>%`

  shinyFiles::shinyFileChoose(
    input = input,
    id = 'LoadRT',
    roots = c(shinyFiles::getVolumes()(),
              Home = Sys.getenv("HOME")),
    defaultRoot = 'Home'
  )

  shinyFiles::shinyFileSave(
    input = input,
    id = 'Sdist',
    roots = c(shinyFiles::getVolumes()(),
              Home = Sys.getenv("HOME")),
    defaultRoot = 'Home',
    allowDirCreate = T
  )

  shinyFiles::shinyFileSave(
    input = input,
    id = 'save_grouped_plot',
    roots = c(shinyFiles::getVolumes()(),
              Home = Sys.getenv("HOME")),
    defaultRoot = 'Home',
    allowDirCreate = T
  )


  #load file
  shiny::observeEvent(input$LoadRT, {
    if (!is.numeric(input$LoadRT)) {
      varialbes$path = shinyFiles::parseFilePaths(
        roots = c(shinyFiles::getVolumes()(),
                  Home = Sys.getenv("HOME")),
        selection = input$LoadRT
      )
      varialbes$path = varialbes$path$datapath

      if (Kronos.scRT::right_format(
        file_path = varialbes$path,
        columns_to_check = c('chr', 'start', 'end', 'RT', 'basename', 'group'),
        logical = T
      )) {
        varialbes$data = rbind(
          varialbes$data,
          readr::read_tsv(file = varialbes$path, col_types = readr::cols())
        )

        varialbes$data = varialbes$data %>%
          unique()

        #if basenames are repeted merege basename and group
        if(length(unique(varialbes$data$group))>length(unique(varialbes$data$basename))){
          varialbes$data = varialbes$data %>%
            tidyr::unite(basename,basename,group,sep = ' - ',remove = F)
        }

        output$RT_loaded = shiny::renderText(NULL)
        #set colors
        varialbes$basename = unique(varialbes$data$basename)
        l = length(varialbes$basename)
        if (l > 8) {
          varialbes$colors = colorRampPalette(RColorBrewer::brewer.pal(name = 'Dark2', n = 8))(l)
        } else{
          varialbes$colors = RColorBrewer::brewer.pal(name = 'Dark2', n = 8)[1:l]
        }
        names(varialbes$colors) = varialbes$basename
        #update color options
        color_module$ui = lapply(1:l, function(x)
          colors_ui(
            varialbes$basename[x],
            color = varialbes$colors[x],
            label = varialbes$basename[x]
          ))
        color_module$server = lapply(1:l, function(x) {
          colors_server(
            id = varialbes$basename[x],
            input = input,
            output = output ,
            session = session
          )
        })
        output$Colors = shiny::renderUI({
          color_module$ui
        })
      } else{
        output$RT_loaded = shiny::renderText(
          '<p style="color:#FF0000">Selected file is not correctly formatted. Files have to be tab spaced named tables with the following columns: chr, start, end, RT, basename and group.</p>'
        )
      }

    } else{
      output$RT_loaded = shiny::renderText(NULL)
    }
  })


  # change colors
  shiny::observeEvent(input$apply_color_changes, {
    if (input$apply_color_changes > 0) {
      varialbes$colors = sapply(color_module$server, function(x)
        x())
      names(varialbes$colors) = varialbes$basename
    }

  })


  #plot distributions
  shiny::observeEvent(c(
    varialbes$data,
    varialbes$colors,
    input$w_dist,
    input$h_dist
  ),
  {
    if (nrow(varialbes$data) > 1) {
      #plot densities
      varialbes$dist_plot = varialbes$data %>%
        ggplot2::ggplot(ggplot2::aes(RT, color = basename)) +
        ggplot2::geom_density() +
        ggplot2::scale_color_manual(values = varialbes$colors)

      output$distribution = shiny::renderPlot({
        varialbes$dist_plot
        },height = function() {
          session$clientData[[paste0('output_distribution_width')]]*0.7
      })
      shinyjs::show('distribution')
      shinyjs::show('SdivB')

    } else{
      output$distribution = shiny::renderPlot(NULL)
    }
  })

  #save distribution
  shiny::observeEvent(input$Sdist, {
    if (!is.numeric(input$Sdist)) {
      varialbes$path = shinyFiles::parseSavePath(
        roots = c(shinyFiles::getVolumes()(),
                  Home = Sys.getenv("HOME")),
        selection = input$Sdist
      )
      varialbes$path = varialbes$path$datapath
      varialbes$sizes=hw_plot_server('Sdist')
      varialbes$sizes=varialbes$sizes()
      ggplot2::ggsave(
        filename = varialbes$path,
        plot = varialbes$dist_plot,
        device = grDevices::cairo_pdf,
        units = varialbes$sizes$unit,
        height =varialbes$sizes$height,
        width = varialbes$sizes$width
      )
    }

  })

  ###RT plots
  shiny::observeEvent(input$Sidebar, {
    if (input$Sidebar == 'Plots' & nrow(varialbes$data) != 0) {
      shiny::updateSelectInput(inputId = 'Chr__RTPlot',
                               choices = unique(varialbes$data$chr))

      # change chrom
      shiny::observeEvent(input$Chr__RTPlot, {
        if (input$Chr__RTPlot != 'Chrom') {
          # update min max range
          Max = varialbes$data %>% dplyr::filter(chr == input$Chr__RTPlot) %>% dplyr::pull(end) %>%
            max() / 10 ^ 6
          Min = varialbes$data %>% dplyr::filter(chr == input$Chr__RTPlot) %>% dplyr::pull(end) %>%
            min() / 10 ^ 6
          Step = abs(varialbes$data[1, 'start'] - varialbes$data[1, 'end']) / 10 ^ 6

          shiny::updateSliderInput(
            inputId = 'range__RTPlot',
            min = Min,
            max = Max,
            value = c(Min, Max),
            step = Step
          )

        }
      })

      #set UI
      Samples_colors_group_plot$ui = lapply(1:length(varialbes$colors), function(x)
        grouped_plots_colors_ui(
          id = names(varialbes$colors)[x],
          color = varialbes$colors[x]
        ))
      #render UI
      output$Samples_colors_group_plot = shiny::renderUI({
        Samples_colors_group_plot$ui
      })

      #recover output server
      Samples_colors_group_plot$server = lapply(1:length(varialbes$colors), function(x)
        grouped_plots_colors_server(id = names(varialbes$colors)[x]))

      shiny::observeEvent(input$Samples_colors_group_apply_change, {
        Samples_colors_group_plot$results = lapply(Samples_colors_group_plot$server, function(x)
          x())
        Samples_colors_group_plot$results = Samples_colors_group_plot$results[!sapply(Samples_colors_group_plot$results, is.null)]

        #reformat output

        Samples_colors_group_plot$Sample = sapply(Samples_colors_group_plot$results, function(x)
          x$Sample)
        Samples_colors_group_plot$Color = sapply(Samples_colors_group_plot$results, function(x)
          x$Color)
        Samples_colors_group_plot$Group = lapply(Samples_colors_group_plot$results, function(x)
          x$Group)

        names(Samples_colors_group_plot$Color) = Samples_colors_group_plot$Sample
        names(Samples_colors_group_plot$Group) = Samples_colors_group_plot$Sample

        #reset option highlight
        Samples_colors_group_plot$Use_color_all_selections = F
      })

      #call plotting function
      shiny::observeEvent(c(
        input$range__RTPlot,
        input$Samples_colors_group_apply_change
      ),
      {
        #plot only if coord are available
        if (!is.null(input$range__RTPlot) &
            all(input$range__RTPlot != 0)) {
          if (!is.null(Samples_colors_group_plot$Group) |
              !is.null(Samples_colors_group_plot$Color)) {
            #if grouping and color selections are present
            Samples_colors_group_plot$plot = Kronos.scRT::Plot_bulkRT(
              varialbes$data %>%
                dplyr::filter(basename %in% Samples_colors_group_plot$Sample),
              Coordinates = list(
                chr = input$Chr__RTPlot,
                start =
                  input$range__RTPlot[1] *
                  10 ^ 6,
                end =
                  input$range__RTPlot[2] *
                  10 ^ 6
              ),
              plotting_groups = Samples_colors_group_plot$Group,
              sample_colors = Samples_colors_group_plot$Color
            )
          } else{
            #if no grouping or color selections is present
            Samples_colors_group_plot$plot = Kronos.scRT::Plot_bulkRT(
              varialbes$data,
              Coordinates = list(
                chr = input$Chr__RTPlot,
                start =
                  input$range__RTPlot[1] *
                  10 ^ 6,
                end =
                  input$range__RTPlot[2] *
                  10 ^ 6
              ),
              sample_colors = varialbes$colors
            )
          }
          #save backup to restore it when dbclick
          Samples_colors_group_plot$plot_backup = Samples_colors_group_plot$plot

          #render plot
          output$Grouped_RT_plot = shiny::renderPlot({
            Samples_colors_group_plot$plot
          },height = function() {
            session$clientData[[paste0('output_Grouped_RT_plot_width')]]/7 * ifelse(input$Samples_colors_group_apply_change==0,
                                                                                    length(varialbes$colors),
                                                                                    length(unique(unlist(Samples_colors_group_plot$Group))))
          })

          #show save
          shinyjs::show('save_grouped_plot_div')
          output$save_grouped_plot = renderUI({
            hw_plot_ui(
              'save_grouped_plot',
              width = 14,
              height = 2 * ifelse(
                input$Samples_colors_group_apply_change == 0,
                length(varialbes$colors),
                length(unique(
                  unlist(Samples_colors_group_plot$Group)
                ))
              )
            )
          })
        }
      })

      #select brush
      shiny::observeEvent(input$Grouped_RT_plot_brush, {
        if (!Samples_colors_group_plot$Use_color_all_selections) {
          #first usage and if the user did not decide to use always the same color
          #pop up widow with colors
          shiny::showModal(
            ui = shiny::modalDialog(
              shinyWidgets::colorPickr(
                inputId = 'Color_modal_Diag',
                label = NULL,
                swatches = c(
                  'set1' = RColorBrewer::brewer.pal(name = 'Set1', n = 7),
                  'set2' = RColorBrewer::brewer.pal(name = 'Set2', n = 8),
                  'dark2' = RColorBrewer::brewer.pal(name = 'Dark2', n = 7)
                ),
                width = '100%',
                interaction = list(
                  hex = F,
                  clear =
                    F,
                  rgba =
                    F,
                  save = F
                ),
                update = 'change',
                theme = 'nano',
                inline = T
              ),
              shiny::fluidRow(column(
                width = 12,
                shiny::checkboxInput(
                  inputId = 'Use_color_all_selections',
                  label = 'Use this color for all further selections',
                  value = F,
                  width = '100%'
                )
              )),
              footer = shiny::tagList(
                shiny::actionButton(inputId = 'modalButton_done', label = 'Done')
              )
            )
          )
          #save xmin xmax
          Samples_colors_group_plot$Min_Max = c(input$Grouped_RT_plot_brush$xmin,
                                                input$Grouped_RT_plot_brush$xmax)
        } else{
          #if after the first time the used decided to use always the same color
          Samples_colors_group_plot$plot = Samples_colors_group_plot$plot +
            ggplot2::annotate(
              'rect',
              xmin = input$Grouped_RT_plot_brush$xmin,
              xmax =  input$Grouped_RT_plot_brush$xmax,
              ymin = -Inf,
              ymax = Inf,
              alpha = 0.2,
              fill = input$Color_modal_Diag
            )
        }
      })

      #if Use_color_all_selections is selected modify variable
      shiny::observeEvent(input$Use_color_all_selections, {
        Samples_colors_group_plot$Use_color_all_selections = input$Use_color_all_selections
      })

      #when pop up window is closed (only if the used selects a new color)
      shiny::observeEvent(input$modalButton_done, {
        if (input$modalButton_done > 0) {
          shiny::removeModal()
          #plot
          Samples_colors_group_plot$plot = Samples_colors_group_plot$plot +
            ggplot2::annotate(
              'rect',
              xmin = Samples_colors_group_plot$Min_Max[1],
              xmax =  Samples_colors_group_plot$Min_Max[2],
              ymin = -Inf,
              ymax = Inf,
              alpha = 0.2,
              fill = input$Color_modal_Diag
            )
        }
      })

      # dbclick resets plot to original state
      shiny::observeEvent(input$Grouped_RT_plot_dbclick, {
        Samples_colors_group_plot$plot = Samples_colors_group_plot$plot_backup
        Samples_colors_group_plot$Use_color_all_selections = F
        output$Grouped_RT_plot = shiny::renderPlot({
          Samples_colors_group_plot$plot
        },height = function() {
          session$clientData[[paste0('output_Grouped_RT_plot_width')]]/7 * ifelse(input$Samples_colors_group_apply_change==0,
                                                                                  length(varialbes$colors),
                                                                                  length(unique(unlist(Samples_colors_group_plot$Group))))
          })
      })

      #save
      shiny::observeEvent(input$save_grouped_plot, {
        if (!is.numeric(input$save_grouped_plot)) {

          sizes=hw_plot_server('save_grouped_plot')
          sizes=sizes()

          path = shinyFiles::parseSavePath(
            roots = c(shinyFiles::getVolumes()(),
                      Home = Sys.getenv("HOME")),
            selection = input$save_grouped_plot
          )
          path = path$datapath
          ggplot2::ggsave(
            filename = path,
            plot = Samples_colors_group_plot$plot,
            device = grDevices::cairo_pdf,
            units = sizes$unit,
            width = sizes$width,
            height = sizes$height
          )

        }
      })
    }
  })

  ### correlation
  shiny::observeEvent(input$Sidebar, {
    if (input$Sidebar == 'Correlation' & nrow(varialbes$data) != 0) {
      Corr = Correlation_server(
        id = 'Correlation',
        RT = varialbes$data %>%
          dplyr::mutate(group = basename) %>%
          unique()
      )

    } else{
      Corr = NULL
    }
  })



  shiny::observeEvent(c(input$Sidebar), {
    if (input$Sidebar == 'Changes' & nrow(varialbes$data) != 0) {
      Changes = RTChanges_server(id = 'Changes',
                                 scRT = varialbes$data)

    } else{
      Changes = NULL
    }
  })


  ### HCluster
  shiny::observeEvent(c(input$Sidebar),
                      {
                        if (input$Sidebar == 'HClust' & nrow(varialbes$data) != 0) {
                          HClust = RTClustering_server(id = 'HClust',
                                                       scRT = varialbes$data)

                        } else{
                          HClust = NULL
                        }
                      },
                      priority = 2)
  #exit
  shiny::observeEvent(input$Sidebar, {
    if (input$Sidebar == 'Exit') {
      shinyjs::js$closeWindow()
      shiny::stopApp()
    }
  })

}

# Run the application
shiny::shinyApp(ui = ui,
                server = server)
