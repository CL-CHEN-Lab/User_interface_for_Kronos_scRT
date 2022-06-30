#options
options(
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE,
  warn = 1,
  scipen = 999
)
#load operator
`%>%` = tidyr::`%>%`
#set ggplot theme
ggplot2::theme_set(new = ggplot2::theme_bw())
#java function to close this app
jscode <- "shinyjs.closeWindow = function() { window.close(); }"
#find max number of cores
maxCores = parallel::detectCores()

# Define UI
ui <- shinydashboard::dashboardPage(
  title = 'Kronos scRT',
  skin = 'green',
  shinydashboard::dashboardHeader(title = shiny::span(
    shiny::img(src = 'KronosLogo.png', width = '100%')
  )),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      id = 'Sidebar',
      shinydashboard::menuItem(text = "Home",
                               tabName = "Home"),
      shinydashboard::menuItem(text = "Dimensionality Reduction",
                               tabName = "DRed"),
      shinydashboard::menuItem(text = "scPlots",
                               tabName = "scPlots"),
      shinydashboard::menuItem(text = "scCN",
                               tabName = "scCN"),
      shinydashboard::menuItem(text = "BinRep",
                               tabName = "BinRep"),
      shinydashboard::menuItem(text = "T-width",
                               tabName = "Twidth"),
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
      color = 'green',
      height = '200px',
      width = '200px'
    ),

    shinydashboard::tabItems(#home
      {
        shinydashboard::tabItem(
          tabName = "Home",
          shinydashboard::box(
            width = 6,
            height = '320',
            title = shiny::div(
              'Input folder(s)',
              bsplus::shiny_iconlink() %>%
                bsplus::bs_embed_popover(title = 'One or more folders created by Kronos scRT processing.' , placement = 'right')
            ),
            solidHeader = T,
            background = 'black',
            align = 'center',
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shinyFiles::shinyDirButton(
                  id = 'Input_dir',
                  label = 'Select',
                  title = 'Input folder(s)',
                  style = 'width:100%;'
                )
              ),
              shiny::column(
                width = 6,
                shiny::selectInput(
                  inputId  = 'RemoveDirectory',
                  label = NULL,
                  choices = 'Remove Folder',
                  selected = 'Remove Folder',
                  width = '100%'
                )
              )
            ),
            shiny::tableOutput('Input_table')
          ),
          shinydashboard::box(
            width = 6,
            height = '320',
            title = shiny::div(
              'Subgroup file(s)',
              bsplus::shiny_iconlink() %>%
                bsplus::bs_embed_popover(title = 'Optional file created by the subgroup option in Dimensionality Reduction. Be aware that this option can be used only if there are no differences in CN between the G1/G2 sub populations. If this is not the case, cells have to be divided accordingly and reprocessed.' , placement = 'right')
            ),
            solidHeader = T,
            background = 'black',
            align = 'center',
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shinyFiles::shinyFilesButton(
                  id = 'SubgroupFile',
                  label = 'Select',
                  title = 'Subgroup file(s)',
                  style = 'width:100%;',
                  multiple = T
                )
              ),
              shiny::column(
                width = 6,
                shiny::selectInput(
                  inputId  = 'RemoveSubgroupFile',
                  label = NULL,
                  choices = 'Remove Subgroup File',
                  selected = 'Remove Subgroup File',
                  width = '100%'
                )
              )
            ),
            shiny::tableOutput('SubgroupFile_out')
          ),
          shinydashboard::box(
            width = 3,
            height = '150',
            title = 'Analysis name',
            solidHeader = T,
            background = 'black',
            shiny::textInput(
              inputId = 'Analysis_Name',
              label = NULL,
              value = 'Analysis',
              width = '100%'
            )
          ),
          shinydashboard::box(
            width = 3,
            height = '150',
            title = 'Output folder',
            solidHeader = T,
            background = 'black',
            align = 'center',
            shinyFiles::shinyDirButton(
              id = 'Output_dir',
              label = 'Select',
              title = 'Output folder',
              style = 'width:100%;'
            ),
            shiny::htmlOutput('Output_dir_out')
          ),
          shinydashboard::box(
            width = 3,
            height = '150',
            title = 'Cores to use for the analysis',
            solidHeader = T,
            background = 'black',
            align = 'center',
            #cores
            shiny::sliderInput(
              width = '100%',
              inputId = 'cores',
              label = NULL,
              value = trunc(maxCores / 2),
              min = 1,
              max = maxCores,
              step = 1,
              ticks = F
            )
          ),
          shinydashboard::box(
            width = 3,
            height = '150',
            title = 'Apply settings',
            solidHeader = T,
            background = 'black',
            align = 'center',
            shiny::actionButton(
              inputId =  'ApplySettings',
              label =
                'Apply settings',
              width = '100%'
            )
          ),
          shinyjs::hidden(shiny::div(id = 'color_option_div',
                                     shiny::fluidRow(
                                       column(
                                         width = 6,
                                         offset = 3,
                                         shinydashboard::box(
                                           width = 12,
                                           title = 'Color Options',
                                           solidHeader = T,
                                           background = 'black',
                                           collapsible = T,
                                           collapsed = F,
                                           align = 'left',
                                           shiny::fluidRow(
                                             shiny::column(width = 10, shiny::h4('Group - Basename')),
                                             shiny::column(width = 2, shiny::h4('Color'))
                                           ),
                                           shiny::uiOutput('General_color'),
                                           shiny::fluidRow(shiny::column(
                                             width = 4,
                                             offset = 4,
                                             shiny::actionButton(
                                               inputId =  'ApplySettingsColors',
                                               label =
                                                 'Apply color settings',
                                               width = '100%'
                                             )
                                           ))
                                         )
                                       )
                                     )))
        )

      },
      #twidth
      {
        shinydashboard::tabItem(
          tabName = "Twidth",
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::radioButtons(
                inputId = 'Regions_tw',
                choices = c('RT categories', 'Customized categories'),
                label = NULL,
                inline = T,
                selected = 'RT categories',
                width = '100%'
              )
            ),
            shiny::column(width = 3,
                          shinyjs::hidden(
                            shiny::div(
                              id = 'div_load_regions_tw',
                              shinyFiles::shinyFilesButton(
                                id = 'load_regions_tw',
                                label = 'Genome Annotation',
                                title = 'Genome Annotation',
                                style = 'width:100%;',
                                multiple = F
                              )
                            )
                          )),
            shiny::column(width = 3,
                          shiny::htmlOutput('GenomeAnnotationFile'))
          ),
          shiny::uiOutput('TW_ui')
        )
      },
      #BinRep
      {
        shinydashboard::tabItem(
          tabName = "BinRep",
          shiny::fluidPage(
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shiny::sliderInput(
                  inputId = 'BinRep_G1_Ploidy',
                  label = 'G1 ploidy quantile',
                  min = 0,
                  max = 1,
                  step = 0.01,
                  dragRange = T,
                  value = c(0.25, 0.75),
                  width = '100%'
                )
              ),
              shiny::column(
                width = 3,
                shiny::sliderInput(
                  inputId = 'BinRep_Early_Cells',
                  label = '% Replication Early cells',
                  min = 0,
                  max = 100,
                  step = 1,
                  dragRange = T,
                  value = c(0, 30),
                  post = '%',
                  width = '100%'
                )
              ),
              shiny::column(
                width = 3,
                shiny::sliderInput(
                  inputId = 'BinRep_Mid_Cells',
                  label = '% Replication Mid cells',
                  min = 0,
                  max = 100,
                  step = 1,
                  dragRange = T,
                  value = c(40, 60),
                  post = '%',
                  width = '100%'
                )
              ),
              shiny::column(
                width = 3,
                shiny::sliderInput(
                  inputId = 'BinRep_Late_Cells',
                  label = '% Replication Late cells',
                  min = 0,
                  max = 100,
                  step = 1,
                  dragRange = T,
                  value = c(70, 100),
                  post = '%',
                  width = '100%'
                )
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shiny::actionButton(
                  inputId = 'Save__BinRep',
                  label = 'Plot',
                  width = '100%'
                )
              ),
              shiny::column(width = 1,
                            shinyjs::hidden(
                              shiny::div(
                                id = 'Save_group__BinRep',
                                hw_plot_ui(
                                  'hw_BinRep',
                                  right = F,
                                  up = F,
                                  height = 7,
                                  width = 14
                                )
                              )
                            )),
              shiny::column(
                offset = 7,
                width = 1,
                shinyWidgets::dropdown(
                  inputId = 'color_dropdown_BinRep',
                  colors_ui(id = 'ES_color_BinRep', "#a7001b", label = 'Early S cells'),
                  colors_ui(id = 'MS_color_BinRep', "#dfbd31", label = 'Mid S cells'),
                  colors_ui(id = 'LS_color_BinRep', "#005095", label = 'Late S cells'),
                  colors_ui(id = 'GG_color_BinRep', "grey", label = 'G1/G2 cells'),
                  shiny::fluidRow(
                    shiny::actionButton(
                      inputId = 'apply_color_changes_BinRep',
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
              )
            ),
            shiny::uiOutput('BinRep_ui')
          )
        )
      },
      #Dred
      {
        shinydashboard::tabItem(tabName = "DRed",
                                shiny::fluidPage(Dim_red_sub_pop_ui('Dred')))
      },
      #scPlots
      {
        shinydashboard::tabItem(
          tabName = 'scPlots',
          shiny::fluidRow(
            shiny::column(
              width = 2,
              shiny::selectInput(
                inputId = 'Chr__scPlot',
                label = 'Chromosome',
                choices =  list('Chrom'),
                selected = 'Chrom',
                multiple = F,
                width = '100%'
              )
            ),
            shiny::column(
              width = 7,
              shiny::sliderInput(
                width = '100%',
                inputId = 'range__scPlot',
                label = 'Coordinates',
                min = 0,
                max = 0,
                value = c(0, 0),
                dragRange = TRUE,
                post = 'Mb'
              )
            ),
            shiny::column(
              width = 3,
              shiny::radioButtons(
                inputId = 'what__scPlot',
                label = 'Filling',
                choices = c('scRT', 'scCN', 'Norm. scCN'),
                selected = 'scRT',
                inline = T,
                width = '100%'
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::actionButton(
                inputId = 'Save__scPlot',
                label = 'Plot',
                width = '100%'
              )
            ),
            shiny::column(width = 1,
                          hw_plot_ui(
                            'HW_scPlot',
                            up = F,
                            height = 15,
                            width = 10
                          )),
            shiny::column(
              offset = 6,
              width = 1,
              shinyWidgets::dropdown(
                inputId = 'color_dropdown_scPlot',
                shiny::uiOutput('Color_UI_scPlot'),
                shiny::fluidRow(
                  shiny::actionButton(
                    inputId = 'apply_color_changes_scPlot',
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
            )

          ),
          shiny::uiOutput('scPlots_UI')
        )
      },
      #scCN
      {
        shinydashboard::tabItem(
          tabName = 'scCN',
          shiny::fluidRow(
            shiny::column(
              width = 2,
              shiny::selectInput(
                inputId = 'Chr__scCN',
                label = 'Chromosome',
                choices =  list('Chrom'),
                selected = 'Chrom',
                multiple = F,
                width = '100%'
              )
            ),
            shiny::column(
              width = 5,
              shiny::sliderInput(
                width = '100%',
                inputId = 'range__scCN',
                label = 'Coordinates',
                min = 0,
                max = 0,
                value = c(0, 0),
                dragRange = TRUE,
                post = 'Mb'
              )
            ),
            shiny::column(
              width = 2,
              shiny::numericInput(
                inputId = 'Levels__scCN',
                step = 1,
                value = 10,
                min = 4,
                label = 'Max CN Levels',
                width = '100%'
              )
            ),
            shiny::column(
              width = 2,
              shiny::actionButton(
                inputId = 'Save__scCN',
                label = 'Save',
                width = '100%'
              )
            )
          ),
          shiny::uiOutput('scCN_UI')
        )
      })
  )
)


server <- function(input, output, session) {
  #variables
  variables = shiny::reactiveValues(
    roots = c(
      shinyFiles::getVolumes()(),
      Home = Sys.getenv("HOME"),
      OutputFolder = file.path(Sys.getenv("HOME"))
    ),
    Save__scPlot = F,
    Save__scCN = F,
    SubgroupFile = dplyr::tibble(),
    Save__BinRep = F,
    Save__BinRep_label = 'Plot',
    colors_BinRep = c(
      "Early S cells" = '#a7001b',
      "Late S cells" = '#005095',
      'Mid S cells' = '#dfbd31',
      'G1/G2 cells' = 'grey'
    )
  )

  #store scCN module info
  scCN_module_ls = shiny::reactiveValues(ui = list(),
                                         server = list())
  #store scPlots module info
  scPlot_module_ls = shiny::reactiveValues(ui = list(),
                                           server = list())

  #store TW module info
  Twidth_module_ls = shiny::reactiveValues(ui = list(),
                                           server = list())
  #store BinRep module info
  BinRep_module_ls = shiny::reactiveValues(ui = list(),
                                           server = list())

  #store data
  data = shiny::reactiveValues(
    S = dplyr::tibble(),
    G = dplyr::tibble(),
    RT = dplyr::tibble(),
    PC = dplyr::tibble(),
    Variability = dplyr::tibble(),
    Reference = dplyr::tibble(),
    input = dplyr::tibble(),
    folder_list = dplyr::tibble(),
    SubgroupFile = dplyr::tibble(),
    variabilityBR = dplyr::tibble(),
    GenomeAnnotation_TW = dplyr::tibble(),
    GenomeAnnotationFile = dplyr::tibble()
  )

  #stop app when the session ends
  session$onSessionEnded(function() {
    shiny::stopApp()
  })

  #home
  {
    shiny::observe({
      shinyFiles::shinyDirChoose(
        input = input,
        id = 'Input_dir',
        session = session,
        roots = variables$roots,
        defaultRoot = 'Home'
      )
      shinyFiles::shinyDirChoose(
        input = input,
        id = 'Output_dir',
        session = session,
        roots = variables$roots,
        defaultRoot = 'Home'
      )
      shiny::observeEvent(input$Analysis_Name, {
        updateTextInput(
          session = session,
          inputId = 'Analysis_Name',
          value = stringr::str_replace_all(
            string = input$Analysis_Name,
            pattern = ' ',
            replacement = '_'
          )
        )
        shinyFiles::shinyFileChoose(
          input = input,
          id = 'SubgroupFile',
          session = session,
          roots = variables$roots,
          defaultRoot = 'Home'
        )
      })

      output$Output_dir_out <-
        renderText(paste(
          '<H4><b>',
          file.path(variables$roots['OutputFolder'], input$Analysis_Name),
          '</H4></b>'
        ))

      if (nrow(data$folder_list) > 0) {
        if (any(data$folder_list$Resolution == 'WrongFolder!')) {
          shinyjs::disable('ApplySettings')
        } else{
          shinyjs::enable('ApplySettings')
        }
      } else{
        shinyjs::disable('ApplySettings')

      }

    })

    #load files
    shiny::observeEvent(input$Input_dir, {
      if (!is.numeric(input$Input_dir)) {
        # id input folder
        variables$Input_folder =  shinyFiles::parseDirPath(roots = variables$roots,
                                                           selection = input$Input_dir)


        #create table with data
        files = list.files(variables$Input_folder, full.names = T)
        data$input =
          rbind(data$input,
                tryCatch(
                  dplyr::tibble(
                    scRT = files[stringr::str_detect(files, pattern = '_calculated_replication_timing')],
                    scCN_S = files[stringr::str_detect(files, pattern = 'G1_G2_single_cells_CNV', negate = T) &
                                     stringr::str_detect(files, pattern = 'single_cells_CNV')],
                    scCN_G = files[stringr::str_detect(files, pattern = 'G1_G2_single_cells_CNV')],
                    Variability = files[stringr::str_detect(files, pattern = '_scRT_variability')],
                    Reference = ifelse(any(
                      stringr::str_detect(files, pattern = '_reference_replication_timing_')
                    ),
                    files[stringr::str_detect(files, pattern = '_reference_replication_timing_')], NA)
                  ),
                  error = function(x)
                    dplyr::tibble()
                )) %>%
          unique()

        data$folder_list = rbind(data$folder_list,
                                 tryCatch(
                                   dplyr::tibble(
                                     Folder = variables$Input_folder,
                                     Resolution = paste0(
                                       readr::read_tsv(file = files[stringr::str_detect(files, pattern = '_calculated_replication_timing')], n_max = 1) %>%
                                         dplyr::mutate(r = end - start) %>%
                                         dplyr::pull(r) / 10 ^ 6,
                                       'Mb'
                                     )

                                   ),
                                   error = function(x)
                                     dplyr::tibble(
                                       Folder = variables$Input_folder ,
                                       Resolution = 'WrongFolder!'
                                     )
                                 )) %>%
          unique()

        shiny::updateSelectInput(
          inputId = 'RemoveDirectory',
          choices = c('Remove Folder', data$folder_list$Folder)
        )

        #update remove and output

        output$Input_table = shiny::renderTable({
          if (nrow(data$folder_list) == 0) {
            NULL
          } else{
            data$folder_list
          }
        })

      }
    })

    shiny::observeEvent(input$RemoveDirectory, {
      # if any folder has been provided it is possible to remove it with RemoveDirectory
      if (input$RemoveDirectory != 'Remove Folder') {
        data$folder_list = data$folder_list %>%
          dplyr::filter(Folder != input$RemoveDirectory)

        data$input = data$input %>% dplyr::filter(
          stringr::str_detect(
            string = scRT,
            pattern = input$RemoveDirectory,
            negate = T
          )
        )

        shiny::updateSelectInput(
          inputId = 'RemoveDirectory',
          choices = ifelse(
            nrow(data$folder_list) != 0,
            c('Remove Folder', data$folder_list$Folder),
            'Remove Folder'
          ),
          selected = 'Remove Folder'
        )
      }
    })

    #load subgroups
    shiny::observeEvent(input$SubgroupFile, {
      if (!is.numeric(input$SubgroupFile)) {
        New_files = shinyFiles::parseFilePaths(roots = variables$roots,
                                               selection = input$SubgroupFile) %>%
          dplyr::mutate(
            State = Kronos.scRT::right_format(
              file_path = datapath,
              columns_to_check = c('Cell', 'basename', 'group', 'subpopulation'),
              delim = '\t',
              wrong_message = 'WrongFormat! Not uploaded',
              rigth_message = 'Uploaded'
            )
          ) %>%
          dplyr::select('Name' = name, State, datapath)

        # id input folder
        variables$SubgroupFile =  rbind(variables$SubgroupFile,
                                        New_files)

        shiny::updateSelectInput(
          inputId = 'RemoveSubgroupFile',
          choices = c('Remove Subgroup File', variables$SubgroupFile$Name)
        )

        #if data have already being loaded
        if (input$ApplySettings > 0) {
          #upload file and add it to the already loaded ones
          temp_subpop = Kronos.scRT::load_multiple_df(New_files$datapath)

          data$SubgroupFile = rbind(data$SubgroupFile, temp_subpop)

          #apply changes to RT and scCN
          temp_variable = Kronos.scRT::extractSubpop(
            scCN = data$S,
            scRT = data$RT,
            scVariability = data$Variability ,
            subpopulation =  temp_subpop,
            RefRT = data$Reference
          )

          data$S = temp_variable$scCN

          data$RT = temp_variable$scRT

          data$Reference = temp_variable$RefRT

          data$Variability = temp_variable$scVariability


          rm('temp_variable')
          rm('temp_subpop')

        }



        #data output
        output$SubgroupFile_out = shiny::renderTable({
          if (nrow(variables$SubgroupFile) == 0) {
            NULL
          } else{
            variables$SubgroupFile[c('Name', 'State')]
          }
        })

      }
    })

    shiny::observeEvent(input$RemoveSubgroupFile, {
      # if any folder has been provided it is possible to remove it with RemoveDirectory
      if (input$RemoveSubgroupFile != 'Remove Subgroup File') {
        temptoRemove = Kronos.scRT::load_multiple_df(
          variables$SubgroupFile %>%
            dplyr::filter(Name == input$RemoveSubgroupFile) %>%
            dplyr::pull(datapath)
        )

        variables$SubgroupFile = variables$SubgroupFile %>%
          dplyr::filter(Name != input$RemoveSubgroupFile)

        shiny::updateSelectInput(
          inputId = 'RemoveSubgroupFile',
          choices = ifelse(
            nrow(variables$SubgroupFile) != 0,
            c('Remove Subgroup File', variables$SubgroupFile$Name),
            'Remove Subgroup File'
          ),
          selected = 'Remove Subgroup File'
        )

        #if data have been applied, restore old setting
        if (input$ApplySettings > 0) {
          #apply changes to RT and scCN
          temp_variable = Kronos.scRT::rejoinSubpop(
            scCN = data$S,
            scRT = data$RT,
            scVariability = data$Variability ,
            subpopulation =  temptoRemove,
            RefRT = data$Reference
          )

          data$S = temp_variable$scCN

          data$RT = temp_variable$scRT

          data$Reference = temp_variable$RefRT

          data$Variability = temp_variable$scVariability

          #upload file and add it to the already loaded ones
          data$SubgroupFile = data$SubgroupFile %>%
            dplyr::left_join(temptoRemove %>% dplyr::mutate(keep = F)) %>%
            dplyr::filter(keep) %>%
            dplyr::select(-keep)

          rm('temptoRemove')
          rm('temp_variable')

        }
      }
    })

    #upload data
    shiny::observeEvent(input$ApplySettings, {
      if (input$ApplySettings > 0) {
        shinyjs::disable('Output_dir')
        shinyjs::disable('Analysis_Name')
        shinyjs::disable('ApplySettings')
        shinyjs::disable('Upload_subgroups')

        data$S = Kronos.scRT::load_multiple_df(data$input$scCN_S)
        data$G = Kronos.scRT::load_multiple_df(data$input$scCN_G)
        data$Variability = Kronos.scRT::load_multiple_df(data$input$Variability)
        data$Reference = Kronos.scRT::load_multiple_df(data$input$Reference)
        data$RT = Kronos.scRT::load_multiple_df(data$input$scRT)

        if (ncol(variables$SubgroupFile) != 0) {
          data$SubgroupFile = Kronos.scRT::load_multiple_df(variables$SubgroupFile$datapath)

          temp_variable = Kronos.scRT::extractSubpop(
            scCN = data$S,
            scRT = data$RT,
            scVariability = data$Variability ,
            subpopulation =  data$SubgroupFile,
            RefRT = data$Reference
          )

          data$S = temp_variable$scCN

          data$RT = temp_variable$scRT

          data$Reference = temp_variable$RefRT

          data$Variability = temp_variable$scVariability

          rm('temp_variable')
        }

        data$samples_names_and_colors =rbind(data$RT %>%
            dplyr::mutate(basename = group) %>%
            dplyr::select(group, basename) %>%
            unique() %>%
            dplyr::mutate(
              type = 'Sample',
              id = paste(group, basename, sep = ' - ')
            ),data$Reference %>%
            dplyr::select(group, basename) %>%
            unique() %>%
            dplyr::mutate(
              type = 'Reference',
              id = paste(group, basename, sep = ' - ')
            ))


        if (nrow(data$samples_names_and_colors%>%dplyr::filter(type == 'Sample')) <= 6) {
          Paired_colors = RColorBrewer::brewer.pal(name = 'Paired', n = 12)
        } else{
          Paired_colors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(name = 'Set1', n =
                                                                                 8))(nrow(data$samples_names_and_colors) * 2)
        }

        Sequence = seq(2, nrow(data$samples_names_and_colors) * 2, 2)
        Paired_colors = lapply(1:nrow(data$samples_names_and_colors%>%dplyr::filter(type == 'Sample')), function(x)
          dplyr::tibble(
            group = data$samples_names_and_colors$group[x],
            Sample = Paired_colors[Sequence[x]],
            Reference = Paired_colors[Sequence[x] - 1]
          ))
        Paired_colors = do.call('rbind', Paired_colors)
        Paired_colors=Paired_colors%>%tidyr::gather(type,color,-group)


        data$samples_names_and_colors=data$samples_names_and_colors%>%dplyr::inner_join(Paired_colors)

        #color options
        output$General_color = shiny::renderUI({
          lapply(1:nrow(data$samples_names_and_colors), function(x)
            colors_ui(
              id = data$samples_names_and_colors$id[x],
              data$samples_names_and_colors$color[x],
              data$samples_names_and_colors$id[x]
            ))
        })


        shinyjs::show('color_option_div')

      }
    })

    #apply color options
    shiny::observeEvent(input$ApplySettingsColors,{

      General_colrs=lapply(1:nrow(data$samples_names_and_colors), function(x)
        colors_server(
          id = data$samples_names_and_colors$id[x]
        ))
      General_colrs=sapply(General_colrs, function(x) x())
      data$samples_names_and_colors$color=General_colrs

    })

    shiny::observeEvent(input$Output_dir, {
      if (!is.numeric(input$Output_dir)) {
        variables$roots = c(
          shinyFiles::getVolumes()(),
          Home = Sys.getenv("HOME"),
          OutputFolder = shinyFiles::parseDirPath(
            roots = variables$roots,
            selection = input$Output_dir
          )
        )
      }
    })
  }
  #scPlots
  {
    shiny::observeEvent(input$Sidebar, {
      if (input$Sidebar == 'scPlots' & ncol(data$RT) != 0) {
        shiny::updateSelectInput(inputId = 'Chr__scPlot',
                                 choices = unique(data$RT$chr))

        data$summary = data$S %>%
          dplyr::ungroup() %>%
          dplyr::summarise(CN_bg = round(stats::quantile(CN_bg, c(0.01, 0.99)), 1),
                           CN = round(stats::quantile(CN, c(0.01, 0.99)), 1))


      }

    })
    # change chrom
    shiny::observeEvent(input$Chr__scPlot, {
      if (input$Chr__scPlot != 'Chrom') {
        # update min max range
        Max = data$RT %>% dplyr::filter(chr == input$Chr__scPlot) %>% dplyr::pull(end) %>%
          max() / 10 ^ 6
        Step = abs(data$RT[1, 'start'] - data$RT[1, 'end']) / 10 ^ 6

        shiny::updateSliderInput(
          inputId = 'range__scPlot',
          min = 0,
          max = Max,
          value = c(0, Max),
          step = Step
        )

      }

    })
    ####create folder plots and set to save
    shiny::observeEvent(input$Save__scPlot, {
      if (!dir.exists(file.path(variables$roots['OutputFolder'],
                                input$Analysis_Name, 'scPlots'))) {
        dir.create(file.path(variables$roots['OutputFolder'],
                             input$Analysis_Name, 'scPlots'),
                   recursive = T)
      }

      variables$Save__scPlot = T
    })


    #reset all if something changes
    shiny::observeEvent(c(
      input$range__scPlot,
      input$what__scPlot,
      input$Save__scPlot
    ),
    {
      if (nrow(data$RT) > 0) {
        G = unique(data$RT$group)

        scPlot_module_ls$ui = lapply(G, function(g)
          scPlots_ui(paste0('scPlots', g), title = g))
        scPlot_module_ls$server = lapply(G, function(g)
          scPlots_server(
            paste0('scPlots', g),
            RTs = rbind(
              data$RT %>%
                dplyr::mutate(basename = group) %>%
                dplyr::select(chr, start, end, group, basename, RT),
              data$Reference
            ) %>%
              dplyr::filter(
                group == g,
                chr == input$Chr__scPlot,
                start > input$range__scPlot[1] *
                  10 ^ 6,
                end < input$range__scPlot[2] *
                  10 ^ 6
              ),
            scCN = data$S %>%
              dplyr::filter(
                group == g,
                chr == input$Chr__scPlot,
                start >= input$range__scPlot[1] *
                  10 ^ 6,
                end <= input$range__scPlot[2] *
                  10 ^ 6
              ),
            filling = input$what__scPlot,
            Extreme_values = data$summary,
            out = file.path(variables$roots['OutputFolder'],
                            input$Analysis_Name, 'scPlots'),
            colors=data$samples_names_and_colors%>%
              dplyr::filter(
                group == g),
            save = variables$Save__scPlot

          ))

        if (variables$Save__scPlot) {
          variables$Save__scPlot = F
        }

        output$scPlots_UI = shiny::renderUI({
          scPlot_module_ls$ui
        })
      }
    })
  }


  # #scCN
  {
    shiny::observeEvent(input$Sidebar, {
      if (input$Sidebar == 'scCN' & nrow(data$G) != 0) {
        shiny::updateSelectInput(inputId = 'Chr__scCN',
                                 choices = unique(data$RT$chr))
      }
    })
    # change chrom
    shiny::observeEvent(input$Chr__scCN, {
      if (input$Chr__scCN != 'Chrom') {
        # update min max range
        Max = data$RT %>% dplyr::filter(chr == input$Chr__scCN) %>% dplyr::pull(end) %>%
          max() / 10 ^ 6
        Step = abs(data$RT[1, 'start'] - data$RT[1, 'end']) / 10 ^ 6

        shiny::updateSliderInput(
          inputId = 'range__scCN',
          min = 0,
          max = Max,
          value = c(0, Max),
          step = Step
        )

      }

    })
    ####create folder plots and set to save
    shiny::observeEvent(input$Save__scCN, {
      if (!dir.exists(file.path(variables$roots['OutputFolder'],
                                input$Analysis_Name, 'scCN'))) {
        dir.create(file.path(variables$roots['OutputFolder'],
                             input$Analysis_Name, 'scCN'),
                   recursive = T)
      }
      variables$Save__scCN = T

    })
    #disable all if something changes
    shiny::observeEvent(c(input$range__scCN,
                          input$Levels__scCN,
                          input$Save__scCN),
                        {
                          if (nrow(data$RT) > 0) {
                            G = unique(data$RT$group)

                            scCN_module_ls$ui = lapply(G, function(g)
                              scCN_ui(paste0('CN', g), title = g))
                            scCN_module_ls$server = lapply(G, function(g)
                              scCN_server(
                                paste0('CN', g),
                                S_Traks = data$S %>%
                                  dplyr::filter(
                                    group == g,
                                    chr == input$Chr__scCN,
                                    start >= input$range__scCN[1] *
                                      10 ^ 6,
                                    end <= input$range__scCN[2] *
                                      10 ^ 6
                                  ),
                                G_Traks = data$G %>%
                                  dplyr::filter(
                                    group == g,
                                    chr == input$Chr__scCN,
                                    start >= input$range__scCN[1] *
                                      10 ^ 6,
                                    end <= input$range__scCN[2] *
                                      10 ^ 6
                                  ),
                                Levels = input$Levels__scCN,
                                out = file.path(variables$roots['OutputFolder'],
                                                input$Analysis_Name, 'scCN'),
                                save = variables$Save__scCN
                              ))
                            if (variables$Save__scCN) {
                              variables$Save__scCN = F
                            }

                            output$scCN_UI = shiny::renderUI({
                              scCN_module_ls$ui
                            })

                          }
                        })
  }

  # # Dimensionality reduction
  shiny::observeEvent(input$Sidebar, {
    if (input$Sidebar == 'DRed' & ncol(data$G) != 0) {
      #call Dim_red_server
      Dred_module = Dim_red_sub_pop_server(
        id = 'Dred',
        scCN = rbind(
          data$S %>% dplyr::ungroup() %>% dplyr::mutate(Phase = 'S'),
          data$G %>% dplyr::ungroup() %>% dplyr::mutate(Phase = 'G1/G2')
        ),
        out =  file.path(variables$roots['OutputFolder'],
                         input$Analysis_Name),
        Inputfolder = data$folder_list,
        cores = input$cores,
        colors=data$samples_names_and_colors%>%
          dplyr::filter(
            type == 'Sample')
      )

    } else {
      Dred_module = NULL
    }
  })
  ###binprobrep start
  {
    #disable save if we are not
    shiny::observeEvent(input$Sidebar, {
      if (nrow(data$G) > 0 & input$Sidebar == 'BinRep') {
        shinyjs::enable('Save__BinRep')
      } else{
        shinyjs::disable('Save__BinRep')
        #store BinRep module info
        BinRep_module_ls$ui = NULL
        BinRep_module_ls$server = NULL
        if (variables$Save__BinRep_label != 'Plot') {
          variables$Save__BinRep_label = 'Plot'
          shiny::updateActionButton(inputId = 'Save__BinRep',
                                    label = variables$Save__BinRep_label)
          shinyjs::hide('Save_group__BinRep')
        }

      }
    })

    #recover binrep colors
    shiny::observeEvent(input$apply_color_changes_BinRep, {
      #recover colors
      ES_color_BinRep = colors_server(id = 'ES_color_BinRep')
      MS_color_BinRep = colors_server(id = 'MS_color_BinRep')
      LS_color_BinRep = colors_server(id = 'LS_color_BinRep')
      GG_color_BinRep = colors_server(id = 'GG_color_BinRep')

      variables$colors_BinRep = c(ES_color_BinRep(),
                                  MS_color_BinRep(),
                                  LS_color_BinRep(),
                                  GG_color_BinRep())
    })


    shiny::observeEvent(input$Save__BinRep, {
      if (input$Save__BinRep > 0) {
        #avoid user clicking twice
        shinyjs::disable('Save__BinRep')

        if (variables$Save__BinRep_label == 'Plot') {
          #calculate variability for G1/G2- and S- phase cells
          data$variabilityBR = rbind(
            Kronos.scRT::Prepare_G1G2_phase_cells_forBinRepProb(
              G1.G2 = data$G,
              RT = data$RT,
              quantile.range = input$BinRep_G1_Ploidy
            ),
            Kronos.scRT::Prepare_S_phase_cells_forBinRepProb(
              S = data$S,
              RT = data$RT,
              Early.cells = input$BinRep_Early_Cells,
              Mid.cells = input$BinRep_Mid_Cells,
              Late.cells = input$BinRep_Late_Cells
            )
          )
          #change label
          variables$Save__BinRep_label = 'Save'
          shiny::updateActionButton(inputId = 'Save__BinRep',
                                    label = variables$Save__BinRep_label)
          shinyjs::show('Save_group__BinRep')
        } else{
          variables$Save__BinRep = T
          if (!dir.exists(file.path(
            variables$roots['OutputFolder'],
            input$Analysis_Name,
            'BinsRepProb'
          ))) {
            dir.create(file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              'BinsRepProb'
            ))
          }
        }
        # how many windows
        Groups = unique(data$variabilityBR$group)
        #call modules
        BinRep_module_ls$ui = lapply(Groups, function(x)
          BinRepProb_ui(paste0('BinRepProb', x), title = x))

        sizes = hw_plot_server('hw_BinRep')
        sizes = sizes()

        if (variables$Save__BinRep) {
          variables$Save__BinRep = F
          BinRep_module_ls$server = lapply(Groups, function(x)
            BinRepProb_server(
              id = paste0('BinRepProb', x),
              variabilityBR = data$variabilityBR,
              out = file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                'BinsRepProb'
              ) ,
              colors = variables$colors_BinRep,
              sizes = sizes,
              file_name = x,
              save = T
            ))

        } else{
          BinRep_module_ls$server = lapply(Groups, function(x)
            BinRepProb_server(
              id = paste0('BinRepProb', x),
              variabilityBR = data$variabilityBR,
              out = file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                'BinsRepProb'
              ) ,
              colors = variables$colors_BinRep,
              sizes = sizes,
              file_name = x,
              save = F
            ))
        }
        output$BinRep_ui <- shiny::renderUI(BinRep_module_ls$ui)



        shinyjs::enable('Save__BinRep')


      }
    })

    shiny::observeEvent(
      c(
        input$BinRep_G1_Ploidy,
        input$BinRep_Early_Cells,
        input$BinRep_Mid_Cells,
        input$BinRep_Late_Cells,
        input$apply_color_changes_BinRep
      ),
      {
        variables$Save__BinRep_label = 'Plot'
        shiny::updateActionButton(inputId = 'Save__BinRep',
                                  label = variables$Save__BinRep_label)
        shinyjs::hide('Save_group__BinRep')

      }
    )

  }
  ### Twidth
  {
    shiny::observeEvent(input$Regions_tw, {
      if (input$Regions_tw == 'Customized categories') {
        shinyjs::show('div_load_regions_tw')
      } else{
        shinyjs::hide('div_load_regions_tw')
        data$GenomeAnnotationFile = dplyr::tibble()
        output$GenomeAnnotationFile = shiny::renderText({
          NULL
        })


      }
    })

    shiny::observe({
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'load_regions_tw',
        session = session,
        roots = variables$roots,
        defaultRoot = 'Home'
      )
    })
    shiny::observeEvent(input$load_regions_tw, {
      if (!is.numeric(input$load_regions_tw)) {
        data$GenomeAnnotationFile = shinyFiles::parseFilePaths(roots = variables$roots,
                                                               selection = input$load_regions_tw) %>%
          dplyr::mutate(
            State = Kronos.scRT::right_format(
              file_path = datapath,
              columns_to_check = c('chr', 'start', 'end', 'annotation'),
              delim = '\t',
              wrong_message = 'Wrong Format!',
              rigth_message = basename(datapath)
            )
          ) %>%
          dplyr::select('Name' = name, State, datapath)

        output$GenomeAnnotationFile = shiny::renderText({
          data$GenomeAnnotationFile$State
        })
      } else{
        data$GenomeAnnotationFile = dplyr::tibble()
        output$GenomeAnnotationFile = shiny::renderText({
          NULL
        })
      }
    })
    shiny::observeEvent(input$Sidebar, {
      if (nrow(data$Variability) > 0 & input$Sidebar == 'Twidth') {
        #load annotation
        if (input$Regions_tw == 'Costumized categores' &
            nrow(data$GenomeAnnotationFile) > 0) {
          if (data$GenomeAnnotationFile$State == 'Wrong Format!') {

          }
        } else if (input$Regions_tw == 'Costumized categores' &
                   nrow(data$GenomeAnnotationFile) == 0) {

        } else if (input$Regions_tw != 'Costumized categores') {

        }
      }
    })

    shiny::observeEvent(input$Sidebar, {
      #if right format upload file
      if (nrow(data$Variability) > 0 & input$Sidebar == 'Twidth') {
        if (nrow(data$GenomeAnnotationFile) > 0) {
          data$GenomeAnnotation_TW = readr::read_tsv(data$GenomeAnnotationFile$datapath)
        } else{
          data$GenomeAnnotation_TW = dplyr::tibble()
        }

        Twidth_module_ls$ui = lapply(unique(data$Variability$group), function(x)
          Twidth_ui(paste0('Twidth', x), title = x))
        output$TW_ui <- shiny::renderUI(Twidth_module_ls$ui)

        Twidth_module_ls$server = lapply(unique(data$Variability$group), function(x)
          Twidth_server(
            id = paste0('Twidth', x),
            file_name = x,
            variability = data$Variability %>% dplyr::filter(group == x),
            out = file.path(variables$roots['OutputFolder'],
                            input$Analysis_Name, 'Twidth'),
            GenomeAnnotation = shiny::isolate(data$GenomeAnnotation_TW),
            cores = input$cores
          ))
      }
    })
  }
  ####exit
  {
    shiny::observeEvent(input$Sidebar, {
      if (input$Sidebar == 'Exit') {
        shinyjs::js$closeWindow()
        shiny::stopApp()
      }
    })
  }

}

# Run the application
shiny::shinyApp(ui = ui,
                server = server)
