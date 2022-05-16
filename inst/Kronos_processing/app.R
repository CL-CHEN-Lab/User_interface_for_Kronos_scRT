#load operators
`%>%` = tidyr::`%>%`
`%dopar%` = foreach::`%dopar%`
`%do%` = foreach::`%do%`

#set theme for plots
ggplot2::theme_set(new = ggplot2::theme_bw())

#funtion to close app
jscode <- "shinyjs.closeWindow = function() { window.close(); }"

#number of available cores
maxCores = parallel::detectCores()

# Define UI for application that draws a histogram
ui <- shinydashboard::dashboardPage(
  title = 'Kronos scRT',
  skin = 'yellow',
  shinydashboard::dashboardHeader(title = shiny::span(
    shiny::img(src = 'KronosLogo.png', width = '100%')
  )),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      id = 'Sidebar',
      shinydashboard::menuItem(text = "Home",
                               tabName = "Home"),
      shinydashboard::menuItem(text = "Upload data",
                               tabName = "Upload"),
      shinydashboard::menuItem(text = "Diagnostic",
                               tabName = 'Diagnostic'),
      shinydashboard::menuItem(text = "Filter Cells",
                               tabName = "FilterCells"),
      shinydashboard::menuItem(
        text = "Variability",
        tabName = "Variability",
        shinydashboard::menuSubItem(text = 'Bin Replication', tabName = 'BinRep'),
        shinydashboard::menuSubItem(text = 'T-width', tabName = 'Twidth')
      ),
      shinydashboard::menuItem(text = "Dimensionality Reduction",
                               tabName = "DRed"),
      shinydashboard::menuItem(text = "scPlots",
                               tabName = "scPlots"),
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
      color = 'orange',
      height = '200px',
      width = '200px'
    ),

    shinydashboard::tabItems(#home
      {
        shinydashboard::tabItem(
          tabName = "Home",
          shinydashboard::box(
            width = 4,
            height = '150',
            title = 'Analysis name',
            solidHeader = T,
            background = 'black',
            shiny::textInput(
              inputId = 'Analysis_Name',
              label = '',
              value = 'Analysis'
            )
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = 'Cores to use for the analysis',
            solidHeader = T,
            background = 'black',
            align = 'center',
            #cores
            sliderInput(
              width = '100%',
              inputId = 'cores',
              label = '',
              value = trunc(maxCores / 2),
              min = 1,
              max = maxCores,
              step = 1,
              ticks = F
            )
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = 'Output folder',
            solidHeader = T,
            background = 'black',
            align = 'center',
            shinyFiles::shinyDirButton(
              id = 'Output_dir',
              label = 'Output folder',
              title = 'Output folder'
            ),
            shiny::htmlOutput('Output_dir_out')
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = shiny::div(
              'Chr prefix to use',
              bsplus::shiny_iconlink() %>%
                bsplus::bs_embed_popover(title = 'Depending on your reference genome, chromosomes or scaffolds could be named differntly. Please, privide the prefix that preceeds the chromosome or scaffold number. If no prefix is present this box has to be emptied.' , placement = 'right')
            ),
            solidHeader = T,
            background = 'black',
            shiny::textInput(
              inputId = 'chr_prefix',
              label = '',
              value = 'chr'
            )
          ),
          shiny::column(
            width = 4,
            height = '150',
            align = 'center',
            shiny::img(src = "Kronos_logo.png", height = 150)
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = shiny::div(
              'Chr to include',
              bsplus::shiny_iconlink() %>%
                bsplus::bs_embed_popover(title = 'Chromosome to include in the analysis. You can define ranges using ":" and individual chromosomes using ",". Don\'t include the chr prefix!' , placement = 'right')
            ),
            solidHeader = T,
            background = 'black',
            shiny::textInput(
              inputId = 'chr_range',
              label = '',
              value = '1:22'
            )
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = shiny::div(
              'Chr size file',
              bsplus::shiny_iconlink() %>%
                bsplus::bs_embed_popover(title = 'For most genome assemblies this file can be downloaded from UCSC Genome Browser website.' , placement = 'right')
            ),
            solidHeader = T,
            background = 'black',
            align = 'center',
            shinyFiles::shinyFilesButton(
              id = 'Chr_size',
              label = 'Chr size file',
              title = 'Chr size file',
              multiple = F
            ),
            shiny::htmlOutput('Chr_size_out', style = 'font-size: 16px;font-weight: bold;')
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = 'Final RT bin size',
            solidHeader = T,
            background = 'black',
            align = 'center',
            #cores
            sliderInput(
              width = '100%',
              inputId = 'binsize',
              label = '',
              value = '0.5',
              min = 0.1,
              max = 2,
              step = 0.05,
              ticks = F,
              post = 'Mb'
            )
          ),
          shinydashboard::box(
            width = 4,
            height = '150',
            title = 'Apply settings',
            solidHeader = T,
            background = 'black',
            align = 'center',
            shiny::actionButton(inputId =  'ApplySettings', label =
                                  'Apply settings')
          )
        )

      },
      #upload
      {
        shinydashboard::tabItem(
          tabName = "Upload",
          shinydashboard::box(
            title = 'Upload samples',
            solidHeader = T,
            status = 'primary',
            collapsible = T,
            width = 12,
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shinyFiles::shinyFilesButton(
                  id = 'PerCell_file',
                  label = shiny::div('PerCell file',
                                     bsplus::shiny_iconlink() %>%
                                       bsplus::bs_embed_popover(title = '_PerCell.csv file created by Kronos scRT processing.' , placement = 'right')
                  ),
                  title = 'PerCell File',
                  multiple = F,
                  style = "width: 100%"
                )
              ),
              shiny::column(width = 9,
                            shiny::htmlOutput(outputId =
                                                'PerCell_file_out'))
            ),
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shinyFiles::shinyFilesButton(
                  id = 'scCN',
                  label = shiny::div('sc tracks file',
                                     bsplus::shiny_iconlink() %>%
                                       bsplus::bs_embed_popover(title = '_scCNV.tsv file created by Kronos scRT processing.' , placement = 'right')
                  ),
                  title = 'sc tracks file',
                  multiple = F,
                  style = "width: 100%"
                )
              ),
              shiny::column(width = 9,
                            shiny::htmlOutput(outputId =
                                                'scCN_out'))
            ),
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shinyFiles::shinyFilesButton(
                  id = 'setting_file',
                  label = shiny::div(
                    'Setting file',
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'Optinal. If this dataset has been used in a previous analysis the produced setting file can be reused.' , placement = 'right')
                  ),
                  title =  "Setting file",
                  multiple = F,
                  style = "width: 100%"
                )
              ),
              shiny::column(width = 9,
                            shiny::htmlOutput(outputId =
                                                'setting_file_out'))
            ),
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shinyFiles::shinyFilesButton(
                  id = 'whoiswho_file',
                  label = shiny::div(
                    "Who's who file",
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'Optinal. A file containing experimental cell staging information. Such a file has to containin cell identifiers under the colums Cell and logical values under the colums S_Phase. (TRUE = a cell is in S phase, FALSE = a cell is in G1 or G2 phase)' , placement = 'right')
                  ),
                  title =  "Who's who file",
                  multiple = F,
                  style = "width: 100%"
                )
              ),
              shiny::column(width = 9,
                            shiny::htmlOutput(outputId = 'whoiswho_file_out'))
            ),
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::textInput(
                  inputId = ('FileName'),
                  value =  'Exp',
                  label = shiny::div(
                    "Basename",
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'This Name identifies each individual experiment' , placement = 'right')
                  ),
                  width = '100%'
                )
              ),
              shiny::column(
                width = 6,
                shiny::textInput(
                  inputId = ('GroupName'),
                  value =  'Exp',
                  label = shiny::div(
                    "Group name",
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'If cells have been sequenced in different experiments, providind the same group name allows to merge them after normalisation.' , placement = 'right')
                  ),
                  width = '100%'
                )
              )
            ),
            shiny::fluidRow(shiny::column(
              width = 3,
              shiny::actionButton(
                inputId = 'Add_sample',
                label = "Add sample",
                width = '100%'
              )
            )),
            shiny::fluidRow(shiny::column(
              width = 12,
              shiny::tableOutput('File_paths')
            )),
            shiny::fluidRow(shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = 'Remove_sample',
                label = "Remove sample",
                choices = 'Sample',
                selected = 'Sample',
                width = '100%'
              )
            ))
          ),
          shinydashboard::box(
            title = shiny::div(
              'Upload references',
              bsplus::shiny_iconlink() %>%
                bsplus::bs_embed_popover(title = 'Bulk RT references are optional. The program allows to upload maximum one reference per group.' , placement = 'right')
            ),
            solidHeader = T,
            status = 'primary',
            collapsible = T,
            collapsed = T,
            width = 12,
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::selectInput(
                  inputId = 'reference_group',
                  label = 'group name',
                  choices = 'Select a group',
                  multiple = F,
                  selected = 'Select a group',
                  width = '100%'
                )
              ),
              shiny::column(
                width = 6,
                shiny::textInput(
                  inputId = 'reference_name',
                  value =  'Reference',
                  label = shiny::div(
                    'Reference name',
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'This name will be visualize on plots.' , placement = 'right')
                  ),
                  width = '100%'
                )
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shinyFiles::shinyFilesButton(
                  id = 'reference_file',
                  label = "Reference RT file",
                  title =  "Reference RT file",
                  multiple = F,
                  style = "width: 100%"
                )
              ),
              shiny::column(width = 9,
                            shiny::htmlOutput(outputId =
                                                'reference_file_out'))
            ),
            shiny::fluidRow(shiny::column(
              width = 3,
              shiny::actionButton(
                inputId = 'Add_reference',
                label = "Add reference",
                width = '100%'
              )
            )),
            shiny::fluidRow(shiny::column(
              width = 12,
              shiny::tableOutput('Reference_paths')
            ))
          ),
          shiny::fluidRow(shiny::column(
            width = 3,
            shinyjs::disabled(
              shiny::actionButton(
                inputId = 'Done_upload',
                label = 'Done',
                width = '100%'
              )
            )
          ))
        )
      },
      #Diagnostic
      {
        shinydashboard::tabItem(tabName = 'Diagnostic',
                                shiny::uiOutput('Diagnostic_UI'),
                                shiny::fluidRow(shiny::column(
                                  width = 3,
                                  shiny::actionButton(
                                    inputId = 'Diagnostic_done',
                                    label = "Start RT analysis",
                                    width = '100%'
                                  )
                                )))

      },
      #FilterCells
      {
        shinydashboard::tabItem(
          tabName = 'FilterCells',
          shiny::fluidRow(shiny::column(
            width = 6,
            offset = 3,
            shiny::sliderInput(
              width = '100%',
              inputId = 'min_correlation_matrix',
              label = 'Min correlation to be kept',
              min = 0,
              max = 1,
              step = 0.01,
              value = 0.25
            )

          )),
          shiny::fluidRow(
            shiny::column(width = 6,
                          plotly::plotlyOutput('plot1_RT')),
            shiny::column(width = 6,
                          plotly::plotlyOutput('plot2_RT'))
          ),

          shiny::fluidRow(shiny::column(
            width = 3,
            shiny::actionButton(
              inputId = 'RT_next',
              label = "Finalise",
              width = '100%'
            )
          )),
          shiny::fluidRow(
            shiny::column(
              width = 12,
              align = 'center',
              shiny::plotOutput('plot3_RT', width = '100%')
            )
          )

        )
      },
      #twidth
      {
        shinydashboard::tabItem(tabName = "Twidth",
                                shiny::uiOutput('TW_ui'))
      },
      #bin prob
      {
        shinydashboard::tabItem(tabName = "BinRep",
                                shiny::uiOutput('BinRep_ui'),
                                shiny::fluidRow(shiny::column(
                                  width = 2,
                                  shiny::actionButton(
                                    inputId = 'Save__BinRep',
                                    label = 'Save',
                                    width = '100%'
                                  )
                                )))
      },
      #Dred
      {
        shinydashboard::tabItem(tabName = "DRed",
                                Dim_red_ui('Dred'))
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
                label = 'Chrmosome',
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
                inputId = 'range__scPlot',
                label = 'Coodrinates',
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
            ),
            shiny::column(
              width = 2,
              shiny::actionButton(
                inputId = 'Save__scPlot',
                label = 'Save',
                width = '100%'
              )
            )
          ),
          shiny::uiOutput('scPlots_UI')
        )
      })
  )
)


server <- function(input, output, session) {
  #variables
  variables = shiny::reactiveValues(
    setting_file = '',
    File_paths = dplyr::tibble(),
    Reference_paths = dplyr::tibble(),
    PerCell_file = '',
    whoiswho_file = '',
    PerCell_file_error = '',
    first_line_Percell = NULL,
    scCN = '',
    scCN_error = '',
    first_line_CNV = NULL,
    reference_file = '',
    groups = '',
    groups_with_a_reference = '',
    Chr_size = '<font color=\"#FF0000\"><b> This is a mandatory input.</b></font>',
    roots = c(
      shinyFiles::getVolumes()(),
      Home = Sys.getenv("HOME"),
      OutputFolder = file.path(Sys.getenv("HOME"))
    ),
    to_keep = NULL,
    Save__scPlots = F,
    Save__BinRep = F

  )
  #store diagnostic module info
  diagnostic_module_ls = shiny::reactiveValues(ui = list(),
                                               menu = list(),
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
    PerCell = dplyr::tibble(),
    scCN = dplyr::tibble(),
    whoswho = dplyr::tibble(),
    setting = dplyr::tibble(),
    reference = dplyr::tibble(),
    scRT = dplyr::tibble(),
    variability = dplyr::tibble(),
    variabilityBR = dplyr::tibble(),
    G1G2 = dplyr::tibble()
  )

  #load
  shiny::observe({
    #add tabitems in diagnostic
    output$Diagnostic_UI = shiny::renderUI({
      diagnostic_module_ls$ui
    })
  })
  #stop app when the session ends
  session$onSessionEnded(function() {
    shiny::stopApp()
  })

  #home
  {
    shiny::observe({
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'Chr_size',
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
      })

      output$Output_dir_out <-
        renderText(paste(
          '<H4><b>',
          file.path(variables$roots['OutputFolder'], input$Analysis_Name),
          '</H4></b>'
        ))

      if (!stringr::str_detect(string = variables$Chr_size,
                               pattern = '<font color=\"#FF0000\">')) {
        shinyjs::enable('Upload_chr_size')
      } else{
        shinyjs::disable('Upload_chr_size')
      }

    })

    #load chrs_size
    shiny::observeEvent(input$Chr_size, {
      if (!is.numeric(input$Chr_size)) {
        #get path
        variables$Chr_size = dplyr::tibble(
          shinyFiles::parseFilePaths(
            roots = variables$roots ,
            selection = input$Chr_size
          )
        ) %>% dplyr::pull(datapath)

        variables$Chr_size = Kronos.scRT::right_format(
          file_path = variables$Chr_size,
          delim = '\t',
          columns_to_check = 2,
          wrong_message = paste0(
            '<b>',
            variables$Chr_size,
            '<font color=\"#FF0000\"> does not have the right format. </b></font>'
          ),
          rigth_message = variables$Chr_size
        )

      } else{
        variables$Chr_size = '<font color=\"#FF0000\"><b> This is a mandatory input.</b></font>'
      }

    })


    shiny::observe(output$Chr_size_out <-
                     renderText(variables$Chr_size))

    shiny::observeEvent(input$Output_dir, {
      if (!is.numeric(input$Output_dir)) {
        variables$roots = c(
          shinyFiles::getVolumes()(),
          Home = Sys.getenv("HOME"),
          OutputFolder = parseDirPath(
            roots = variables$roots,
            selection = input$Output_dir
          )
        )
      }
    })

    shiny::observe({
      if (!is.numeric(input$Chr_size) & input$ApplySettings == 0) {
        shinyjs::enable('ApplySettings')
      } else{
        shinyjs::disable('ApplySettings')
      }
    })

    shiny::observeEvent(input$ApplySettings, {
      if (input$ApplySettings == 1) {
        shinyjs::disable('ApplySettings')
        shinyjs::disable('binsize')
        shinyjs::disable('Chr_size')
        shinyjs::disable('chr_range')
        shinyjs::disable('chr_prefix')
        shinyjs::disable('Output_dir')
        shinyjs::disable('cores')
        shinyjs::disable('Analysis_Name')

        #move to next tab
        shinydashboard::updateTabItems(session = session,
                                       inputId = 'Sidebar',
                                       selected = 'Upload')

      }
    })
  }
  #load data
  {
    #samples
    {
      shiny::observe({
        shinyFiles::shinyFileChoose(
          input = input,
          id = 'PerCell_file',
          roots = variables$roots,
          session = session,
          defaultRoot =  'Home'
        )
        shinyFiles::shinyFileChoose(
          input = input,
          id = 'scCN',
          roots = variables$roots,
          session = session,
          defaultRoot = 'Home'
        )
        shinyFiles::shinyFileChoose(
          input = input,
          id = 'whoiswho_file',
          roots = variables$roots,
          session = session,
          defaultRoot =  'Home'
        )
        shinyFiles::shinyFileChoose(
          input = input,
          id = 'setting_file',
          roots = variables$roots,
          session = session,
          defaultRoot = 'Home'
        )
      })

      #activate deactivate Add input butto
      shiny::observe({
        if (variables$PerCell_file == '' | variables$scCN == '' |
            stringr::str_detect(string = variables$PerCell_file, pattern = ' does not have the right format') |
            stringr::str_detect(string = variables$scCN, pattern = ' does not have the right format') |
            stringr::str_detect(string = variables$whoiswho_file, pattern = ' does not have the right format') |
            input$ApplySettings != 1) {
          shinyjs::disable(id = 'Add_sample')
          output$PerCell_file_out = shiny::renderText(variables$PerCell_file)
          output$scCN_out = shiny::renderText(variables$scCN)
        } else{
          if (all(!c('basename', 'group') %in% names(variables$first_line_Percell)) &
              all(!c('basename', 'group') %in% names(variables$first_line_CNV))) {
            #rest outputs if they had changed
            output$PerCell_file_out = shiny::renderText(variables$PerCell_file)
            output$scCN_out = shiny::renderText(variables$scCN)

            shinyjs::enable(id = 'Add_sample')

          } else if (all(c('basename', 'group') %in% names(variables$first_line_Percell)) &
                     all(c('basename', 'group') %in% names(variables$first_line_CNV))) {
            if (variables$first_line_Percell$basename == variables$first_line_CNV$basename &
                variables$first_line_Percell$group == variables$first_line_CNV$group) {
              #update file name and group name if they match

              updateTextInput(
                session = session,
                inputId = 'FileName',
                value = variables$first_line_Percell$basename
              )
              updateTextInput(
                session = session,
                inputId = 'GroupName',
                value = variables$first_line_Percell$group
              )

              #rest outputs if they had changed
              output$PerCell_file_out = shiny::renderText(variables$PerCell_file)
              output$scCN_out = shiny::renderText(variables$scCN)

              shinyjs::enable(id = 'Add_sample')

            } else{
              #sebd worning

              output$PerCell_file_out = shiny::renderText(
                '<font color=\"#FF0000\"> Basename and Group do not match between PerCell and sc tracks files!!! </b></font>'
              )
              output$scCN_out = shiny::renderText(
                '<font color=\"#FF0000\"> Basename and Group do not match between PerCell and sc tracks files!!! </b></font>'
              )
              shinyjs::disable(id = 'Add_sample')
            }
          } else{
            output$PerCell_file_out = shiny::renderText(
              '<font color=\"#FF0000\"> Basename and Group do not match between PerCell and sc tracks files!!! </b></font>'
            )
            output$scCN_out = shiny::renderText(
              '<font color=\"#FF0000\"> Basename and Group do not match between PerCell and sc tracks files!!! </b></font>'
            )
            shinyjs::disable(id = 'Add_sample')
          }
        }
        #activate button to srat RT analysis
        if (nrow(variables$File_paths) > 0) {
          shinyjs::enable(id = 'Diagnostic_done')
          shinyjs::enable(id = 'Done_upload')
        } else{
          shinyjs::disable(id = 'Diagnostic_done')
          shinyjs::disable(id = 'Done_upload')
        }

      })

      #move to next page
      shiny::observeEvent(input$Done_upload, {
        shinydashboard::updateTabItems(session = session,
                                       inputId = 'Sidebar',
                                       selected = 'Diagnostic')

      })


      #replace ' ' with '_'
      shiny::observeEvent(input$FileName, {
        updateTextInput(
          session = session,
          inputId = 'FileName',
          value = stringr::str_replace_all(
            string = input$FileName,
            pattern = ' ',
            replacement = '_'
          )
        )

      })
      shiny::observeEvent(input$GroupName, {
        updateTextInput(
          session = session,
          inputId = 'GroupName',
          value = stringr::str_replace_all(
            string = input$GroupName,
            pattern = ' ',
            replacement = '_'
          )
        )
      })

      # get PerCell_file path and check if it has the right format
      shiny::observeEvent(input$PerCell_file, {
        if (!is.numeric(input$PerCell_file)) {
          #get path
          variables$PerCell_file = dplyr::tibble(
            shinyFiles::parseFilePaths(
              roots = variables$roots ,
              selection = input$PerCell_file
            )
          ) %>% dplyr::pull(datapath)


          #check if it is the right file
          variables$PerCell_file = Kronos.scRT::right_format(
            file_path = variables$PerCell_file,
            delim = ',',
            columns_to_check = c(
              "Cell",
              "normalized_dimapd",
              "mean_ploidy",
              "ploidy_confidence",
              "is_high_dimapd",
              "is_noisy",
              "coverage_per_1Mbp"
            ),
            wrong_message = paste0(
              '<b>',
              variables$PerCell_file,
              '<font color=\"#FF0000\"> does not have the right format. </b></font>'
            ),
            rigth_message = variables$PerCell_file
          )
          #read line one and look for group and basename
          if (stringr::str_detect(
            string = variables$PerCell_file,
            pattern = ' does not have the right format',
            negate = T
          )) {
            variables$first_line_Percell = readr::read_csv(variables$PerCell_file, n_max = 1)
          }

        } else{
          variables$PerCell_file = ''
        }
      })

      # get scCN path and check if it has the right format
      shiny::observeEvent(input$scCN, {
        if (!is.numeric(input$scCN)) {
          #get path
          variables$scCN = dplyr::tibble(
            shinyFiles::parseFilePaths(roots = variables$roots ,
                                       selection = input$scCN)
          ) %>% dplyr::pull(datapath)

          #check if it is the right file
          variables$scCN = Kronos.scRT::right_format(
            file_path = variables$scCN,
            delim = '\t',
            columns_to_check = c('chr',
                                 'start',
                                 'end' ,
                                 'copy_number',
                                 'reads' ,
                                 'Cell'),
            wrong_message = paste0(
              '<b>',
              variables$scCN,
              '<font color=\"#FF0000\"> does not have the right format. </b></font>'
            ),
            rigth_message = variables$scCN
          )

          #read line one and look for group and basename
          if (stringr::str_detect(
            string = variables$scCN,
            pattern = ' does not have the right format',
            negate = T
          )) {
            variables$first_line_CNV = readr::read_tsv(variables$scCN, n_max = 1)
          }

        } else{
          variables$scCN = ''
        }
      })

      # get whoiswho path and check if it has the right format
      shiny::observeEvent(input$whoiswho_file, {
        if (!is.numeric(input$whoiswho_file)) {
          #get path
          variables$whoiswho_file = dplyr::tibble(
            shinyFiles::parseFilePaths(
              roots = variables$roots ,
              selection = input$whoiswho_file
            )
          ) %>% dplyr::pull(datapath)

          #check if it is the right file
          variables$whoiswho_file = Kronos.scRT::right_format(
            file_path = variables$whoiswho_file,
            delim = '\t',
            columns_to_check = c('Cell', 'Phase'),
            wrong_message = paste0(
              '<b>',
              variables$whoiswho_file,
              '<font color=\"#FF0000\"> does not have the right format. </b></font>'
            ),
            rigth_message = variables$whoiswho_file
          )

          output$whoiswho_file_out = shiny::renderText(variables$whoiswho_file)

        } else{
          variables$whoiswho_file = ''
          output$whoiswho_file_out = shiny::renderText(variables$whoiswho_file)

        }
      })

      # get setting path and check if it has the right format
      shiny::observeEvent(input$setting_file, {
        if (!is.numeric(input$setting_file)) {
          #get path
          variables$setting_file = dplyr::tibble(
            shinyFiles::parseFilePaths(
              roots = variables$roots ,
              selection = input$setting_file
            )
          ) %>% dplyr::pull(datapath)

          #check if it is the right file
          variables$setting_file = Kronos.scRT::right_format(
            file_path = variables$setting_file,
            delim = '\t',
            columns_to_check = c(
              'threshold_Sphase' ,
              'threshold_G1G2phase' ,
              'Sphase_first_part',
              'Sphase_second_part' ,
              'RPMPH_TH',
              'RPM_TH'
            ),
            wrong_message = paste0(
              '<b>',
              variables$setting_file,
              '<font color=\"#FF0000\"> does not have the right format. </b></font>'
            ),
            rigth_message = variables$setting_file
          )

          output$setting_file_out = shiny::renderText(variables$setting_file)

        } else{
          variables$setting_file = ''
          output$setting_file_out = shiny::renderText(variables$setting_file)

        }
      })
      #add samples to dataframe
      shiny::observeEvent(input$Add_sample, {
        shinyjs::disable(id = 'Add_sample')
        shinyjs::disable(id = 'Remove_sample')


        #create directory if it does not exist
        if (!dir.exists(file.path(variables$roots['OutputFolder'], input$Analysis_Name))) {
          dir.create(file.path(variables$roots['OutputFolder'], input$Analysis_Name),
                     recursive = T)
        }

        #check if it is a new group
        if (!(
          input$GroupName %in% variables$groups |
          input$GroupName %in% variables$groups_with_a_reference
        )) {
          variables$groups = c(variables$groups, input$GroupName)
        }
        #new line
        newline = dplyr::tibble(
          PerCellFile = basename(variables$PerCell_file),
          scCNV = basename(variables$scCN),
          basename = input$FileName,
          group = input$GroupName,
          WhoSWho = ifelse(
            !is.numeric(input$whoiswho_file),
            basename(variables$whoiswho_file),
            ''
          ),
          Setting = ifelse(
            !is.numeric(input$setting_file),
            basename(variables$setting_file),
            ''
          )
        )

        #load data
        #if whoiswho_file's been provided
        if (!is.numeric(input$whoiswho_file)) {
          data$PerCell = rbind(
            data$PerCell,
            readr::read_csv(variables$PerCell_file, col_types = 'cnnnlln') %>%
              dplyr::mutate(
                basename = input$FileName,
                group = input$GroupName
              )
          ) %>%
            dplyr::inner_join(readr::read_tsv(variables$whoiswho_file, col_types = 'cc'),
                              by = c("Cell")) %>%
            dplyr::mutate(
              is_high_dimapd = ifelse(Phase == 'S', T, F),
              is_noisy = ifelse(is_high_dimapd, T, is_noisy)
            ) %>%
            dplyr::select(-Phase)

        } else{
          data$PerCell = rbind(
            data$PerCell,
            readr::read_csv(variables$PerCell_file, col_types = 'cnnnlln') %>%
              dplyr::mutate(
                basename = input$FileName,
                group = input$GroupName
              )
          )
        }
        data$scCN = rbind(
          data$scCN,
          readr::read_tsv(variables$scCN, col_types = 'ccnnnn') %>%
            dplyr::mutate(
              basename = input$FileName,
              group = input$GroupName
            )
        )

        if (!is.numeric(input$whoiswho_file) &
            variables$whoiswho_file != '') {
          data$whoswho = rbind(
            data$whoswho,
            readr::read_tsv(variables$whoiswho_file, col_types = 'cc') %>%
              dplyr::mutate(
                basename = input$FileName,
                group = input$GroupName
              )
          )
        }
        #if setting_file's been provided
        if (!is.numeric(input$setting_file) &
            variables$setting_file != '') {
          data$setting = rbind(
            data$setting,
            readr::read_tsv(variables$setting_file, col_types = 'nnnnnnncc') %>%
              dplyr::mutate(
                basename = input$FileName,
                group = input$GroupName
              )
          )
        }

        #add line to dataframe
        variables$File_paths = rbind(variables$File_paths,
                                     newline)

        #call diagnosti module
        diagnostic_module_ls$ui[[paste0(input$FileName)]] = diagnostic_ui(id = input$FileName)
        diagnostic_module_ls$server[[paste0(input$FileName)]] = diagnostic_server(
          id = input$FileName,
          PerCell = data$PerCell %>%
            dplyr::filter(basename ==
                            input$FileName),
          #if setting is not there return embpyt tibble
          Setting = tryCatch(
            data$setting %>%
              dplyr::filter(basename ==
                              input$FileName),
            error =
              function(x)
                return(dplyr::tibble())
          ),
          Cores = input$cores,
          OutFolder = file.path(variables$roots['OutputFolder'], input$Analysis_Name)
        )

        #upadte Remove_sample
        shiny::updateSelectInput(
          inputId = 'Remove_sample',
          choices = c('Sample', variables$File_paths$basename)
        )

        # reset inputs
        shiny::updateTextInput(
          session = session,
          inputId = 'FileName',
          value = paste('Exp', input$Add_sample, sep = '_')
        )
        shiny::updateTextInput(
          session = session,
          inputId = 'GroupName',
          value = paste('Exp', input$Add_sample, sep = '_')
        )
        variables$PerCell_file = ''
        variables$scCN = ''
        variables$whoiswho_file = ''
        variables$setting_file = ''
        output$scCN_out = shiny::renderText('')
        output$PerCell_file_out = shiny::renderText('')
        output$whoiswho_file_out = shiny::renderText('')
        output$setting_file_out = shiny::renderText('')

        output$File_paths = shiny::renderTable(variables$File_paths)

        #update groups for references
        shiny::updateSelectInput(
          inputId = 'reference_group',
          choices = c('Select a group',
                      variables$groups),
          selected = 'Select a group'
        )

        #reactivate remove sample
        shinyjs::enable('Remove_sample')
      })

      #remove samples
      shiny::observeEvent(input$Remove_sample, {
        if (input$Remove_sample != 'Sample') {
          #remove
          variables$File_paths = variables$File_paths %>%
            dplyr::filter(basename != input$Remove_sample)
          # stop module
          diagnostic_module_ls$ui[[input$Remove_sample]] = NULL
          diagnostic_module_ls$server[[input$Remove_sample]] = NULL

          #upadte Remove_sample
          shiny::updateSelectInput(
            inputId = 'Remove_sample',
            choices = c('Sample', variables$File_paths$basename),
            selected = 'Sample'
          )

          #update groups for Reference
          variables$groups = variables$groups[variables$groups %in% unique(variables$File_paths$group)]

          shiny::updateSelectInput(
            inputId = 'reference_group',
            choices = c('Select a group',
                        variables$groups),
            selected = 'Select a group'
          )

        }

      })
      shiny::observeEvent(input$Diagnostic_done, {
        #disable button
        shinyjs::disable('Diagnostic_done')
        #find all the setting files created by diagnostic
        files_to_load = list.files(
          path = file.path(variables$roots['OutputFolder'], input$Analysis_Name),
          pattern = '_settings.txt',
          full.names = T
        )
        # select those with the right format
        files_to_load = files_to_load[unlist(
          lapply(
            files_to_load,
            Kronos.scRT::right_format,
            delim = '\t',
            columns_to_check = c(
              'threshold_Sphase' ,
              'threshold_G1G2phase' ,
              'Sphase_first_part',
              'Sphase_second_part' ,
              'Ploidy',
              'RPMPH_TH',
              'RPM_TH',
              'basename',
              'group'
            ),
            wrong_message = F,
            rigth_message = T
          )
        )]

        #upload new setting file
        data$setting = Kronos.scRT::load_multiple_df(dirs = files_to_load,
                                                     delim = '\t',
                                                     col_types = 'nnnnnnncc')

        ########RT Start first part
        {
          variables$chr_list = paste0(input$chr_prefix,
                                      unlist(
                                        Kronos.scRT::String_to_Range(x = stringr::str_split(input$chr_range, ',')[[1]])
                                      ))

          #upload chr sizes
          data$Chr_size = readr::read_tsv(
            variables$Chr_size,
            col_names = c('chr', 'size'),
            col_types = 'cn'
          ) %>%
            dplyr::filter(chr %in% variables$chr_list) %>%
            dplyr::mutate(chr = factor(x =  chr, levels = variables$chr_list)) %>%
            tidyr::drop_na()

          variables$chr_list = variables$chr_list[variables$chr_list %in% unique(data$Chr_size$chr)]

          #Apply Settings to PerCell

          data$PerCell = Kronos.scRT::AdjustPerCell(
            PerCell = data$PerCell ,
            Settings = data$setting,
            Basename_leves = variables$File_paths$basename ,
            Group_leves = unique(variables$File_paths$group)
          )


          #Correct scCN based on PerCell

          data$scCN = Kronos.scRT::AdjustCN(
            scCN = data$scCN,
            PerCell =  data$PerCell ,
            Basename_leves = variables$File_paths$basename ,
            Group_leves = unique(variables$File_paths$group),
            Chr_filter = variables$chr_list
          )

          #create bins
          bins = Kronos.scRT::GenomeBinning(
            Chr_size = data$Chr_size,
            size = input$binsize * 10 ^ 6,
            Chr_filter = variables$chr_list,
            Cores = input$cores
          )
          #rebin data
          data$signal_smoothed = Kronos.scRT::Rebin(
            PerCell = data$PerCell,
            scCN = data$scCN,
            Bins = bins,
            Sphase = T
          )

          G1G2_smoothed = Kronos.scRT::Rebin(
            PerCell = data$PerCell,
            scCN = data$scCN,
            Bins = bins,
            Sphase = F
          )

          #calculate background
          median_G1G2_profile = Kronos.scRT::BackGround(G1_scCN = G1G2_smoothed)

          # create single G1/G2 single cell file
          G1G2_smoothed = Kronos.scRT::Replication_state(
            Samples = G1G2_smoothed,
            background = median_G1G2_profile,
            Chr_filter = variables$chr_list,
            cores = input$cores
          )

          #reshape and save
          G1G2_smoothed = G1G2_smoothed %>%
            dplyr::select(
              chr,
              start,
              end,
              CN,
              background,
              CN_bg,
              th,
              Rep,
              PercentageReplication,
              Cell,
              basename,
              group,
              newIndex
            )

          G1G2_smoothed %>%
            readr::write_tsv(
              file = file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                paste0(
                  input$Analysis_Name,
                  '_G1_G2_single_cells_CNV_',
                  input$binsize,
                  'Mb.tsv'
                )
              ),
              col_names = T
            )

          data$G1G2 = G1G2_smoothed

          if (nrow(data$reference) != 0) {
            # rebin reference RT
            data$reference = Kronos.scRT::RebinRT(
              RT = data$reference,
              Bins = bins,
              Chr_filter = variables$chr_list
            )
            #write output
            readr::write_tsv(
              x = data$reference,
              file = file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                paste0(
                  input$Analysis_Name,
                  '_reference_replication_timing_',
                  input$binsize,
                  'Mb.tsv'
                )
              ),
              col_names = T
            )
          }

          #Calculate replication state S
          data$signal_smoothed = Kronos.scRT::Replication_state(
            Samples = data$signal_smoothed,
            background = median_G1G2_profile,
            Chr_filter = variables$chr_list,
            cores = input$cores
          )

          # remove control track
          rm('median_G1G2_profile')

          #matrix for the correlation
          data$signal_smoothed = data$signal_smoothed %>%
            dplyr::ungroup() %>%
            dplyr::arrange(group, newIndex) %>%
            tidyr::unite(index, c(group, newIndex), sep = ' _ ') %>%
            dplyr::mutate(index = factor(index, levels = unique(index)))
          mat = data$signal_smoothed %>%
            tidyr::unite(pos, c(chr, start), sep = ':') %>%
            dplyr::mutate(Rep = as.numeric(Rep)) %>%
            dplyr::select(pos, index, Rep) %>%
            tidyr::spread(key = index, value = Rep) %>%
            tibble::column_to_rownames('pos') %>%
            dplyr::filter(complete.cases(.)) %>%
            as.matrix()

          #index
          variables$Index = colnames(mat)
          #correlation similarity distance
          variables$results = 1 - as.matrix(ade4::dist.binary(
            t(mat),
            method = 2,
            diag = T,
            upper = T
          ))
          variables$basenames = dplyr::tibble(Group = stringr::str_remove(colnames(mat), ' _ [0-9]{1,10}$'))

          #write matrix and plot heatmap before filtering
          saveRDS(
            object = variables$results,
            file = file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              paste0(
                input$Analysis_Name,
                '_correlation_per_cell_before_filtering.rds'
              )
            )
          )

          #prepare color patterns
          variables$color = grDevices::colorRampPalette(
            colors = c(
              "#00204DFF",
              "#233E6CFF",
              "#575C6DFF",
              "#7C7B78FF",
              "#A69D75FF",
              "#D3C164FF",
              "#FFEA46FF"
            )
          )


          output$plot1_RT <- plotly::renderPlotly(
            heatmaply::heatmaply(
              x = variables$results,
              colors = variables$color,
              dendrogram = F,
              showticklabels = F,
              row_side_colors = variables$basenames,
              col_side_colors = variables$basenames,
              limits = c(0, 1)
            ) %>%
              plotly::layout(
                showlegend = FALSE,
                legend = FALSE,
                annotations = list(visible = FALSE)
              )
          )

          variables$to_keep = foreach::foreach(i = unique(variables$basenames$Group)) %do% {
            sub_mat = variables$results[variables$basenames$Group == i, variables$basenames$Group == i]
            diag(sub_mat) = 0
            ! matrixStats::rowQuantiles(x = sub_mat,
                                        probs = 0.60,
                                        na.rm = T) <= input$min_correlation_matrix
          }

          variables$to_keep = unlist(variables$to_keep)
          variables$results_after_filtering = variables$results[variables$to_keep, variables$to_keep]
          variables$basenames_after_filtering = variables$basenames[variables$to_keep,]

          output$plot2_RT <- plotly::renderPlotly(
            heatmaply::heatmaply(
              x = variables$results_after_filtering,
              colors = variables$color,
              dendrogram = F,
              showticklabels = F,
              row_side_colors = variables$basenames_after_filtering,
              col_side_colors = variables$basenames_after_filtering,
              limits = c(0, 1)
            ) %>%
              plotly::layout(
                showlegend = FALSE,
                legend = FALSE,
                annotations = list(visible = FALSE)
              )
          )
          #bring the user to FilterCells
          shinydashboard::updateTabItems(session = session,
                                         inputId = 'Sidebar',
                                         selected = 'FilterCells')


        }
        ########RT end first part
      })

      ###filter
      shiny::observeEvent(input$min_correlation_matrix, {
        if (!is.null(variables$to_keep)) {
          variables$to_keep = foreach::foreach(i = unique(variables$basenames$Group)) %do% {
            sub_mat = variables$results[variables$basenames$Group == i, variables$basenames$Group == i]
            diag(sub_mat) = 0
            ! matrixStats::rowQuantiles(x = sub_mat,
                                        probs = 0.60,
                                        na.rm = T) <= input$min_correlation_matrix
          }

          variables$to_keep = unlist(variables$to_keep)
          variables$results_after_filtering = variables$results[variables$to_keep, variables$to_keep]
          variables$basenames_after_filtering = variables$basenames[variables$to_keep,]

          output$plot2_RT <- plotly::renderPlotly(
            heatmaply::heatmaply(
              x = variables$results_after_filtering,
              colors = variables$color,
              dendrogram = F,
              showticklabels = F,
              row_side_colors = variables$basenames_after_filtering,
              col_side_colors = variables$basenames_after_filtering,
              limits = c(0, 1)
            ) %>%
              plotly::layout(
                showlegend = FALSE,
                legend = FALSE,
                annotations = list(visible = FALSE)
              )
          )
        }


      })

      ########RT second part start
      shiny::observeEvent(input$RT_next, {
        if (input$RT_next == 1) {
          shinyjs::disable('RT_next')
          #write matrix and plot heatmap after filtering
          saveRDS(
            object = variables$results_after_filtering,
            file = file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              paste0(
                input$Analysis_Name,
                '_correlation_per_cell_after_filtering.rds'
              )
            )
          )

          variables$Index = variables$Index[variables$to_keep]

          #filter out samples that don't correlate and save
          data$signal_smoothed = data$signal_smoothed %>%
            dplyr::filter(index %in% variables$Index) %>%
            tidyr::separate(index, c('group', 'index'), sep = ' _ ') %>%
            dplyr::mutate(group = factor(group, level = unique(variables$File_paths$group)))

          rep_percentage = data$signal_smoothed %>%
            dplyr::group_by(Cell, basename, group, index) %>%
            dplyr::summarise(Rep_percentage = mean(Rep))

          #new index
          new_index_list = rep_percentage %>%
            dplyr::ungroup() %>%
            dplyr::arrange(Rep_percentage) %>%
            dplyr::group_by(group) %>%
            dplyr::mutate(newIndex = 1:dplyr::n()) %>%
            dplyr::arrange(group, newIndex) %>%
            dplyr::select(oldIndex = index, newIndex, Cell, basename, group)

          data$signal_smoothed = data$signal_smoothed %>%
            dplyr::ungroup() %>%
            dplyr::inner_join(new_index_list,
                              by = c('Cell', 'index' = 'oldIndex', 'basename', 'group')) %>%
            dplyr::select(
              chr,
              start,
              end,
              CN,
              background,
              CN_bg,
              th,
              Rep,
              PercentageReplication,
              Cell,
              basename,
              group,
              newIndex
            )


          readr::write_tsv(
            x = data$signal_smoothed,
            file = file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              paste0(
                input$Analysis_Name,
                '_single_cells_CNV_',
                input$binsize,
                'Mb.tsv'
              )
            ),
            col_names = T
          )

          #select used data and save the new per cell files
          used_cells =  rbind(
            new_index_list %>%
              dplyr::select(Cell, basename, group) %>%
              dplyr::ungroup(),
            data$PerCell %>%
              dplyr::filter(Type == 'G1/G2-phase cells') %>%
              dplyr::select(Cell, basename, group) %>%
              dplyr::ungroup()
          )

          data$PerCell = data$PerCell %>%
            dplyr::inner_join(used_cells,
                              Joining,
                              by = c("Cell", "basename", "group"))

          #create folder
          if (!dir.exists(
            file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              'Cells_used_in_the_analysis_info'
            )
          )) {
            dir.create(
              file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                'Cells_used_in_the_analysis_info'
              ),
              recursive = T
            )
          }

          #isolate for foreach
          PerCell = shiny::isolate(data$PerCell)

          bs = foreach::foreach(i = unique(PerCell$basename)) %do% {
            PerCell %>%
              dplyr::filter(basename == i) %>%
              dplyr::select(
                Cell,
                normalized_dimapd,
                mean_ploidy,
                ploidy_confidence,
                is_high_dimapd,
                is_noisy,
                coverage_per_1Mbp
              ) %>%
              readr::write_csv(
                file.path(
                  variables$roots['OutputFolder'],
                  input$Analysis_Name,
                  'Cells_used_in_the_analysis_info',
                  paste0(
                    input$Analysis_Name,
                    i,
                    '_per_Cell_summary_metrics.csv'
                  )
                )
              )
            i
          }

          rm('bs')
          rm('PerCell')
          rm('new_index_list')

          #calculate psudobulk
          data$scRT = Kronos.scRT::pseudoBulkRT(data$signal_smoothed)

          readr::write_tsv(
            x = data$scRT,
            file = file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              paste0(
                input$Analysis_Name,
                '_calculated_replication_timing_',
                input$binsize,
                'Mb.tsv'
              )
            ),
            col_names = T
          )

          #merege reference and scRT
          data$RTs = rbind(data$scRT %>% dplyr::ungroup(),
                           data$reference %>% dplyr::ungroup())
          ###correlationplot
          if (length(unique(data$RTs$basename)) > 1) {
            #correlation plot
            p_pairs = Kronos.scRT::KCorr_plot(data$RTs)

            output$plot3_RT <- shiny::renderPlot({
              p_pairs
            },
            height = function() {
              session$clientData[[paste0('output_plot3_RT_out_width')]]
            })

            ggplot2::ggsave(
              plot = p_pairs,
              filename = file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                paste0(
                  input$Analysis_Name,
                  '_pair_scatter_plot_RTs.pdf'
                )
              ),
              device = grDevices::cairo_pdf
            )
          }
          #prepare for twidth
          data$variability = Kronos.scRT::Variability(S_scCN = data$signal_smoothed,
                                                      scRT = data$scRT)
          readr::write_tsv(
            x = data$variability,
            file = file.path(
              variables$roots['OutputFolder'],
              input$Analysis_Name,
              paste0(input$Analysis_Name,
                     '_scRT_variability.tsv')
            ),
            col_names = T
          )

        }
      })
      ########RT second part end

    }
    ##### twidth start
    {
      shiny::observe({
        if (nrow(data$variability) > 0 & input$Sidebar == 'Twidth') {
          Twidth_module_ls$ui = lapply(unique(data$variability$group), function(x)
            Twidth_ui(x))
          output$TW_ui <- shiny::renderUI(Twidth_module_ls$ui)

          Twidth_module_ls$server = lapply(unique(data$variability$group), function(x)
            Twidth_server(
              id = x,
              variability = data$variability %>% dplyr::filter(group == x),
              out = file.path(variables$roots['OutputFolder'],
                              input$Analysis_Name, 'Twidth'),
              cores = input$cores
            ))
        }
      })

    }
    ##### twidth end

    ###binprobrep start
    {
      shiny::observe({
        if (nrow(data$G1G2) > 0 & input$Sidebar == 'BinRep') {
          data$variabilityBR = rbind(
            Kronos.scRT::Prepare_G1G2_phase_cells_forBinRepProb(G1.G2 = data$G,
                                                                RT = data$RT),
            Kronos.scRT::Prepare_S_phase_cells_forBinRepProb(S = data$S,
                                                             RT = data$RT)
          )
        } else{
          shinyjs::disable('Save__BinRep')
        }
      })

      shiny::observeEvent(input$Save__BinRep, {
        if (input$Save__BinRep > 0) {
          variables$Save__BinRep = T
          shinyjs::disable('Save__BinRep')
        }
      })

      shiny::observe({
        if (nrow(data$variabilityBR) > 0) {
          Groups = unique(data$variabilityBR$group)

          BinRep_module_ls$ui = lapply(Groups, function(x)
            BinRepProb_ui(x))
          BinRep_module_ls$server = lapply(Groups, function(x)
            BinRepProb_server(
              x,
              variabilityBR = data$variabilityBR,
              out = file.path(
                variables$roots['OutputFolder'],
                input$Analysis_Name,
                'BinsRepProb'
              ) ,
              save = variables$Save__BinRep
            ))
          output$BinRep_ui <- shiny::renderUI(BinRep_module_ls$ui)
          shinyjs::enable('Save__BinRep')
          if (variables$Save__BinRep) {
            variables$Save__BinRep = F
          }
        }



      })

    }
    ### binprobrep end


    #reference
    {
      shiny::observe({
        shinyFiles::shinyFileChoose(
          input = input,
          id = 'reference_file',
          roots = variables$roots,
          session = session,
          defaultRoot = 'Home'
        )

        if (file.exists(variables$reference_file) &
            input$reference_group != 'Select a group') {
          shinyjs::enable('Add_reference')
        } else{
          shinyjs::disable('Add_reference')

        }
      })

      #select and check reference file
      shiny::observeEvent(input$reference_file, {
        if (!is.numeric(input$reference_file)) {
          #get path
          variables$reference_file = dplyr::tibble(
            shinyFiles::parseFilePaths(
              roots = variables$roots ,
              selection = input$reference_file
            )
          ) %>% dplyr::pull(datapath)

          #check if it is the right file
          variables$reference_file = Kronos.scRT::right_format(
            file_path = variables$reference_file,
            delim = '\t',
            columns_to_check = 4,
            wrong_message = paste0(
              '<b>',
              variables$reference_file,
              '<font color=\"#FF0000\"> does not have the right format. </b></font>'
            ),
            rigth_message = variables$reference_file
          )

          output$reference_file_out = shiny::renderText(variables$reference_file)

        } else{
          variables$reference_file = ''
          output$reference_file_out = shiny::renderText(variables$reference_file)

        }
      })

      #replace ' ' with '_' for reference name
      shiny::observeEvent(input$reference_name, {
        updateTextInput(
          session = session,
          inputId = 'reference_name',
          value = stringr::str_replace_all(
            string = input$reference_name,
            pattern = ' ',
            replacement = '_'
          )
        )
      })

      #add reference
      shiny::observeEvent(input$Add_reference, {
        shinyjs::disable('Add_reference')
        data$reference = rbind(
          data$reference,
          readr::read_tsv(
            variables$reference_file,
            col_names = c('chr', 'start', 'end', 'RT'),
            col_types = 'cnnn'
          ) %>%
            dplyr::mutate(
              basename = input$reference_name,
              group = input$reference_group
            )
        )

        variables$Reference_paths = rbind(
          variables$Reference_paths,
          dplyr::tibble(
            File = basename(variables$reference_file),
            basename = input$reference_name,
            group = input$reference_group
          )
        )

        output$Reference_paths = shiny::renderTable(variables$Reference_paths)


        variables$reference_file = ''
        output$reference_file_out = shiny::renderText('')
        variables$groups_with_a_reference = c(variables$groups_with_a_reference,
                                              input$reference_group)
        variables$groups = variables$groups[!variables$groups %in% variables$groups_with_a_reference]

        #update groups for Reference
        shiny::updateSelectInput(
          inputId = 'reference_group',
          choices = c('Select a group',
                      variables$groups),
          selected = 'Select a group'
        )
        updateTextInput(inputId = 'reference_name',
                        value = paste0('Reference_', input$Add_reference))
      })

    }
  }
  #scPlots

  shiny::observe({
    if (input$Sidebar == 'scPlots' & ncol(data$scRT) != 0) {
      shiny::updateSelectInput(inputId = 'Chr__scPlot',
                               choices = unique(data$scRT$chr))

      data$summary = data$signal_smoothed %>%
        dplyr::ungroup() %>%
        dplyr::summarise(CN_bg = round(stats::quantile(CN_bg, c(0.01, 0.99)), 1),
                         CN = round(stats::quantile(CN, c(0.01, 0.99)), 1))


      # change chrom
      shiny::observeEvent(input$Chr__scPlot, {
        if (input$Chr__scPlot != 'Chrom') {
          # update min max range
          Max = data$scRT %>% dplyr::filter(chr == input$Chr__scPlot) %>% dplyr::pull(end) %>%
            max() / 10 ^ 6
          Step = abs(data$scRT[1, 'start'] - data$scRT[1, 'end']) / 10 ^ 6

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
        shinyjs::disable('Save__scPlots')
        variables$Save__scPlots = T
      })


      shiny::observeEvent(c(
        input$range__scPlot,
        input$what__scPlot,
        input$Save__scPlots
      ),
      {
        for (g in unique(data$RTs$group)) {
          scPlot_module_ls$ui[[g]] = scPlots_ui(g)
          scPlot_module_ls$server[[g]] = scPlots_server(
            g,
            RTs = data$RTs %>%
              dplyr::filter(
                group == g,
                chr == input$Chr__scPlot,
                start > input$range__scPlot[1] *
                  10 ^ 6,
                end < input$range__scPlot[2] *
                  10 ^ 6
              ),
            scCN = data$signal_smoothed %>%
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
            save = variables$Save__scPlots
          )
          if (variables$Save__scPlots) {
            variables$Save__scPlots = F
            shinyjs::enable('Save__scPlots')

          }


        }
        output$scPlots_UI = shiny::renderUI({
          scPlot_module_ls$ui
        })
      })


    } else{
      shinyjs::disable('Save__scPlots')

    }
  })


  # Dimensionality reduction
  shiny::observe({
    if (input$Sidebar == 'DRed' & ncol(data$scRT) != 0) {
      #call Dim_red_server
      Dred = Dim_red_server(
        id = 'Dred',
        scCN = data$signal_smoothed,
        out =  file.path(variables$roots['OutputFolder'],
                         input$Analysis_Name),
        cores = input$cores
      )
    }
  })

  ####exit
  {
    shiny::observe({
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
