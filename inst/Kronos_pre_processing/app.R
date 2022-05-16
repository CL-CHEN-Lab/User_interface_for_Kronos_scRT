#load operator
`%>%` = tidyr::`%>%`
#set ggplot theme
ggplot2::theme_set(new = ggplot2::theme_bw())
#java function to close this app
jscode <- "shinyjs.closeWindow = function() { window.close(); }"
#find max numebr of cores
maxCores = parallel::detectCores()

# Define UI for application that draws a histogram
ui <- shinydashboard::dashboardPage(
  title = 'Kronos scRT',
  skin = 'red',
  shinydashboard::dashboardHeader(title = shiny::span(
    shiny::img(src = 'KronosLogo.png', width = '100%')
  )),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      id = 'Sidebar',
      shinydashboard::menuItem(text = "FastqToBam",
                               tabName = "Home"),
      shinydashboard::menuItem(text = "FastqToBam Metrics",
                               tabName = "Metrics"),
      shinydashboard::menuItem(text = "Binning",
                               tabName = "Binning"),
      shinydashboard::menuItem(text = "Copy Number Calling",
                               tabName = "CN"),
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
      color = 'red',
      height = '200px',
      width = '200px'
    ),

    shinydashboard::tabItems(#####home
      {
        shinydashboard::tabItem(tabName = "Home",
                                shiny::fluidPage(
                                  shiny::fluidRow(
                                    shiny::fluidRow(
                                      shiny::column(
                                        width = 3,

                                        shinyFiles::shinyDirButton(
                                          id = 'Fastq',
                                          label = 'Fastq files folder',
                                          title = 'Fastq files folder',
                                          style = 'width:100%;'
                                        )
                                      ),
                                      shiny::column(width = 9, shiny::htmlOutput('Fastq_out'))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(
                                        width = 3,
                                        shinyFiles::shinyFilesButton(
                                          id = 'index',
                                          label = shiny::div(
                                            'Bowtie2 index',
                                            bsplus::shiny_iconlink() %>%
                                              bsplus::bs_embed_popover(title = "Select one of the index files. If you don't have a compiled index, you can download one from http://bowtie-bio.sourceforge.net/bowtie2."  , placement = 'right')
                                          ),
                                          title = 'Bowtie2 index',
                                          style = 'width:100%;',
                                          multiple = F
                                        )
                                      ),
                                      shiny::column(width = 9, shiny::htmlOutput('index_out'))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(
                                        width = 3,
                                        shinyFiles::shinyFilesButton(
                                          id = 'Adapters',
                                          label = shiny::div(
                                            'Adapters list',
                                            bsplus::shiny_iconlink() %>%
                                              bsplus::bs_embed_popover(title = 'A text files with one adapter per line with no header. If not provided Kronos will look for standard illumina adapters' , placement = 'right')
                                          ),
                                          title = 'Adapters list',
                                          style = 'width:100%;',
                                          multiple = F
                                        )
                                      ),
                                      shiny::column(width = 9, shiny::htmlOutput('Adapters_out'))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(
                                        width = 3,
                                        shinyFiles::shinyDirButton(
                                          id = 'out',
                                          label = 'Output files folder',
                                          title = 'Output files folder',
                                          style = 'width:100%;'
                                        )
                                      ),
                                      shiny::column(width = 9, shiny::htmlOutput('Output_out'))
                                    ),
                                    shiny::fluidRow(
                                      shiny::column(
                                        width = 3,
                                        shiny::radioButtons(
                                          inputId = 'reads_type',
                                          label = 'Sequencing',
                                          choices = c('PE', 'SE', 'PE to treat as SE'),
                                          inline = T,
                                          selected = 'PE',
                                          width = '100%'
                                        )
                                      ),
                                      shiny::column(
                                        width = 3,
                                        shiny::numericInput(
                                          inputId = 'min_size',
                                          label = 'Read min size after trimming',
                                          value = 25,
                                          min = 15,
                                          step = 1,
                                          width = '100%'
                                        )
                                      ),
                                      shiny::column(
                                        width = 3,
                                        shiny::radioButtons(
                                          inputId = 'phred33',
                                          label = 'Quality score',
                                          choiceValues = c(T, F) ,
                                          choiceNames = c('phred33', 'phred64'),
                                          selected = T,
                                          inline = T,
                                          width = '100%'
                                        )
                                      )
                                    ),
                                    shiny::uiOutput('patter_sample'),
                                    shiny::fluidRow(shiny::column(
                                      width = 3,
                                      shiny::sliderInput(
                                        inputId = 'cores',
                                        label = 'Cores to use',
                                        min = 1,
                                        max = maxCores,
                                        value = trunc(maxCores / 2),
                                        step = 1
                                      )
                                    )),
                                    shiny::fluidRow(shiny::column(
                                      width = 3,
                                      shinyjs::disabled(
                                        shiny::actionButton(
                                          inputId = 'run',
                                          label = 'Run',
                                          width = '100%'
                                        )
                                      )
                                    )),
                                    shiny::fluidRow(shiny::column(width = 12,
                                                                  shiny::htmlOutput('final')))
                                  )
                                ))
      },
      ######metrics
      {
        shinydashboard::tabItem(
          tabName = "Metrics",
          shiny::fluidRow(
            shiny::column(
              width = 4,
              align = 'center',
              shinyFiles::shinyFilesButton(
                id = 'Load_metrics',
                label = 'Load old Metrics',
                title = 'Load old Metrics',
                multiple = F,
                style = 'width:75%;'
              )
            ),
            shiny::column(
              width = 4,
              align = 'center',
              shinyjs::disabled(
                shinyFiles::shinySaveButton(
                  id = 'Save_plots',
                  label = 'Save diagnostic Plot',
                  title = 'Save diagnostic Plot',
                  filetype = '.pdf',
                  style = 'width:75%;'
                )
              )
            ),
            shiny::column(
              width = 4,
              align = 'center',
              shinyjs::disabled(
                shinyFiles::shinySaveButton(
                  id = 'Save_table',
                  label = 'Save table',
                  title = 'Save table',
                  filetype = '.tsv',
                  style = 'width:75%;'
                )
              )
            )
          ),
          shiny::fluidRow(shiny::column(
            width = 12,
            shiny::plotOutput(
              'final_plot1',
              # click = shiny::clickOpts(id = 'final_plot_click1'),
              brush = shiny::brushOpts(id = 'final_plot_brush1',
                                       resetOnNew = T),
              dblclick = shiny::dblclickOpts(id = 'final_plot_dclick1', clip = T, delay = 200)
            )
          )),

          shiny::fluidRow(shiny::column(
            width = 12,
            shinyjs::hidden(
              shiny::radioButtons(
                inputId = 'final_table_selector',
                label = 'Parameter',
                choices = '',
                inline = F,
                width = '100%'
              )
            )
          )),
          shiny::fluidRow(shiny::column(
            width = 12,
            shiny::tableOutput('final_table')
          ))
        )
      },
      #####binning
      {
        shinydashboard::tabItem(
          tabName = "Binning",
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyFilesButton(
                id = 'RefGenome_binning',
                label = shiny::div(
                  'Reference genome fasta file',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = "For most genomes, fasta file can be downloaded from UCSC Genome Browser. N.B. you need to download the file containing all chromosomes." , placement = 'right')
                ),
                title = 'Reference genome fasta file',
                multiple = F,
                style = 'width:100%;'
              )
            ),
            shiny::column(width = 9,
                          shiny::htmlOutput('RefGenome_binning_out'))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyFilesButton(
                id = 'index_binning',
                label = shiny::div(
                  'Bowtie2 index',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = "Select one of the index files. If you don't have a compiled index, you can download one from http://bowtie-bio.sourceforge.net/bowtie2." , placement = 'right')
                ),
                title = 'Bowtie2 index',
                style = 'width:100%;',
                multiple = F
              )
            ),
            shiny::column(width = 9, shiny::htmlOutput('index_binning_out'))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyDirButton(
                id = 'out_binning',
                label = 'Output files folder',
                title = 'Output files folder',
                style = 'width:100%;'
              )
            ),
            shiny::column(width = 9, shiny::htmlOutput('Output_binning_out'))
          ),
          shiny::fluidRow(
            shiny::column(
            width = 3,
            shiny::radioButtons(
              inputId = 'Param_binning',
              label = 'Parameters estimation',
              choices = c('Auto', 'Manual'),
              selected = 'Auto',
              inline = T,
              width = '100%'
            )
          )),
          shiny::uiOutput('Paramui_binning_out1'),
          shiny::uiOutput('Paramui_binning_out2'),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::sliderInput(
                inputId = 'Coverage_binning',
                label = shiny::div(
                  'Coverage',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = 'Simulated dept of sequencing to calculate bin mappability.' , placement = 'right')
                ),
                min = 1,
                max = 10,
                value = 1,
                step = 1,
                post = 'X',
                width = '100%'
              )
            ),
            shiny::column(
              width = 3,
              shiny::sliderInput(
                inputId = 'BinSize_binning',
                label = shiny::div(
                  'BinSize',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = 'Bin size that will be used to call CNV. This size does not correspond to the final scRT resolution. The bin size should increased if the average cell coverage is very low.' , placement = 'right')
                ),
                min = 5,
                max = 100,
                value = 20,
                step = 1,
                post = 'kb',
                width = '100%'
              )
            ),
            shiny::column(
              width = 3,
              shiny::sliderInput(
                inputId = 'Cores_binning',
                label = 'Cores',
                min = 1,
                max = maxCores,
                value = trunc(maxCores / 2),
                step = 1,
                width = '100%'
              )
            )
          ),
          shinydashboard::box(
            title = 'Advance Options',
            status = 'primary',
            solidHeader = T,
            width = '100%',
            collapsible = T,
            collapsed = T,
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shiny::sliderInput(
                  inputId = 'Mappability_binning',
                  label = shiny::div(
                    'Acceptable mappability range',
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'Bins with mappability higher or lower than the acceptable mappability range will be excluded from the analysis.' , placement = 'right')
                  ),
                  min = 0.5,
                  max = 3,
                  value = c(0.8, 1.5),
                  step = 0.1,
                  width = '100%'
                )
              ),
              shiny::column(
                width = 3,
                shiny::sliderInput(
                  inputId = 'ErrorRate_binning',
                  label = shiny::div(
                    'Sequencer Error Rate',
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'Is the error rate of the sequencer you used to sequence your data. This parameter is used to simulates reads and calculate the bin mappability.' , placement = 'right')
                  ),
                  min = 0,
                  max = 0.3,
                  value = 0.1,
                  post = '%',
                  step = 0.001,
                  width = '100%'
                )
              )
            ),
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shinyFiles::shinyFilesButton(
                  id = 'Blacklist_binning',
                  label = shiny::div(
                    'Blacklist',
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'Because of their repetitive nature, some regions of the genome are difficult to map. The user can hide these regions form the analysis providing a blacklist file. Such lists are available in literature. Some useful ones can be found at https://github.com/Boyle-Lab/Blacklist' , placement = 'right')
                  ),
                  title = 'Blacklist',
                  multiple = F
                )
              ),
              shiny::column(width = 9,
                            shiny::htmlOutput('Blacklist_binning_out'))

            )
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::textInput(
                inputId = 'chr_prefix_binning',
                label = shiny::div(
                  'Chromosome prefix',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = 'Depending on your reference genome, chromosomes or scaffolds could be named differntly. Please, privide the prefix that preceeds the chromosome or scaffold number. If no prefix is present this box has to be emptied.' , placement = 'right')),
                  value = 'chr',
                  width = '100%'
                )
              ),
              shiny::column(
                width = 3,
                shiny::textInput(
                  inputId = 'chr_range_binning',
                  label = shiny::div(
                    'Chromosomes to consider',
                    bsplus::shiny_iconlink() %>%
                      bsplus::bs_embed_popover(title = 'Chromosome to include in the analysis. You can define ranges using ":" and individual chromosomes using ",". Don\'t include the chr prefix!' , placement = 'right')
                  ),
                  value = '1:22',
                  width = '100%'
                )
              )
            ),
            shiny::fluidRow(shiny::column(
              width = 3,
              shinyjs::disabled(
                shiny::actionButton(
                  inputId = 'run_binning',
                  label = 'Run',
                  width = '100%'
                )
              )
            ))
          )
      },
      #####CN
      {
        shinydashboard::tabItem(
          tabName = "CN",
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyFilesButton(
                id = 'Bins_CN',
                label = div('Bins file',
                            bsplus::shiny_iconlink() %>%
                              bsplus::bs_embed_popover(title = 'Bins file created during the binning step.' , placement = 'right')),
                title = 'Bins file',
                style = 'width:100%;',
                multiple = F
              )
            ),
            shiny::column(width = 9,
                          shiny::htmlOutput('Bins_CN'))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyFilesButton(
                id = 'ChrSize_CN',
                label = shiny::div(
                  'Chrom size file',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = 'For most genome assemblies this file can be downloaded from UCSC Genome Browser website.' , placement = 'right')),
                title = 'Chrom size file',
                style = 'width:100%;',
                multiple = F
              )
            ),
            shiny::column(width = 9,
                          shiny::htmlOutput('ChrSize_CN'))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyDirButton(
                id = 'Output_CN',
                label = 'Results folder',
                title = 'Results folder',
                style = 'width:100%;'
              )
            ),
            shiny::column(width = 9,
                          shiny::htmlOutput('Output_CN'))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::textInput(
                inputId = 'CHR_CN',
                label = shiny::div('Chromosome prefix',
                                   bsplus::shiny_iconlink() %>%
                                     bsplus::bs_embed_popover(title = 'Depending on your reference genome, chromosomes or scaffolds could be named differntly. Please, privide the prefix that preceeds the chromosome or scaffold number. If no prefix is present this box has to be emptied.' , placement = 'right')
                ),
                value = 'chr',
                width = '100%'
              )
            ),
            shiny::column(
              width = 3,
              shiny::textInput(
                inputId = 'CHR_RANGE_CN',
                label = shiny::div(
                  'Chromosomes to consider',
                  bsplus::shiny_iconlink() %>%
                    bsplus::bs_embed_popover(title = 'Chromosome to include in the analysis. You can define ranges using ":" and individual chromosomes using ",". Don\'t include the chr prefix!' , placement = 'right')
                ),
                value = '1:22',
                width = '100%'
              )
            ),
            shiny::column(
              width = 3,
              shiny::sliderInput(
                inputId = 'Cores_CN',
                label = 'Number of cores',
                min = 1,
                max = maxCores,
                value = trunc(maxCores / 2),
                step = 1
              )
            )
          ),
          shinydashboard::box(
            title = 'Advance Options',
            status = 'primary',
            solidHeader = T,
            width = 12,
            collapsible = T,
            collapsed = T,
            shiny::fluidRow(shiny::column(
              width = 3,
              shiny::numericInput(
                inputId = 'min_reads_CN',
                label = shiny::div('Min number of reads per cell',
                                   bsplus::shiny_iconlink() %>%
                                     bsplus::bs_embed_popover(title = 'This threshold depends on the genome size of your cells and on the chromosomes you want to keep in the analysis.' , placement = 'right')
                ),
                value = 2 * 10 ^ 5,
                step = 10 ^ 4,
                width = '100%'
              )
            )),
            shiny::fluidRow(shiny::column(
              width = 12,
              shiny::radioButtons(
                inputId = 'Ploidy_Restrictions_CN',
                label = div('Ploidy Restrictions',
                            bsplus::shiny_iconlink() %>%
                              bsplus::bs_embed_popover(title = 'If Fix Ploidy Range is selected, the software selects the most appropriated ploidy in the ragne and calculates the cofidence associated. If Fix Average Ploidy is selected, cells are imposed the closest possible ploidy to the average one and does not calculate the confidence of this calling.' , placement = 'right')
                ),
                choices = c('Fix Average Ploidy', 'Fix Ploidy Range'),
                inline = T,
                selected = 'Fix Ploidy Range',
                width = '100%'
              )
            )),
            shiny::uiOutput('FixAveragePloidy'),
            shiny::uiOutput('FixPloidyRange')
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shinyFiles::shinyDirButton(
                id = 'Bamdir_CN',
                label = 'Bam files folder',
                title = 'Bam files folder',
                style = 'width:100%;'
              )
            ),
            column(width = 9,
                   shiny::htmlOutput('Bamdir_CN'))
          ),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::textInput(
                inputId = 'basename_CN',
                label = shiny::div('Sample basename',
                                   bsplus::shiny_iconlink() %>%
                                     bsplus::bs_embed_popover(title = 'This Name identifies each individual experiment. It can be modified later on in the analysis.' , placement = 'right')),
                value = 'Exp',
                width = '100%'
              )
            ),
            shiny::column(
              width = 3,
              shiny::textInput(
                inputId = 'group_CN',
                label = shiny::div('Sample group',
                                   bsplus::shiny_iconlink() %>%
                                     bsplus::bs_embed_popover(title = 'If cells have been sequenced in different experiments, providind the same group name allows to merge them after normalisation. It can be modified later on in the analysis.' , placement = 'right')),
                value = 'Exp',
                width = '100%'
              )
            )
          ),
          shiny::fluidRow(shiny::column(
            width = 3,
            shinyjs::disabled(
              shiny::actionButton(
                inputId = 'AddSample_CN',
                label = 'Add Sample',
                width = '100%'
              )
            )
          )),
          shiny::fluidRow(shiny::column(
            width = 9,
            shiny::tableOutput('Samples_CN')
          )),
          shiny::fluidRow(shiny::column(
            width = 3,
            shinyjs::disabled(
              shiny::actionButton(
                inputId = 'Run_CN',
                label = 'Run',
                width = '100%'
              )
            )
          ))
        )

      })
    )
  )

  server <- function(input, output, session) {
    #define operators
    `%>%` = tidyr::`%>%`
    `%do%` = foreach::`%do%`
    `%dopar%` = foreach::`%dopar%`

    #find max number of cores
    maxCores = parallel::detectCores()

    #set roots
    roots = c(shinyFiles::getVolumes()(),
              Home = Sys.getenv("HOME"))

    #set shinyFiles buttos
    shinyFiles::shinyDirChoose(
      input = input,
      id = 'Fastq',
      session = session,
      roots = roots,
      defaultRoot = 'Home'
    )
    shinyFiles::shinyDirChoose(
      input = input,
      id = 'out',
      session = session,
      roots = roots,
      defaultRoot = 'Home'
    )
    shinyFiles::shinyFileChoose(
      input = input,
      id = 'index',
      session = session,
      roots = roots,
      defaultRoot = 'Home'
    )
    shinyFiles::shinyFileChoose(
      input = input,
      id = 'Adapters',
      session = session,
      roots = roots,
      defaultRoot = 'Home'
    )
    shinyFiles::shinyFileChoose(
      input = input,
      id = 'fasta',
      session = session,
      roots = roots,
      defaultRoot = 'Home'
    )
    #define variables
    variables = reactiveValues(
      File1 = NULL,
      File2 = NULL,
      Adapters = NULL,
      out = NULL,
      fastq_dir = NULL,
      index = NULL,
      activate_start1 = F,
      activate_start2 = F,
      activate_start3 = F,
      proceed_run = F,
      info = dplyr::tibble()
    )

    variables_binning = reactiveValues(
      fasta = NULL,
      index = NULL,
      out = NULL,
      Bam_dir = NULL,
      blacklist = NULL,
      activate_start1 = F,
      activate_start2 = F,
      activate_start3 = F,
      activate_start4 = F,
      activate_start5 = T,
      read_size = 40,
      frag_size = 200,
      read_type = 'PE'
    )

    ##fastqToBam
    {
      #stop app when the session ends
      session$onSessionEnded(function() {
        shiny::stopApp()
      })

      #loadr tex input if PE format
      shiny::observeEvent(input$reads_type, {
        if (input$reads_type != 'SE') {
          output$patter_sample = shiny::renderUI({
            shiny::fluidRow(column(
              width = 3,
              shiny::textInput(
                inputId = 'R1',
                label = 'Pattern to identiry R1 fastq',
                value = '_R1'
              )
            ),
            column(
              width = 3,
              shiny::textInput(
                inputId = 'R2',
                label = 'Pattern to identiry R2 fastq',
                value = '_R2'
              )
            ))
          })
        } else{
          output$patter_sample = shiny::renderUI({
            NULL
          })
        }
      })


      ###select fastq folder
      shiny::observeEvent(input$Fastq, {
        output$final <- shiny::renderText(NULL)

        if (!is.numeric(input$Fastq)) {
          #recover fastq directory
          variables$fastq_dir = shinyFiles::parseDirPath(selection = input$Fastq, roots = roots)

          #check for fastq files
          if (length(list.files(path =  variables$fastq_dir, pattern = 'fastq|fq')) >
              0) {
            output$Fastq_out <- shiny::renderText(variables$fastq_dir)
            variables$activate_start1 = T
          } else{
            output$Fastq_out <-
              shiny::renderText(
                paste(
                  '<b><p style="color:#FF0000";> selected directory does not contain fastq files</b></p>'
                )
              )
            variables$activate_start1 = F
          }


        } else{
          variables$activate_start1 = F
          output$Fastq_out <- shiny::renderText(NULL)
        }
      })

      ###select index
      shiny::observeEvent(input$index, {
        output$final <- shiny::renderText(NULL)

        if (!is.numeric(input$index)) {
          #recover index directory
          variables$index = shinyFiles::parseFilePaths(selection = input$index, roots = roots)
          variables$index = variables$index$datapath

          #check that the selected file is actually part of an index
          if (stringr::str_detect(string = variables$index, pattern = '.bt2')) {
            variables$index = stringr::str_remove(variables$index, pattern = '.[0-9]{1,2}.bt2|.rev.[0-9]{1,2}.bt2')

            output$index_out <- shiny::renderText(variables$index)
            variables$activate_start3 = T
          } else{
            output$index_out <-
              shiny::renderText(
                paste(
                  '<b><p style="color:#FF0000";> selected file does not belong to a bowtie2 index.</b></p>'
                )
              )
            variables$activate_start3 = F
          }
        } else{
          variables$activate_start3 = F
          output$index_out <-
            shiny::renderText(NULL)
        }
      })

      ###select adapter
      shiny::observeEvent(input$Adapters,
                          {
                            if (!is.numeric(input$Adapters))
                            {
                              #recover index directory
                              variables$Adapters = shinyFiles::parseFilePaths(selection = input$Adapters, roots = roots)
                              variables$Adapters = variables$Adapters$datapath

                              #check that the selected file is actually part of an index
                              if (Kronos.scRT::right_format(
                                file_path = variables$Adapters,
                                columns_to_check = 1,
                                logical = T
                              )) {
                                variables$Adapters = readr::read_lines(variables$Adapters)

                                output$Adapters_out <-
                                  shiny::renderText(variables$Adapters)
                              } else
                              {
                                output$Adapters_out <-
                                  shiny::renderText(
                                    paste(
                                      '<b><p style="color:#FF0000";> Adapters file has to be a text file with one adapter per line, No header!</b></p>'
                                    )
                                  )
                              }
                            }
                          })

      ###select out folder
      shiny::observeEvent(input$out,
                          {
                            output$final <- shiny::renderText(NULL)
                            if (!is.numeric(input$out))
                            {
                              #recover fastq directory
                              variables$out = shinyFiles::parseDirPath(selection = input$out, roots = roots)

                              #check for fastq files
                              if (length(list.files(path =  variables$out)) == 0)
                              {
                                output$Output_out <- shiny::renderText(variables$out)
                                variables$activate_start2 = T
                              } else{
                                output$Output_out <-
                                  shiny::renderText(
                                    paste(
                                      '<b><p style="color:#FF0000";> selected directory is not empty</b></p>'
                                    )
                                  )
                                variables$activate_start2 = F
                              }
                            } else{
                              variables$activate_start2 = F
                              output$Output_out <-
                                shiny::renderText(NULL)
                            }
                          })

      ##activate start
      shiny::observeEvent(
        c(
          variables$activate_start1,
          variables$activate_start2,
          variables$activate_start3
        ),
        {
          if (variables$activate_start1 &
              variables$activate_start2 &
              variables$activate_start3)
          {
            shinyjs::enable('run')
          } else
          {
            shinyjs::disable('run')

          }

        }
      )

      ###start
      shiny::observeEvent(input$run,
                          {
                            if (input$run > 0)
                            {
                              #disable all buttons till the end of the run
                              shinyjs::disable('run')
                              shinyjs::disable('Fastq')
                              shinyjs::disable('out')
                              shinyjs::disable('reads_type')
                              shinyjs::disable('cores')
                              shinyjs::disable('Adapters')
                              shinyjs::disable('index')

                              all_files = list.files(
                                path = variables$fastq_dir,
                                pattern = 'fastq|fq',
                                ignore.case = T,
                                full.names = T
                              )
                              if (input$reads_type == 'PE')
                              {
                                R1 = stringr::str_detect(string = basename(all_files),
                                                         pattern = input$R1)
                                R2 = stringr::str_detect(string = basename(all_files),
                                                         pattern = input$R2)

                                variables$File1 = all_files[R1]
                                variables$File2 = all_files[R2]
                                #check if file names order is matching
                                if (all(
                                  stringr::str_split(
                                    string = variables$File1,
                                    pattern = input$R1,
                                    simplify = T
                                  )
                                  [, 1] == stringr::str_split(
                                    string = variables$File2,
                                    pattern = input$R2,
                                    simplify = T
                                  )[, 1]
                                )) {
                                  variables$proceed_run = T
                                } else
                                {
                                  variables$proceed_run = F

                                }

                              } else
                                if (input$reads_type == 'SE')
                                {
                                  variables$File1 = all_files
                                  variables$proceed_run = T

                                } else
                                {
                                  R1 = stringr::str_detect(string = all_files, pattern = input$R1)
                                  R2 = stringr::str_detect(string = all_files, pattern = input$R2)

                                  variables$File1 = all_files[R1]
                                  variables$File2 = all_files[R2]

                                  #check if file names order is matching
                                  if (all(
                                    stringr::str_split(
                                      string = variables$File1,
                                      pattern = input$R1,
                                      simplify = T
                                    )
                                    [, 1] == stringr::str_split(
                                      string = variables$File2,
                                      pattern = input$R2,
                                      simplify = T
                                    )[, 1]
                                  )) {
                                    variables$proceed_run = T

                                    cl = snow::makeCluster(input$cores)
                                    doSNOW::registerDoSNOW(cl)
                                    shiny::isolate({
                                      file1 = variables$File1
                                      file2 = variables$File2
                                      out = variables$out
                                    })
                                    variables$File1 = foreach::foreach(i = 1:length(variables$File1)) %dopar%
                                    {
                                      Kronos.scRT::mergefastq(
                                        FileList = c(file1[i], file2[i]),
                                        output = file.path(out, 'MergedFastq')
                                      )
                                    }
                                    snow::stopCluster(cl)
                                    variables$File2 = NULL
                                    variables$File1 = unlist(variables$File1)
                                    variables$proceed_run = T
                                  } else{
                                    variables$proceed_run = F
                                  }
                                }

                              if (variables$proceed_run) {
                                if (input$reads_type == 'PE') {
                                  variables$info = Kronos.scRT::FastqToBam(
                                    bowtie2_index = variables$index,
                                    File1 = variables$File1,
                                    File2 = variables$File2,
                                    outputdir = variables$out,
                                    adapters_list = variables$Adapters,
                                    cores = input$cores,
                                    Read_min_size_after_trimming = input$min_size,
                                    phred33 = input$phred33
                                  )
                                } else{
                                  variables$info = Kronos.scRT::FastqToBam(
                                    bowtie2_index = variables$index,
                                    File1 = variables$File1,
                                    outputdir = variables$out,
                                    adapters_list = variables$Adapters,
                                    cores = input$cores,
                                    Read_min_size_after_trimming = input$min_size,
                                    phred33 = input$phred33
                                  )
                                }

                                #save info
                                variables$info %>% readr::write_tsv(file = file.path(variables$out, 'FastqToBamMetrics.tsv'))

                                shinydashboard::updateTabItems(session = session,
                                                               inputId = 'Sidebar',
                                                               selected = 'Metrics')

                              } else {
                                output$final_plot1 <-
                                  shiny::renderPlot({
                                    NULL
                                  })
                                output$final_table <-
                                  shiny::renderTable({
                                    NULL
                                  })
                                output$final <-
                                  shiny::renderText(
                                    '<b><p style="color:#FF0000";> files do not match, check your fastq folder.</b></p>'
                                  )
                              }

                              #enable all buttons but run at the end of the run
                              shinyjs::enable('Fastq')
                              shinyjs::enable('out')
                              shinyjs::enable('reads_type')
                              shinyjs::enable('cores')
                              shinyjs::enable('Adapters')
                              shinyjs::enable('index')
                            }

                          })

      ####exit
      {
        shiny::observeEvent(input$Sidebar,
                            {
                              if (input$Sidebar == 'Exit') {
                                shinyjs::js$closeWindow()
                                shiny::stopApp()
                              }
                            })
      }

    }

    ####metrics
    {
      #shinyFiles buttons
      shinyFiles::shinyFileSave(
        input = input,
        id = 'Save_plots',
        session = session,
        defaultRoot = 'Home',
        roots = roots
      )
      shinyFiles::shinyFileSave(
        input = input,
        id = 'Save_table',
        session = session,
        defaultRoot = 'Home',
        roots = roots
      )
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'Load_metrics',
        session = session,
        defaultRoot = 'Home',
        roots = roots
      )


      #upload data
      shiny::observeEvent(input$Load_metrics, {
        if (!is.numeric(input$Load_metrics)) {
          #recover path and upload data
          path = shinyFiles::parseFilePaths(roots = roots, selection = input$Load_metrics)

          #check type
          if (Kronos.scRT::right_format(
            file_path = path$datapath,
            columns_to_check = c('File', 'Parameter', 'Values', 'Category'),
            logical = T
          )) {
            variables$Metrics = readr::read_tsv(path$datapath) %>%
              dplyr::mutate(
                Parameter = factor(Parameter, levels = unique(Parameter)),
                Category = factor(Category, levels = unique(Category)),
                x = ifelse(
                  Category %in% c('Total Reads', 'Unmapped Fraction', 'Proper Pair'),
                  1,
                  as.numeric(Category)
                ) + sample(
                  seq(-0.2, 0.2, 0.001),
                  size = dplyr::n(),
                  replace = T
                )
              ) %>%
              dplyr::ungroup()
            #initialise Plot_selection
            variables$Metrics$Color = 'red'
            #active buttons
            shinyjs::enable('Save_plots')
            shinyjs::enable('Save_table')

            #update selec parameters
            shiny::updateRadioButtons(
              session = session,
              inputId = 'final_table_selector',
              choices = unique(variables$Metrics$Parameter),
              inline = T
            )
            shinyjs::show('final_table_selector')

          } else{
            output$final_plot1 <- shiny::renderPlot(NULL)
            output$final_table <- shiny::renderText('Wrong File!')
            shinyjs::disable('Save_plots')
            shinyjs::disable('Save_table')
            shinyjs::hide('final_table_selector')

          }
        }
      })


      shiny::observeEvent(variables$info, {
        if (nrow(variables$info) > 0) {
          #add x variability for plotting
          variables$Metrics = variables$info %>%
            dplyr::mutate(
              Category = factor(Category, levels = unique(Category)),
              x = ifelse(
                Category %in% c('Total Reads', 'Unmapped Fraction', 'Proper Pair'),
                1,
                as.numeric(Category)
              ) + sample(
                seq(-0.2, 0.2, 0.001),
                size = dplyr::n(),
                replace = T
              )
            ) %>%
            dplyr::ungroup()
          #initialise Plot_selection
          variables$Metrics$Color = 'red'
          #active buttons
          shinyjs::enable('Save_plots')
          shinyjs::enable('Save_table')
          shiny::updateRadioButtons(
            session = session,
            inputId = 'final_table_selector',
            choices = unique(variables$info$Parameter),
            inline = T
          )
          shinyjs::show('final_table_selector')
        }
      })

      #brush plot 1
      shiny::observeEvent(input$final_plot_brush1, {
        if (!is.null(input$final_plot_brush1)) {
          #selected files
          Files = variables$Metrics %>%
            dplyr::filter(
              Parameter == input$final_plot_brush1$panelvar1,
              Values >= input$final_plot_brush1$ymin,
              Values <= input$final_plot_brush1$ymax,
              x >= input$final_plot_brush1$xmin,
              x <= input$final_plot_brush1$xmax
            ) %>% dplyr::pull(File)

          variables$Metrics = variables$Metrics %>%
            dplyr::mutate(Color = ifelse(File %in%  Files, 'red', 'black'))
        }
      })


      #reset with dbclick
      shiny::observeEvent(input$final_plot_dclick1,
                          {
                            #initialise reset
                            variables$Metrics = variables$Metrics %>%
                              dplyr::mutate(Color = 'red')
                          })

      #if something has been selected
      shiny::observeEvent(variables$Metrics, {
        variables$plot1 = variables$Metrics %>%
          dplyr::arrange(Color) %>%
          ggplot2::ggplot() +
          ggplot2::geom_boxplot(
            ggplot2::aes(Category, Values),
            outlier.size =  0,
            outlier.colour = NA,
            width = 0.8
          ) +
          ggplot2::geom_point(ggplot2::aes(x, Values, color = Color), alpha =
                                0.5) +
          ggplot2::theme_bw() +
          ggplot2::labs(x = '', y = '')  +
          ggplot2::theme(
            legend.position = 'none',
            axis.text.x = ggplot2::element_text(
              angle = 45,
              hjust = 1,
              vjust = 1
            )
          ) +
          ggplot2::scale_y_continuous(
            labels = function(x)
              sprintf(fmt = '%.1f', x)
          ) +
          ggforce::facet_row(
            ~ Parameter,
            scales = 'free',
            space = 'free',
            strip.position = 'right'
          ) +
          ggplot2::scale_color_manual(values = c('red' = 'red', 'black', 'black'))

        output$final_plot1 <-
          shiny::renderPlot({
            variables$plot1
          })
      })


      shiny::observeEvent(c(variables$Metrics, input$final_table_selector), {
        if (input$final_table_selector != '') {
          output$final_table <-
            shiny::renderTable(
              variables$Metrics %>%
                dplyr::filter(Parameter == input$final_table_selector,
                              Color == 'red') %>%
                dplyr::select(-x, -Color) %>%

                tidyr::spread(Category, Values)
            )
        }
      })

      ##save
      shiny::observeEvent(input$Save_plots, {
        if (!is.numeric(input$Save_plots)) {
          path = shinyFiles::parseSavePath(roots = roots, selection = input$Save_plots)

          ggplot2::ggsave(
            filename = path$datapath,
            plot = variables$plot1,
            device = grDevices::cairo_pdf,
            width = 12,
            height = 5
          )

        }
      })
      shiny::observeEvent(input$Save_table, {
        if (!is.numeric(input$Save_table)) {
          path = shinyFiles::parseSavePath(roots = roots, selection = input$Save_table)

          variables$Metrics %>%
            dplyr::select(-x, -Color) %>%
            readr::write_tsv(file = path$datapath)

        }
      })

    }
    ### Binning
    {
      #ui for auto selection of parameters
      output$Paramui_binning_out1 = shiny::renderUI({
        shiny::req(input$Param_binning == 'Auto')

        shiny::fluidRow(
          shiny::column(
            width = 3,
            shinyFiles::shinyDirButton(
              id = 'Bam_binning',
              label = 'Bam Files directory',
              title = 'Bam Files directory',
              style = 'width:100%;'
            )
          ),
          shiny::column(width = 6,
                        shiny::htmlOutput('Bam_binning_out'))
        )

      })

      #ui for auto selection of parameters

      output$Paramui_binning_out2 = shiny::renderUI({
        shiny::req(input$Param_binning == 'Manual')


        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::radioButtons(
              inputId = 'Reads_type_binning',
              label = 'Reads type',
              choices = c('PE', 'SE'),
              selected = variables_binning$read_type,
              inline = T,
              width = '100%'
            )
          ),
          shiny::column(
            width = 3,
            shiny::sliderInput(
              inputId = 'Read_size_binning',
              label = 'Read size',
              min = 20,
              max = 300,
              value = variables_binning$read_size,
              step = 1,
              width = '100%'
            )
          ),
          shiny::column(width = 3,
                        shiny::uiOutput('Reads_type_binning'))
        )

      })

      output$Reads_type_binning = shiny::renderUI({
        shiny::req(input$Reads_type_binning == 'PE')


        shiny::sliderInput(
          inputId = 'Fragment_size_binning',
          label = 'Fragment size',
          min = 0,
          max = 1000,
          value = variables_binning$frag_size,
          step = 1,
          width = '100%'
        )
      })

      #set shinyfiles
      shinyFiles::shinyDirChoose(
        input = input,
        id = 'Bam_binning',
        session = session,
        roots = roots,
        defaultRoot = 'Home'
      )
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'RefGenome_binning',
        session = session,
        roots = roots,
        defaultRoot = 'Home'
      )
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'index_binning',
        session = session,
        roots = roots,
        defaultRoot = 'Home'
      )
      shinyFiles::shinyDirChoose(
        input = input,
        id = 'out_binning',
        session = session,
        roots = roots,
        defaultRoot = 'Home'
      )
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'Blacklist_binning',
        session = session,
        roots = roots,
        defaultRoot = 'Home'
      )

      # if manual setting is put on SE hide fragment size
      shiny::observeEvent(input$Reads_type_binning, {
        if (input$Reads_type_binning == 'SE') {
          shinyjs::hide('Fragment_size_binning')
        } else{
          shinyjs::show('Fragment_size_binning')

        }
      })

      #select fasta file
      shiny::observeEvent(input$RefGenome_binning, {
        if (!is.numeric(input$RefGenome_binning)) {
          #recover path
          variables_binning$fasta = shinyFiles::parseFilePaths(roots = roots,
                                                               selection = input$RefGenome_binning)
          variables_binning$fasta = variables_binning$fasta$datapath
          output$RefGenome_binning_out = shiny::renderText({
            variables_binning$fasta
          })
          # first step to activate run
          variables_binning$activate_start1 = T
        } else{
          variables_binning$fasta = NULL
          output$RefGenome_binning_out = shiny::renderText({
            variables_binning$fasta
          })
          variables_binning$activate_start1 = F
        }
      })



      ###select index
      shiny::observeEvent(input$index_binning, {
        if (!is.numeric(input$index_binning)) {
          #recover index directory
          variables_binning$index = shinyFiles::parseFilePaths(selection = input$index_binning,
                                                               roots = roots)
          variables_binning$index = variables_binning$index$datapath

          #check that the selected file is actually part of an index
          if (stringr::str_detect(string = variables_binning$index, pattern = '.bt2')) {
            variables_binning$index = stringr::str_remove(variables_binning$index, pattern = '.[0-9]{1,2}.bt2|.rev.[0-9]{1,2}.bt2')

            output$index_binning_out <-
              shiny::renderText(variables_binning$index)
            variables_binning$activate_start2 = T
          } else{
            output$index_binning_out <-
              shiny::renderText(
                paste(
                  '<b><p style="color:#FF0000";> selected file does not belong to a bowtie2 index.</b></p>'
                )
              )
            variables$activate_start2 = F
          }
        } else{
          variables$activate_start2 = F
          output$index_binning_out <-
            shiny::renderText(NULL)
        }
      })

      ###select out folder
      shiny::observeEvent(input$out_binning,
                          {
                            if (!is.numeric(input$out_binning)) {
                              #recover out directory
                              variables_binning$out = shinyFiles::parseDirPath(selection = input$out_binning, roots = roots)

                              #render out dir
                              output$Output_binning_out <-
                                shiny::renderText(variables_binning$out)

                              variables_binning$activate_start3 = T
                            } else{
                              variables_binning$activate_start3 = F
                              output$Output_binning_out <-
                                shiny::renderText(NULL)

                            }
                          })

      #select bamfile folder
      shiny::observeEvent(input$Bam_binning, {
        if (!is.numeric(input$Bam_binning)) {
          #recover fastq directory
          variables_binning$Bam_dir = shinyFiles::parseDirPath(selection = input$Bam_binning, roots = roots)

          #check for fastq files
          if (length(list.files(path =  variables_binning$Bam_dir, pattern = '.bam')) >
              0) {
            output$Bam_binning_out <-
              shiny::renderText(variables_binning$Bam_dir)
            variables_binning$activate_start4 = T
          } else{
            output$Bam_binning_out <-
              shiny::renderText(
                paste(
                  '<b><p style="color:#FF0000";> selected directory does not contain fastq files</b></p>'
                )
              )
            variables_binning$activate_start4 = F
          }
        } else{
          variables_binning$activate_start4 = F
          output$Bam_binning_out <-
            shiny::renderText(NULL)
          variables_binning$Bam_dir = NULL
        }
      })

      #activate start4 if Manual parameters are selected
      shiny::observeEvent(input$Param_binning, {
        if (input$Param_binning == 'Auto') {
          variables_binning$activate_start4 = F
        } else{
          variables_binning$activate_start4 = T
        }
      })


      #select blacklist
      shiny::observeEvent(input$Blacklist_binning, {
        if (!is.numeric(input$Blacklist_binning)) {
          #recover index directory
          variables_binning$blacklist = shinyFiles::parseFilePaths(selection = input$Blacklist_binning,
                                                                   roots = roots)
          variables_binning$blacklist = variables_binning$blacklist$datapath

          #check that the selected file is actually part of an index
          if (Kronos.scRT::right_format(
            file_path = variables_binning$blacklist,
            columns_to_check = 3,
            logical = T
          )) {
            output$Blacklist_binning_out <-
              shiny::renderText(variables_binning$blacklist)
            variables_binning$activate_start5 = T
          } else{
            output$Blacklist_binning_out <-
              shiny::renderText(
                paste(
                  '<b><p style="color:#FF0000";> selected file does not have the required format.</b></p>'
                )
              )
            variables_binning$activate_start5 = F
          }
        } else{
          output$Blacklist_binning_out <-
            shiny::renderText(NULL)
          variables_binning$activate_start5 = T
          variables_binning$blacklist = NULL
        }

      })

      ##activate start
      shiny::observeEvent(
        c(
          variables_binning$activate_start1,
          variables_binning$activate_start2,
          variables_binning$activate_start3,
          variables_binning$activate_start4,
          variables_binning$activate_start5
        ),
        {
          if (variables_binning$activate_start1 &
              variables_binning$activate_start2 &
              variables_binning$activate_start3 &
              variables_binning$activate_start4 &
              variables_binning$activate_start5)
          {
            shinyjs::enable('run_binning')
          } else
          {
            shinyjs::disable('run_binning')

          }

        }
      )

      #manual parameters change
      shiny::observeEvent(input$Reads_type_binning, {
        variables_binning$read_type = input$Reads_type_binning
      })
      shiny::observeEvent(input$Read_size_binning, {
        variables_binning$read_size = input$Read_size_binning
      })
      shiny::observeEvent(input$Fragment_size_binning, {
        variables_binning$frag_size = input$Fragment_size_binning
      })

      #run
      shiny::observeEvent(input$run_binning, {
        if (input$run_binning > 0) {
          #disable buttons
          shinyjs::disable('RefGenome_binning')
          shinyjs::disable('index_binning')
          shinyjs::disable('out_binning')
          shinyjs::disable('Param_binning')
          shinyjs::disable('Bam_binning')
          shinyjs::disable('Coverage_binning')
          shinyjs::disable('Reads_type_binning')
          shinyjs::disable('Read_size_binning')
          shinyjs::disable('Fragment_size_binning')
          shinyjs::disable('Cores_binning')
          shinyjs::disable('Mappability_binning')
          shinyjs::disable('ErrorRate_binning')
          shinyjs::disable('Blacklist_binning')
          shinyjs::disable('chr_prefix_binning')
          shinyjs::disable('chr_range_binning')
          shinyjs::disable('BinSize_binning')
          shinyjs::disable('run_binning')

          #start binning
          bins = Kronos.scRT::binning(
            RefGenome = variables_binning$fasta,
            bowtie2_index =  variables_binning$index,
            bin_size =  input$BinSize_binning * 1000,
            read_size = variables_binning$read_size,
            fragment_size = variables_binning$frag_size,
            paired_ends = variables_binning$read_type == 'PE',
            directory_to_bamfiles = variables_binning$Bam_dir,
            tmp_dir = file.path(variables_binning$out, 'tmp'),
            upper_mappability_th = max(input$Mappability_binning),
            lower_mappability_th = min(input$Mappability_binning),
            black_list = variables_binning$blacklist,
            coverage = input$Coverage_binning,
            errorRate = input$ErrorRate_binning / 100,
            cores = input$Cores_binning,
            chr_prefix = ifelse(
              input$chr_prefix_binning == '',
              NULL,
              input$chr_prefix_binning
            ),
            chr_range = input$chr_range_binning,
            return_estimated_param = T
          )

          #if directory_to_bamfiles is given return parameters
          param = bins$param
          bins = bins$bins

          #write file
          readr::write_tsv(bins, file.path(
            variables_binning$out,
            paste0(
              basename(variables_binning$index),
              '_',
              input$BinSize_binning,
              'kb_bins_coverage_',
              input$Coverage_binning,
              'X_',
              paste0(
                'reads_',
                param$read_size,
                'bp_',
                ifelse(
                  param$paired_ends,
                  paste0('PE_fragmentSize_', param$fragment_size, 'bp_')
                  ,
                  'SE_'
                )
              ),
              ifelse(
                !is.null(variables_binning$blacklist),
                'blacklisted_',
                ''
              ),
              paste0(
                'error_rate_',
                input$ErrorRate_binning,
                '_min_mappability_',
                min(input$Mappability_binning),
                '_max_mappability_',
                max(input$Mappability_binning)
              ),
              '.tsv'
            )
          ))
          #enable
          shinyjs::enable('RefGenome_binning')
          shinyjs::enable('index_binning')
          shinyjs::enable('out_binning')
          shinyjs::enable('Param_binning')
          shinyjs::enable('Bam_binning')
          shinyjs::enable('Coverage_binning')
          shinyjs::enable('Reads_type_binning')
          shinyjs::enable('Read_size_binning')
          shinyjs::enable('Fragment_size_binning')
          shinyjs::enable('Cores_binning')
          shinyjs::enable('Mappability_binning')
          shinyjs::enable('ErrorRate_binning')
          shinyjs::enable('Blacklist_binning')
          shinyjs::enable('chr_prefix_binning')
          shinyjs::enable('chr_range_binning')
          shinyjs::enable('run_binning')
          shinyjs::enable('BinSize_binning')
        }
      })
    }
    #### CN
    {
      #set shiny Files buttons
      shinyFiles::shinyDirChoose(
        input = input,
        id = 'Bamdir_CN',
        session = session,
        roots = roots,
        defaultRoot =  'Home'
      )
      shinyFiles::shinyDirChoose(
        input = input,
        id = 'Output_CN',
        session = session,
        roots = roots,
        defaultRoot =  'Home'
      )
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'Bins_CN',
        session = session,
        roots = roots,
        defaultRoot =  'Home'
      )
      shinyFiles::shinyFileChoose(
        input = input,
        id = 'ChrSize_CN',
        session = session,
        roots = roots,
        defaultRoot =  'Home'
      )

      #set variables
      variables_cn = reactiveValues(
        Bamdir_CN = NULL,
        Output_CN = NULL,
        Bins_CN = NULL,
        ChrSize_CN = NULL,
        Samples = dplyr::tibble(),
        run_ok1 = F,
        run_ok2 = F,
        sample_counter = 0,
        group_counter = 0
      )

      # ui if user wants fix ploidy average
      output$FixAveragePloidy = shiny::renderUI({
        shiny::req(input$Ploidy_Restrictions_CN == 'Fix Average Ploidy')
        shiny::fluidRow(shiny::column(
          width = 3,
          shiny::numericInput(
            inputId = 'meanPloidy_CN',
            label = 'Average Ploidy',
            value = 2,
            min = 1,
            step = 0.1,
            width = '100%'
          )
        ))

      })

      # ui if user wants fix ploidy average
      output$FixPloidyRange = shiny::renderUI({
        shiny::req(input$Ploidy_Restrictions_CN == 'Fix Ploidy Range')
        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::numericInput(
              inputId = 'minPloidy_CN',
              label = 'Min accepted Ploidy',
              value = 2,
              min = 0,
              step = 0.1,
              width = '100%'
            )
          ),
          shiny::column(
            width = 3,
            shiny::numericInput(
              inputId = 'maxPloidy_CN',
              label = 'Max accepted ploidy',
              value = 8,
              min = 0,
              step = 0.1,
              width = '100%'
            )
          )
        )

      })

      #set output folder
      shiny::observeEvent(input$Output_CN, {
        if (!is.numeric(input$Output_C)) {
          variables_cn$Output_CN = shinyFiles::parseDirPath(roots = roots, selection = input$Output_CN)
        } else{
          variables_cn$Output_CN = NULL
        }
        output$Output_CN = shiny::renderText(variables_cn$Output_CN)

      })

      #set BinsFIle
      shiny::observeEvent(input$Bins_CN, {
        if (!is.numeric(input$Bins_CN)) {
          variables_cn$Bins_CN = shinyFiles::parseFilePaths(roots = roots, selection = input$Bins_CN)
          variables_cn$Bins_CN = variables_cn$Bins_CN$datapath
          variables_cn$Bins_CN = Kronos.scRT::right_format(
            file_path = variables_cn$Bins_CN,
            columns_to_check = c(
              'chr',
              'start',
              'end',
              'mappability',
              'mappability_th',
              'gc_frequency',
              'type'
            ),
            wrong_message = paste(
              '<b><p  style="color:#FF0000";>',
              variables_cn$Bins_CN,
              'does not have the right format!</p></b>'
            ),
            rigth_message = variables_cn$Bins_CN
          )
        } else{
          variables_cn$Bins_CN = NULL
        }
        output$Bins_CN = shiny::renderText(variables_cn$Bins_CN)

      })

      #authorize run if binsfile is good
      shiny::observeEvent(variables_cn$Bins_CN, {
        if (!is.null(variables_cn$Bins_CN) &
            stringr::str_detect(string = variables_cn$Bins_CN,
                                pattern = 'does not have the right format!</p></b>',
                                negate = T)) {
          variables_cn$run_ok1 = T
        } else{
          variables_cn$run_ok1 = F
        }
      })

      #set Chrsize
      shiny::observeEvent(input$ChrSize_CN, {
        if (!is.numeric(input$ChrSize_CN)) {
          variables_cn$ChrSize_CN = shinyFiles::parseFilePaths(roots = roots, selection = input$ChrSize_CN)
          variables_cn$ChrSize_CN = variables_cn$ChrSize_CN$datapath
          variables_cn$ChrSize_CN = Kronos.scRT::right_format(
            file_path = variables_cn$ChrSize_CN,
            columns_to_check = c('chr',
                                 'size'),
            wrong_message = paste(
              '<b><p  style="color:#FF0000";>',
              variables_cn$ChrSize_CN,
              'does not have the right format!</p></b>'
            ),
            rigth_message = variables_cn$ChrSize_CN
          )
        } else{
          variables_cn$ChrSize_CN = NULL
        }
        output$ChrSize_CN = shiny::renderText(variables_cn$ChrSize_CN)

      })

      #authorize run if Chrsize is good
      shiny::observeEvent(variables_cn$ChrSize_CN, {
        if (!is.null(variables_cn$ChrSize_CN) &
            stringr::str_detect(
              string = variables_cn$ChrSize_CN,
              pattern = 'does not have the right format!</p></b>',
              negate = T
            )) {
          variables_cn$run_ok2 = T
        } else{
          variables_cn$run_ok2 = F
        }
      })

      #set BAM folder
      shiny::observeEvent(input$Bamdir_CN, {
        if (!is.numeric(input$Bamdir_CN)) {
          variables_cn$Bamdir_CN = shinyFiles::parseDirPath(roots = roots, selection = input$Bamdir_CN)
          #check if file contains BAMS
          if (!any(stringr::str_detect(
            string = list.files(variables_cn$Bamdir_CN),
            pattern = '.bam$'
          ))) {
            variables_cn$Bamdir_CN = paste(
              '<b><p  style="color:#FF0000";>',
              variables_cn$Bamdir_CN,
              'does not contain bam files</p></b>'
            )
            shinyjs::disable('AddSample_CN')
          } else{
            shinyjs::enable('AddSample_CN')
          }
        } else{
          variables_cn$Bamdir_CN = NULL
          shinyjs::disable('AddSample_CN')
        }
        output$Bamdir_CN = shiny::renderText(variables_cn$Bamdir_CN)

      })

      #add samples
      shiny::observeEvent(input$AddSample_CN, {
        if (input$AddSample_CN > 0) {
          shinyjs::disable('AddSample_CN')
          variables_cn$Samples = rbind(
            variables_cn$Samples,
            dplyr::tibble(
              Basename = input$basename_CN,
              Group = input$group_CN,
              Bamfiles = variables_cn$Bamdir_CN
            )
          )

          output$Samples_CN = shiny::renderTable(variables_cn$Samples)

          #reset interface
          if (any(
            stringr::str_detect(
              string = variables_cn$Samples$Basename,
              pattern = '^Exp_{0,1}[0-9]{0,3}$'
            )
          )) {
            variables_cn$sample_counter = variables_cn$sample_counter + 1
            shiny::updateTextInput(
              inputId = 'basename_CN',
              value = paste0('Exp_', variables_cn$sample_counter)
            )
          } else{
            shinyjs::reset(input$basename_CN)
          }
          if (any(
            stringr::str_detect(string = variables_cn$Samples$Group, pattern = '^Exp_{0,1}[0-9]{0,3}$')
          )) {
            variables_cn$group_counter = variables_cn$group_counter + 1
            shiny::updateTextInput(
              inputId = 'group_CN',
              value = paste0('Exp_', variables_cn$group_counter)
            )
          } else{
            shinyjs::reset(input$group_CN)
          }

        }
      })

      #activate run
      shiny::observeEvent(
        c(
          variables_cn$Samples,
          variables_cn$run_ok1,
          variables_cn$run_ok2,
          variables_cn$Output_CN
        ),
        {
          if (variables_cn$run_ok1 &
              variables_cn$run_ok2 &
              nrow(variables_cn$Samples) > 0 &
              !is.null(variables_cn$Output_CN)) {
            shinyjs::enable('Run_CN')
          } else{
            shinyjs::disable('Run_CN')
          }
        }
      )

      shiny::observeEvent(input$Run_CN, {
        if (input$Run_CN > 0) {
          #load df
          Bins_CN = readr::read_tsv(variables_cn$Bins_CN)
          ChrSize_CN = readr::read_tsv(variables_cn$ChrSize_CN)
          shinyjs::disable('Run_CN')
          for (i in 1:nrow(variables_cn$Samples)) {
            if (input$Ploidy_Restrictions_CN != 'Fix Average Ploidy') {
              results = Kronos.scRT::CallCNV(
                directory = variables_cn$Samples$Bamfiles[i],
                bins = Bins_CN,
                chrom_size =  ChrSize_CN,
                basename = variables_cn$Samples$Basename[i],
                group = variables_cn$Samples$Group[i],
                tmp_dir = file.path(variables_cn$Output_CN, 'tmp'),
                min_n_reads = input$min_reads_CN,
                mim_mean_CN_accepted = ifelse(
                  is.null(input$minPloidy_CN),
                  2,
                  input$minPloidy_CN
                ),
                max_mean_CN_accepted = ifelse(
                  is.null(input$maxPloidy_CN),
                  8,
                  input$maxPloidy_CN
                ),
                chr_prefix = input$CHR_CN,
                chr_range = input$CHR_RANGE_CN,
                cores = input$Cores_CN
              )
            } else{
              results = Kronos.scRT::CallCNV(
                directory = variables_cn$Samples$Bamfiles[i],
                bins = Bins_CN,
                chrom_size =  ChrSize_CN,
                basename = variables_cn$Samples$Basename[i],
                group = variables_cn$Samples$Group[i],
                tmp_dir = file.path(variables_cn$Output_CN, 'tmp'),
                min_n_reads = input$min_reads_CN,
                ploidy = ifelse(
                  is.null(input$meanPloidy_CN),
                  2,
                  input$meanPloidy_CN
                ),
                chr_prefix = input$CHR_CN,
                chr_range = input$CHR_RANGE_CN,
                cores = input$Cores_CN
              )
            }
            dir = file.path(variables_cn$Output_CN,
                            variables_cn$Samples$Basename[i])
            if (!dir.exists(dir)) {
              dir.create(dir)
            }
            results$PerCell %>%
              readr::write_csv(file = file.path(
                dir,
                paste0(variables_cn$Samples$Basename[i], '_PerCell.csv')
              ))
            results$CNV %>%
              readr::write_tsv(file = file.path(
                dir,
                paste0(variables_cn$Samples$Basename[i], '_scCNV.tsv')
              ))
          }
        }
      })

    }


  }

  # Run the application
  shiny::shinyApp(ui = ui,
                  server = server)
