Dim_red_sub_pop_ui = function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::fluidRow(
      align = 'center',
      shiny::column(
        width = 3,
        offset = 1,
        shiny::radioButtons(
          inputId = ns('setting'),
          label = 'Dimensionality reduction type',
          choices = c('UMAP', 'T-SNE'),
          inline = T,
          selected = 'UMAP',
          width = '100%'
        )
      ),
      shiny::column(
        width = 3,
        shiny::numericInput(
          inputId = ns('seed'),
          value = as.integer(Sys.Date()),
          min = 1,
          max = 10 ^ 10,
          step = 1,
          label = 'Seed',
          width = '100%'
        )
      ),
      shiny::column(
        width = 3,
        shiny::radioButtons(
          inputId = ns('what'),
          label = 'Parameter',
          choices = c('scRep', 'scCN'),
          inline = T,
          selected = 'scRep',
          width = '100%'
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 3,
        offset = 1,
        shiny::actionButton(
          inputId = ns('random'),
          label = 'Random seed',
          width = '100%'
        )
      ),
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
          inputId = ns('save'),
          label = 'Save',
          width = '100%'
        )
      ),
      shiny::column(width = 1,
                    hw_plot_ui(
                      ns('save_First'),
                      up = F,
                      height = 7,
                      width = 10
                    ))
    ),
    shiny::fluidRow(align = 'center',
                    shiny::plotOutput(ns('dred_out'),width = '50%',height = 'auto')),
    shiny::fluidRow(
      shiny::column(
        width = 3,
        shiny::radioButtons(
          inputId = ns('color'),
          label = 'Color',
          choices = c('Sample', '%Rep'),
          width = '100%',
          inline = T
        )
      ),
      shiny::column(width = 3,
                    shinyjs::hidden(
                      shiny::checkboxGroupInput(
                        inputId = ns("Sample"),
                        label = "Samples to keep",
                        choices = '',
                        width = '100%',
                        inline = T
                      )
                    )),
      shiny::column(width = 3,
                    shinyjs::hidden(
                      shiny::checkboxGroupInput(
                        inputId = ns("phase"),
                        label = "Phase",
                        choices = '',
                        width = '100%',
                        inline = T
                      )
                    ))
    ),
    shiny::fluidRow(shiny::column(
      width =  3,
      shinyWidgets::switchInput(
        inputId = ns('start_subpop'),
        label = 'Sub-pop mode',
        value = F,
        inline = T,
        labelWidth = '100%',
        disabled = T,
        onStatus = "success",
        offStatus = "danger"
      )
    )),
    shiny::uiOutput(ns('SubPopUI')),
    shiny::uiOutput(ns('SubPopUI2'))
  )

}


Dim_red_sub_pop_server = function(id, scCN, out, Inputfolder,colors, cores = 3) {
  shiny::moduleServer(id,
                      function(input,
                               output,
                               session,
                               scCNV = scCN,
                               Cores = cores,
                               Out = out,
                               col=colors,
                               IN = Inputfolder,
                               ID=id) {
                        #load required operators
                        `%>%` = tidyr::`%>%`

                        #recover colors
                        Colors=col$color
                        names(Colors)=col$basename

                        #initialize
                        Data_to_plot = reactiveValues(
                          selected_data = dplyr::tibble(),
                          scCNV = dplyr::tibble(),
                          plot = ggplot2::ggplot()
                        )
                        shinyjs::disable('save')
                        shinyjs::disable('run')

                        if (length(unique(scCNV$basename)) > 1) {
                          shinyjs::show('Sample')
                        }

                        shiny::updateCheckboxGroupInput(
                          inputId = 'Sample',
                          choices = unique(scCNV$basename),
                          selected = unique(scCNV$basename),
                          inline = T
                        )
                        #depending on Dim red mode make different options available
                        shiny::req('what')
                        shiny::observe({
                          if (input$what == 'scCN') {
                            shiny::updateCheckboxGroupInput(
                              inputId = 'phase',
                              choices = unique(scCNV$Phase),
                              selected = unique(scCNV$Phase),
                              inline = T
                            )
                            shinyjs::show('phase')
                          } else{
                            shiny::updateCheckboxGroupInput(
                              inputId = 'phase',
                              choices = 'S',
                              selected = 'S',
                              inline = T
                            )
                            shinyjs::hide('phase')
                          }
                        })
                        #activate run
                        shinyjs::enable('run')

                        #random seed
                        shiny::observeEvent(input$random, {
                          shiny::updateNumericInput(inputId = 'seed',
                                                    value = round(stats::runif(1), 9) * 10 ^ 9)
                        })

                        shiny::observeEvent(input$run, {
                          if (input$run > 0 &
                              length(input$Sample) > 0 &
                              length(input$phase) > 0) {
                            shinyjs::disable('run')
                            shinyjs::disable('save')

                            shinyWidgets::updateSwitchInput(
                              inputId = 'start_subpop',
                              value = F,
                              session = session,
                              disabled = T
                            )
                            #select columns of interest and reformat
                            Data_to_plot$selected_data = scCNV %>%
                              dplyr::mutate(pos = paste0(chr, ':', start, '-', end)) %>%
                              dplyr::filter(basename %in% input$Sample,
                                            Phase %in% input$phase) %>%
                              dplyr::select(
                                'pos',
                                'Cell',
                                'PercentageReplication',
                                'basename',
                                'group',
                                'data' = dplyr::case_when(input$what ==
                                                            'scRep' ~ 'Rep',
                                                          input$what == 'scCN' ~ 'CN')
                              ) %>%
                              dplyr::mutate(data = as.numeric(data)) %>%
                              tidyr::spread(pos, data)
                            #remove col with na
                            Data_to_plot$selected_data = Data_to_plot$selected_data[, colSums(is.na(Data_to_plot$selected_data)) == 0]



                            #create matrix data and save extra info
                            mat = Data_to_plot$selected_data[, -c(1:4)]
                            Data_to_plot$selected_data = Data_to_plot$selected_data[, c(1:4)]


                            if (input$what == 'scRep') {
                              #calculate simple matching coef
                              results = as.matrix(ade4::dist.binary(df = mat, method = 2))
                            } else{
                              #convert to matrix
                              results = as.matrix(mat)
                            }

                            rm('mat')

                            #if random set new seed
                            shiny::req(input$seed)
                            if (input$setting == 'T-SNE') {
                              Perplex = ceiling(nrow(results) / 50)
                              tsne <-
                                Kronos.scRT::TSNE(
                                  X = results,
                                  seed = input$seed,
                                  dims = 2,
                                  perplexity = ifelse(Perplex < 10, 10, Perplex),
                                  check_duplicates = F,
                                  theta = 0.25,
                                  is_distance = dplyr::case_when(input$what ==
                                                                   'scRep' ~ T,
                                                                 input$what == 'scCN' ~ F),
                                  verbose = F,
                                  max_iter = 5000,
                                  num_threads = Cores,
                                  partial_pca = T
                                )

                              Data_to_plot$scCNV = Data_to_plot$selected_data %>% dplyr::mutate(x = tsne$Y[, 1],
                                                                                                y = tsne$Y[, 2])

                              shinyjs::enable('save')
                              shinyjs::enable('run')
                              shinyWidgets::updateSwitchInput(inputId = 'start_subpop', disabled = F)
                            } else if (input$setting == 'UMAP') {
                              umap <-
                                umap::umap(
                                  d = results,
                                  input = dplyr::case_when(
                                    input$what == 'scRep' ~ 'dist',
                                    input$what == 'scCN' ~ 'data'
                                  ),
                                  random_state = input$seed
                                )
                              Data_to_plot$scCNV = Data_to_plot$selected_data %>% dplyr::mutate(x =
                                                                                                  umap$layout[, 1],
                                                                                                y = umap$layout[, 2])
                            }

                            shinyjs::enable('save')
                            shinyjs::enable('run')
                            shinyWidgets::updateSwitchInput(inputId = 'start_subpop', disabled = F)
                          }

                          shiny::observeEvent(c(Data_to_plot$scCNV, input$color), {
                            #plot
                            Data_to_plot$plot = Data_to_plot$scCNV %>%
                              dplyr::select(
                                'x',
                                'y',
                                'color' = dplyr::case_when(
                                  input$color == 'Sample' ~ 'basename',
                                  input$color == '%Rep' ~ 'PercentageReplication'
                                )
                              ) %>%
                              ggplot2::ggplot() +
                              ggplot2::geom_point(
                                ggplot2::aes(x, y, color = color),
                                alpha = 0.4,
                                size = 2
                              ) +
                              ggplot2::labs(
                                x = paste(input$setting, '1', sep = '-'),
                                y = paste(input$setting, '2', sep = '-'),
                                color = input$color
                              )

                            if (input$color == '%Rep') {
                              Data_to_plot$plot = Data_to_plot$plot +
                                ggplot2::scale_color_gradient2(
                                  low = "#FFEA46FF",
                                  mid = "#7C7B78FF",
                                  high = "#00204DFF",
                                  lim = c(0, 1),
                                  midpoint = 0.5
                                )
                            }else{
                              Data_to_plot$plot = Data_to_plot$plot +
                                ggplot2::scale_color_manual(values = Colors)
                            }

                            output$dred_out = shiny::renderPlot(Data_to_plot$plot,
                                                                height = function() {
                                                                  session$clientData[[paste0('output_', ID, '-dred_out_width')]]*0.8})
                          })
                        })

                        shiny::observeEvent(input$save, {
                          #create dir
                          if (!dir.exists(file.path(Out, input$setting))) {
                            dir.create(file.path(Out, input$setting), recursive = T)
                          }
                          #plot size

                          First_plot_dim=hw_plot_server('save_First')
                          First_plot_dim=First_plot_dim()
                          #plot
                          ggplot2::ggsave(
                            plot = Data_to_plot$plot,
                            device = grDevices::cairo_pdf,
                            width = First_plot_dim$width,
                            height = First_plot_dim$height,
                            units = First_plot_dim$unit,
                            filename = file.path(
                              Out,
                              input$setting,
                              paste0(
                                basename(Out),
                                '_',
                                input$setting,
                                '_on_',
                                dplyr::case_when(input$what == 'scRep' ~ 'Rep',
                                                 input$what == 'scCN' ~ 'CN'),
                                '_seed_',
                                input$seed,
                                '_Samples_',
                                paste(input$Sample, collapse = '_'),
                                '_',
                                paste(
                                  stringr::str_replace_all(
                                    input$phase,
                                    pattern = '/',
                                    replacement = '_'
                                  ),
                                  collapse = '_and_'
                                ),
                                '_phase',
                                '_color_',
                                ifelse(input$color == '%Rep', 'RepPerc', input$color),
                                '.pdf'
                              )
                            )
                          )

                          #if file does not exist save plot info
                          if (!file.exists(file.path(
                            Out,
                            input$setting,
                            paste0(
                              basename(Out),
                              '_',
                              input$setting,
                              '_seed_',
                              input$seed,
                              '_Samples_',
                              paste(input$Sample, collapse = '_'),
                              '_',
                              paste(
                                stringr::str_replace_all(
                                  input$phase,
                                  pattern = '/',
                                  replacement = '_'
                                ),
                                collapse = '_and_'
                              ),
                              '_phase.tsv'
                            )
                          ))) {
                            readr::write_tsv(x = Data_to_plot$scCNV,
                                             file = file.path(
                                               Out,
                                               input$setting,
                                               paste0(
                                                 basename(Out),
                                                 '_',
                                                 input$setting,
                                                 '_seed_',
                                                 input$seed,
                                                 '_Samples_',
                                                 paste(input$Sample, collapse = '_'),
                                                 '_',
                                                 paste(
                                                   stringr::str_replace_all(
                                                     input$phase,
                                                     pattern = '/',
                                                     replacement = '_'
                                                   ),
                                                   collapse = '_and_'
                                                 ),
                                                 '_phase.tsv'
                                               )
                                             ))
                          }

                        })

                        shiny::observeEvent(c(
                          input$what,
                          input$Sample,
                          input$phase,
                          input$Parameter
                        ),
                        {
                          shinyjs::disable('save')
                          shinyWidgets::updateSwitchInput(inputId = 'start_subpop', disabled = T,value = FALSE)
                        })

                        shiny::observeEvent(input$start_subpop, {
                          if (input$start_subpop) {
                            # create a reactive variable
                            data_subpop = shiny::reactiveValues(
                              counter = 1,
                              G1 = NULL,
                              S = NULL,
                              subpopulation = NULL,
                              final = dplyr::tibble(a = 1),
                              click = 1,
                              dbclick = 1,
                              poligon = dplyr::tibble(),
                              summary = dplyr::tibble(phase = NA),
                              found_matches = F
                            )
                            ns <- session$ns
                            #subpop UI
                            output$SubPopUI <-
                              shiny::renderUI({
                                shiny::fluidPage(
                                  shinyjs::useShinyjs(),
                                  title = 'Define Sub population',
                                  width = 12,
                                  height = '100%',
                                  shiny::tags$div(
                                    class = "header",
                                    checked = NA,
                                    shiny::tags$h4(
                                      "Select your sub-population on this first plot. Click on two different point to create the first side of you gate and then continue till you are satisfied. Use a double click once you have done and you want to close a gate."
                                    ),
                                    style = 'padding:30px;'
                                  ),
                                  #upper row with two plots
                                  shiny::fluidRow(
                                    shiny::column(
                                      width = 6,
                                      shiny::plotOutput(
                                        ns("plot1_subpop"),
                                        click = ns('click_subpop'),
                                        dblclick = ns('dbclick_subpop'),
                                        width = '100%',
                                        height = 'auto'
                                      ),
                                      align = "center"
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::plotOutput(ns(
                                                    "plot2_subpop"
                                                  ),
                                                  width = '100%',
                                                  height = 'auto'))
                                  ),
                                  shiny::fluidRow(
                                    shiny::column(
                                      width = 3,
                                      align = 'center',
                                      shinyjs::disabled(
                                        shiny::actionButton(
                                          inputId = ns('SaveGroups_subpop'),
                                          label = 'Save sub-groups',
                                          width = '100%'
                                        )
                                      )
                                    ),
                                    shiny::column(width = 1,
                                                  hw_plot_ui(
                                                    ns('save_second'),
                                                    up = F,
                                                    height = 7,
                                                    width = 10
                                                  )),
                                    shiny::column(
                                      width = 3,
                                      align = 'center',
                                      shinyjs::disabled(
                                        shiny::actionButton(
                                          inputId = ns('RmLast_subpop'),
                                          label = 'Remove last group',
                                          width = '100%'
                                        )
                                      )
                                    ),
                                    shiny::column(
                                      width = 3,
                                      align = 'center',
                                      shinyjs::disabled(
                                        shiny::actionButton(
                                          inputId = ns('Cancel_subpop'),
                                          label = 'Cancel',
                                          width = '100%'
                                        )
                                      )
                                    ),
                                    shiny::column(
                                      width =  3,
                                      shinyWidgets::switchInput(
                                        inputId = ns('AssignGroups_subpop'),
                                        label = 'Assign G to S-phase',
                                        value = F,
                                        inline = T,
                                        labelWidth = '100%',
                                        disabled = T,
                                        onStatus = "success",
                                        offStatus = "danger"
                                      )
                                    )
                                  ),
                                  shiny::fluidRow(shiny::tableOutput(ns(
                                    'group_phase_subpop'
                                  ))),
                                  shiny::fluidRow()
                                )
                              })
                            #initialize
                            output$group_phase_subpop <-
                              shiny::renderTable(NULL)
                            output$plot2_subpop <-
                              shiny::renderPlot({
                                NULL
                              })
                            #subpop functions
                            {
                              # copy TSNE or UMAP data
                              data_subpop$data <- shiny::isolate({
                                Data_to_plot$scCNV
                              }) %>%
                                dplyr::mutate(phase = NA, subgroup = NA)

                              data_subpop$plot = ggplot2::ggplot(Data_to_plot$scCNV,
                                                                 ggplot2::aes(x, y, color = PercentageReplication)) +
                                ggplot2::geom_point() +
                                ggplot2::scale_color_gradient2(
                                  low = "#FFEA46FF",
                                  mid = "#7C7B78FF",
                                  high = "#00204DFF",
                                  lim = c(0, 1),
                                  midpoint = 0.5
                                ) +
                                ggplot2::theme(legend.position = 'top') +
                                ggplot2::labs(
                                  x = paste(input$setting, 1, sep = '-'),
                                  y = paste(input$setting, 2, sep = '-')
                                )

                              #plot TSNE or UMAP data using % Rep as color
                              output$plot1_subpop <-
                                shiny::renderPlot({
                                  if (nrow(data_subpop$poligon) == 1) {
                                    data_subpop$plot +
                                      ggplot2::geom_point(
                                        data = data_subpop$poligon,
                                        ggplot2::aes(x, y, group = group),
                                        color = 'red'
                                      )
                                  } else if (nrow(data_subpop$poligon) >
                                             1) {
                                    data_subpop$plot +
                                      ggplot2::geom_path(
                                        data = data_subpop$poligon,
                                        ggplot2::aes(x, y, group = group),
                                        color = 'red'
                                      ) +
                                      ggplot2::geom_point(
                                        data = data_subpop$poligon,
                                        ggplot2::aes(x, y, group = group),
                                        color = 'red'
                                      )
                                  } else{
                                    data_subpop$plot
                                  }
                                },
                                height = function() {
                                  session$clientData[[paste0('output_', ID, '-plot1_subpop_width')]]*0.8})

                              #disable save, cancel and remove last if no subgrous are defined
                              shiny::observe({
                                if (nrow(data_subpop$summary) > 1) {
                                  shinyjs::enable('SaveGroups_subpop')
                                } else{
                                  shinyjs::disable('SaveGroups_subpop')
                                }

                                if (nrow(data_subpop$summary) > 0) {
                                  shinyjs::enable('RmLast_subpop')
                                  shinyjs::enable('Cancel_subpop')

                                } else{
                                  shinyjs::disable('RmLast_subpop')
                                  shinyjs::disable('Cancel_subpop')

                                }
                              })


                              #on click
                              shiny::observeEvent(input$click_subpop, {
                                if (data_subpop$click == 1) {
                                  data_subpop$startX = input$click_subpop$x
                                  data_subpop$startY = input$click_subpop$y
                                }

                                data_subpop$poligon = rbind(
                                  data_subpop$poligon,
                                  dplyr::tibble(
                                    x = input$click_subpop$x,
                                    y = input$click_subpop$y,
                                    group = paste0('Sub-group', data_subpop$dbclick)
                                  )
                                )
                                data_subpop$click = data_subpop$click +
                                  1

                              })

                              #on double click

                              shiny::observeEvent(input$dbclick_subpop, {
                                data_subpop$poligon = rbind(
                                  data_subpop$poligon,
                                  dplyr::tibble(
                                    x = data_subpop$startX,
                                    y = data_subpop$startY,
                                    group = paste0('Sub-group', data_subpop$dbclick)
                                  )
                                )

                                selected_list = as.logical(
                                  sp::point.in.polygon(
                                    point.x = data_subpop$data$x,
                                    point.y = data_subpop$data$y,
                                    pol.x = data_subpop$poligon$x[data_subpop$poligon$group ==
                                                                    paste0('Sub-group', data_subpop$dbclick)],
                                    pol.y = data_subpop$poligon$y[data_subpop$poligon$group ==
                                                                    paste0('Sub-group', data_subpop$dbclick)]
                                  )
                                )
                                data_subpop$data$subgroup[selected_list] = paste0('Sub-group', data_subpop$dbclick)
                                data_subpop$data$phase[selected_list] <-
                                  ifelse(
                                    stats::median(data_subpop$data$PercentageReplication[selected_list]) < .05,
                                    'G1/G2',
                                    'S'
                                  )

                                data_subpop$dbclick = data_subpop$dbclick + 1
                                data_subpop$click = 1
                                #if there is nothing to plot, don't plot
                                output$plot2_subpop <-
                                  shiny::renderPlot({
                                    if (nrow(data_subpop$data %>%
                                             tidyr::drop_na()) != 0) {
                                      ggplot2::ggplot(data_subpop$data,
                                                      ggplot2::aes(x, y, color = subgroup)) +
                                        ggplot2::geom_point(na.rm = TRUE) +
                                        ggplot2::theme(legend.position = 'top') +
                                        ggplot2::labs(
                                          x = paste(input$setting, 1, sep = '-'),
                                          y = paste(input$setting, 2, sep = '-')
                                        )
                                    } else{
                                      NULL
                                    }
                                  },
                                  height = function() {
                                    session$clientData[[paste0('output_', ID, '-plot2_subpop_width')]]*0.8})

                                data_subpop$summary = data_subpop$data %>%
                                  dplyr::group_by(subgroup, phase) %>%
                                  tidyr::drop_na() %>%
                                  dplyr::summarise(n = dplyr::n()) %>%
                                  dplyr::ungroup()%>%
                                  dplyr::arrange(subgroup) %>%
                                  `colnames<-`(c('Sub-group', 'Phase', 'Number of cells'))

                                output$group_phase_subpop <-
                                  shiny::renderTable(data_subpop$summary)

                              })

                              #remove last
                              shiny::observeEvent(input$RmLast_subpop, {
                                if (input$RmLast_subpop > 0) {
                                  #decrease dbclick of 1
                                  data_subpop$dbclick = data_subpop$dbclick - 1

                                  #filter out all the df
                                  data_subpop$data$subgroup[data_subpop$data$subgroup == paste0('Sub-group', data_subpop$dbclick)] = NA

                                  data_subpop$poligon = data_subpop$poligon %>%
                                    dplyr::filter(group !=  paste0('Sub-group', data_subpop$dbclick))

                                  data_subpop$summary = data_subpop$summary  %>%
                                    dplyr::filter(`Sub-group` !=  paste0('Sub-group', data_subpop$dbclick))

                                }
                              })
                              # clean selections
                              shiny::observeEvent(input$Cancel_subpop, {
                                if (input$Cancel_subpop > 0) {
                                  data_subpop$data$phase = NA
                                  data_subpop$data$subgroup = NA
                                  data_subpop
                                  data_subpop$counter = 1
                                  data_subpop$dbclick =  1
                                  data_subpop$click = 1
                                  data_subpop$poligon = dplyr::tibble()
                                  data_subpop$summary = dplyr::tibble(Phase= character())

                                  shinyWidgets::updateSwitchInput(
                                    inputId = 'AssignGroups_subpop',
                                    value = F,
                                    disabled = T
                                  )
                                }
                              })

                              shiny::observeEvent(input$SaveGroups_subpop, {
                                if (input$SaveGroups_subpop > 0) {
                                  if (!dir.exists(file.path(Out, input$setting, 'SubPop'))) {
                                    dir.create(file.path(Out, input$setting, 'SubPop'),
                                               recursive = T)
                                  }

                                  if (input$AssignGroups_subpop &
                                      data_subpop$found_matches) {
                                    data_subpop$data %>%
                                      dplyr::select(Cell, basename, group, subpopulation) %>%
                                      readr::write_tsv(file.path(
                                        Out,
                                        input$setting,
                                        'SubPop',
                                        paste0(paste(
                                          unique(data_subpop$data$group),
                                          collapse = '_'
                                        ),
                                        '_subpop.tsv')
                                      ))

                                    Second_plot_dim=hw_plot_server('save_second')
                                    Second_plot_dim=Second_plot_dim()

                                    ggplot2::ggsave(
                                      plot = data_subpop$plot1(),
                                      filename = file.path(
                                        Out,
                                        input$setting,
                                        'SubPop',
                                        paste0(
                                          paste(unique(data_subpop$data$group),
                                                collapse = '_'),
                                          '_subgroups.pdf'
                                        )
                                      ),
                                      device = grDevices::cairo_pdf,
                                      width = Second_plot_dim$width,
                                      height = Second_plot_dim$height,
                                      units = Second_plot_dim$unit,
                                    )
                                  } else{

                                    data_subpop$data %>%
                                      dplyr::filter(!is.na(subgroup)) %>%
                                      dplyr::select(Cell, basename, group, "subpopulation" =
                                                      subgroup) %>%
                                      readr::write_tsv(file.path(
                                        Out,
                                        input$setting,
                                        'SubPop',
                                        paste0(paste(
                                          unique(data_subpop$data$group),
                                          collapse = '_'
                                        ),
                                        #depending on the mode save different ending
                                        ifelse(input$what == 'scCN',
                                        '_unmatched_subpop.tsv', '_S_phase_subpop.tsv'))
                                      ))
                                    }
                                }
                              })

                              shiny::observeEvent(input$AssignGroups_subpop, {
                                #if assign groups is on
                                if (input$AssignGroups_subpop) {
                                  shinyjs::disable('SaveGroups_subpop')
                                  output$SubPopUI2 <-
                                    shiny::renderUI({
                                      shiny::fluidPage(
                                        title = 'Results assignment',
                                        width = 12,
                                        height = '100%',
                                        class = "header",
                                        checked = NA,
                                        #plot grouping results
                                        shiny::fluidRow(
                                          shiny::htmlOutput(ns('results_subpop')),
                                          shiny::plotOutput(ns("plot3_subpop"),width = '50%',height = 'auto')
                                        )
                                      )
                                    })

                                  data_subpop$scCNV = shiny::isolate(
                                    dplyr::inner_join(
                                    scCNV %>%
                                      dplyr::select(
                                        chr,
                                        start,
                                        end,
                                        CN,
                                        PercentageReplication,
                                        Cell,
                                        basename,
                                        group
                                      ),
                                    data_subpop$data %>%
                                      dplyr::select(
                                        Cell,
                                        PercentageReplication,
                                        basename,
                                        group,
                                        subgroup,
                                        phase
                                      ) %>%
                                      tidyr::drop_na(),
                                    by = c(
                                      "PercentageReplication",
                                      "Cell",
                                      "basename",
                                      "group"
                                    )
                                  ))

                                  data_subpop$G1 = data_subpop$scCNV %>%
                                    dplyr::filter(phase == 'G1/G2') %>%
                                    dplyr::group_by(chr, start, end, subgroup) %>%
                                    dplyr::summarise(CN = stats::median(CN)) %>%
                                    dplyr::ungroup() %>%
                                    tidyr::spread(key = subgroup,
                                                  value = CN,
                                                  fill = 0) %>%
                                    tidyr::gather(key = 'subgroup',
                                                  value = 'CN',
                                                  -chr,
                                                  -start,
                                                  -end) %>%
                                    dplyr::group_by(chr, start, end) %>%
                                    dplyr::filter(CN != stats::median(CN)) %>%
                                    dplyr::ungroup()

                                  data_subpop$S = data_subpop$scCNV %>%
                                    dplyr::filter(phase == 'S') %>%
                                    dplyr::group_by(chr, start, end, subgroup) %>%
                                    dplyr::summarise(CN = list(stats::quantile(CN, seq(
                                      0.1, 0.9, 0.01
                                    )))) %>%
                                    dplyr::ungroup() %>%
                                    dplyr::mutate(percentile = list(seq(0.1, 0.9, 0.01))) %>%
                                    tidyr::unnest(c('CN', 'percentile')) %>%
                                    dplyr::ungroup()

                                  #keep only variable bins
                                  bins_to_keep = data_subpop$G1 %>%
                                    dplyr::group_by(chr, start, end) %>%
                                    dplyr::filter(CN != stats::median(CN)) %>%
                                    dplyr::select(chr, start, end) %>%
                                    unique() %>%
                                    dplyr::ungroup()

                                  if (nrow(bins_to_keep) == 0) {
                                    output$results_subpop = shiny::renderText("<b> G1 cells don't have differet CNVs. </b>")
                                    data_subpop$found_matches = F
                                  } else{
                                    # associate a G1/G2 to an S phase using sum(abs(log2(x+1/y+1)))
                                    data_subpop$subpopulation = data_subpop$G1 %>%
                                      dplyr::inner_join(bins_to_keep,
                                                        by = c("chr", "start", "end")) %>%
                                      dplyr::rename('subgroup_G1' = subgroup,
                                                    'CN_G1' = CN) %>%
                                      dplyr::left_join(
                                        data_subpop$S,
                                        by = c("chr", "start", "end"),
                                        fill = 0
                                      ) %>%
                                      dplyr::group_by(subgroup_G1,
                                                      subgroup,
                                                      percentile) %>%
                                      dplyr::summarise(dist = sum(abs(log2((CN_G1 + 1) / (CN + 1)
                                      )))) %>%
                                      dplyr::group_by(subgroup_G1, percentile) %>%
                                      dplyr::filter(dist == min(dist)) %>%
                                      dplyr::group_by(subgroup_G1, subgroup) %>%
                                      dplyr::summarise(n = dplyr::n(),
                                                       min_dist = min(dist)) %>%
                                      dplyr::group_by(subgroup_G1) %>%
                                      dplyr::filter(n == max(n)) %>%
                                      #deals with ties
                                      dplyr::group_by(subgroup_G1) %>%
                                      dplyr::mutate(n_uG1 = dplyr::n()) %>%
                                      dplyr::filter(n_uG1 == 1 |
                                                      min_dist == min(min_dist))

                                    #check if there is only one result per group
                                    if (length(unique(data_subpop$subpopulation$subgroup_G1)) == length(unique(data_subpop$subpopulation$subgroup))) {
                                      # select smallest distance
                                      data_subpop$subpopulation = data_subpop$subpopulation %>%
                                        dplyr::ungroup() %>%
                                        dplyr::mutate(subpopulation = paste0('Sub-population_', 1:dplyr::n())) %>%
                                        dplyr::select(-min_dist, -n_uG1, -n) %>%
                                        dplyr::rename('G1/G2' = subgroup_G1,
                                                      'S' = subgroup) %>%
                                        tidyr::gather(phase,
                                                      subgroup, -subpopulation)


                                      data_subpop$data = data_subpop$data %>%
                                        dplyr::inner_join(data_subpop$subpopulation,
                                                          by = c("phase", "subgroup"))
                                      #plot resutls
                                      data_subpop$plot1 <- reactive(
                                        ggplot2::ggplot(
                                          data_subpop$data,
                                          ggplot2::aes(x,
                                                       y,
                                                       color = subpopulation)
                                        ) +
                                          ggplot2::geom_point() +
                                          ggplot2::theme(legend.position = 'top') +
                                          ggplot2::labs(
                                            x = paste(input$setting, 1, sep = '-'),
                                            y = paste(input$setting, 2, sep = '-')
                                          )
                                      )

                                      output$plot3_subpop <-
                                        renderPlot({
                                          data_subpop$plot1()
                                        },
                                        height = function() {
                                          session$clientData[[paste0('output_', ID, '-plot3_subpop_width')]]*0.8
                                        } )
                                      output$results_subpop = shiny::renderText(
                                        '<h2><b><p style="color:red", align="center">Results</p></b></h2>'
                                      )
                                      data_subpop$found_matches = T

                                    } else{
                                      output$results_subpop = renderText('<b> It was not possible to match G1/G2- and S-phase groups. </b>')
                                      data_subpop$found_matches = F
                                    }
                                  }
                                  shinyjs::enable('SaveGroups_subpop')
                                  shinyjs::enable('save_subpop')

                                } else{
                                  output$SubPopUI2 <- shiny::renderUI({
                                    NULL
                                  })
                                }
                              })
                              }

                            shiny::observeEvent(data_subpop$summary, {
                              if (all(c('S', 'G1/G2') %in% input$phase) &
                                  nrow(data_subpop$summary) > 3 &
                                  sum(data_subpop$summary$Phase == 'S') == sum(data_subpop$summary$Phase == 'G1/G2') &
                                  input$what == 'scCN') {
                                shinyWidgets::updateSwitchInput(inputId = 'AssignGroups_subpop', disabled = F)
                              } else{
                                shinyWidgets::updateSwitchInput(inputId = 'AssignGroups_subpop', disabled = T)
                              }

                            })

                          } else{
                            data_subpop = shiny::reactiveValues(
                              data=dplyr::tibble(Phase = character()),
                              counter = 1,
                              G1 = NULL,
                              S = NULL,
                              subpopulation = NULL,
                              final = dplyr::tibble(a = 1),
                              click = 1,
                              dbclick = 1,
                              poligon = dplyr::tibble(),
                              summary = dplyr::tibble(Phase = character())
                            )

                            shinyWidgets::updateSwitchInput(inputId = 'AssignGroups_subpop',
                                                            value = F,
                                                            disabled = T)

                            output$group_phase_subpop <-
                              shiny::renderTable(NULL)

                            output$SubPopUI <- shiny::renderUI({
                              NULL
                            })
                            output$SubPopUI2 <-
                              shiny::renderUI({
                                NULL
                              })
                            output

                          }
                        })

                      })

}
