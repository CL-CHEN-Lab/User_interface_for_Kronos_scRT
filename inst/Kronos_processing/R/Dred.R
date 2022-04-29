Dim_red_ui = function(id) {
  ns <- shiny::NS(paste0('Dred', id))

  shiny::fluidPage(
    shiny::fluidRow(
      align = 'center',
      shiny::column(
        width = 3,
        offset = 3,
        shiny::radioButtons(
          inputId = ns('setting'),
          label = 'Chose dimentsionality reduction',
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
      )
    ),
    shiny::fluidRow(
      shiny::column(width = 10,
                    shiny::plotOutput(ns('dred_out'))),
      shiny::column(
        width = 2,
        shiny::radioButtons(
          inputId = ns('color'),
          label = 'Color',
          choices = c('Group', 'Sample', '%Rep'),
          width = '100%'
        ),
        shiny::radioButtons(
          inputId = ns('shape'),
          label = 'Shape',
          choices = c('Group', 'Sample'),
          width = '100%'
        )
      )

    )
  )

}

Dim_red_server = function(id, scCN, out, cores = 3) {
  shiny::moduleServer(paste0('Dred', id),
                      function(input,
                               output,
                               session,
                               scCNV = scCN,
                               Cores = cores,
                               Out = out) {

                        #load required operators
                        `%>%` = tidyr::`%>%`
                        #initialise
                        shinyjs::disable('save')
                        shinyjs::disable('run')

                        Data_to_plot = reactiveValues(scCNV = tibble(), plot = ggplot2::ggplot())

                        #select colums of interest and reformat
                        scCNV = scCNV %>%
                          dplyr::mutate(pos = paste0(chr, ':', start, '-', end)) %>%
                          dplyr::select('pos',
                                        'Cell',
                                        'PercentageReplication',
                                        'basename',
                                        'group',
                                        'Rep') %>%
                          dplyr::mutate(Rep = as.numeric(Rep)) %>%
                          tidyr::spread(pos, Rep)
                        #remove col with na
                        scCNV = scCNV[, colSums(is.na(scCNV)) == 0]



                        #create matrix data and save extra info
                        mat = scCNV[, -c(1:4)]
                        scCNV = scCNV[, c(1:4)]

                        #calculate simple matching coef
                        results = as.matrix(ade4::dist.binary(df = mat, method = 2))
                        rm('mat')

                        #activate run
                        shinyjs::enable('run')

                        shiny::observeEvent(input$random, {
                          shiny::updateNumericInput(inputId = 'seed',
                                                    value = round(stats::runif(1), 9) * 10 ^ 9)
                        })

                        shiny::observeEvent(input$run, {

                          #if random set new seed
                          shiny::req(input$seed)
                          if (input$run > 0) {

                            #disable whilre running
                            shinyjs::disable('run')

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
                                  is_distance = T,
                                  verbose = F,
                                  max_iter = 5000,
                                  num_threads = cores,
                                  partial_pca = T
                                )

                              Data_to_plot$scCNV = scCNV %>% dplyr::mutate(x = tsne$Y[, 1],
                                                                           y = tsne$Y[, 2])



                            } else if (input$setting == 'UMAP') {
                              umap <-
                                umap::umap(d = results,
                                           input = 'dist',
                                           random_state = input$seed)
                              Data_to_plot$scCNV = scCNV %>% dplyr::mutate(x =
                                                                             umap$layout[, 1],
                                                                           y = umap$layout[, 2])
                            }

                            #activate run and save
                            shinyjs::enable('save')
                            shinyjs::enable('run')
                          }

                          shiny::observeEvent(c(Data_to_plot$scCNV, input$color, input$shape), {
                            #plot
                            Data_to_plot$plot = Data_to_plot$scCNV %>%
                              dplyr::select(
                                'x',
                                'y',
                                'color' = dplyr::case_when(
                                  input$color == 'Group' ~ 'group',
                                  input$color == 'Sample' ~ 'basename',
                                  input$color == '%Rep' ~ 'PercentageReplication'
                                ),
                                'shape' = dplyr::case_when(
                                  input$shape == 'Group' ~ 'group',
                                  input$shape == 'Sample' ~ 'basename'
                                )
                              ) %>%
                              ggplot2::ggplot() +
                              ggplot2::geom_point(ggplot2::aes(x, y, color = color, shape = shape),
                                                  alpha = 0.4,
                                                  size = 2) +
                              ggplot2::labs(
                                x = paste(input$setting, '1', sep = '-'),
                                y = paste(input$setting, '2', sep = '-'),
                                color = input$color,
                                shape = input$shape
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
                            }

                            output$dred_out = shiny::renderPlot(Data_to_plot$plot)
                          })
                        })

                        shiny::observeEvent(input$save, {
                          if (!dir.exists(file.path(Out, input$setting))) {
                            dir.create(file.path(Out, input$setting),recursive = T)
                          }

                          ggplot2::ggsave(
                            plot = Data_to_plot$plot,
                            filename = file.path(
                              Out,
                              input$setting,
                              paste0(
                                basename(Out),
                                '_',
                                input$setting,
                                '_seed_',
                                input$seed,
                                '_color_',
                                ifelse(input$color == '%Rep', 'RepPerc', input$color),
                                '_shape_',
                                input$shape,
                                '.pdf'
                              )
                            ),
                            device = grDevices::cairo_pdf
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
                              '.tsv'
                            )
                          ))) {
                            readr::write_tsv(
                              x = Data_to_plot$scCNV,
                              file = file.path(
                                Out,
                                input$setting,
                                paste0(
                                  basename(Out),
                                  '_',
                                  input$setting,
                                  '_seed_',
                                  input$seed,
                                  '.tsv'
                                )
                              )
                            )
                          }

                        })

                      })
}
