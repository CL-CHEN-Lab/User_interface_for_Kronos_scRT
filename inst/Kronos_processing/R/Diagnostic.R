diagnostic_ui <- function(id) {
  ns <- shiny::NS(id)

  shinydashboard::box(
    id = ns('box'),
    status = 'primary',
    title = id,
    width = 12,
    solidHeader = T,
    collapsible = T,
    shinyjs::useShinyjs(),
    shinydashboard::tabBox(
      id = ns('Diagnostic_tab'),
      width = 12,
      shiny::tabPanel(
        title = 'Setting RPMH and define G1/G2- S-phase populations',
        align = 'center',
        shiny::fluidRow(
          shiny::column(width = 6,
                        plotOutput(ns(
                          'plot__diagnostic'
                        ))),
          shiny::column(width = 6,
                        plotOutput(
                          ns('plot2__diagnostic'),
                          brush = shiny::brushOpts(
                            id = ns("plot2_brush__diagnostic"),
                            direction = 'y'
                          )
                        ))
        ),
        shiny::fluidRow(
          shiny::sliderInput(
            inputId = ns('min_n_reads__diagnostic'),
            label = 'Reads per Mb per Haplotype',
            min = 50,
            max = 300,
            value = 160,
            width = '100%'

          ),
          shiny::tags$div(
            class = "header",
            checked = NA,
            tags$h4(
              "To modify the S- and G1/G2- phase definition draw a line or a box between the two on the right plot."
            )
          ),
          radioButtons(
            inputId = ns('Auto__diagnostic'),
            choices = c('Auto', 'Manual'),
            label = 'S-phase correction',
            inline = T,
            selected = 'Auto',
            width = '25%'
          ),
          shiny::actionButton(
            inputId = ns('CleanStaging__diagnostic'),
            label = 'Clean',
            width = '25%'
          ),
          shiny::actionButton(
            inputId = ns('Apply__diagnostic'),
            label = 'Apply',
            width = '25%'
          )
        )
      ),

      shiny::tabPanel(
        title = 'Fix S-phase progression',
        align = 'center',
        shiny::fluidRow(
          shiny::column(width = 6,
                        plotOutput(ns(
                          'plot3__diagnostic'
                        ))),
          shiny::column(width = 6,
                        plotOutput(ns(
                          'plot4__diagnostic'
                        )))
        ),
        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::htmlOutput(ns('parameters_text__diagnostic')),
            shiny::sliderInput(
              inputId = ns('First_p_S__diagnostic'),
              label = 'Parameter fist part of the S phase',
              min = 0.3,
              max = 1.7,
              step = 0.001,
              value = 1,
              width = '100%'
            ),
            shiny::sliderInput(
              inputId = ns('Second_p_S__diagnostic'),
              label = 'Parameter second part of the S phase',
              min = 0.3,
              max = 1.7,
              step = 0.001,
              value = 1,
              width = '100%'
            ),
            shiny::actionButton(
              inputId = ns('Reset_Correction__diagnostic'),
              'Reset',
              width = '50%'
            ),
            shiny::actionButton(
              inputId = ns('Apply_Correction__diagnostic'),
              'Apply correction',
              width = '50%'
            )
          ),
          shiny::column(width = 6,
                        plotOutput(ns(
                          'plot5__diagnostic'
                        )))
        )
      )
    )
  )
}

diagnostic_server <-
  function(id,
           PerCell,
           Sample,
           Setting = tibble(),
           Cores = 3,
           OutFolder = '~/') {
    shiny::moduleServer(id,
                        function(input,
                                 output,
                                 session,
                                 sample = id,
                                 PC_df = PerCell,
                                 S_df = Setting,
                                 cores = Cores,
                                 out = OutFolder) {
                          #load required operators
                          `%>%` = tidyr::`%>%`
                          `%dopar%` = foreach::`%dopar%`
                          `%do%` = foreach::`%do%`

                          #initialize
                          R_df = shiny::reactiveValues(result = dplyr::tibble(),
                                                       PC_df = PC_df)
                          plots = shiny::reactiveValues()

                          #if setting file is provided change df
                          if (nrow(S_df) == 1) {
                            if (!(is.na(S_df$threshold_Sphase) |
                                  is.na(S_df$threshold_G1G2phase))) {
                              R_df$PC_df = R_df$PC_df %>%
                                dplyr::filter(
                                  normalized_dimapd > S_df$threshold_Sphase |
                                    normalized_dimapd < S_df$threshold_G1G2phase
                                ) %>%
                                dplyr::mutate(
                                  is_high_dimapd = dplyr::case_when(
                                    normalized_dimapd > S_df$threshold_Sphase ~ T,
                                    normalized_dimapd < S_df$threshold_G1G2phase ~ F,
                                    T ~ is_high_dimapd
                                  ),
                                  is_noisy = dplyr::case_when(
                                    normalized_dimapd > S_df$threshold_Sphase ~ T,
                                    normalized_dimapd < S_df$threshold_G1G2phase ~ F,
                                    T ~ is_noisy
                                  )
                                )
                            }

                          }
                          shiny::observe({
                            #calculate median ploidy of the suppose G1
                            median_ploidy_not_noisy__diagnostic =  R_df$PC_df %>% dplyr::filter(is_noisy == F) %>%
                              dplyr::pull(mean_ploidy) %>% stats::median()


                            #disable buttons
                            shinyjs::disable('First_p_S__diagnostic')
                            shinyjs::disable('Second_p_S__diagnostic')
                            shinyjs::disable('Reset_Correction__diagnostic')
                            shinyjs::disable('Apply_Correction__diagnostic')

                            if (is.numeric(input$min_n_reads__diagnostic)) {
                              #if setting file is provided update widgets
                              if (nrow(S_df) == 1) {
                                shiny::updateSliderInput(
                                  session = session,
                                  inputId = 'min_n_reads__diagnostic',
                                  value = S_df$RPMPH_TH
                                )
                                shiny::updateSliderInput(
                                  session = session,
                                  inputId = 'First_p_S__diagnostic',
                                  value = S_df$Sphase_first_part
                                )
                                shiny::updateSliderInput(
                                  session = session,
                                  inputId = 'Second_p_S__diagnostic',
                                  value = S_df$Sphase_second_part
                                )
                                shiny::updateRadioButtons(session = session,
                                                          inputId = 'Auto__diagnostic',
                                                          selected = 'Manual')
                              }

                              #prepare df
                              R_df$PC_df = R_df$PC_df %>%
                                dplyr::mutate(
                                  Type = dplyr::case_when(
                                    as.logical(is_high_dimapd) == T &
                                      as.logical(is_noisy) == T ~ 'S-phase cells',
                                    as.logical(is_high_dimapd) == F &
                                      as.logical(is_noisy) == T ~ 'unknown cells',
                                    T ~ 'G1/G2-phase cells'
                                  ),
                                  Type = dplyr::case_when(
                                    coverage_per_1Mbp < input$min_n_reads__diagnostic * median_ploidy_not_noisy__diagnostic ~ 'Low Coverage',
                                    ploidy_confidence < 2 &
                                      ploidy_confidence != -100 ~ 'Low ploidy confidence',
                                    mean_ploidy < median_ploidy_not_noisy__diagnostic / 1.5 ~ 'Too low ploidy compared to G1/G2-phase pool',
                                    mean_ploidy > median_ploidy_not_noisy__diagnostic * 2 ~ 'Too high ploidy compared to G1/G2-phase pool',
                                    T ~ Type
                                  )
                                )

                              #select potentially suitable cells
                              R_df$PC_df_subset = R_df$PC_df %>%
                                dplyr::filter(
                                  Type %in% c(
                                    'S-phase cells',
                                    'G1/G2-phase cells',
                                    'Too low ploidy compared to G1/G2-phase pool',
                                    'Too high ploidy compared to G1/G2-phase pool'
                                  )
                                ) %>%
                                dplyr::mutate(Type_mem = Type)
                            }
                            # first plot
                            plots$plot__diagnostic = shiny::reactive({
                              R_df$PC_df %>%
                                ggplot2::ggplot(ggplot2::aes(mean_ploidy,
                                                             normalized_dimapd,
                                                             color = Type)) +
                                ggplot2::geom_point(alpha = 0.3) +
                                ggplot2::scale_color_manual(
                                  values = c(
                                    'Low Coverage' = "#ff7949",
                                    'Low ploidy confidence' = "#70001e",
                                    'Too low ploidy compared to G1/G2-phase pool' = "#01e7ab",
                                    'Too high ploidy compared to G1/G2-phase pool' = "#a7001b",
                                    'G1/G2-phase cells' = "#005095",
                                    'S-phase cells' = "#78bd3e",
                                    'unknown cells' = "#dfbd31"
                                  )
                                ) +
                                ggplot2::theme(legend.position = 'top',
                                               legend.title = ggplot2::element_blank()) +
                                ggplot2::xlab('Ploidy') + ggplot2::ylab('Variability')
                            })
                            output$plot__diagnostic <-
                              shiny::renderPlot({
                                plots$plot__diagnostic()
                              })
                          })
                          #calculate new median ploidy after brush selection
                          shiny::observeEvent(input$plot2_brush__diagnostic, {
                            if (!(
                              is.null(input$plot2_brush__diagnostic) |
                              is.numeric(input$plot2_brush__diagnostic)
                            )) {
                              R_df$median_ploidy_not_noisy_2__diagnostic = stats::median(
                                R_df$PC_df_subset %>% dplyr::filter(Type ==
                                                                      'G1/G2-phase cells') %>% dplyr::pull(mean_ploidy)
                              )

                              #reassign cells based on new median ploidy
                              R_df$PC_df_subset = R_df$PC_df_subset %>%
                                dplyr::mutate(
                                  Type = dplyr::case_when(
                                    mean_ploidy < R_df$median_ploidy_not_noisy_2__diagnostic / 1.5 ~ 'Too low ploidy compared to G1/G2-phase pool',
                                    mean_ploidy > R_df$median_ploidy_not_noisy_2__diagnostic * 2 ~ 'Too high ploidy compared to G1/G2-phase pool',
                                    normalized_dimapd > input$plot2_brush__diagnostic$ymax ~ 'S-phase cells',
                                    normalized_dimapd < input$plot2_brush__diagnostic$ymin & !is_noisy ~
                                      'G1/G2-phase cells',
                                    T ~ 'unknown cells'
                                  )
                                )

                              #save thresholds
                              S_th = input$plot2_brush__diagnostic$ymax
                              G1G2_th = input$plot2_brush__diagnostic$ymin
                            }

                          })
                          #reset S/G1G2 th
                          shiny::observeEvent(input$CleanStaging__diagnostic, {
                            R_df$PC_df_subset = R_df$PC_df_subset %>%
                              mutate(Type = Type_mem)


                          })

                          plots$plot2__diagnostic = shiny::reactive({
                            R_df$PC_df_subset %>%
                              ggplot2::ggplot(ggplot2::aes(mean_ploidy,
                                                           normalized_dimapd,
                                                           color = Type)) +
                              ggplot2::geom_point(alpha = 0.3) +
                              ggplot2::scale_color_manual(
                                values = c(
                                  'G1/G2-phase cells' = "#005095",
                                  'S-phase cells' = "#78bd3e",
                                  'unknown cells' = "#dfbd31",
                                  'Too low ploidy compared to G1/G2-phase pool' = "#01e7ab",
                                  'Too high ploidy compared to G1/G2-phase pool' = "#a7001b"
                                )
                              ) +
                              ggplot2::theme(legend.position = 'top',
                                             legend.title = ggplot2::element_blank()) +
                              ggplot2::xlab('Ploidy') + ggplot2::ylab('Variability')
                          })

                          output$plot2__diagnostic <-
                            shiny::renderPlot({
                              plots$plot2__diagnostic()
                            })


                          shiny::observeEvent(input$Apply__diagnostic, {
                            shinyjs::disable('min_n_reads__diagnostic')
                            shinyjs::disable('Auto__diagnostic')
                            shinyjs::disable('CleanStaging__diagnostic')
                            shinyjs::disable('Apply__diagnostic')
                            shinyjs::enable('First_p_S__diagnostic')
                            shinyjs::enable('Second_p_S__diagnostic')
                            shinyjs::enable('Reset_Correction__diagnostic')
                            shinyjs::enable('Apply_Correction__diagnostic')


                            #start spinner
                            shinybusy::show_spinner()
                            #adjust S-phase
                            R_df$PC_df_subset_2 = R_df$PC_df_subset %>%
                              dplyr::mutate(
                                is_high_dimapd = ifelse(Type == 'S-phase cells', T, F),
                                is_noisy = ifelse(is_high_dimapd, T, is_noisy)
                              ) %>%
                              dplyr::filter(Type %in% c('G1/G2-phase cells', 'S-phase cells'))

                            R_df$median_ploidy_not_noisy_3__diagnostic = stats::median(
                              R_df$PC_df_subset_2 %>% dplyr::filter(Type ==
                                                                      'G1/G2-phase cells') %>% dplyr::pull(mean_ploidy)
                            )

                            plots$plot3__diagnostic <-
                              shiny::reactive({
                                R_df$PC_df_subset_2 %>%
                                  ggplot2::ggplot(ggplot2::aes(mean_ploidy, normalized_dimapd, color = Type)) +
                                  ggplot2::geom_point(alpha = 0.3) +
                                  ggplot2::scale_color_manual(
                                    values = c(
                                      'Low Coverage' = "#ff7949",
                                      'Low ploidy confidence' = "#70001e",
                                      'Too low ploidy compared to G1/G2' = "#01e7ab",
                                      'Too high ploidy compared to G1/G2' = "#a7001b",
                                      'G1/G2-phase cells' = "#005095",
                                      'S-phase cells' = "#78bd3e",
                                      'unknown cells' = "#dfbd31"
                                    )
                                  ) +
                                  ggplot2::geom_vline(xintercept = R_df$median_ploidy_not_noisy_3__diagnostic) +
                                  ggplot2::theme(legend.position = 'top',
                                                 legend.title = ggplot2::element_blank()) +
                                  ggplot2::xlab('Ploidy') + ggplot2::ylab('Variability')
                              })

                            output$plot3__diagnostic <-
                              shiny::renderPlot({
                                plots$plot3__diagnostic()
                              })

                            if (input$Auto__diagnostic == 'Auto') {
                              shiny::isolate({
                                isolated_PC_df_subset_2 <- R_df$PC_df_subset_2
                                isolated_median_ploidy_not_noisy_3__diagnostic =
                                  R_df$median_ploidy_not_noisy_3__diagnostic
                              })


                              # correct mean ploidy
                              cl = snow::makeCluster(cores)
                              doSNOW::registerDoSNOW(cl)

                              # test multiple parameters to correct S phase
                              R_df$distributions = foreach::foreach(
                                a = seq(0.95, 1, by = 0.001),
                                .combine = 'rbind',
                                .packages = c('tidyverse',
                                              'LaplacesDemon',
                                              'foreach',
                                              'shiny'),
                                .export = c("session")
                              ) %dopar% {
                                dist = foreach::foreach(
                                  b = seq(0.5, 0.55, by = 0.001),
                                  .combine = 'rbind',
                                  .packages = c('tidyverse',
                                                'LaplacesDemon',
                                                'shiny'),
                                  .export = c("session")
                                ) %do% {
                                  #adjust S-phase based on a and b
                                  x = isolated_PC_df_subset_2 %>%
                                    dplyr::filter(Type == 'S-phase cells') %>%
                                    dplyr::mutate(
                                      corrected_mean_ploidy = ifelse(
                                        mean_ploidy >= isolated_median_ploidy_not_noisy_3__diagnostic,
                                        mean_ploidy / a,
                                        mean_ploidy / b
                                      )
                                    ) %>%
                                    dplyr::select(corrected_mean_ploidy) %>%
                                    dplyr::pull()

                                  # are the data unimodal?
                                  if (LaplacesDemon::is.unimodal(x)) {
                                    # d is the distance between theoretical center of the S phase (G1 median ploidy *1.5)
                                    #and the average ploidy of the corrected S-phase
                                    # d is divided by sd(x) in order to select parameters that keep the distribution as wide as possible
                                    dplyr::tibble(
                                      A = a,
                                      B = b,
                                      d = 1 / stats::sd(x),
                                      unimodal = T
                                    )

                                  } else{
                                    dplyr::tibble(
                                      A = a,
                                      B = b,
                                      d = stats::sd(x),
                                      unimodal = F
                                    )
                                  }

                                }
                                dist
                              }

                              snow::stopCluster(cl)
                            } else{
                              R_df$distributions = dplyr::tibble()
                            }


                            #select the minimum value of d
                            if (nrow(R_df$distributions) == 0) {
                              output$parameters_text__diagnosc = renderText('Please provide manual ones')
                              R_df$PC_df_corrected = R_df$PC_df_subset_2 %>% dplyr::mutate(
                                mean_ploidy_corrected = mean_ploidy,
                                Type = dplyr::case_when(
                                  Type == 'S-phase cells' &
                                    mean_ploidy < R_df$median_ploidy_not_noisy_3__diagnostic ~ 'Second-part-S-phase cells',
                                  Type == 'S-phase cells' &
                                    mean_ploidy > R_df$median_ploidy_not_noisy_3__diagnostic ~ 'First-part-S-phase cells',
                                  T ~ Type
                                ),
                                mean_ploidy_corrected_bu = mean_ploidy_corrected
                              )
                            } else{
                              R_df$distributions = R_df$distributions %>%
                                dplyr::filter(unimodal == ifelse(any(unimodal == T), T, F)) %>%
                                dplyr::filter(d == min(d))

                              if (nrow(R_df$distributions) > 1) {
                                R_df$PC_df_corrected  = R_df$PC_df_subset_2 %>%
                                  dplyr::mutate(
                                    mean_ploidy_corrected = mean_ploidy,
                                    Type = dplyr::case_when(
                                      Type == 'S-phase cells' &
                                        mean_ploidy < R_df$median_ploidy_not_noisy_3__diagnostic ~ 'Second-part-S-phase cells',
                                      Type == 'S-phase cells' &
                                        mean_ploidy > R_df$median_ploidy_not_noisy_3__diagnostic ~ 'First-part-S-phase cells',
                                      T ~ Type
                                    ),
                                    mean_ploidy_corrected_bu = mean_ploidy_corrected
                                  )
                                output$parameters_text__diagnostic = renderText(
                                  'S phase correction parameters could not be established, please provide manual ones.'
                                )

                              } else{
                                output$parameters_text__diagnostic = renderText(
                                  paste(
                                    'S phase correction parameters have been estimated.\n',
                                    'Parameter first part S phase: ',
                                    R_df$distributions$A,
                                    '\n',
                                    'Parameter Second part S phase: ',
                                    R_df$distributions$B
                                  )
                                )
                                #save parameters
                                shiny::updateSliderInput(
                                  session,
                                  "First_p_S__diagnostic",
                                  value = R_df$distributions$A,
                                  min = 0.3,
                                  max = 1.7,
                                  step = 0.001
                                )

                                shiny::updateSliderInput(
                                  session,
                                  "Second_p_S__diagnostic",
                                  value = R_df$distributions$B,
                                  min = 0.3,
                                  max = 1.7,
                                  step = 0.001
                                )

                                R_df$firstpart_S = R_df$distributions$A
                                R_df$secondpart_S = R_df$distributions$B

                                #correct
                                R_df$PC_df_corrected = R_df$PC_df_subset_2 %>%
                                  dplyr::mutate(
                                    mean_ploidy_corrected = dplyr::case_when(
                                      Type == 'S-phase cells' &
                                        mean_ploidy < R_df$median_ploidy_not_noisy_3__diagnostic ~ mean_ploidy / R_df$distributions$B,
                                      Type == 'S-phase cells' &
                                        mean_ploidy > R_df$median_ploidy_not_noisy_3__diagnostic ~ mean_ploidy /  R_df$distributions$A,
                                      T ~ mean_ploidy
                                    ),
                                    Type = dplyr::case_when(
                                      Type == 'S-phase cells' &
                                        mean_ploidy < R_df$median_ploidy_not_noisy_3__diagnostic ~ 'Second-part-S-phase cells',
                                      Type == 'S-phase cells' &
                                        mean_ploidy > R_df$median_ploidy_not_noisy_3__diagnostic ~ 'First-part-S-phase cells',
                                      T ~ Type

                                    ),
                                    mean_ploidy_corrected_bu = mean_ploidy_corrected
                                  )


                              }
                            }
                            #stop spinner
                            shinybusy::hide_spinner()

                            #change page
                            shiny::updateTabsetPanel(session = session,
                                                     "Diagnostic_tab",
                                                     selected = 'Fix S-phase progression')


                            #sliders
                            # manual correction
                            shiny::observeEvent(input$First_p_S__diagnostic,
                                                {
                                                  R_df$PC_df_corrected = R_df$PC_df_corrected %>% dplyr::mutate(
                                                    mean_ploidy_corrected = dplyr::case_when(
                                                      Type == 'First-part-S-phase cells' ~ mean_ploidy / input$First_p_S__diagnostic,
                                                      T ~ mean_ploidy_corrected
                                                    )

                                                  )

                                                  R_df$firstpart_S = input$First_p_S__diagnostic

                                                })
                            shiny::observeEvent(input$Second_p_S__diagnostic,
                                                {
                                                  R_df$PC_df_corrected = R_df$PC_df_corrected %>% dplyr::mutate(
                                                    mean_ploidy_corrected = dplyr::case_when(
                                                      Type == 'Second-part-S-phase cells' ~ mean_ploidy / as.numeric(input$Second_p_S__diagnostic),
                                                      T ~ mean_ploidy_corrected
                                                    )

                                                  )

                                                  R_df$secondpart_S = input$Second_p_S__diagnostic

                                                })

                            #plot corrected data
                            plots$plot4__diagnostic <-
                              shiny::reactive({
                                R_df$PC_df_corrected %>%
                                  ggplot2::ggplot(ggplot2::aes(mean_ploidy_corrected,
                                                               normalized_dimapd,
                                                               color = Type)) +
                                  ggplot2::geom_point(alpha = 0.3) +
                                  ggplot2::scale_color_manual(
                                    values = c(
                                      'G1/G2-phase cells' = '#005095',
                                      'First-part-S-phase cells' = '#78bd3e',
                                      'Second-part-S-phase cells' = '#83007e',
                                      'unknown cells' = '#dfbd31'
                                    )
                                  ) +
                                  ggplot2::geom_vline(xintercept = R_df$median_ploidy_not_noisy_3__diagnostic) +
                                  ggplot2::theme(legend.position = 'top',
                                                 legend.title = ggplot2::element_blank()) +
                                  ggplot2::xlab('Ploidy') + ggplot2::ylab('Variability')
                              })
                            output$plot4__diagnostic <-
                              shiny::renderPlot({
                                plots$plot4__diagnostic()
                              })

                            #plot density
                            plots$plot5__diagnostic <-
                              shiny::reactive({
                                R_df$PC_df_corrected %>% dplyr::filter(Type %in% c(
                                  'Second-part-S-phase cells',
                                  'First-part-S-phase cells'
                                )) %>%
                                  ggplot2::ggplot() +
                                  ggplot2::geom_density(ggplot2::aes(x = mean_ploidy_corrected, y = (..density..)),
                                                        color = "black") +
                                  ggplot2::xlab('Ploidy') + ggplot2::ylab('Variability')
                              })
                            output$plot5__diagnostic <-
                              shiny::renderPlot({
                                plots$plot5__diagnostic()
                              })


                          })

                          shiny::observeEvent(input$Reset_Correction__diagnostic,
                                              {
                                                R_df$PC_df_corrected = R_df$PC_df_corrected %>%
                                                  dplyr::mutate(mean_ploidy_corrected = mean_ploidy_corrected_bu)

                                                if (nrow(S_df) == 1) {
                                                  R_df$firstpart_S = S_df$Sphase_first_part
                                                  R_df$secondpart_S = S_df$Sphase_second_part
                                                } else if (nrow(R_df$distributions) == 1) {
                                                  R_df$firstpart_S = R_df$distributions$A
                                                  R_df$secondpart_S =  R_df$distributions$B
                                                } else{
                                                  R_df$firstpart_S = 1
                                                  R_df$secondpart_S = 1
                                                }


                                                shiny::updateSliderInput(
                                                  session,
                                                  "First_p_S__diagnostic",
                                                  value = R_df$firstpart_S,
                                                  min = 0.3,
                                                  max = 1.7,
                                                  step = 0.001
                                                )

                                                shiny::updateSliderInput(
                                                  session,
                                                  "Second_p_S__diagnostic",
                                                  value = R_df$secondpart_S,
                                                  min = 0.3,
                                                  max = 1.7,
                                                  step = 0.001
                                                )


                                              })

                          shiny::observeEvent(input$Apply_Correction__diagnostic, {
                            #disable buttons
                            shinyjs::disable('First_p_S__diagnostic')
                            shinyjs::disable('Second_p_S__diagnostic')
                            shinyjs::disable('Reset_Correction__diagnostic')
                            shinyjs::disable('Apply_Correction__diagnostic')
                            #save setting file

                            if (nrow(S_df) == 1) {
                              S_th = S_df$threshold_Sphase
                              G_th = S_df$threshold_G1G2phase
                            } else if (!is.numeric(input$plot2_brush__diagnostic)) {
                              S_th = input$plot2_brush__diagnostic$ymax
                              G_th = input$plot2_brush__diagnostic$ymin
                            } else{
                              S_th = NULL
                              G_th = NULL
                            }

                            R_df$setting = dplyr::tibble(
                              threshold_Sphase = ifelse(is.null(S_th), NA, S_th),
                              threshold_G1G2phase = ifelse(is.null(G_th), NA, G_th),
                              Sphase_first_part = R_df$firstpart_S,
                              Sphase_second_part = R_df$secondpart_S,
                              Ploidy = R_df$median_ploidy_not_noisy_3__diagnostic,
                              RPMPH_TH = input$min_n_reads__diagnostic,
                              RPM_TH = round(
                                input$min_n_reads__diagnostic * R_df$median_ploidy_not_noisy_3__diagnostic
                              ),
                              basename = unique(R_df$PC_df_corrected$basename),
                              group = unique(R_df$PC_df_corrected$group)
                            )

                            R_df$setting %>%
                              readr::write_tsv(file = file.path(out, paste0(sample, '_settings.txt')))

                            #save plots

                            ggplot2::ggsave(
                              filename = file.path(out,
                                                   paste0(sample,
                                                          '_all_cells.pdf')),
                              plot = plots$plot__diagnostic(),
                              device = grDevices::cairo_pdf
                            )
                            ggplot2::ggsave(
                              filename = file.path(out,
                                                   paste0(sample,
                                                          '_first_filtering.pdf')),
                              plot = plots$plot2__diagnostic(),
                              device = grDevices::cairo_pdf
                            )
                            ggplot2::ggsave(
                              filename = file.path(out,
                                                   paste0(
                                                     sample,
                                                     '_selected_G1_S_cells.pdf'
                                                   )),
                              plot = plots$plot3__diagnostic(),
                              device = grDevices::cairo_pdf
                            )
                            ggplot2::ggsave(
                              filename = file.path(out,
                                                   paste0(sample,
                                                          '_fixed_S_phase.pdf')),
                              plot = plots$plot4__diagnostic(),
                              device = grDevices::cairo_pdf
                            )
                            ggplot2::ggsave(
                              filename = file.path(out,
                                                   paste0(
                                                     sample,
                                                     '_S_phase_cell_distribution.pdf'
                                                   )),
                              plot = plots$plot5__diagnostic(),
                              device = grDevices::cairo_pdf
                            )
                          })
                        })

  }
