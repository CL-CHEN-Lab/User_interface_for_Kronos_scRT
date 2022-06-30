scPlots_ui <- function(id,title=NULL) {
  ns <- shiny::NS(id)
  shinydashboard::box(
    id = ns('box'),
    status = 'primary',
    title = title,
    width = 12,
    solidHeader = T,
    collapsible = T,
    shiny::fluidRow(plotOutput(ns('plot__scPlot'), width = '100%',height = '100%'))
  )
}


scPlots_server <-
  function(id,
           RTs,
           scCN,
           filling,
           Extreme_values,
           out,
           colors,
           save) {
    shiny::moduleServer(id,
                        function(input,
                                 output,
                                 session,
                                 scRT = RTs,
                                 scTracks = scCN,
                                 what = filling,
                                 extremes = Extreme_values,
                                 Out = out,
                                 Save = save,
                                 col=colors,
                                 ID=id) {
                          #load required operators
                          `%>%` = tidyr::`%>%`

                          #recover colors
                          Colors=col$color
                          names(Colors)=col$basename


                          if (nrow(scTracks) != 0) {
                            #depending on selection change aesthetics of the plot
                            if (what == 'scRT') {
                              scTracks = scTracks %>%
                                dplyr::rename(Value = Rep) %>%
                                dplyr::mutate(Value = ifelse(Value, "Replicated", "Unreplicated"))

                              ggEXTRA = ggplot2::ggplot() +
                                ggplot2::scale_fill_manual(values = c("Replicated" = '#a7001b',
                                                                      "Unreplicated" = '#005095')) +
                                ggplot2::labs(fill = 'State', color = 'Sample')

                            } else if (what == 'scCN') {
                              scTracks = scTracks %>%
                                dplyr::rename(Value = CN)

                              ggEXTRA = ggplot2::ggplot() +
                                ggplot2::scale_fill_gradient(
                                  low = '#c56700',
                                  high = '#6d0042',
                                  limits = c(extremes$CN[1], extremes$CN[2]),
                                  oob = scales::squish
                                ) +
                                ggplot2::labs(fill = 'CNV ', color = 'Sample')

                            } else if (what == 'Norm. scCN') {
                              scTracks = scTracks %>%
                                dplyr::rename(Value = CN_bg)

                              ggEXTRA = ggplot2::ggplot() +
                                ggplot2::scale_fill_gradient(
                                  low = '#00c5a4',
                                  high = '#e58225',
                                  limits = c(extremes$CN_bg[1], extremes$CN_bg[2]),
                                  oob = scales::squish
                                ) +
                                ggplot2::labs(fill = expression(over(S[CNV[i]],
                                                                     bar(G1 / G2[CNV]))),
                                              color = 'Sample')

                            }

                            Chr = unique(scRT$chr)
                            Maxi = max(scTracks$newIndex)
                            n = stringr::str_count(Maxi)

                            plot =
                              ggEXTRA +
                              ggplot2::geom_path(
                                data = scRT %>%
                                  dplyr::mutate(
                                    mid = (start + end) / 2,
                                    RT = RT * Maxi / 6 + Maxi / 20
                                  ) %>%
                                  tidyr::gather(what, pos, start, end) %>%
                                  dplyr::arrange(mid, pos),
                                ggplot2::aes(pos, RT, color = basename)
                              ) +
                              ggplot2::annotate(
                                "rect",
                                xmin = -Inf,
                                xmax = Inf,
                                ymin = Maxi / 100,
                                ymax = Maxi / 20 - Maxi / 100,
                                fill = 'white'
                              ) +
                              ggplot2::geom_rect(
                                data = scTracks ,
                                ggplot2::aes(
                                  xmin = start,
                                  xmax = end,
                                  ymin = -newIndex,
                                  ymax = -newIndex - 1,
                                  fill = Value
                                )
                              ) +
                              ggplot2::scale_y_continuous(
                                breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                                labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                                name = 'RT',
                                sec.axis = ggplot2::sec_axis(
                                  ~ .,
                                  breaks = c(-seq(1, round(
                                    Maxi / 10 ^ (n - 1)
                                  ), 1) * 10 ^ (n - 1)) - 0.5,
                                  labels = as.character(c(seq(
                                    1, round(Maxi / 10 ^ (n - 1)), 1
                                  ) * 10 ^ (n - 1))),
                                  name = 'Single Cell tracks ordered by S-phase progression'
                                )
                              ) +
                              ggplot2::scale_x_continuous(
                                labels = function(x)
                                  paste(x / 10 ^ 6, 'Mb', sep = ' ')
                              ) + ggplot2::theme(
                                legend.position = 'right',
                                legend.direction = "vertical",
                                axis.text.x = ggplot2::element_text(
                                  angle = 45,
                                  hjust = 1,
                                  vjust = 1
                                ),
                                axis.title.y.right = ggplot2::element_text(hjust = 0.6, vjust =
                                                                             2),
                                axis.title.y.left  = ggplot2::element_text(hjust = 0.92, vjust =
                                                                             2)

                              ) + ggplot2::xlab(Chr) +
                              ggplot2::scale_color_manual(values = Colors)+
                              ggplot2::theme(aspect.ratio = 2)


                          } else{
                            plot = ggplot2::ggplot()

                          }



                          #save
                          if (Save) {
                            start = min(scRT$start) / 10 ^ 6
                            end = max(scRT$end) / 10 ^ 6
                            group = unique(scRT$group)
                            Chr = unique(scRT$chr)

                            ggplot2::ggsave(
                              plot = plot,
                              filename = file.path(
                                Out,
                                paste0(group, '_', Chr, '_', start, 'Mb_', end, 'Mb.pdf')
                              ),
                              device = grDevices::cairo_pdf
                            )
                          }
                          #render plot
                          output$plot__scPlot <- renderPlot({
                            plot
                          },
                          height = function() {
                            session$clientData[[paste0('output_', ID, '-plot__scPlot_width')]]*1.5
                          })
                        })

  }
