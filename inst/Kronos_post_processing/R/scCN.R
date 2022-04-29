scCN_ui <- function(id) {
  ns <- shiny::NS(paste0('scCN',id))
  shinydashboard::box(
    id = ns('box'),
    status = 'primary',
    title = id,
    width = 12,
    solidHeader = T,
    collapsible = T,
    shiny::fluidRow(plotOutput(ns('plot__scCN'), height = "100%",width = '100%'))
  )
}


scCN_server <-
  function(id, S_Traks,G_Traks,Levels=10,out,save) {
    shiny::moduleServer(paste0('scCN',id),
                        function(input,
                                 output,
                                 session,
                                 S = S_Traks,
                                 G = G_Traks,
                                 L =Levels,
                                 Out=out,
                                 Save=save,
                                 ID=paste0('scCN',id)) {
                          #load required operators
                          `%>%` = tidyr::`%>%`

                          tracks=tryCatch(
                            rbind(
                              G %>% dplyr::select(chr,start,end,CN,newIndex,group)%>%
                                dplyr::mutate(phase='G1G2-phase')%>%
                                dplyr::ungroup(),
                              S %>% dplyr::select(chr,start,end,CN,newIndex,group)%>%
                                dplyr::mutate(phase='S-phase')%>%
                                dplyr::ungroup()
                            ),
                              error= function(x) dplyr::tibble()

                          )

                          if (nrow(tracks) != 0) {

                            tracks= tracks%>%dplyr::mutate(CN=round(CN),
                                             MaxAllowed=min(CN,na.rm = T)+L,
                                             CN= ifelse(CN > MaxAllowed,MaxAllowed,CN))%>%
                              dplyr::arrange(CN)%>%
                              dplyr::mutate(CN=factor(CN))%>%
                              tidyr::drop_na()%>%
                              dplyr::arrange(chr,start)

                            Chr = unique(tracks$chr)

                            plot =
                              ggplot2::ggplot() +
                              ggplot2::geom_rect(
                                data = tracks ,
                                ggplot2::aes(
                                  xmin = start,
                                  xmax = end,
                                  ymin = -newIndex,
                                  ymax = -newIndex - 1,
                                  fill = CN
                                )
                              )+ggplot2::facet_grid(phase ~ .,scales = 'free_y',space =  'free_y') +
                              ggplot2::scale_x_continuous(
                                labels = function(x)
                                  paste(x / 10 ^ 6, 'Mb', sep = ' ')
                              ) + ggplot2::theme(
                                legend.position = 'right',
                                legend.direction = "vertical",
                                axis.text.y = ggplot2::element_blank(),
                                axis.ticks.y=ggplot2::element_blank(),
                                axis.text.x = ggplot2::element_text(
                                  angle = 45,
                                  hjust = 1,
                                  vjust = 1
                                )

                              ) +
                              ggplot2::xlab(Chr) +
                              ggplot2::scale_fill_manual(values = viridis::viridis(L+1))

                          } else{
                            plot = ggplot2::ggplot()

                          }
                          output$plot__scCN <- renderPlot({plot},
                                                          height = function() {
                                                            session$clientData[[paste0('output_', ID, '-plot__scCN_width')]]*1.5
                                                          })
                          #save
                          if(Save){

                            start=min(tracks$start)/10^6
                            end=max(tracks$end)/10^6
                            group=unique(tracks$group)
                            Chr = unique(tracks$chr)

                            ggplot2::ggsave(plot = plot,
                                            filename = file.path(
                                              Out,
                                              paste0(group, '_', Chr, '_', start, 'Mb_', end, 'Mb.pdf')
                                            ),device = grDevices::cairo_pdf)
                          }

                        })
  }

