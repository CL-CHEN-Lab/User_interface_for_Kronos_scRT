#'  RT changes in function of ΔRT
#'
#' @return ggplot element
#'
#' @importFrom  dplyr as_label select tibble mutate pull
#' @importFrom  tidyr %>% drop_na spread
#' @importFrom  ggplot2 aes annotate coord_cartesian element_blank element_text geom_density geom_hex geom_polygon ggplot scale_fill_gradient2 scale_fill_gradient  scale_y_continuous theme theme_bw theme_set
#' @importFrom scales percent_format
#' @importFrom foreach %:% %do% foreach
#' @import GGally
#'
#' @param ..., one or more dataframes containing bulk and pseudo-bulk RTs and the following columns:chr, start, end, group, basename, RT
#' @param deltaRT, if deltaRT is set, the function behavior changes and it shows RT changes at that deltaRT threshold
#' @param colors, plot colors
#' @export
#'
RT_changes_plot = function(...,
                           deltaRT = NULL,
                           colors=NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`
  `%:%` = foreach::`%:%`
  `%do%` = foreach::`%do%`

  #set theme
  ggplot2::theme_set(new = ggplot2::theme_bw())

  #set colors
  if(is.null(colors)){
    colors= c(
      'Late to Early' = 'darkred',
      'to Earlier' = 'red',
      'to Later' = 'blue',
      'Early to Late' = 'darkblue'
    )
  }else{
    colors=colorRampPalette(colors)(4)
    names(colors)=c('Late to Early', 'to Earlier', 'to Later', 'Early to Late')
  }

  df = do.call('rbind', list(...))

  df = df %>%
    dplyr::select(-group) %>%
    tidyr::spread(basename, RT) %>%
    tidyr::drop_na() %>%
    dplyr::select(-chr, -start, -end)
  Samples = names(df)
  len = nrow(df)
  if (is.null(deltaRT)) {
    data = foreach::foreach(a = Samples, .combine = 'rbind') %:%
      foreach::foreach(b = Samples, .combine = 'rbind') %do% {
        if (a == b) {
          dplyr::tibble()
        } else{
          sub1 = df[a] %>% dplyr::pull()
          sub2 = df[b] %>% dplyr::pull()

          new_data = lapply(seq(
            from = 1,
            to = 0,
            by = -0.1
          ), function(x) {
            diff = sub1 - sub2

            earlier = 100 * sum(abs(diff) > x & diff < 0) / len
            later = 100 * sum(abs(diff) > x & diff > 0) / len
            crossing_earlier = 100 * sum(abs(diff) > x &
                                           diff < 0  &
                                           sub1 < 0.5 & sub2 > 0.5) / len
            crossing_later = 100 * sum(abs(diff) > x &
                                         diff > 0 &
                                         sub1 > 0.5 & sub2 < 0.5) / len
            earlier = earlier - crossing_earlier
            later = later - crossing_later

            dplyr::tibble(
              Group = paste('ΔRT =', a, '-', b),
              deltaRT = x,
              change = c(
                'to Earlier',
                'to Later',
                'Late to Early',
                'Early to Late'
              ),
              percent = c(earlier, later, crossing_earlier, crossing_later)
            )
          })

          do.call('rbind', new_data)

        }
      }

    Max = data %>%
      dplyr::filter(percent > 0) %>%
      dplyr::pull(deltaRT) %>%
      max()

    plot = data %>%
      dplyr::mutate(change = factor(
        change,
        c('Late to Early', 'to Earlier', 'to Later', 'Early to Late')
      )) %>%
      dplyr::filter(deltaRT <= Max) %>%
      ggplot2::ggplot(ggplot2::aes(deltaRT, percent, color = change)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ Group, ncol = (length(Samples) - 1)) +
      ggplot2::labs(x = 'ΔRT', y = '% of changing regions', color = 'Direction') +
      ggplot2::theme(legend.position = 'top') +
      ggplot2::scale_color_manual(
        values =colors
      )
  } else if (deltaRT <= 1 & deltaRT >= 0) {
    data = foreach::foreach(a = Samples, .combine = 'rbind') %:%
      foreach::foreach(b = Samples, .combine = 'rbind') %do% {
        if (a == b) {
          dplyr::tibble()
        } else{
          sub1 = df[a] %>% dplyr::pull()
          sub2 = df[b] %>% dplyr::pull()

          diff = sub1 - sub2
          earlier = 100 * sum(abs(diff) > deltaRT & diff < 0) / len
          later = 100 * sum(abs(diff) > deltaRT & diff > 0) / len
          crossing_earlier = 100 * sum(abs(diff) > deltaRT &
                                         diff < 0  &
                                         sub1 < 0.5 & sub2 > 0.5) / len
          crossing_later = 100 * sum(abs(diff) > deltaRT &
                                       diff > 0 &
                                       sub1 > 0.5 & sub2 < 0.5) / len
          earlier = earlier - crossing_earlier
          later = later - crossing_later

          dplyr::tibble(
            Group = paste('ΔRT =', a, '-', b),
            deltaRT = deltaRT,
            change = c('to Earlier', 'to Later', 'Late to Early', 'Early to Late'),
            percent = c(earlier, later, crossing_earlier, crossing_later)
          )


        }
      }

    data = data %>%
      dplyr::mutate(change = factor(
        change,
        c('Late to Early', 'to Earlier', 'to Later', 'Early to Late')
      ))

    plot = data %>%
      dplyr::mutate(change = factor(
        change,
        c('Late to Early', 'to Earlier', 'to Later', 'Early to Late')
      )) %>%
      ggplot2::ggplot(ggplot2::aes(change, percent, fill = change)) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~ Group, ncol = (length(Samples) - 1)) +
      ggplot2::labs(x = paste('ΔRT =', deltaRT),
                    y = '% of changing regions',
                    color = 'Direction') +
      ggplot2::theme(legend.position = 'top') +
      ggplot2::scale_fill_manual(
        values = colors
      )

  } else{
    stop('deltaRT is supposed to be a value between 0 and 1')
  }

  return(plot)

}
