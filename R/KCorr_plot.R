#'  Customized GGally::ggparis function
#'
#' @return ggplot element
#'
#' @importFrom  dplyr as_label select tibble mutate
#' @importFrom  tidyr %>% drop_na spread
#' @importFrom  ggplot2 aes annotate coord_cartesian element_blank element_text geom_density geom_hex geom_polygon ggplot scale_fill_gradient2 scale_fill_gradientn scale_x_continuous scale_y_continuous theme theme_bw theme_set
#' @importFrom stringr str_to_title
#' @importFrom GGally ggpairs
#'
#' @param df, a dataframe containing bulk and pseudo-bulk RTs and the following columns:chr, start, end, group, basename, RT
#'
#' @export
#'
KCorr_plot = function(df, method = 'spearman') {
  #load required operators
  `%>%` = tidyr::`%>%`
  #set theme
  ggplot2::theme_set(new = ggplot2::theme_bw())



  plot = GGally::ggpairs(
    title = paste(stringr::str_to_title(method), 'correlation'),
    df %>%
      dplyr::select(-group) %>%
      tidyr::spread(basename, RT) %>%
      tidyr::drop_na() %>%
      dplyr::select(-chr, -start, -end) ,
    diag = list(
      continuous = function(data, mapping) {
        p <- ggplot2::ggplot(data, mapping) +
          ggplot2::geom_density(ggplot2::aes(y = ..density.. / max(..density..)),
                                fill = 'grey') +
          ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1)) +
          ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1))
        return(p)
      }
    ),
    upper = list(
      continuous = function(data, mapping) {
        Correlation = stats::cor(data[dplyr::as_label(mapping$x)], data[dplyr::as_label(mapping$y)], method = method)
        data = dplyr::tibble(x = seq(0, 2 * pi, length.out = 200),
                             Corr = Correlation) %>%
          dplyr::mutate(y = 0.5 + Corr / 2 * sin(x),
                        x = 0.5 + Corr / 2 * cos(x))
        p <-
          ggplot2::ggplot(data, ggplot2::aes(x = x,
                                             y = y,
                                             fill = Corr)) +
          ggplot2::geom_polygon() +
          ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank()
          ) +
          ggplot2::annotate(
            'text',
            x = 0.5,
            y = 0.5,
            label = paste("Corr:", round(unique(Correlation), 3), sep = '\n'),
            color = 'red'
          ) +
          ggplot2::scale_fill_gradient2(
            low = '#BCAF6FFF',
            high = '#00204DFF',
            mid = '#7C7B78FF',
            midpoint = 0,
            limits = c(-1, 1)
          ) + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

        return(p)
      }
    ),
    lower = list(
      continuous = function(data, mapping, ...) {
        p <- ggplot2::ggplot(data = data, mapping = mapping) +
          ggplot2::geom_hex(bins = 50, ggplot2::aes(fill = ..ndensity..)) +
          ggplot2::scale_fill_gradientn(
            'Density',
            colours = c(
              "#FFEA46FF",
              "#D3C164FF",
              "#A69D75FF",
              "#7C7B78FF",
              "#575C6DFF",
              "#233E6CFF",
              "#00204DFF"
            )
          ) +
          ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
          ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1)) +
          ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1))

        return(p)
      }
    ),
    legend = c(2, 1)
  ) + ggplot2::theme(
    legend.position = "right",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

  return(plot)
}
