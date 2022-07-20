#'  Customized GGally::ggparis function
#'
#' @return ggplot element
#'
#' @importFrom  dplyr as_label select tibble mutate
#' @importFrom  tidyr %>% drop_na spread unite
#' @importFrom  ggplot2 aes annotate coord_cartesian element_blank element_text geom_density geom_hex geom_polygon ggplot scale_fill_gradient2 scale_fill_gradientn scale_x_continuous scale_y_continuous theme theme_bw theme_set
#' @importFrom stringr str_to_title
#' @importFrom viridis viridis inferno magma plasma cividis
#' @import GGally
#'
#' @param ..., one or multiple dataframes containing bulk and pseudo-bulk RTs and the following columns:chr, start, end, group, basename, RT
#' @param method, correlation method
#' @param density_fill, diagonal plot fill color
#' @param correlation_gradient, min, midpoint and max colors for the upper triangle
#' @param correlation_text_color, text color upper triangle correlation
#' @param correlation_text_size, text size upper triangle correlation
#' @param hex_color_palette, color palette for the lower triangle density gradient
#' @param hex_fill, multiple manual colors to used for the lower triangle density gradient, if provided hex_color_palette is ignored
#' @export
#'
KCorr_plot = function(..., method = 'spearman',density_fill=NULL,correlation_gradient=NULL,correlation_text_color='red',correlation_text_size=15,hex_color_palette=c('viridis','inferno','magma','plasma','cividis'),hex_fill=NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`
  #set theme
  ggplot2::theme_set(new = ggplot2::theme_bw())
  #adjust font size
  correlation_text_size = correlation_text_size * 0.35
  #set colors
  if(is.null(density_fill)){
    density_fill='grey'
  }

  if(is.null(correlation_gradient)){
    correlation_gradient=c('#BCAF6FFF','#7C7B78FF','#00204DFF')
  }else if(length(correlation_gradient)<3){
    correlation_gradient=colorRampPalette(correlation_gradient)(3)
  }

  if(is.null(hex_fill)){
      if(hex_color_palette[1]=='inferno'){
        hex_fill=viridis::inferno(10)
      }else if(hex_color_palette[1]=='magma'){
        hex_fill=viridis::magma(10)
      }else if(hex_color_palette[1]=='plasma'){
        hex_fill=viridis::plasma(10)
      }else if(hex_color_palette[1]=='cividis'){
        hex_fill=viridis::cividis(10)
      }else{
          hex_fill=viridis::viridis(10)
      }}

  df=do.call('rbind',list(...))

  plot = GGally::ggpairs(
    title = paste(stringr::str_to_title(method), 'correlation'),
    tryCatch(expr = df %>%
      dplyr::select(-group) %>%
      tidyr::spread(basename, RT) %>%
      tidyr::drop_na() %>%
      dplyr::select(-chr, -start, -end),
      error=function(x) df%>%
      tidyr::unite(basename,basename,group,sep = ' - ') %>%
        tidyr::spread(basename, RT) %>%
        tidyr::drop_na() %>%
        dplyr::select(-chr, -start, -end)),
    diag = list(
      continuous = function(data, mapping) {
        p <- ggplot2::ggplot(data, mapping) +
          ggplot2::geom_density(ggplot2::aes(y = ..density.. / max(..density..)),
                                fill = density_fill) +
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
            color = correlation_text_color,
            size=correlation_text_size
          ) +
          ggplot2::scale_fill_gradient2(
            low = correlation_gradient[1],
            high = correlation_gradient[3],
            mid = correlation_gradient[2],
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
            colours = hex_fill
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
