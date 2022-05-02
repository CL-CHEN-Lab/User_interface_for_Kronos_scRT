#'  scRT plot of a region of interest
#'
#' @return plot
#'
#' @importFrom  dplyr tibble filter rename arrange mutate summarise ungroup
#' @importFrom  tidyr %>% gather
#' @importFrom  ggplot2 aes annotate element_text facet_grid geom_path geom_rect ggplot labs scale_color_manual scale_fill_gradient scale_fill_manual scale_x_continuous scale_y_continuous sec_axis theme theme_bw theme_set xlab
#' @importFrom ggrastr geom_tile_rast
#' @importFrom stringr str_count
#'
#' @param pseudoBulkRT, pseudoBulkRT that can be rbinded with a refrence RT (generatede bu RebinRT)
#' @param S_scCN, S-phase scCN dataframe created by Replication_state (optionally filtered)
#' @param Coordinates, named list containing chr,start and end position
#' @param rasterized_heatmap, logic: whether or not the heatmap should be rasterized
#'
#' @export
#'

scRTplot = function(pseudoBulkRT,
                    S_scCN,
                    Coordinates,
                    Plot = 'scRT',
                    rasterized_heatmap = F) {
  #load required operators
  `%>%` = tidyr::`%>%`

  #set theme
  ggplot2::theme_set(new = ggplot2::theme_bw())


  #check Plot

  if (!Plot %in% c('scRT', 'S_scCN', 'Norm. S_scCN')) {
    stop('Plot can be either scRT, S_scCN or Norm. S_scCN')
  }


  #calculate exremes
  extremes = S_scCN %>%
    dplyr::ungroup() %>%
    dplyr::summarise(CN_bg = round(stats::quantile(CN_bg, c(0.01, 0.99)), 1),
                     CN = round(stats::quantile(CN, c(0.01, 0.99)), 1))
  #filter data
  pseudoBulkRT = pseudoBulkRT %>%
    dplyr::filter(chr == Coordinates$chr,
                  start >= Coordinates$start,
                  end <= Coordinates$end)
  S_scCN = S_scCN %>%
    dplyr::filter(chr == Coordinates$chr,
                  start >= Coordinates$start,
                  end <= Coordinates$end)

  #depending on selection change aestetics of the plot
  if (Plot == 'scRT') {
    S_scCN = S_scCN %>%
      dplyr::rename(Value = Rep) %>%
      dplyr::mutate(Value = ifelse(Value, "Replicated", "Unreplicated"))

    ggEXTRA = ggplot2::ggplot() +
      ggplot2::scale_fill_manual(values = c("Replicated" = '#a7001b',
                                            "Unreplicated" = '#005095')) +
      ggplot2::labs(fill = 'State', color = 'Sample')

  } else if (Plot == 'S_scCN') {
    S_scCN = S_scCN %>%
      dplyr::rename(Value = CN)

    ggEXTRA = ggplot2::ggplot() +
      ggplot2::scale_fill_gradient(
        low = '#c56700',
        high = '#6d0042',
        limits = c(extremes$CN[1], extremes$CN[2]),
        oob = scales::squish
      ) +
      ggplot2::labs(fill = 'CNV ', color = 'Sample')

  } else if (Plot == 'Norm. scCN') {
    S_scCN = S_scCN %>%
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

  Chr = Coordinates$chr
  Maxi = max(S_scCN$newIndex)
  n = stringr::str_count(Maxi)

  plot =
    ggEXTRA +
    ggplot2::geom_path(
      data = pseudoBulkRT %>%
        dplyr::mutate(mid = (start + end) / 2,
                      RT = RT * Maxi / 6 + Maxi / 20) %>%
        tidyr::gather(Plot, pos, start, end) %>%
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
    )
  if (rasterized_heatmap) {
    plot = plot + ggrastr::geom_tile_rast(data = S_scCN ,
                                          ggplot2::aes(
                                            x = start,
                                            y = -newIndex,
                                            fill = Value
                                          ))


  } else{
    plot = plot +
      ggplot2::geom_rect(
        data = S_scCN ,
        ggplot2::aes(
          xmin = start,
          xmax = end,
          ymin = -newIndex,
          ymax = -newIndex - 1,
          fill = Value
        )
      )
  }

  plot = plot +
    ggplot2::scale_y_continuous(
      breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
      labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
      name = 'RT',
      sec.axis = ggplot2::sec_axis(
        ~ .,
        breaks = c(-seq(1, round(Maxi / 10 ^ (
          n - 1
        )), 1) * 10 ^ (n - 1)) - 0.5,
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
    ggplot2::scale_color_manual(values = c("#0262b5",
                                           "#f15bba")) +
    ggplot2::theme(aspect.ratio = 2) +
    ggplot2::facet_grid(~ group)
  return(plot)

}


#'  scCN plot
#'
#' @return plot
#'
#' @importFrom  dplyr filter select arrange mutate ungroup
#' @importFrom  tidyr %>% drop_na
#' @importFrom  ggplot2 aes element_blank element_text facet_grid geom_rect ggplot scale_fill_manual scale_x_continuous theme theme_bw theme_set
#' @importFrom  foreach %do% foreach
#' @importFrom ggrastr geom_tile_rast
#' @importFrom viridis viridis
#'
#' @param S_scCN, S-phase scCN dataframe created by Replication_state (optionally filtered)
#' @param G_scCN, G1/G2-phase scCN dataframe created by Replication_state (optionally filtered)
#' @param CN_limit, Max number of different level to visulize from the lowest one.
#' @param Coordinates, data frame containing chr,start and end position (multiple regions allowed, facultative)
#' @param rasterized_heatmap, logic: whether or not the heatmap should be rasterized
#'
#' @export
#'
#'
scCNplot = function(S_scCN,
                    G_scCN,
                    CN_limit = 4,
                    Coordinates = NULL,
                    rasterized_heatmap = F) {
  #load required operators
  `%>%` = tidyr::`%>%`
  `do` = foreach::`%do%`

  #set theme
  ggplot2::theme_set(new = ggplot2::theme_bw())

  #combine data
  tracks =
    rbind(
      G_scCN %>% dplyr::select(chr, start, end, CN, newIndex, group) %>%
        dplyr::mutate(phase = 'G1G2-phase') %>%
        dplyr::ungroup(),
      S_scCN %>% dplyr::select(chr, start, end, CN, newIndex, group) %>%
        dplyr::mutate(phase = 'S-phase') %>%
        dplyr::ungroup()
    )

  tracks = tracks %>% dplyr::mutate(
    CN = round(CN),
    MaxAllowed = min(CN, na.rm = T) + CN_limit,
    CN = ifelse(CN > MaxAllowed, MaxAllowed, CN)
  ) %>%
    dplyr::arrange(CN) %>%
    dplyr::mutate(CN = factor(CN)) %>%
    tidyr::drop_na() %>%
    dplyr::arrange(chr, start)

  #filter data if provided
  if (!is.null(Coordinates)) {
    tracks = foreach::foreach(i = 1:nrow(Coordinates), .combine = 'rbind') %do%
    {
      tracks %>%
        dplyr::filter(chr == Coordinates$chr[i],
                      start >= Coordinates$start[i],
                      end <= Coordinates$end[i])
    }
  }

  if (rasterized_heatmap) {
    plot =
      ggplot2::ggplot() +
      ggrastr::geom_tile_rast(data = tracks ,
                              ggplot2::aes(
                                x = start,
                                y = -newIndex,
                                fill = CN
                              ))

  } else{
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
      )
  }

  plot = plot + ggplot2::facet_grid(phase ~ chr, scales = 'free', space =  'free') +
    ggplot2::scale_x_continuous(
      labels = function(x)
        paste(x / 10 ^ 6, 'Mb', sep = ' ')
    ) + ggplot2::theme(
      legend.position = 'right',
      legend.direction = "vertical",
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )

    ) +
    ggplot2::scale_fill_manual(values = viridis::viridis(CN_limit + 1))

  return(plot)
}
