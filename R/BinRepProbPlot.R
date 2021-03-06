#'  Bin Replication probability plot
#'
#' @importFrom  cowplot plot_grid
#' @importFrom  dplyr group_by mutate summarise ungroup
#' @importFrom  tidyr %>%
#' @importFrom  ggplot2 aes annotate coord_cartesian element_blank facet_grid geom_boxplot ggplot labs scale_fill_manual scale_y_continuous theme unit xlab ylab
#'
#' @param Variability, a row bind combination of dataframes produce by Prepare_G1G2_phase_cells_forBinRepProb and Prepare_S_phase_cells_forBinRepProb or ether of the 2
#' @param upper_range, limits of the upper panel
#' @param lower_range, limits of the lower panel
#' @param relative_heights, relative heights between the 3 panels
#' @param colors, boxplot fill colors to pass to ggplot
#'
#' @export
#'

BinRepProbPlot = function(Variability,
                          upper_range = c(0.95, 1),
                          lower_range = c(0, 0.05),
                          relative_heights = c(1, 2, 1),
                          colors=NULL) {
  #load operator
  `%>%` = tidyr::`%>%`

  #number of groups
  plot_groups=unique(Variability$type)
  #define colors
  if(is.null(colors)){
    colors=c(
      "Early S cells" = '#a7001b',
      "Late S cells" = '#005095',
      'Mid S cells' = '#dfbd31',
      'G1/G2 cells' = 'grey'
    )
  }else if(length(colors)< length(plot_groups) & is.null(names(colors))){
    #if the given colors are less than the number of types and the vector is not named
    default_colors=c(
      "Early S cells" = '#a7001b',
      "Late S cells" = '#005095',
      'Mid S cells' = '#dfbd31',
      'G1/G2 cells' = 'grey'
    )
    #select default colors
    default_colors=default_colors[plot_groups]

    # substitute default colors with provided ones
    colors=sapply(1:length(plot_groups), function(x) {
      ifelse(is.na(colors[x]),default_colors[x],colors[x])
    } )

    #name the vector
    names(colors)=plot_groups

  }else if(length(colors)<  length(plot_groups) & !is.null(names(colors))){

    #if not enough colors have been provided and the list is named
    default_colors=c(
      "Early S cells" = '#a7001b',
      "Late S cells" = '#005095',
      'Mid S cells' = '#dfbd31',
      'G1/G2 cells' = 'grey'
    )
    # replace default colors with provided ones
    colors=sapply(plot_groups, function(x){
      ifelse(is.na(colors[x]),default_colors[x],colors[x])
    })
    #reassign names
    names(colors)=plot_groups
  }

  # select
  p = Variability %>%
    ggplot2::ggplot(ggplot2::aes(RT, BinRepProb, fill = type)) +
    ggplot2::geom_boxplot() +
    ggplot2::ylab('Bin probability of replication') +
    ggplot2::scale_fill_manual(
      values = colors
    ) +
    ggplot2::labs(fill = NULL) +
    ggplot2::facet_grid(~ type) +
    ggplot2::theme(legend.position = 'none')

  q = cowplot::plot_grid(
    p +
      ggplot2::coord_cartesian(ylim = upper_range) +
      ggplot2::labs(fill = NULL, y = '', x = '') +
      ggplot2::scale_y_continuous(breaks = upper_range) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        legend.position = 'none',
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0.5, 0, 0, 0), "cm")
      ),
    p +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = 'none',
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::xlab('') +
      ggplot2::annotate(
        "rect",
        xmin = 0.5,
        xmax = 11.5,
        ymin = lower_range[1],
        ymax = lower_range[2],
        fill = NA,
        color = 'black',
        linetype = 'dashed'
      ) +
      ggplot2::annotate(
        "rect",
        xmin = 0.5,
        xmax = 11.5,
        ymin = upper_range[1],
        ymax = upper_range[2],
        fill = NA,
        color = 'black',
        linetype = 'dashed'
      ) + ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")),
    p + ggplot2::coord_cartesian(ylim = lower_range) +
      ggplot2::labs(fill = NULL, y = '', x = 'RT') +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        legend.position = 'none',
        plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
      ) +
      ggplot2::scale_y_continuous(breaks = lower_range),
    ncol = 1,
    rel_heights = relative_heights
  )

  return(q)

}


#'  Prepare G1/G2 phase cells for BinRepProb plots
#'
#' @importFrom  dplyr filter group_by inner_join mutate pull select summarise ungroup
#' @importFrom  tidyr %>%
#'
#' @param G1.G2, G1/G2 single cell CNV dataframe
#' @param RT, pseudo bulk replication timing dataframe
#' @param quantile.range, Ploidy quantile range to keep a G1/G2 cell
#'
#'
#' @export
#'

Prepare_G1G2_phase_cells_forBinRepProb = function(G1.G2, RT, quantile.range =
                                                    c(0.25, 0.75)) {
  #load operator
  `%>%` = tidyr::`%>%`

  #identify G1 cells to carry on
  MP = G1.G2 %>%
    dplyr::group_by(Cell,basename, group)  %>%
    dplyr::summarise(Mean_ploidy = mean(CN)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      a = stats::quantile(Mean_ploidy, min(quantile.range)),
      b = stats::quantile(Mean_ploidy, max(quantile.range))
    ) %>%
    dplyr::filter(Mean_ploidy < b, Mean_ploidy >
                    a) %>% dplyr::pull(Cell)

  #filter cell and calculate variability
  G1.G2 %>%
    dplyr::filter(Cell %in% MP) %>%
    dplyr::mutate(type = 'G1/G2 cells') %>%
    dplyr::select(chr, start, end, Rep, newIndex, type, Cell,basename, group) %>%
    dplyr::inner_join(RT %>% dplyr::select(chr, start, end, RT, group),
                      by = c("chr", "start", "end", "group")) %>%
    dplyr::mutate(type = factor(
      type,
      levels = c('Early S cells',
                 'Mid S cells',
                 'Late S cells',
                 'G1/G2 cells')
    )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(RT = round(RT, 1)) %>%
    dplyr::group_by(chr, start, end, RT, type, basename,group) %>%
    dplyr::summarise(BinRepProb = mean(Rep, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(RT = factor(RT, levels = seq(1, 0, -0.1))) %>%
    return()

}


#'  Prepare S phase cells for BinRepProb plots
#'
#' @importFrom  dplyr filter case_when inner_join mutate select ungroup
#' @importFrom  tidyr %>%
#'
#' @param S, S single cell CNV dataframe
#' @param RT, pseudobulk replication timing dataframe
#' @param Early.cells, percentage of replication range defying Early replicating cells (0-100)
#' @param Mid.cells, percentage of replication range defying Mid replicating cells (0-100)
#' @param Late.cells, percentage of replication range defying Late replicating cells (0-100)
#' @export
#'

Prepare_S_phase_cells_forBinRepProb = function(S,
                                               RT,
                                               Early.cells = c(0, 30),
                                               Mid.cells = c(40, 60),
                                               Late.cells = c(70, 100)) {
  Early.cells = Early.cells / 100
  Mid.cells = Mid.cells / 100
  Late.cells = Late.cells / 100
  #load operator
  `%>%` = tidyr::`%>%`

  #select cells within replication ranges
  S %>%
    dplyr::filter((
      PercentageReplication <= max(Early.cells) &
        PercentageReplication >= min(Early.cells)
    ) |
      (
        PercentageReplication <= max(Mid.cells) &
          PercentageReplication >= min(Mid.cells)
      ) |
      (
        PercentageReplication <= max(Late.cells) &
          PercentageReplication >= min(Late.cells)
      )
    ) %>%
    dplyr::mutate(type = dplyr::case_when((
      PercentageReplication <= max(Early.cells) &
        PercentageReplication >= min(Early.cells)
    ) ~ 'Early S cells',
    (
      PercentageReplication <= max(Late.cells) &
        PercentageReplication >= min(Late.cells)
    ) ~ 'Late S cells',
    (
      PercentageReplication <= max(Mid.cells) &
        PercentageReplication >= min(Mid.cells)
    ) ~ 'Mid S cells'
    )) %>%
    dplyr::select(chr, start, end, Rep, newIndex, type, Cell,basename, group) %>%
    dplyr::inner_join(RT %>% dplyr::select(chr, start, end, RT, group),
                      by = c("chr", "start", "end", "group")) %>%
    dplyr::mutate(type = factor(
      type,
      levels = c('Early S cells',
                 'Mid S cells',
                 'Late S cells',
                 'G1/G2 cells')
    )) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(RT = round(RT, 1)) %>%
    dplyr::group_by(chr, start, end, RT, type, basename,group) %>%
    dplyr::summarise(BinRepProb = mean(Rep, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(RT = factor(RT, levels = seq(1, 0, -0.1)))  %>%
    return()

}
