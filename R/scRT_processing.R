#'  Bins the genome in equally sized bins
#'
#' @return GenomicRanges object
#'
#' @importFrom  dplyr tibble mutate
#' @importFrom  tidyr %>% drop_na
#' @importFrom  foreach foreach %do% %dopar%
#' @importFrom  doSNOW registerDoSNOW
#' @importFrom snow makeCluster stopCluster
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#'
#' @param Chr_size, a dataframe (chr=chrom name, size=chrom size)
#' @param size, bin size
#' @param Cores, number of threads to parallelise
#'
#' @export
#'

GenomeBinning <-
  function(Chr_size,
           size = 5 * 10 ^ 6,
           Chr_filter = NULL,
           Cores = 1) {
    #load required operators
    `%dopar%` = foreach::`%dopar%`
    `%do%` = foreach::`%do%`
    `%>%` = tidyr::`%>%`

    #filter
    if (!is.null(Chr_filter)) {
      Chr_size = Chr_size %>%
        dplyr::filter(chr %in% Chr_filter)
    }


    #parallelise
    cl <- snow::makeCluster(Cores)
    doSNOW::registerDoSNOW(cl)
    on.exit(snow::stopCluster(cl))

    #create bins
    bins = foreach::foreach(
      chr = 1:nrow(Chr_size),
      .combine = 'rbind'
    ) %dopar% {
      bins = seq(from = size,
                 to = Chr_size$size[chr] ,
                 by =  size)
      bins = dplyr::tibble(chr =  Chr_size$chr[chr],
                           start = bins - size,
                           end = bins)
      bins
    }
    bins = bins %>%
      #convert bins into granges
      GenomicRanges::makeGRangesFromDataFrame(
        seqnames.field  = 'chr',
        start.field = 'start' ,
        end.field = 'end'
      )

    return(bins)
  }

#'  Applies diagnostic settings to PerCell dataframe
#'
#' @return tibble
#'
#' @importFrom  dplyr case_when inner_join mutate
#' @importFrom  tidyr %>%
#'
#' @param PerCell, PerCell dataframe created by CallCNV
#' @param Settings, Setting dataframe created by diagnostic
#' @param Basename_leves, Sample basenames order (optional)
#' @param Group_leves, Group basenames order (optional)
#'
#' @export
#'

AdjustPerCell <-
  function(PerCell,
           Settings,
           Basename_leves = NULL,
           Group_leves = NULL) {
    #load required operators
    `%>%` = tidyr::`%>%`

    #merge PerCell with settings
    PerCell = PerCell %>%
      dplyr::inner_join(Settings, by = c('basename', 'group')) %>%
      #filter out cells that won't be used
      dplyr::filter(coverage_per_1Mbp >= RPM_TH)

    #convert basename and group to factors

    if (is.null(Basename_leves)) {
      PerCell = PerCell %>%
        dplyr::mutate(basename = factor(basename,
                                        levels = sort(unique(basename))))
    } else{
      PerCell = PerCell %>%
        dplyr::mutate(basename = factor(basename,
                                        levels = Basename_leves))
    }

    if (is.null(Group_leves)) {
      PerCell = PerCell %>%
        dplyr::mutate(group = factor(group,
                                     levels = sort(unique(group))))
    } else{
      PerCell = PerCell %>%
        dplyr::mutate(group = factor(group,
                                     levels = Group_leves))
    }

    PerCell = PerCell %>%
      dplyr::mutate(
        #use setting file to change is_noisy and is_high_dimap
        is_noisy = as.logical(is_noisy),
        is_high_dimapd = ifelse(
          is.na(threshold_Sphase),
          is_high_dimapd,
          ifelse(normalized_dimapd > threshold_Sphase, T, F)
        ),
        is_noisy = ifelse(is_high_dimapd, T, is_noisy)
      ) %>%
      #filter out cells based on G1G2 median Ploidy (from settings) and ploidy confidence
      dplyr::filter(
        mean_ploidy > Ploidy / 1.50 ,
        mean_ploidy < Ploidy * 2,
        !ploidy_confidence <= 2 |
          ploidy_confidence == -100
      ) %>%
      # add staging
      dplyr::mutate(
        Type = ifelse(
          # if is_high_dimapd and is_noisy a cell is in S-phase
          as.logical(is_high_dimapd) == T &
            as.logical(is_noisy) == T,
          'S-phase cells',
          ifelse(
            #if a G1G2 threshold was provided include it into the formula
            ifelse(
              is.na(threshold_G1G2phase),
              # if there is not a G1G2 th G1/G2 cells are is_high_dimapd ==F and is noisy==F
              as.logical(is_high_dimapd) == F &
                as.logical(is_noisy) == F,
              # if there is  a G1G2 th G1/G2 cells are is_high_dimapd ==F and is noisy==F and normalized_dimapd < threshold_G1G2phase
              as.logical(is_high_dimapd) == F &
                as.logical(is_noisy) == F &
                normalized_dimapd < threshold_G1G2phase
            ),
            'G1/G2-phase cells',
            # other cells are unknown
            'unknown cells'
          )
        ),
        #correct cell mean ploidy if cell is in S-phase
        mean_ploidy_corrected = dplyr::case_when(
          Type == 'S-phase cells' &
            mean_ploidy < Ploidy ~ mean_ploidy / Sphase_second_part,
          Type == 'S-phase cells' &
            mean_ploidy > Ploidy ~ mean_ploidy / Sphase_first_part,
          T ~ mean_ploidy
        )
      )
    return(PerCell)
  }

#'  Corrects CN and filters cells based on AdjustPerCell output
#'
#' @return tibble
#'
#' @importFrom  dplyr inner_join mutate select
#' @importFrom  tidyr %>% drop_na
#'
#' @param scCN, scCN dataframe created by CallCNV
#' @param PerCell, PerCell dataframe created by AdjustPerCell
#' @param Basename_leves, Sample basenames order (optional)
#' @param Group_leves, Group basenames order (optional)
#' @param Chr_filter, filters chromosomes to keep in the analysis (optional)
#'
#' @export
#'

AdjustCN <-
  function(scCN,
           PerCell,
           Basename_leves = NULL,
           Group_leves = NULL,
           Chr_filter = NULL) {
    #load required operators
    `%>%` = tidyr::`%>%`
    #assign levels
    if (is.null(Basename_leves)) {
      scCN = scCN %>%
        dplyr::mutate(basename = factor(basename,
                                        levels = sort(unique(basename))))
    } else{
      scCN = scCN %>%
        dplyr::mutate(basename = factor(basename,
                                        levels = Basename_leves))
    }
    if (is.null(Group_leves)) {
      scCN = scCN %>%
        dplyr::mutate(group = factor(group,
                                     levels = sort(unique(group))))
    } else{
      scCN = scCN %>%
        dplyr::mutate(group = factor(group,
                                     levels = Group_leves))
    }
    if (is.null(Chr_filter)) {
      scCN = scCN %>%
        dplyr::mutate(chr = factor(chr,
                                   levels = sort(unique(chr))))
    } else{
      scCN = scCN %>%
        dplyr::mutate(chr = factor(chr,
                                   levels = Chr_filter))
    }

    scCN = scCN %>%
      #drop na
      tidyr::drop_na() %>%
      #join with Cell and basename columns from Per cell to filter cells to use
      dplyr::inner_join(
        PerCell %>% dplyr::select(
          Cell,
          basename,
          mean_ploidy,
          mean_ploidy,
          mean_ploidy_corrected
        ),
        by = c("Cell", "basename")
      ) %>%
      #correct CN if needed
      dplyr::mutate(
        copy_number_corrected = ifelse(
          mean_ploidy == mean_ploidy_corrected,
          copy_number,
          copy_number * mean_ploidy_corrected / mean_ploidy
        )
      )

    return(scCN)
  }

#'  rebins kronos CN data (preprocessing output)
#'
#' @return tibble
#'
#' @importFrom  dplyr n filter select as_tibble arrange group_by inner_join mutate summarise ungroup
#' @importFrom  tidyr %>%
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom IRanges findOverlaps
#' @importFrom matrixStats weightedMedian rowQuantiles
#' @importFrom S4Vectors subjectHits queryHits
#'
#' @param PerCell, PerCell dataframe created by CallCNV
#' @param scCN, scCN dataframe created by CallCNV
#' @param Bins, a GenomicRanges object with genomic binning
#' @param Sphase, either True or False
#'
#' @export
#'

Rebin <- function(PerCell, scCN, Bins, Sphase = NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`

  if (!is.logical(Sphase)) {
    stop('Sphase must be set as True or False')
  }

  #filter data
  if (Sphase) {
    selected_data = PerCell %>%
      dplyr::filter(Type == 'S-phase cells') %>%
      dplyr::arrange(mean_ploidy_corrected) %>%
      dplyr::mutate(index = 1:dplyr::n()) %>%
      dplyr::select(index,
                    Cell,
                    basename,
                    group)
  } else{
    selected_data = PerCell %>%
      dplyr::filter(Type == 'G1/G2-phase cells') %>%
      dplyr::mutate(index = as.numeric(factor(Cell))) %>%
      dplyr::select(index,
                    Cell,
                    basename,
                    group)

  }

  #convert into Granges
  scCN = scCN %>%
    dplyr::inner_join(selected_data, by = c("Cell", "basename", "group")) %>%
    GenomicRanges::makeGRangesFromDataFrame(
      seqnames.field  = 'chr',
      start.field = 'start' ,
      end.field = 'end',
      keep.extra.columns = T
    )

  #calculate median CN across a bin
  hits = IRanges::findOverlaps(Bins, scCN)

  #recover info
  scCN = cbind(
    dplyr::as_tibble(Bins[S4Vectors::queryHits(hits)]) %>% dplyr::select('chr' =
                                                                           seqnames, start, end),
    dplyr::as_tibble(scCN[S4Vectors::subjectHits(hits)]) %>% dplyr::select(-seqnames, -start, -end)
  ) %>%
    dplyr::group_by(Cell, index, basename, group, chr, start, end) %>%
    dplyr::summarise(CN = matrixStats::weightedMedian(x = copy_number_corrected,
                                                      w = width,
                                                      na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::select(chr,
                  start,
                  end,
                  CN,
                  Cell,
                  basename,
                  group,
                  index)

  return(scCN)
}

#'  rebins RT data
#'
#' @return tibble
#'
#' @importFrom  dplyr select as_tibble group_by mutate summarise ungroup
#' @importFrom  tidyr %>% drop_na
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom IRanges findOverlaps
#' @importFrom matrixStats weightedMedian rowQuantiles
#' @importFrom S4Vectors subjectHits queryHits
#'
#' @param RT, RT dataframe (chr,start,end,RT)
#' @param Bins, a GenomicRanges object with genomic binning
#' @param Chr_filter, filters chromosomes to keep in the analysis
#'
#' @export
#'

RebinRT <-
  function(RT,
           Bins,
           Chr_filter = NULL,
           Basename = NULL,
           Group = NULL) {
    #load required operators
    `%>%` = tidyr::`%>%`

    if (!is.null(Basename)) {
      RT = RT %>%
        dplyr::mutate(basename = Basename,
                      group = ifelse(is.null(Group), Basename, Group))
    }
    if (is.null(Basename) & !is.null(Group)) {
      stop('Basename is mandatory when Group is provided')
    }


    RT = RT %>%
      #make Granges
      GenomicRanges::makeGRangesFromDataFrame(
        seqnames.field = 'chr',
        start.field = 'start',
        end.field = 'end',
        keep.extra.columns = T
      )

    hits = IRanges::findOverlaps(Bins, RT)
    # calculate weighted median RT on new bins
    RT = cbind(
      dplyr::as_tibble(Bins[S4Vectors::queryHits(hits)]) %>% dplyr::select('chr' =
                                                                             seqnames, start, end) ,
      dplyr::as_tibble(RT[S4Vectors::subjectHits(hits)]) %>% dplyr::select(-seqnames, -start, -end)
    ) %>%
      dplyr::group_by(chr, start, end, group, basename) %>%
      dplyr::summarise(RT = matrixStats::weightedMedian(x = RT,
                                                        w =
                                                          width,
                                                        na.rm = T)) %>%
      dplyr::ungroup()



    #Assign level
    if (is.null(Chr_filter)) {
      RT = RT %>%
        dplyr::mutate(chr = factor(chr,
                                   levels = sort(unique(chr))))
    } else{
      RT = RT %>%
        dplyr::mutate(chr = factor(chr,
                                   levels = Chr_filter))
    }

    #resize RT
    RT = RT %>%
      dplyr::mutate(RT = (RT - min(RT)) / (max(RT) - min(RT))) %>%
      dplyr::ungroup() %>%
      tidyr::drop_na()

    return(RT)
  }

#'  Calculates median G1/G2 profile
#'
#' @return tibble
#'
#' @importFrom  dplyr group_by summarise
#' @importFrom  tidyr %>%
#'
#' @param G1_cells, scCN dataframe from Rebin function
#'
#' @export
#'

BackGround <- function(G1_scCN) {
  #load required operators
  `%>%` = tidyr::`%>%`

  #calculate background
  bg = G1_scCN %>%
    dplyr::group_by(basename, group, chr, start, end) %>%
    dplyr::summarise(background = median(CN))

  return(bg)
}

#'  binarizes data into replicated or not replicated bins
#'
#' @return tibble
#'
#' @importFrom  dplyr n filter select tibble arrange group_by inner_join mutate summarise ungroup
#' @importFrom  tidyr %>% drop_na
#' @importFrom  foreach %do% %dopar% foreach
#' @importFrom snow makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#'
#' @param Samples, a dataframe containing scCN info with the following columns:chr, start, end, group, basename, CN, index
#' @param background, a dataframe containing the median scCN info of the G1/G2 population with the following columns:chr, start, end, group, basename, background
#' @param Chr_filter, filters chromosomes to keep in the analysis
#' @param cores, number of threads to use while to parallelise
#'
#' @export
#'

Replication_state = function(Samples,
                             background,
                             Chr_filter = NULL,
                             cores = 3) {
  #load required operators
  `%dopar%` = foreach::`%dopar%`
  `%do%` = foreach::`%do%`
  `%>%` = tidyr::`%>%`

  #declare cluster
  cl <- snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))

  #assign factors
  if (is.null(Chr_filter)) {
    Samples = Samples %>%
      dplyr::ungroup() %>%
      dplyr::mutate(chr = factor(x =  chr,
                                 levels =  sort(unique(chr))))
  } else{
    Samples = Samples %>%
      dplyr::ungroup() %>%
      dplyr::mutate(chr = factor(x =  chr,
                                 levels =  Chr_filter))
  }

  #merge signal and bg and calculate their ratio
  Samples = Samples %>%
    dplyr::inner_join(background,
                      by = c("chr", "start", "end", "basename", 'group')) %>%
    dplyr::mutate(CN_bg = log2(CN / background)) %>%
    tidyr::drop_na() %>%
    dplyr::filter(is.finite(CN_bg))

  # identify threshold that minimizes the difference of the real data with a binary state (1 or 2)
  selecte_th = foreach::foreach(line = unique(Samples$basename),
                                .combine = 'rbind') %dopar% {
                                  #load required operators
                                  `%do%` = foreach::`%do%`
                                  `%>%` = tidyr::`%>%`

                                  sub_sig = Samples %>%
                                    dplyr::filter(basename == line)

                                  #identify range within looking for a CNV threshold to define replicated and not replicated values
                                  range = seq(0, 1, 0.01)

                                  th_temp = foreach::foreach(
                                    i = range,
                                    .combine = 'rbind'
                                  ) %do% {
                                    summary = sub_sig %>%
                                      dplyr::mutate(Rep = ifelse(CN_bg >= i, T, F),
                                                    Error = (Rep - CN_bg) ^ 2) %>%
                                      dplyr::group_by(Cell, index, basename, group)  %>%
                                      dplyr::summarise(summary = sum(Error))

                                    dplyr::tibble(
                                      th = i,
                                      basename = summary$basename,
                                      index = summary$index,
                                      sum_error = summary$summary
                                    )
                                  }
                                  th_temp
                                }

  selecte_th = selecte_th %>%
    dplyr::group_by(index, basename) %>%
    dplyr::filter(sum_error == min(sum_error)) %>%
    dplyr::summarise(th = min(th)) %>%
    dplyr::ungroup()

  # mark replicated bins
  Samples = Samples %>%
    dplyr::inner_join(selecte_th, by = c("index", "basename")) %>%
    dplyr::mutate(Rep = ifelse(CN_bg >= th, T, F))

  #identify new distribution in the S phase based the amount of replicated bins
  new_index_list = Samples %>%
    dplyr::group_by(Cell, index, basename, group) %>%
    dplyr::summarise(PercentageReplication = mean(Rep)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(PercentageReplication) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(newIndex = 1:dplyr::n()) %>%
    dplyr::select(oldIndex = index,
                  newIndex,
                  Cell,
                  basename,
                  group,
                  PercentageReplication)

  Samples = Samples %>%
    dplyr::inner_join(new_index_list,
                      by = c('Cell', 'index' = 'oldIndex', 'basename', 'group'))


  return(Samples)
}

#'  filter cells
#'
#' @return list
#'
#' @importFrom  dplyr n filter select tibble inner_join
#' @importFrom  tidyr %>% separate spread unite
#' @importFrom  foreach foreach %do%
#' @importFrom ade4 dist.binary
#' @importFrom heatmaply heatmaply
#' @importFrom grDevices colorRampPalette
#' @importFrom matrixStats rowQuantiles
#' @importFrom stringr str_remove
#' @importFrom tibble column_to_rownames
#' @importFrom gplots heatmap.2
#'
#' @param scCN, scCN dataframe from Replication_state function
#' @param min_cor, min correlation to be kept
#' @param save_plot_name, if provided the function saves plots instead of returning them
#' @export
#'

FilterCells <- function(scCN,
                        min_cor = 0.25,
                        save_plot_basename = NULL) {
  #load required operators
  `%do%` = foreach::`%do%`
  `%>%` = tidyr::`%>%`

  scCN = scCN %>%
    dplyr::ungroup() %>%
    dplyr::arrange(group, newIndex) %>%
    tidyr::unite(index, c(group, newIndex), sep = ' _ ') %>%
    dplyr::mutate(index = factor(index, levels = unique(index)))
  #build matrix
  mat = scCN %>%
    tidyr::unite(pos, c(chr, start), sep = ':') %>%
    dplyr::mutate(Rep = as.numeric(Rep)) %>%
    dplyr::select(pos, index, Rep) %>%
    tidyr::spread(key = index, value = Rep) %>%
    tibble::column_to_rownames('pos') %>%
    dplyr::filter(complete.cases(.)) %>%
    as.matrix()

  #index
  Index = colnames(mat)
  #correlation similarity distance
  results = 1 - as.matrix(ade4::dist.binary(
    t(mat),
    method = 2,
    diag = T,
    upper = T
  ))

  basenames = dplyr::tibble(Group = stringr::str_remove(colnames(mat), ' _ [0-9]{1,10}$'))%>%
    dplyr::mutate(Group=factor(Group))

  #prepare color patterns
  color = grDevices::colorRampPalette(
    colors = c(
      "#00204DFF",
      "#233E6CFF",
      "#575C6DFF",
      "#7C7B78FF",
      "#A69D75FF",
      "#D3C164FF",
      "#FFEA46FF"
    )
  )

  if (is.null(save_plot_basename)) {
    Plot1 <-
      heatmaply::heatmaply(
        x = results,
        colors = color,
        dendrogram = F,
        showticklabels = F,
        row_side_colors = basenames,
        col_side_colors = basenames,
        limits = c(0, 1)
      ) %>%
      plotly::layout(
        showlegend = FALSE,
        legend = FALSE,
        annotations = list(visible = FALSE)
      )
  } else{
    selcol <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))
    color_basenames = selcol(length(unique(basenames$Group)))

    jpeg(filename = paste0(
      save_plot_basename,
      '_correlation_plot_before_filter.jpeg'
    ),width = 800,height = 800)
    gplots::heatmap.2(
      results,
      trace = "none",
      dendrogram = 'none',
      Colv = F,
      Rowv = F,
      breaks = seq(0, 1, length.out = 101),
      col = color(100),
      density.info =  'density',
      keysize = 1,
      key.title = 'Simple matching coefficient',
      RowSideColors = color_basenames[as.numeric(basenames$Group)],
      ColSideColors = color_basenames[as.numeric(basenames$Group)],
      labRow = FALSE,
      labCol = FALSE
    )
    dev.off()
  }

  to_keep = foreach::foreach(i = unique(basenames$Group)) %do% {
    sub_mat = results[basenames$Group == i, basenames$Group == i]
    diag(sub_mat) = 0
    ! matrixStats::rowQuantiles(x = sub_mat,
                                probs = 0.60,
                                na.rm = T) <= min_cor
  }

  to_keep = unlist(to_keep)
  results_after_filtering = results[to_keep, to_keep]
  basenames_after_filtering = basenames[to_keep, ]

  if (is.null(save_plot_basename)) {
    Plot2 <-
      heatmaply::heatmaply(
        x = results_after_filtering,
        colors = color,
        dendrogram = F,
        showticklabels = F,
        row_side_colors = basenames_after_filtering$Group,
        col_side_colors = basenames_after_filtering$Group,
        limits = c(0, 1)
      ) %>%
      plotly::layout(
        showlegend = FALSE,
        legend = FALSE,
        annotations = list(visible = FALSE)
      )
  } else{
    jpeg(filename = paste0(
      save_plot_basename,
      '_correlation_plot_after_filter.jpeg'
    ),width = 800,height = 800)
    gplots::heatmap.2(
      results_after_filtering,
      trace = "none",
      dendrogram = 'none',
      Colv = F,
      Rowv = F,
      breaks = seq(0, 1, length.out = 101),
      col = color(100),
      density.info =  'density',
      keysize = 1,
      key.title = 'Simple matching coefficient',
      RowSideColors = color_basenames[as.numeric(basenames_after_filtering$Group)],
      ColSideColors = color_basenames[as.numeric(basenames_after_filtering$Group)],
      labRow = FALSE,
      labCol = FALSE
    )
    dev.off()
  }


  #filter out samples that don't correlate and save
  scCN = scCN %>%
    dplyr::filter(index %in% Index) %>%
    tidyr::separate(index, c('group', 'index'), sep = ' _ ')

  rep_percentage = scCN %>%
    dplyr::group_by(Cell, basename, group, index) %>%
    dplyr::summarise(Rep_percentage = mean(Rep))

  #new index
  new_index_list = rep_percentage %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Rep_percentage) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(newIndex = 1:dplyr::n()) %>%
    dplyr::arrange(group, newIndex) %>%
    dplyr::select(oldIndex = index, newIndex, Cell, basename, group)

  scCN = scCN %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(new_index_list,
                      by = c('Cell', 'index' = 'oldIndex', 'basename', 'group')) %>%
    dplyr::select(
      chr,
      start,
      end,
      CN,
      background,
      CN_bg,
      th,
      Rep,
      PercentageReplication,
      Cell,
      basename,
      group,
      newIndex
    )
  if (is.null(save_plot_basename)) {
    return(list = c(
      Before_filtering = list(Plot1),
      After_filter = list(Plot2),
      FilteredData = list(scCN)
    ))
  } else{
    return(scCN)
  }
}


#'  calculates pseudobulk RT
#'
#' @return tibble
#'
#' @importFrom  dplyr filter arrange group_by mutate summarise ungroup select tibble inner_join pull
#' @importFrom  tidyr %>%
#' @importFrom  foreach %do% %:% foreach
#'
#' @param S_scCN, scCN dataframe from Rebin function (after filtering,optional)
#'
#' @export
#'
#'
pseudoBulkRT <- function(S_scCN) {
  #load required operators
  `%:%` = foreach::`%:%`
  `%do%` = foreach::`%do%`
  `%>%` = tidyr::`%>%`

  #calculate replication percentage
  rep_percentage = S_scCN %>%
    dplyr::group_by(Cell, basename, group, newIndex) %>%
    dplyr::summarise(Rep_percentage = mean(Rep))

  #calculate replication timing normalizing each bin by the number of cells in each bin and then calculating the average of the average
  # select symmetrically distributed cells.

  rep_percentage = rep_percentage %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      min_perc = 1 - max(Rep_percentage),
      max_perc = 1 - min(Rep_percentage)
    ) %>%
    dplyr::filter(Rep_percentage >= round(min_perc, 2),
                  Rep_percentage <= round(max_perc, 2)) %>%
    dplyr::mutate(
      min_perc = 1 - max(Rep_percentage),
      max_perc = 1 - min(Rep_percentage)
    ) %>%
    dplyr::filter(Rep_percentage >= round(min_perc, 2),
                  Rep_percentage <= round(max_perc, 2))

  # bin the cells based on their percentage of replication in order have continuous bins with at least one cell.

  RT_binning = foreach::foreach(Group = unique(rep_percentage$group),
                                .combine = 'rbind') %:%
    foreach::foreach(bins = 1:30, .combine = 'rbind') %do% {
      rep_group = rep_percentage %>%
        dplyr::ungroup() %>%
        dplyr::filter(group == Group) %>%
        dplyr::mutate(rep_group = ceiling(100 * Rep_percentage / bins)) %>%
        dplyr::arrange(rep_group) %>%
        dplyr::pull(rep_group) %>%
        unique()

      cont_rep_group = min(rep_group):max(rep_group)

      if (length(cont_rep_group) == length(rep_group)) {
        if (all(rep_group == cont_rep_group)) {
          dplyr::tibble(group = Group,
                        Binning_step = bins)
        } else{
          dplyr::tibble()
        }
      } else{
        dplyr::tibble()

      }
    }

  RT_binning = RT_binning %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(Binning_step = min(Binning_step))

  scRT = S_scCN %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      min_perc = 1 - max(PercentageReplication),
      max_perc = 1 - min(PercentageReplication)
    ) %>%
    dplyr::filter(PercentageReplication >= min_perc,
                  PercentageReplication <= max_perc) %>%
    dplyr::mutate(
      min_perc = 1 - max(PercentageReplication),
      max_perc = 1 - min(PercentageReplication)
    ) %>%
    dplyr::filter(PercentageReplication >= min_perc,
                  PercentageReplication <= max_perc) %>%
    dplyr::inner_join(RT_binning, by = "group") %>%
    dplyr::mutate(RepGroup = (ceiling(
      100 * PercentageReplication / Binning_step
    ))) %>%
    dplyr::group_by(chr, start, end, RepGroup, group) %>%
    dplyr::summarise(Rep = mean(Rep)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(chr, start, end, group) %>%
    dplyr::summarise(RT = mean(Rep)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(RT = (RT - min(RT)) / (max(RT) - min(RT)),
                  basename = group) %>%
    dplyr::select(chr, start, end, group, basename, RT) %>%
    dplyr::ungroup()%>%
    dplyr::arrange(basename)

  return(scRT)
}

#'  prepares variability file
#'
#' @return tibble
#'
#' @importFrom  dplyr group_by mutate summarise ungroup select inner_join
#' @importFrom  tidyr %>%
#'
#' @param S_scCN, scCN dataframe from Replication_state function (after filtering,optional)
#' @param scRT, pseudo-bulk-RT dataframe from pseudoBulkRT
#'
#' @export
#'
#'
Variability <- function(S_scCN, scRT) {
  #load required operators
  `%>%` = tidyr::`%>%`

  variability = S_scCN %>%
    dplyr::select(group, chr, start, end, Rep, PercentageReplication) %>%
    dplyr::inner_join(scRT, by = c("group", "chr", "start", "end"))  %>%
    dplyr::mutate(time = round(10 * (1 - RT - PercentageReplication), 1)) %>%
    dplyr::group_by(group, time, RT, chr, start, end) %>%
    dplyr::summarise(percentage = mean(Rep)) %>%
    dplyr::ungroup()

  return(variability)

}

#'  reverts extractSubpop
#'
#' @return list
#'
#' @importFrom  dplyr filter arrange group_by mutate group_split inner_join left_join
#' @importFrom  tidyr %>%
#' @importFrom stringr str_remove_all
#'
#' @param scCN, S-phase cells scCN dataframe from extractSubpop
#' @param scRT, pseudo-bulk-RT dataframe from extractSubpop
#' @param scVariability, single cell variability dataframe created by extractSubpop
#' @param subpopulation, dataframe containing sub-population info (Cell,basename,group,subpopulation)
#' @param RefRT, Reference RT of a group created by extractSubpop
#'
#' @export
#'
#'

rejoinSubpop = function(scCN,
                        scRT,
                        scVariability,
                        subpopulation,
                        RefRT = NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`

  subpopulation = subpopulation %>%
    dplyr::mutate(group = paste(group, subpopulation, sep = '_'))

  #select group that is changing
  scCN = scCN %>%
    dplyr::left_join(subpopulation, by = c("Cell", "basename", "group"))

  changing = scCN %>%
    dplyr::filter(!is.na(subpopulation))

  scCN = scCN %>%
    dplyr::filter(is.na(subpopulation)) %>%
    dplyr::select(-subpopulation)

  #filter out RTs to change
  scRT = scRT %>%
    dplyr::filter(!group %in% unique(changing$group))

  #filter out variability to change
  scVariability = scVariability %>%
    dplyr::filter(!group %in% unique(changing$group))

  # dedup Ref
  if (!is.null(RefRT)) {
    RefRT = RefRT %>%
      dplyr::left_join(subpopulation %>%
                         dplyr::select(-Cell) %>%
                         unique(),
                       by = c("group"))
    # RT to be deduplicated
    dedupRefRT = RefRT %>%
      dplyr::filter(!is.na(subpopulation)) %>%
      dplyr::mutate(group = stringr::str_remove_all(
        string = group,
        pattern = paste0('_', subpopulation)
      )) %>%
      dplyr::select(-subpopulation) %>%
      unique()

    RefRT = RefRT %>%
      dplyr::filter(is.na(subpopulation)) %>%
      dplyr::select(-subpopulation) %>%
      rbind(dedupRefRT)

  }
  #put back old group
  changing = changing  %>%
    dplyr::mutate(group = stringr::str_remove_all(string = group,
                                                  pattern = paste0('_', subpopulation))) %>%
    dplyr::select(-subpopulation)
  #changing index
  Index = changing %>%
    dplyr::select(Cell, basename, group, PercentageReplication) %>%
    unique() %>%
    dplyr::arrange(PercentageReplication) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(Ind = 1:dplyr::n())
  #assign new index
  changing = changing %>%
    dplyr::inner_join(Index,
                      by = c("PercentageReplication", "Cell", "basename", "group")) %>%
    dplyr::mutate(newIndex = Ind) %>%
    dplyr::select(-Ind)

  #reconstitute df
  scCN = rbind(scCN,
               changing)
  #calculate scRT changing samples
  changingRT = do.call(
    'rbind',
    lapply(
      FUN = Kronos.scRT::pseudoBulkRT,
      X = changing %>%
        dplyr::group_by(group)
      %>% dplyr::group_split()
    )
  )

  scRT = rbind(scRT, changingRT)

  #calculate new variability

  Changing_Var = Kronos.scRT::Variability(S_scCN = changing, scRT = changingRT)

  scVariability = rbind(scVariability, Changing_Var)

  if (is.null(RefRT)) {
    return(list(
      scCN = scCN,
      scRT = scRT,
      scVariability = scVariability
    ))
  } else{
    return(list(
      scCN = scCN,
      scRT = scRT,
      scVariability = scVariability,
      RefRT = RefRT
    ))
  }

}

#'  separates S phase cells into sub-groups and calculates relative pseudo-bulk RT
#'
#' @return list
#'
#' @importFrom  dplyr n filter select group_split inner_join arrange group_by mutate ungroup
#' @importFrom  tidyr %>%
#'
#' @param scCN, S-phase cells scCN dataframe from Rebin function (after filtering,optional)
#' @param scRT, pseudo-bulk-RT dataframe
#' @param scVariability, single cell variability dataframe the Variability function
#' @param subpopulation, dataframe containing sub-population info (Cell,basename,group,subpopulation)
#' @param RefRT, Reference RT dataframe from RebinRT
#'
#' @export
#'

extractSubpop = function(scCN,
                         scRT,
                         scVariability,
                         subpopulation,
                         RefRT = NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`

  #select groups that are changing
  changing = scCN %>%
    dplyr::filter(
      group %in% unique(subpopulation$group) &
        basename %in% unique(subpopulation$basename)
    ) %>%
    dplyr::inner_join(subpopulation, by = c('Cell', 'group', 'basename'))

  scCN = scCN %>%
    dplyr::filter(!(
      group %in% unique(subpopulation$group) &
        basename %in% unique(subpopulation$basename)
    ))
  #filter out RTs to change
  scRT = scRT %>%
    dplyr::filter(!group %in% unique(subpopulation$group))

  #put back old group
  changing = changing %>%
    dplyr::mutate(group = paste(group, subpopulation, sep = '_')) %>%
    dplyr::select(-subpopulation)
  #newIndex
  Index = changing %>%
    dplyr::select(Cell, basename, group, PercentageReplication) %>%
    unique() %>%
    dplyr::arrange(PercentageReplication) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(Ind = 1:dplyr::n())
  #assign new index
  changing = changing %>%
    dplyr::inner_join(Index,
                      by = c("PercentageReplication", "Cell", "basename", "group")) %>%
    dplyr::mutate(newIndex = Ind) %>%
    dplyr::select(-Ind) %>%
    dplyr::ungroup()

  #reconstitute df
  scCN = rbind(scCN,
               changing)
  #calculate scRT changing samples

  changingRT = do.call(
    'rbind',
    lapply(
      FUN = Kronos.scRT::pseudoBulkRT,
      X = changing %>%
        dplyr::group_by(group)
      %>% dplyr::group_split()
    )
  )

  scRT = rbind(scRT,
               changingRT)

  #variability

  scVariability = scVariability %>%
    dplyr::filter(!group %in% unique(subpopulation$group))

  Changing_Var = Kronos.scRT::Variability(S_scCN = changing, scRT = changingRT)

  scVariability = rbind(scVariability, Changing_Var)



  if (is.null(RefRT)) {
    return(list(
      scCN = scCN,
      scRT = scRT,
      scVariability = scVariability
    ))
  } else{
    dup_RefRT = RefRT %>%
      dplyr::inner_join(subpopulation %>%
                          dplyr::select(-Cell, -basename) %>%
                          unique(),
                        by = c("group")) %>%
      dplyr::mutate(group = paste0(group, '_', subpopulation)) %>%
      dplyr::select(-subpopulation)

    RefRT = RefRT %>%
      dplyr::filter(!group %in% unique(subpopulation$group))

    RefRT = rbind(RefRT, dup_RefRT)


    return(list(
      scCN = scCN,
      scRT = scRT,
      scVariability = scVariability,
      RefRT = RefRT
    ))
  }
}
