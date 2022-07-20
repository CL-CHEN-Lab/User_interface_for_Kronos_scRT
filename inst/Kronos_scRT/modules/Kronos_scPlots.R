#options
options(dplyr.summarise.inform = FALSE,
        scipen = 999)

#parse input
option_list = list(
  optparse::make_option(
    c("-L", "--List"),
    type = "character",
    default = NULL,
    help = "A Tab separated file containing in each colum scRT_Tracks and scCNV files paths. Alternative to -R and -C options.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-R", "--scRT_Tracks"),
    type = "character",
    default = NULL,
    help = "*calculated_replication_timing* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.  Alternative to -L option.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-C", "--scCNV"),
    type = "character",
    default = NULL,
    help = "*single_cells_CNV* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.  Alternative to -L option.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-E", "--extra_RT_tracks"),
    type = "character",
    default = NULL,
    help = "A reference RT track.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-s", "--order"),
    type = "character",
    default = NULL,
    help = "Groups separated by a comma in the desired order for plotting.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--CNV_values"),
    type = "character",
    default = "scRT",
    help = "What type of date to plot for the single cell traks: ('scRT'=Binarized, 'S_scCN'=Copy number variation, 'Norm.S_scCN'=log2(CNV_Cell/CNV_mean_G1/G2_cells)) [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("-r", "--region"),
    type = "character",
    default = NULL,
    help = "Region to plot  chr:start-end (multiple regins can be separated by a comma) or provided as a bed file",
    metavar = "character"
  ),
  optparse::make_option(
    c("-o", "--out"),
    type = "character",
    default = "output",
    help = "Output directory [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("--rasterazi_heatmap"),
    type = "logical",
    default = F,
    action = 'store_true',
    help = "Base name for the output file [default= %default]",
    metavar = "logical"
  ),
  optparse::make_option(
    c("--HeatmapColors"),
    type = "character",
    default = NULL,
    help = "min, midpoint and max colors separated by a comma",
    metavar = "character"
  )
)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#set plotting theme
ggplot2::theme_set(ggplot2::theme_bw())

#load operators
`%>%` = tidyr::`%>%`
`%do%` = foreach::`%do%`
`%:%` = foreach::`%:%`

#load files
if ('List' %in% names(opt)) {
  opt$List = tryCatch(
    expr = readr::read_tsv(
      opt$List,
      col_names = F,
      col_types = readr::cols()
    ),
    error = function(e) {
      stop('Settings file does not exitst. See script usage (--help)')
    }
  )

  opt$scRT_Tracks = opt$List$X1
  opt$scCNV = opt$List$X2

} else{
  if (!'scRT_Tracks' %in% names(opt)) {
    stop("scRT_Tracks file must be provided. See script usage (--help)")
  }
  if (!'scCNV' %in% names(opt)) {
    stop("scCNV file must be provided. See script usage (--help)")
  }

  opt$scRT_Tracks = stringr::str_split(opt$scRT_Tracks, ',')[[1]]
  opt$scCNV = stringr::str_split(opt$scCNV, ',')[[1]]

}

if (!'region' %in% names(opt)) {
  stop("No region of interest has been selected. See script usage (--help)")
}



#load files
scRT = Kronos.scRT::load_multiple_df(opt$scRT_Tracks, col_types = readr::cols(chr =
                                                                                'c'))
scCNV = Kronos.scRT::load_multiple_df(opt$scCNV, col_types = readr::cols(chr =
                                                                           'c'))

# assign order if any
if ('order' %in% names(opt)) {
  opt$order = str_split(opt$order, ',')[[1]]

  scRT = scRT %>%
    dplyr::mutate(group = factor(group, levels = opt$order))

  scCNV = scCNV %>%
    dplyr::mutate(group = factor(group, levels = opt$order))

}

#create directory
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = T)
}

#load reference if provided
if ('extra_RT_tracks' %in% names(opt)) {
  opt$extra_RT_tracks = stringr::str_split(opt$extra_RT_tracks, ',')[[1]]

  extra_RT_tracks = Kronos.scRT::load_multiple_df(opt$extra_RT_tracks, col_types = readr::cols(chr =
                                                                                                 'c'))

}


if (file.exists(opt$region)) {
  #load bed file if exist
  opt$region = read_tsv(
    opt$region,
    col_names = c('chr', 'start', 'end'),
    col_types = cols(chr = 'c')
  ) %>%
    mutate(
      n_0_start = str_length(str_extract(start, '0{1,10}$')),
      n_0_end = str_length(str_extract(end, '0{1,10}$')),
      unit = min(
        factor(
          case_when(
            is.na(n_0_start) ~ 'bp',
            n_0_start < 3 ~ 'bp',
            n_0_start < 6 ~ 'Kb',
            n_0_start >= 6 ~ 'Mp'
          ),
          levels = c('bp', 'Kb', 'Mp'),
          ordered = TRUE
        ),
        factor(
          case_when(
            is.na(n_0_end) ~ 'bp',
            n_0_end < 3 ~ 'bp',
            n_0_end < 6 ~ 'Kb',
            n_0_end >= 6 ~ 'Mp'
          ),
          levels = c('bp', 'Kb', 'Mp'),
          ordered = TRUE
        )
      ),
      name_reg = paste(
        chr,
        case_when(
          unit == 'bp' ~ paste0(start, unit, '_', end, unit),
          unit == 'Kb' ~ paste0(start / 10 ^ 3, unit, '_', end /
                                  10 ^ 3, unit),
          unit == 'Mb' ~ paste0(start / 10 ^ 6, unit, '_', end /
                                  10 ^ 6, unit)
        ),
        sep = '_'
      )
    ) %>%
    dplyr::select(-unit, -n_0_start, -n_0_end)
} else{
  #reshape regions
  opt$region = dplyr::tibble(coord = stringr::str_split(opt$region, pattern = ',')[[1]]) %>%
    dplyr::mutate(name_reg = stringr::str_replace_all(coord,
                                                      pattern = '[-:]',
                                                      replacement = '_')) %>%
    tidyr::separate(coord, c('chr', 'pos'), ':') %>%
    tidyr::separate(pos, c('start', 'end'), '-') %>%
    dplyr::mutate(
      start_unit = stringr::str_extract(start, pattern = '.{2}$'),
      start = as.numeric(stringr::str_remove(start, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]")) * dplyr::case_when(
        grepl(x = start_unit, pattern =  '[Kk][Bb]') ~ 1000,
        grepl(x = start_unit, pattern = '[Mm][Bb]') ~ 1000000,
        grepl(x = start_unit, pattern = '[Bp][Pp]') ~ 1,
        grepl(x = start_unit, pattern =  '[0-9][0-9]') ~ 1,
        is.na(start_unit) ~ 1
      ),
      end_unit = stringr::str_extract(end, pattern = '.{2}$'),

      end = as.numeric(stringr::str_remove(end, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]")) * dplyr::case_when(
        grepl(x = end_unit, pattern =  '[Kk][Bb]') ~ 1000,
        grepl(x = end_unit, pattern = '[Mm][Bb]') ~ 1000000,
        grepl(x = end_unit, pattern = '[Bp][Pp]') ~ 1,
        grepl(x = end_unit, pattern =  '[0-9][0-9]') ~ 1,
        is.na(end_unit) ~ 1
      )
    ) %>%
    dplyr::select(-start_unit, -end_unit)
}



if (!opt$CNV_values %in% c('scRT', 'S_scCN', 'Norm.S_scCN')) {
  opt$CNV_values = 'scRT'
  warning('unrecognized CNV value to plot. Binarized tracks restored')
}

x = foreach::foreach(i = 1:nrow(opt$region)) %:%
  foreach::foreach(bn = unique(scRT$basename)) %do% {
    Coordinates = list(
      chr = opt$region$chr[i],
      start = opt$region$start[i],
      end = opt$region$end[i]
    )

    gorup = unique(scRT$group[scRT$basename == bn])


    if (exists('extra_RT_tracks')) {
      pseudoBulkRT = rbind(
        scRT %>% dplyr::filter(basename == bn),
        extra_RT_tracks %>% dplyr::filter(group == gorup)
      )
    } else{
      pseudoBulkRT = scRT %>% dplyr::filter(basename == bn)
    }


    plot = Kronos.scRT::scRTplot(
      pseudoBulkRT = pseudoBulkRT,
      S_scCN = scCNV,
      Coordinates = Coordinates,
      Plot = opt$CNV_values,
      rasterized_heatmap = opt$rasterazi_heatmap ,
      heatmap_colors = opt$HeatmapColors
    )

    suppressMessages(ggplot2::ggsave(
      plot = plot,
      filename = paste0(
        file.path(opt$out,
                  bn),
        '_scPlot_',
        opt$region$name_reg[i],
        '.pdf'
      )
    ))

  }
print('done')
