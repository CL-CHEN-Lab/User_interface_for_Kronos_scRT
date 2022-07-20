#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

#parse input
option_list = list(
  optparse::make_option(
    c("-C", "--scCNV"),
    type = "character",
    default = NULL,
    help = "*single_cells_CNV* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--CNV_values"),
    type = "character",
    default = "B",
    help = "What type of date to plot for the single cell traks: ('B'=Binarized, 'CNV'=Copy number variation, 'log2'=log2(CNV_Cell/CNV_mean_G1/G2_cells)) [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("--per_Chr"),
    type = 'logical',
    default = F,
    action = 'store_TRUE',
    help = "Calculate TSNE/UMAP on each chromosome",
    metavar = "logical"
  ),
  optparse::make_option(
    c("-o", "--out"),
    type = "character",
    default = "output",
    help = "Output directory [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("-f", "--output_file_base_name"),
    type = "character",
    default = "out",
    help = "Base name for the output file [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("-c", "--cores"),
    type = "integer",
    default = 3,
    help = "Number of cores to use [default= %default]",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = as.integer(Sys.Date()),
    help = "Set seed for reproducibility (optional).",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-U", "--UMAP"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Skip t-SNE, only plot UMAP.",
    metavar = "logical"
  ),
  optparse::make_option(
    c("-T", "--TSNE"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Skip UMAP, only plot t-SNE.",
    metavar = "logical"
  ),
  optparse::make_option(
    c("--chr_prefix"),
    type = "character",
    action = 'store',
    help = "Chromosome prefix, if there is no prefix use none [default= %default]",
    default = "chr",
    metavar = "character"
  ),
  optparse::make_option(
    c("--chr_range"),
    type = "character",
    action = 'store',
    help = "Chromosomes to consider in the analysis (example 1:5,8,15:18,X) [default= %default]",
    default = "1:22",
    metavar = "character"
  )

)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#define operators
`%>%` = tidyr::`%>%`
`%dopar%` = foreach::`%dopar%`
`%do%` = foreach::`%do%`


#output dir
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = T)
}

#set plotting theme
ggplot2::theme_set(ggplot2::theme_bw())

#Set seed for reproducibility
set.seed(opt$seed)

#check files
if (!'scCNV' %in% names(opt)) {
  stop("scCNV file must be provided. See script usage (--help)")
}

opt$scCNV = stringr::str_split(opt$scCNV, ',')[[1]]

if ('order' %in% names(opt)) {
  opt$order = stringr::str_split(opt$order, ',')[[1]]

}

#check format
if (!all(sapply(opt$scCNV, function(x)
  Kronos.scRT::right_format(
    file_path = x,
    logical = T,
    columns_to_check = c(
      'chr',
      'start',
      'end',
      'CN',
      'background',
      'CN_bg',
      'th',
      'Rep',
      'PercentageReplication',
      'Cell',
      'basename',
      'group',
      'newIndex'
    )
  )))) {
  stop(
    'Provided scCNV file is not a tab delemeted table with the following columns: chr, start, end, CN, background, CN_bg, th, Rep, PercentageReplication, Cell, basename, group, newIndex'
  )
}

#load CNV files
scCNV = Kronos.scRT::load_multiple_df(opt$scCNV,
                                      delim = '\t',
                                      col_types = readr::cols(chr = 'c'))

if ('order' %in% names(opt)) {
  scCNV = scCNV %>%
    dplyr::mutate(group = factor(group, levels = opt$order))
}


# select chrs of interest

chr_list = paste0(
  ifelse(opt$chr_prefix == 'none', '', opt$chr_prefix),
  unlist(
    Kronos.scRT::String_to_Range(x = stringr::str_split(opt$chr_range, ',')[[1]], div = ':')
  )
)

#filter chr and convert into factor
scCNV = scCNV %>%
  dplyr::filter(chr %in% chr_list) %>%
  dplyr::mutate(chr = factor(x =  chr, levels = chr_list)) %>%
  tidyr::drop_na()


if (!opt$CNV_values %in% c('B', 'CNV', 'Log2', 'all')) {
  opt$CNV_values = 'B'
  warning('unrecognized CNV value to plot. Binarized tracks restored')
}

#store chr info
chr = unique(scCNV$chr)

scCNV = scCNV %>%
  dplyr::mutate(pos = paste0(chr, ':', start, '-', end)) %>%
  dplyr::select(
    'pos',
    'Cell',
    'PercentageReplication',
    'basename',
    'group',
    'data' = dplyr::case_when(
      opt$CNV_values == 'B' ~ 'Rep',
      opt$CNV_values == 'CNV' ~ 'CN',
      opt$CNV_values == 'Log2' ~ 'CN_bg'
    )
  ) %>%
  dplyr::mutate(data = as.numeric(data)) %>%
  tidyr::spread(pos, data)

#remove NA columns
scCNV = scCNV[, colSums(is.na(scCNV)) == 0]

mat = scCNV[, -c(1:4)]
scCNV = scCNV[, c(1:4)]

#start cluster if per chr option is selected
if (opt$per_Chr) {
  cl = snow::makeCluster(opt$cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))
}

if (opt$CNV_values == 'B') {
  if (opt$per_Chr) {
    # calculate distance per each chr
    results = foreach::foreach(C = chr, .inorder = T) %dopar% {
      as.matrix(ade4::dist.binary(mat[, grepl(pattern = paste0(C, ':'), x = colnames(mat))], method = 2))
    }

    names(results) = chr

  } else{
    results = list()
    results[['all Chr']] = as.matrix(ade4::dist.binary(mat, method = 2))
  }
} else{
  if (opt$per_Chr) {
    results = foreach::foreach(C = chr, .inorder = T) %dopar% {
      as.matrix(mat[, grepl(pattern = paste0(C, ':'), x = colnames(mat))])
    }
    names(results) = chr

  } else{
    results = list()
    results[['all Chr']] = as.matrix(mat)
  }

}

#free mem
rm('mat')

scCNV_umap = scCNV
# TSNE
if (!opt$UMAP) {
  scCNV = foreach::foreach(C = names(results), .combine = 'rbind') %do% {
    Perplex = ceiling(nrow(results[[C]]) / 50)
    tsne <-
      Rtsne::Rtsne(
        X = results[[C]],
        dims = 2,
        perplexity = ifelse(Perplex < 10, 10, Perplex),
        check_duplicates = F,
        theta = 0.25,
        is_distance = opt$CNV_values == 'B',
        verbose = F,
        max_iter = 5000,
        num_threads = opt$cores,
        partial_pca = T
      )

    scCNV %>% dplyr::mutate(Chr = C,
                            x = tsne$Y[, 1],
                            y = tsne$Y[, 2])
  }
  readr::write_tsv(scCNV, paste0(file.path(opt$out,
                                           opt$output_file_base_name),
                                 '_tsne.txt'))
  X = foreach::foreach(C = names(results), .combine = 'rbind') %do% {
    plot = scCNV %>%
      dplyr::filter(Chr == C) %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x, y, color = group, shape = group),
                          alpha = 0.4,
                          size = 2) + ggplot2::xlab('TSNE1') + ggplot2::ylab('TSNE2') + ggplot2::facet_wrap(~ Chr)

    suppressMessages(ggplot2::ggsave(
      plot = plot,
      dpi = 300,
      filename = paste0(
        file.path(opt$out,
                  opt$output_file_base_name),
        '_',
        stringr::str_replace(string = C,pattern = ' ',replacement = '_'),
        '_tsne_color_by_group.pdf'
      )
    ))

    plot = scCNV %>%
      dplyr::filter(Chr == C) %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(x, y, color = basename, shape = group),
        alpha = 0.4,
        size = 2
      ) + ggplot2::xlab('TSNE1') + ggplot2::ylab('TSNE2') + ggplot2::facet_wrap( ~ Chr)

    suppressMessages(ggplot2::ggsave(
      plot = plot,
      dpi = 300,
      filename = paste0(
        file.path(opt$out,
                  opt$output_file_base_name),
        '_',
        stringr::str_replace(string = C,pattern = ' ',replacement = '_'),
        '_tsne_color_by_basename.pdf'
      )
    ))

    plot = scCNV %>%
      dplyr::filter(Chr == C) %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(x, y, color = PercentageReplication, shape = group),
        alpha = 0.4,
        size = 2
      ) + ggplot2::scale_color_gradient2(
        low = "#FFEA46FF",
        mid = "#7C7B78FF",
        high = "#00204DFF",
        lim = c(0, 1),
        midpoint = 0.5
      ) + ggplot2::xlab('TSNE1') + ggplot2::ylab('TSNE2') + ggplot2::facet_wrap( ~ Chr)
    suppressMessages(ggplot2::ggsave(
      plot = plot,
      dpi = 300,
      filename = paste0(
        file.path(opt$out,
                  opt$output_file_base_name),
        '_',
        stringr::str_replace(string = C,pattern = ' ',replacement = '_'),
        '_tsne_color_by_rep_percentage.pdf'
      )
    ))
    C
  }
}

# UMAP
if (!opt$TSNE) {
  scCNV_umap = foreach::foreach(C = names(results),
                                .combine = 'rbind') %do% {
                                  if (opt$CNV_values == 'B') {
                                    input_mat = 'dist'
                                  } else{
                                    input_mat = "data"
                                  }
                                  umap <-
                                    umap::umap(d = results[[C]],
                                               input = input_mat,
                                               random_state = opt$seed)
                                  scCNV_umap %>% dplyr::mutate(Chr = C,
                                                               x = umap$layout[, 1],
                                                               y = umap$layout[, 2])
                                }
  readr::write_tsv(scCNV_umap,
                   paste0(file.path(opt$out,
                                    opt$output_file_base_name),
                          '_umap.txt'))

  X_umap = foreach::foreach(C = names(results),
                            .combine = 'rbind') %do% {
                              plot = scCNV_umap %>%
                                dplyr::filter(Chr == C) %>%
                                ggplot2::ggplot() +
                                ggplot2::geom_point(ggplot2::aes(x, y, color = group, shape = group),
                                                    alpha = 0.4,
                                                    size = 2) + ggplot2::xlab('UMAP1') + ggplot2::ylab('UMAP2') + ggplot2::facet_wrap(~ Chr)

                              suppressMessages(ggplot2::ggsave(
                                plot = plot,
                                dpi = 300,
                                filename = paste0(
                                  file.path(opt$out,
                                            opt$output_file_base_name),
                                  '_',
                                  stringr::str_replace(string = C,pattern = ' ',replacement = '_'),
                                  '_umap_color_by_group.pdf'
                                )
                              ))

                              plot = scCNV_umap %>%
                                dplyr::filter(Chr == C) %>%
                                ggplot2::ggplot() +
                                ggplot2::geom_point(
                                  ggplot2::aes(x, y, color = basename, shape = group),
                                  alpha = 0.4,
                                  size = 2
                                ) + ggplot2::xlab('UMAP1') + ggplot2::ylab('UMAP2') + ggplot2::facet_wrap( ~ Chr)

                              suppressMessages(ggplot2::ggsave(
                                plot = plot,
                                dpi = 300,
                                filename = paste0(
                                  file.path(opt$out,
                                            opt$output_file_base_name),
                                  '_',
                                  stringr::str_replace(string = C,pattern = ' ',replacement = '_'),
                                  '_umap_color_by_basename.pdf'
                                )
                              ))

                              plot = scCNV_umap %>%
                                dplyr::filter(Chr == C) %>%
                                ggplot2::ggplot() +
                                ggplot2::geom_point(
                                  ggplot2::aes(x, y, color = PercentageReplication, shape = group),
                                  alpha = 0.4,
                                  size = 2
                                ) + ggplot2::scale_color_gradient2(
                                  low = "#FFEA46FF",
                                  mid = "#7C7B78FF",
                                  high = "#00204DFF",
                                  lim = c(0, 1),
                                  midpoint = 0.5
                                ) + ggplot2::xlab('UMAP1') + ggplot2::ylab('UMAP2') + ggplot2::facet_wrap(~ Chr)

                              suppressMessages(ggplot2::ggsave(
                                plot = plot,
                                dpi = 300,
                                filename = paste0(
                                  file.path(opt$out,
                                            opt$output_file_base_name),
                                  '_',
                                  stringr::str_replace(string = C,pattern = ' ',replacement = '_'),
                                  '_umap_color_by_rep_percentage.pdf'
                                )
                              ))
                              C
                            }
}

print('done')
