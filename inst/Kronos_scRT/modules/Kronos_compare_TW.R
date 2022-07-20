#options
options(dplyr.summarise.inform = FALSE,
        scipen = 999)

#parse input
option_list = list(
  optparse::make_option(
    c("-F", "--file"),
    type = "character",
    default = NULL,
    help = "Variability file with groups produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
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
    c("-f", "--output_file_base_name"),
    type = "character",
    default = "out",
    help = "Base name for the output file [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("-A", "--Annotation_file"),
    type = "character",
    help = "Annotatation file",
    metavar = "character"
  ),
  optparse::make_option(
    c("-G", "--RT_Groups"),
    type = "numeric",
    default = 3,
    help = "if an Annotatation file is not provided, TW is calculated over RT groups. The numeber of groups can be set here: 1,2,3 or 5 [default= %default]",
    metavar = "numeric"
  ),
  optparse::make_option(
    c("-C", "--cores"),
    default = 3,
    help = "Number of for bootstrapping [default= %default]",
    metavar = "numeric",
    type = "numeric"
  ),
  optparse::make_option(
    c("-c", "--color"),
    default = 'red',
    help = "Color Plot, multiple colors can be given usins a comma [default= %default]",
    metavar = "character",
    type = "character"
  ),
  optparse::make_option(
    c("-p", "--pval"),
    type = "logical",
    action = "store_true",
    default = F,
    help = "Bootstrap pValue for the difference in TW between groups (it works only with one annotation) [default= %default]",
    metavar = "logical"
  ),
  optparse::make_option(
    c("-T", "--pairs_to_test"),
    type = "character",
    help = "Pairs of groups for which to calculate the pvalue. Groups in a pair have to be separated by a comma while pairs are separated by a semicolon e.g. A,B;A,C",
    metavar = "character"
  ),
  optparse::make_option(
    c("-H", "--pval_alternative_hypotesis"),
    default = 'two.sided',
    help = "greater, lower, two.sided. This option is active only if the option G is provided [default= %default]",
    metavar = "character",
    type = "character"
  ),
  optparse::make_option(
    c("-Μ", "--adjust_methods"),
    default = 'none',
    help = " correction method: holm, hochberg, hommel, bonferroni, BY (Benjamini & Yekutieli ),fdr (false discovery rate), none. [default= %default]",
    metavar = "character",
    type = "character"
  ),
  optparse::make_option(
    c("-i", "--number_of_iterations"),
    default = 10000,
    help = "Number of iterations to calculate bootstrap [default= %default]",
    metavar = "numeric",
    type = "numeric"
  )
)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#set plotting theme
ggplot2::theme_set(ggplot2::theme_bw())

#check inputs
if (!'file' %in% names(opt)) {
  stop("Variability file must be provided. See script usage (--help)")
}

if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

opt$file = stringr::str_split(opt$file, ',')[[1]]
#check file types
if ('Annotation_file' %in% names(opt)) {
  if (!Kronos.scRT::right_format(
    file_path = opt$Annotation_file,
    columns_to_check = c('chr', 'start', 'end', 'annotation'),
    delim = '\t',
    logical = T
  )) {
    stop(
      'Provided annotation file is not a tab delemeted file containing the following columns:chr, start, end, annotation'
    )
  } else{
    Annotation = readr::read_tsv(file = opt$Annotation_file)
  }
}

.tmp = lapply(opt$file, function(x)
  if (!Kronos.scRT::right_format(
    file_path = x,
    columns_to_check = c('group', 'time', 'RT', 'chr', 'start', 'end', 'percentage'),
    delim = '\t',
    logical = T
  )) {
    stop(
      'Provided annotation file is not a tab delemeted file containing the following columns:group, time, RT, chr, start,end, percentage'
    )
  })


#load files
data <-
  Kronos.scRT::load_multiple_df(
    dirs = opt$file,
    delim = '\t',
    col_types = readr::cols(chr = 'c')
  )


#prepare variability file
if ('Annotation_file' %in% names(opt)) {
  variability = Kronos.scRT::TW_GenomeAnnotation(Variability = data,
                                                 GenomeAnnotation = Annotation)

} else{
  variability = Kronos.scRT::TW_RTAnnotation(Variability = data,
                                             RT_Groups = opt$RT_Groups)
}

MaxNcors = length(unique(variability$category)) * length(unique(variability$group))

#fit data
twidth_fitted_data = Kronos.scRT::Twidth_fit_data(df = variability,
                                                  ncores = ifelse(MaxNcors > opt$cores,
                                                                  opt$cores,
                                                                  MaxNcors))
#calculate TW
twidth = Kronos.scRT::Twidth(twidth_fitted_data)

groups = unique(variability$group)

#recover colors
opt$color = stringr::str_split(opt$color, ',', simplify = T)

if (length(opt$color) < length(groups)) {
  opt$color = rep(opt$color[1], length(groups))
}
#extended plot
plot1 = lapply(1:length(groups), function(x)
  Kronos.scRT::Twidth_extended_plot(
    Variability = variability[variability$group == groups[x], ],
    Fitted_data = twidth_fitted_data[twidth_fitted_data$group == groups[x], ],
    Twidth = twidth[twidth$group == groups[x], ],
    Color = opt$color[x]
  ))

#barplot
plot2 = lapply(1:length(groups), function(x)
  Kronos.scRT::Twidth_barplot(
    Variability = variability[variability$group == groups[x], ],
    Twidth = twidth[twidth$group == groups[x], ],
    Color = opt$color[x]
  ))



.tmp = lapply(1:length(groups), function(x)
  suppressMessages(
    ggplot2::ggsave(
      plot =  plot1[[x]],
      filename = paste0(file.path(opt$out,
                                  groups[x]),
                        '_extended_Twidths.pdf'),
      width = 3 * length(unique(variability$category)),
      height = 5,
      device = grDevices::cairo_pdf
    )
  ))
.tmp = lapply(1:length(groups), function(x)
  suppressMessages(
    ggplot2::ggsave(
      plot2[[x]],
      filename = paste0(file.path(opt$out,
                                  groups[x]),
                        '_barplot_Twidths.pdf'),
      width = 1.4 * length(unique(variability$category)),
      height = 6
    )
  ))

readr::write_tsv(x = twidth, file = paste0(file.path(opt$out,
                                                     opt$output_file_base_name),
                                           '_Twidth.tsv'))

# pval option
if (opt$pval) {
  #if provided reshape pairs_to_test
  if (!is.null(opt$pairs_to_test)) {
    opt$pairs_to_test = stringr::str_split(
      string =
        stringr::str_split(
          string = opt$pairs_to_test,
          pattern = ';',
          simplify = T
        )
      ,
      pattern = ',',
      simplify = T
    )

    opt$pairs_to_test = dplyr::tibble(Category1 = opt$pairs_to_test[, 1],
                                      Category2 = opt$pairs_to_test[, 2])

  }

  #calculate pval
  pval = Kronos.scRT::Twidth_pval(
    variability = variability,
    twidth = twidth,
    pairs_to_test = opt$pairs_to_test,
    adjust.methods = opt$adjust_methods,
    alternative = opt$pval_alternative_hypotesis,
    nIterations = opt$number_of_iterations,
    ncores = opt$cores
  )

  #create output folder
  if (!dir.exists(file.path(opt$out, 'pval'))) {
    dir.create(file.path(opt$out, 'pval'))
  }

  readr::write_tsv(x = pval,
                   file = file.path(opt$out, 'pval', 'Pvalues.tsv'))

  plot3 = lapply(1:length(groups), function(x)
    Kronos.scRT::Twidth_barplot(
      Variability = variability[variability$group == groups[x],],
      Twidth = twidth[twidth$group == groups[x],],
      Color = opt$color[x],
      pval = pval[pval$group == groups[x],]
    ))

  .tmp = lapply(1:length(groups), function(x)
    suppressMessages(
      ggplot2::ggsave(
        plot3[[x]],
        filename = paste0(
          file.path(opt$out,
                    'pval',
                    groups[x]),
          '_barplot_Twidths_with_pval.pdf'
        ),
        width = 1.4 * length(unique(variability$category)),
        height = 6
      )
    ))


}

print('done')
