#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)
#parse input
option_list = list(
  optparse::make_option(
    c("-K", "--Kronos_conf_file"),
    type = "character",
    default = NULL,
    help = "Kronos setting file. If provided -R,-N,-G are ignored. Tab file containing: Reference <TAB> ReferenceBasename <TAB> ReferenceGroup",
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
    c("-R", "--referenceRT"),
    type = "character",
    default = NULL,
    help = "Reference RT min=Late, max=Early. To provide multiple files separate",
    metavar = "character"
  ),
  optparse::make_option(
    c("-N", "--ref_name"),
    type = "character",
    default = "Reference",
    help = "Name for the reference track. If multiple RT were provided, relative groups can be separated using a comma. [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c('-G', "--ref_group"),
    type = "character",
    help = "Group to which this a Reference RT belongs. If multiple RT were provided, relative groups can be separated using a comma.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-B", "--bins"),
    type = "character",
    default = NULL,
    help = "bins RDS file created by Kronos RT",
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
    c("-c", "--cores"),
    type = "integer",
    default = 1,
    help = "Number of parallel jobs to run [default= %default]",
    metavar = "integer"
  )
)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#set plotting theme
ggplot2::theme_set(ggplot2::theme_bw())

#define operators
`%do%` = foreach::`%do%`
`%>%` = tidyr::`%>%`
`%dopar%` = foreach::`%dopar%`

#check inputs
if ('Kronos_conf_file' %in% names(opt)) {
  if (!file.exists(opt$Kronos_conf_file)) {
    stop('Provided setting file does not exist')
  }

  settings = tryCatch(
    expr = readr::read_tsv(
      opt$Kronos_conf_file,
      col_names = c('Reference', 'basename', 'groups'),
      col_types = readr::cols()
    ) %>%
      dplyr::mutate(
        basename = ifelse(
          is.na(basename),
          paste0('exp', dplyr::row_number()),
          basename
        ),
        groups = ifelse(is.na(groups), basename, groups)
      ),
    error = function(e) {
      stop('Settings file does not exitst. See script usage (--help)')
    },
    warning = function(w) {
      tmp = suppressWarnings(
        readr::read_tsv(
          opt$Kronos_conf_file,
          col_names = c('Reference',  'basename', 'groups'),
          col_types = readr::cols()
        ) %>%
          tidyr::mutate(
            basename = ifelse(
              is.na(basename),
              paste0('exp', tidyr::row_number()),
              basename
            ),
            groups = ifelse(is.na(groups), basename, groups)
          )
      )
      warning('missing basenames and groups were replaced with default parameters')
      return(tmp)
    }
  )

  #reformat files
  opt$referenceRT = settings$Reference

  opt$ref_name = settings$basename

  opt$ref_group = settings$groups

} else{

  if (!'referenceRT' %in% names(opt)) {
    stop("Reference file or Kronos setting file must be provided. See script usage (--help)")
  }

  if (!'ref_group' %in% names(opt)) {
    stop(
      "ref_group is a mandatory input. See script usage (--help)"
    )
  }

  #reformat files
  opt$referenceRT = stringr::str_split(opt$referenceRT, ',')[[1]]

  opt$ref_name = stringr::str_split(opt$ref_name, ',')[[1]]

  opt$ref_group = stringr::str_split(opt$ref_group, ',')[[1]]

}

#check uboyts
if(length(opt$ref_group)!=length(opt$referenceRT)) {
  stop('The number of reference files does not math the number of groups.')
}

if (length(opt$ref_name) != length(opt$referenceRT)) {
  if (length(opt$ref_name) == 1) {
    warning(paste(opt$ref_name, 'was used as basename for all the samples'))
    opt$ref_name = rep(opt$ref_name, length(opt$referenceRT))
  } else{
    stop('The number of reference files does not math the number of basenames')

  }
}

#check reference file
  results = paste(
    Kronos.scRT::right_format(
      file_path = opt$referenceRT,
      delim = '\t',
      columns_to_check = 4,
      wrong_message = paste(
        opt$referenceRT,
        ',provided as a reference RT file, does not have the right format'
      )
    ),
    collapse = ''
  )

#create directory
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = T)
}
  #load bins
  bins=readRDS(opt$bins)

  cl=snow::makeCluster(opt$cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))

  # load and rebin reference RT
  Reference_RT <-foreach::foreach(x=1:length(opt$referenceRT),.combine = 'rbind')%dopar%{
      Kronos.scRT::RebinRT(
        RT = readr::read_delim(
          opt$referenceRT[x],
          delim = '\t',
          col_names = c('chr', 'start', 'end', 'RT'),
          col_types = readr::cols(chr = 'c')
        ) %>%
          tidyr::drop_na() ,
        Bins = bins$bins,
        Basename = opt$ref_name[x],
        Group = opt$ref_group[x]
      )
}

  #write output
  readr::write_delim(
    x = Reference_RT ,
    file = paste0(
      file.path(
      opt$out,
      opt$output_file_base_name),
      '_reference_replication_timing_',
      bins$bs,
      '.tsv'
    ),
    delim = '\t',
    col_names = T
  )



#plot profile binning
plot = Reference_RT %>%
  ggplot2::ggplot(ggplot2::aes(RT, color = group)) +
  ggplot2::geom_density(ggplot2::aes(y = ..scaled..)) +
  ggplot2::scale_x_continuous(labels = scales::percent) +
  ggplot2::xlab('Reference Replication Timing') +
  ggplot2::ylab('density') + ggplot2::coord_cartesian(xlim = c(0, 1))

suppressMessages(ggplot2::ggsave(plot,
                                 filename = file.path(
                                   opt$out,
                                   paste0(
                                     opt$output_file_base_name,
                                     '_reference_RT_distribution.pdf'
                                   )
                                 )))

print('done')
