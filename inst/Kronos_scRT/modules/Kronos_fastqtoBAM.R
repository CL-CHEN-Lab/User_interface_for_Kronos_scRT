#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

#parse input
option_list = list(
  optparse::make_option(
    c("-O", "--one"),
    type = "character",
    default = NULL,
    help = "Paths to fastq files (R1 files for PE) separated by a comma",
    metavar = "character"
  ),
  optparse::make_option(
    c("-T", "--two"),
    type = "character",
    default = NULL,
    help = "Paths to R2 fastq files for PE separated by a comma.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-i", "--index"),
    type = "character",
    action = 'store',
    help = "Bowtie 2 index",
    metavar = "character"
  ),
  optparse::make_option(
    c("-c", "--cores"),
    type = "integer",
    default = 1,
    action = 'store',
    help = "Number of cores to use. [default= %default]",
    metavar = "integer"
  ),
  optparse::make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = 'output/',
    action = 'store',
    help = "Output folder. [default= %default]",
    metavar = "character"
  ),
  optparse::make_option(
    c("--adapters_list"),
    type = "character",
    default = NULL,
    action = 'store',
    help = "list of adapters to trim, if not provided general illumina adaters will be used",
    metavar = "character"
  ),
  optparse::make_option(
    c("--dont_trim"),
    type = "logical",
    default = FALSE,
    action = 'store_true',
    help = "if selected reads won't be trimmed",
    metavar = "logical"
  ),
  optparse::make_option(
    c("--phred"),
    type = 'numeric',
    default = 33,
    action = 'store',
    help = "either 33 or 64 [default: %default]",
    metavar = "numeric"
  ),
  optparse::make_option(
    c("--min"),
    type = "numeric",
    default = 25,
    action = 'store',
    help = "read minimum size to be kept [default: %default]",
    metavar = "numeric"
  )
)

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list),
                 convert_hyphens_to_underscores = T)

#separate
opt$one = stringr::str_split(string = opt$one, pattern = ',')[[1]]

if ('two' %in% names(opt)) {
  opt$two = stringr::str_split(string = opt$two, pattern = ',')[[1]]
}

#load needed packages
info=Kronos.scRT::FastqToBam(
  bowtie2_index = opt$index,
  File1 = opt$one,
  File2 = opt$two,
  adapters_list = opt$adapters_list,
  outputdir = opt$output_dir,
  cores = opt$cores ,
  trim = !opt$dont_trim,
  phred33 = opt$phred==33 ,
  Read_min_size_after_trimming = opt$min
)

plot= ggplot2::ggplot(data = info) +
  ggplot2::geom_boxplot(
    ggplot2::aes(Category, Values),
    outlier.size =  0,
    outlier.colour = NA,
    width = 0.8
  ) +
  ggplot2::geom_jitter(ggplot2::aes(Category, Values), alpha =
                        0.2,color='red') +
  ggplot2::theme_bw() +
  ggplot2::labs(x = '', y = '')  +
  ggplot2::theme(
    legend.position = 'none',
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x)
      sprintf(fmt = '%.1f', x)
  ) +
  ggforce::facet_row(
    ~ Parameter,
    scales = 'free',
    space = 'free',
    strip.position = 'right'
  )

#create metrics folder
Mfolder=file.path(opt$output_dir,'FastqToBamMetrics')
if(!dir.exists(Mfolder)){
  dir.create(Mfolder,recursive = T)
}

readr::write_tsv(x = info,file = file.path(
    Mfolder,
    'FastqToBamMetrics.tsv'
  ))

suppressMessages(ggplot2::ggsave(plot = plot,filename = file.path(
  Mfolder,
  'FastqToBamMetrics.pdf'
),device = grDevices::cairo_pdf))

print('done')
