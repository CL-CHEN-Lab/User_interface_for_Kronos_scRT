#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)
#parse input
option_list = list(
    optparse::make_option(
        c("-D", "--directory"),
        type = "character",
        default = NULL,
        help = "Single cell Bamfiles directory",
        metavar = "character"
    ),
    optparse::make_option(
        c("-B", "--bins"),
        type = "character",
        default = NULL,
        help = "File with bins produced by Kronos binning",
        metavar = "character"
    ),
    optparse::make_option(
      c("-C", "--chrSizes"),
      type = "character",
      default = NULL,
      help = "Chromosome size file",
      metavar = "character"
    ),
    optparse::make_option(
        c("-n", "--min_n_reads"),
        type = "double",
        default = 200000,
        action = 'store',
        help = "Min n of reads to keep a cell in the analysis [default= %default]",
        metavar = "double"
    ),
    optparse::make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 3,
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
        c("-e", "--ExpName"),
        type = "character",
        default = 'Exp',
        action = 'store',
        help = "Experiment name. [default= %default]",
        metavar = "character"
    ),
    optparse::make_option(
      c("-g", "--ExpGroup"),
      type = "character",
      default = 'Exp',
      action = 'store',
      help = "Group to which this experimetn belongs to. [default= %default]",
      metavar = "character"
    ),
    optparse::make_option(
        c("-p", "--ploidy"),
        type = "numeric",
        help = "user extimated ploidy",
        metavar = "numeric"
    ),
    optparse::make_option(
        c("-m", "--mim_mean_CN_accepted"),
        type = "numeric",
        help = "Min mean CN accepted as result. [default= %default]",
        default = 2,
        metavar = "numeric"
    ),
    optparse::make_option(
        c("-M", "--max_mean_CN_accepted"),
        type = "numeric",
        help = "Max mean CN accepted as result. [default= %default]",
        default = 8,
        metavar = "numeric"
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

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list),
                 convert_hyphens_to_underscores = T)

#declare operator
`%>%`=tidyr::`%>%`

if(!dir.exists(opt$output_dir)){
  dir.create(opt$output_dir,recursive = T)
  }

#check inputs
if (!'bins' %in% names(opt)) {
  stop("Bins with gc percentage not provided. See script usage (--help)")

} else{
  if (file.exists(opt$bins)) {
    if (!Kronos.scRT::right_format(
      file_path = opt$bins,
      columns_to_check = c(
        'chr',
        'start',
        'end',
        'mappability',
        'mappability_th',
        'gc_frequency',
        'type'
      ),
      logical = T
    )) {
      stop(paste(opt$bins,
                 ', does not have the right format'))
    }
  } else{
    stop(paste(opt$bins, 'does not exist'))
  }

}

if (!"directory" %in% names(opt)) {
    stop("Directory to bam files not provided. See script usage (--help)")
} else{
    if (dir.exists(opt$directory)) {
        if (length(list.files(opt$directory, pattern = '.bam$')) == 0) {
            stop(paste0(opt$directory, " does not contain bam files."))
        }
    } else{
        stop(paste0(opt$directory, " does not exist."))
    }
}


if (!'chrSizes' %in% names(opt)) {
  warning("Chromosome sizes file was not provide, sizes will be calculated from the bins file")
  #chr info
  genome.Chromsizes <-  bins %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(size = max(end))
}else{
  if(Kronos.scRT::right_format(file_path = opt$chrSizes,delim = '\t',columns_to_check=2,logical = T)){
    genome.Chromsizes <- readr::read_tsv(opt$chrSizes,col_names = c('chr','size'), col_types = readr::cols(chr = 'c'))
  }else{
    warning("Provided chromosome sizes file does not reflect the required format, sizes will be calculated from the bins file")
    #chr info
    genome.Chromsizes <-  bins %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(size = max(end))
  }
}

# load bins and gc percentage
bins = readr::read_tsv(opt$bins, col_types = readr::cols(chr = 'c')) %>%
    dplyr::arrange(chr, start)

### run CNV
results=Kronos.scRT::CallCNV(directory = opt$directory,
                     bins = bins,
                     chrom_size = genome.Chromsizes,
                     basename =  opt$ExpName,
                     group = opt$ExpGroup,
                     tmp_dir = file.path(opt$output_dir,'tmp'),
                     min_n_reads = opt$min_n_reads,
                     mim_mean_CN_accepted =opt$mim_mean_CN_accepted,
                     max_mean_CN_accepted = opt$max_mean_CN_accepted,
                     ploidy = opt$ploidy,
                     chr_prefix =ifelse(opt$chr_prefix == 'none', '', opt$chr_prefix),
                     chr_range = opt$chr_range,
                     cores = opt$cores
                     )

results$PerCell %>%
  readr::write_csv(file = file.path(
    opt$output_dir,
    paste0(opt$ExpName, '_PerCell.csv')
  ))
results$CNV %>%
  readr::write_tsv(file = file.path(
    opt$output_dir,
    paste0(opt$ExpName, '_scCNV.tsv')
  ))

print('done')
