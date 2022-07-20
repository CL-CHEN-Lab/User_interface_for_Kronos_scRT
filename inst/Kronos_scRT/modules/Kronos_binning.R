#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)
#parse input
option_list = list(
  optparse::make_option(
    c("-R", "--RefGenome"),
    type = "character",
    default = NULL,
    help = "Fasta file of genome of interst",
    metavar = "character"
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
    c("-s", "--reads_size"),
    type = "integer",
    default = 40,
    action = 'store',
    help = "Length of the simulated reads. [default= %default bp]",
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
    c("-i", "--index"),
    type = "character",
    action = 'store',
    help = "Bowtie 2 index",
    metavar = "character"
  ),
  optparse::make_option(
    c("--paired_ends"),
    type = "logical",
    action = 'store_true',
    help = "Generates paired ends reads [default: %default]",
    metavar = "logical",
    default = F
  ),
  optparse::make_option(
    c("--fragment_size"),
    type = "integer",
    action = 'store',
    help = "Fragment size if paired end option is used. [default: %default]",
    metavar = "integer",
    default = '200'
  ),
  optparse::make_option(
    c("--bin_size"),
    type = "character",
    default = '20Kb',
    action = 'store',
    help = "Bins size. [default= %default ]",
    metavar = "character"
  ),
  optparse::make_option(
    c("-d", "--dir_indexed_bam"),
    type = "character",
    action = 'store',
    help = "If provided parameters will be automatically estimated from the data.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-u", "--upper_mappability_th"),
    type = "double",
    default = 1.5,
    action = 'store',
    help = "Maximum mappability for a bin to be considered in the analisys  [default= %default]",
    metavar = "double"
  ),
  optparse::make_option(
    c("-l", "--lower_mappability_th"),
    type = "double",
    action = 'store',
    default = 0.8,
    help = "Minimum mappability for a bin to be considered in the analisys  [default= %default]",
    metavar = "double"
  ),
  optparse::make_option(
    c("-B", "--black_list"),
    type = "character",
    action = 'store',
    help = "Regions to ignore",
    metavar = "character"
  ),
  optparse::make_option(
    c("-x", "--coverage"),
    type = "character",
    action = 'store',
    help = "Coverage for simulated genome. [default= %default]",
    default = '1x',
    metavar = "character"
  ),
  optparse::make_option(
    c("-e", "--errorRate"),
    type = "character",
    action = 'store',
    help = "Simulated sequencing error rate (%) [default= %default]",
    default = "0.1%",
    metavar = "character"
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

opt = optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  convert_hyphens_to_underscores = T
)

# create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir,recursive = T)
}

# check inputs
if (!"RefGenome" %in% names(opt)) {
  stop("Fasta file not provided. See script usage (--help)")
} else{
  if (!file.exists(opt$RefGenome)) {
    stop("Provided Fasta file does not exist")
  }
}

if (!"index" %in% names(opt)) {
  stop("Bowtie2 indexed genome not provided. See script usage (--help)")
}

opt$coverage = tryCatch(
  expr = as.numeric(stringr::str_remove(opt$coverage, '[Xx]')),
  error = function(x)
    tryCatch(
      expr = as.numeric(opt$coverage),
      error = stop('Wrong coverage format')
    )
)

converted_error_rate = tryCatch(
  expr = as.numeric(stringr::str_remove(opt$errorRate, '[%]'))/100,
  error = function(x)
    tryCatch(
      expr = as.numeric(opt$errorRate)/100,
      error = stop('Wrong errorRate format')
    )
)

# convert binsize to numeric
extract_unit = stringr::str_extract(opt$bin_size, pattern = '.{2}$')

if (grepl(x = extract_unit, pattern =  '[0-9][0-9]')) {
  n_of_0 = stringr::str_length(str_extract(opt$bin_size, '0{1,10}$'))
  BS = dplyr::case_when(
    is.na(n_of_0) ~ paste0(opt$bin_size, 'bp'),
    n_of_0 < 3 ~ paste0(opt$bin_size, 'bp'),
    n_of_0 < 6 ~ paste0(str_remove(opt$bin_size, '0{3}$'), 'Kb'),
    n_of_0 >= 6 ~ paste0(str_remove(opt$bin_size, '0{6}$'), 'Mp')
  )
} else{
  BS = opt$bin_size
}

opt$bin_size = as.numeric(stringr::str_remove(opt$bin_size, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]")) * dplyr::case_when(
  grepl(x = extract_unit, pattern =  '[Kk][Bb]') ~ 1000,
  grepl(x = extract_unit, pattern = '[Mm][Bb]') ~ 1000000,
  grepl(x = extract_unit, pattern = '[Bp][Pp]') ~ 1,
  grepl(x = extract_unit, pattern =  '[0-9][0-9]') ~ 1
)

if (is.na(opt$bin_size)) {
  stop('binsize have an incorrect format')
}

if ('black_list' %in% names(opt)) {
  #check if the file exists
  if (!file.exists(opt$black_list)) {
    stop('Blacklist file does not exit')
  } else if (Kronos.scRT::right_format(
    file_path = opt$black_list,
    columns_to_check = 3,
    delim = '\t',
    logical = T
  )) {
    stop('Blacklist file does not have the right format')
  }
}

bins = Kronos.scRT::binning(
  RefGenome = opt$RefGenome,
  cores = opt$cores,
  #if not provided return NULL
  directory_to_bamfiles = opt$dir_indexed_bam,
  bowtie2_index = opt$index,
  bin_size = opt$bin_size,
  read_size = opt$reads_size,
  fragment_size = opt$fragment_size,
  paired_ends = opt$paired_ends,
  tmp_dir = file.path(opt$output_dir, 'tmp'),
  upper_mappability_th = opt$upper_mappability_th,
  lower_mappability_th = opt$lower_mappability_th,
  black_list = opt$black_list,
  coverage = opt$coverage,
  errorRate = converted_error_rate,
  chr_prefix = ifelse(opt$chr_prefix == 'none', '', opt$chr_prefix),
  chr_range = opt$chr_range,
  return_estimated_param = T
)


#separate bins from parameters
param = bins$param
bins = bins$bins

#write bins with info
readr::write_tsv(x = bins, file =  file.path(
  opt$output_dir,
  paste0(
    basename(opt$index),
    '_',
    BS,
    '_bins_coverage_',
    opt$coverage,
    'X_',
    paste0(
      'reads_',
      param$read_size,
      'bp_',
      ifelse(
        param$paired_ends,
        paste0('PE_fragmentSize_',
               param$fragment_size, 'bp_'),
        'SE_'
      )
    ),
    ifelse('black_list' %in% names(opt), 'blacklisted_', ''),
    paste0(
      'error_rate_',
      opt$errorRate,
      '_min_mappability_',
      opt$lower_mappability_th,
      '_max_mappability_',
      opt$upper_mappability_th
    ),
    '.tsv'
  )
))

print('done')
