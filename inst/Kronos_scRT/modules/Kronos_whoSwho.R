#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

#parse input
option_list = list(
  optparse::make_option(
    c("-F", "--file"),
    type = "character",
    default = NULL,
    help = "PerCell stat file path",
    metavar = "character"
  ),
  optparse::make_option(
    c("-W", "--whoSwho"),
    type = "character",
    default = NULL,
    help = "Who's who file path ( tsv file with header: Cell \t S_Phase)",
    metavar = "character"
  ),
  optparse::make_option(
    c("-o", "--out"),
    type = "character",
    default = "./output",
    help = "Output directory [default= %default]",
    metavar = "character"
  )
)

opt = optparse::parse_args(optparse::OptionParser(option_list = option_list))

#check inputs
if (!'file' %in% names(opt)) {
  stop("PerCell stat file must be provided. See script usage (--help)")
  if (!file.exists(opt$file)) {
    stop('Provided PerCell stat file does not exist')
  } else if (!Kronos.scRT::right_format(
    file_path = opt$file,
    columns_to_check = c(
      'Cell',
      'normalized_dimapd',
      'mean_ploidy',
      'ploidy_confidence',
      'is_high_dimapd',
      'is_noisy',
      'coverage_per_1Mbp'
    ),
    delim = ',',
    logical = T
  )) {
    stop('Provided PerCell stat file does not have the right format')

  }
}

if (!'whoSwho' %in% names(opt)) {
  stop("who's who file must be provided. See script usage (--help)")
  if (!file.exists(opt$file)) {
    stop("Provided who's who file does not exist")
  } else if (!Kronos.scRT::right_format(
    file_path = opt$whoSwho,
    columns_to_check = c('Cell',
                         'S_Phase'),
    logical = T
  )) {
    stop("Provided who's who file does not have the right format")
  }
}


#create output directory
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = T)
}

#run whoiswho
data = Kronos.scRT::WhoIsWho(
  PerCell = readr::read_csv(opt$file,
                            col_types = cols()),
  WhoIsWho = readr::read_tsv(opt$whoSwho,
                             col_types = cols())
)
# write data
readr::write_csv(x = data, file = file.path(opt$out, paste0('phased_', basename(opt$file))))

print('done')
