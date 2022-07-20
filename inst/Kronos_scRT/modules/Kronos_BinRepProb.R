#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

#parse input
option_list = list(
  optparse::make_option(
    c("-S", "--S_file"),
    type = "character",
    default = NULL,
    help = "S-phase scCNV file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
    metavar = "character"
  ),
  optparse::make_option(
    c("-G", "--G1G2_file"),
    type = "character",
    help = "Optional: G1/G2-phase scCNV file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
    metavar = "character"
  ),
  optparse::make_option(
    c("-R", "--RT_file"),
    type = "character",
    help = "BulkRT file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
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
    c("-e", "--early_range"),
    type = "character",
    default = "0,30",
    help = "percentage of replication range defying Early replicating cells (0-100) [default= %default]",
    metavar = "range"
  ),
  optparse::make_option(
    c("-m", "--mid_range"),
    type = "character",
    default = "40,60",
    help = "percentage of replication range defying Mid replicating cells (0-100) [default= %default]",
    metavar = "range"
  ),
  optparse::make_option(
    c("-l", "--late_range"),
    type = "character",
    default = "70,100",
    help = "percentage of replication range defying Late replicating cells (0-100) [default= %default]",
    metavar = "range"
  ),
  optparse::make_option(
    c("-q", "--G1G2_quantile_range"),
    type = "character",
    default = "0.25,0.75",
    help = "Ploidy quantile range to keep a G1/G2 cell [default= %default]",
    metavar = "range"
  ),
  optparse::make_option(
    c("--upper_range"),
    type = "character",
    default = "0.95,1.00",
    help = "limits of the upper panel [default= %default]",
    metavar = "range"
  ),
  optparse::make_option(
    c("--lower_range"),
    type = "character",
    default = "0.00,0.05",
    help = "limits of the lower panel [default= %default]",
    metavar = "range"
  ),
  optparse::make_option(
    c("-c", "--colors"),
    help = "Color Plot, multiple colors can be given usins a comma",
    metavar = "character",
    type = "character"
  )
)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))
#load operators
`%>%`=tidyr::`%>%`
#set plotting theme
ggplot2::theme_set(ggplot2::theme_bw())

#check inputs
if (!'S_file' %in% names(opt)) {
    stop("S-phase scCNV file must be provided. See script usage (--help)")
}

if (!'RT_file' %in% names(opt)) {
  stop("bulk/pseudobulk RT file must be provided. See script usage (--help)")
}

opt$S_file = stringr::str_split(opt$S_file, ',')[[1]]
opt$RT_file = stringr::str_split(opt$RT_file, ',')[[1]]

if(!all(sapply(opt$S_file, function (x)
  Kronos.scRT::right_format(
    file_path = x ,
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
    ),
    logical = T
  )))) {
  stop(
    'Provided S-phase scCNV file is not a tab spaced file with the following columns: chr, start, end, CN, background, CN_bg, th, Rep, PercentageReplication, Cell, basename, group and newIndex'
  )
}

if(!all(sapply(opt$RT_file, function (x)
  Kronos.scRT::right_format(
    file_path = x ,
    columns_to_check = c(
      'chr',
      'start',
      'end',
      'basename',
      'group',
      'RT'
    ),
    logical = T
  )))) {
  stop(
    'Provided bulk RT file is not a tab spaced file with the following columns: chr, start, end, basename, group and RT'
  )
}
#load df
Sphase=Kronos.scRT::load_multiple_df(dirs = opt$S_file,col_types=readr::cols(chr='c'))
RT=Kronos.scRT::load_multiple_df(dirs = opt$RT_file,col_types=readr::cols(chr='c'))

# if G1G2 are provided check files and load dfs
if ('G1G2_file' %in% names(opt)) {

  opt$G1G2_file = stringr::str_split(opt$G1G2_file, ',')[[1]]


  if (!all(sapply(opt$G1G2_file, function (x)
    Kronos.scRT::right_format(
      file_path = x ,
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
      ),
      logical = T
    )))) {
    stop(
      'Provided G1G2-phase scCNV file is not a tab spaced file with the following columns: chr, start, end, CN, background, CN_bg, th, Rep, PercentageReplication, Cell, basename, group and newIndex'
    )
  }

  G1G2_file=Kronos.scRT::load_multiple_df(dirs = opt$G1G2_file,col_types=readr::cols(chr='c'))
}

#create output
if(!dir.exists(opt$out)){
  dir.create(opt$out)
}

#set ranges
opt$early_range = as.numeric(stringr::str_split(
  opt$early_range,
  pattern = ',',
  n = 2,
  simplify = T
))
opt$mid_range = as.numeric(stringr::str_split(
  opt$mid_range,
  pattern = ',',
  n = 2,
  simplify = T
))
opt$late_range = as.numeric(stringr::str_split(
  opt$late_range,
  pattern = ',',
  n = 2,
  simplify = T
))
opt$G1G2_quantile_range = as.numeric(stringr::str_split(
  opt$G1G2_quantile_range,
  pattern = ',',
  n = 2,
  simplify = T
))
opt$lower_range = as.numeric(stringr::str_split(
  opt$lower_range,
  pattern = ',',
  n = 2,
  simplify = T
))
opt$upper_range = as.numeric(stringr::str_split(
  opt$upper_range,
  pattern = ',',
  n = 2,
  simplify = T
))

#set colors
if('colors' %in% names(opt)){
  opt$colors=stringr::str_split(
    opt$colors,
    pattern = ',',
    simplify = T
  )
}

#Prepare files
BinProb = Kronos.scRT::Prepare_S_phase_cells_forBinRepProb(
  S = Sphase,
  RT = RT,
  Early.cells = opt$early_range,
  Mid.cells = opt$mid_range ,
  Late.cells = opt$late_range
)
if ('G1G2_file' %in% names(opt)) {
  BinProb = rbind(
    BinProb,
    Kronos.scRT::Prepare_G1G2_phase_cells_forBinRepProb(
      G1.G2 = G1G2_file,
      RT = RT,
      quantile.range = opt$G1G2_quantile_range
    )
  )
}

#run bin prob-plot
BN=unique(BinProb$basename)
plots = lapply(BN, function (x)
  Kronos.scRT::BinRepProbPlot(
    Variability = BinProb  %>%
      dplyr::filter(basename  ==  x),
    colors = opt$colors,
    upper_range = opt$upper_range,
    lower_range = opt$lower_range
  ))

#save bin prob
BinProb%>%
  readr::write_tsv(paste0(file.path(
    opt$out,
    opt$output_file_base_name
  ),'_BinRepProb.tsv'))


#save plots
.tmp=lapply(1:length(BN), function(x)
suppressMessages(
  ggplot2::ggsave(
    plots[[x]],
  filename = paste0(file.path(opt$out,
                              BN[x]),
                    '_BinRepProb.pdf'),width = 28,height = 14,units = 'cm',device = grDevices::cairo_pdf
)))


print('done')
