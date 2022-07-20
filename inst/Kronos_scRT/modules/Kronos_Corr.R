#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)
#parse input
option_list = list(
    optparse::make_option(
        c("-F", "--File"),
        type = "character",
        default = NULL,
        help = "Replication timing files separated by a comma. Format: chr <TAB> start <TAB> end <TAB> group",
        metavar = "character"
    ),
    optparse::make_option(
        c("-s", "--sort"),
        type = "character",
        default = NULL,
        help = "Group names orders",
        metavar = "character"
    ),
    optparse::make_option(
        c("-m", "--method"),
        type = "character",
        default = 'spearman',
        help = "Correlation method: spearman, pearson or kendall [default= %default]",
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
    )
)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#define operators
`%do%`=foreach::`%do%`
`%>%`= tidyr::`%>%`

if (!'File' %in% names(opt)) {
    stop("RT files were not provided. See script usage (--help)")

} else{
    opt$File = stringr::str_split(opt$File, ',')[[1]]
}
if ('sort' %in% names(opt)) {
    opt$sort = stringr::str_split(opt$sort, ',')[[1]]

}

#create directory
if(!dir.exists(opt$out)){
    dir.create(opt$out,recursive = T)
}


scRT = Kronos.scRT::load_multiple_df(dirs = opt$File, delim = '\t',col_types=readr::cols(chr='c'))

if ('sort' %in% names(opt)) {
    scRT=scRT %>%
        dplyr::mutate(group = factor(group, levels = opt$sort))
}

plot = Kronos.scRT::KCorr_plot(df = scRT,method = opt$method)

suppressMessages(ggplot2::ggsave(
    plot = plot,
    filename = file.path(
        opt$out,
        paste(
        opt$output_file_base_name,
        opt$method,'correlation.pdf',sep = '_'
    )
)))

print('done')

