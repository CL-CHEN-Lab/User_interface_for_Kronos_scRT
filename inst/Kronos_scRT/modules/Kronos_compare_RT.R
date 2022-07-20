#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

#parse
option_list = list(
  optparse::make_option(
    c("-R", "--RTs"),
    type = "character",
    default = NULL,
    help = "RT files with same binning",
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
    c('-D', "--deltaRT_threshold"),
    type = "double",
    action = "store",
    help = "DeltaRT threshold to define changes between 0 and 1",
    metavar = "double"
  ),
  optparse::make_option(
    c("-f", "--group_filter"),
    type = "character",
    default = NULL,
    help = "Filter out unwanted samples for RT files",
    metavar = "character"
  ),
  optparse::make_option(
    c("-c", "--colors"),
    type = "character",
    default = NULL,
    help = "min, midpoint and max colors separated by a comma",
    metavar = "character"
  )
)

opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#operators
`%>%`=tidyr::`%>%`
`%do%`=foreach::`%do%`
`%dopar%`=foreach::`%dopar%`

#set plotting theme
ggplot2::theme_set(ggplot2::theme_bw())

#create output directory
if(!dir.exists(opt$out)){
  dir.create(opt$out,recursive = T)
}


if (!'RTs' %in% names(opt)) {
    stop("No RT file has been provided. See script usage (--help)")
}
#check and load files
opt$RTs = stringr::str_split(opt$RTs, ',', simplify = F)[[1]]

if(!all(sapply(opt$RTs , function(x)
  Kronos.scRT::right_format(
    file_path = x,
    columns_to_check = c('chr', 'start', 'end', 'group', 'basename', 'RT'),
    logical = T
  ), simplify = T))) {
  stop(
    'One of the provided RTs it is not a tab separated file with the following columns: chr, start, end, group, basename and RT.'
  )
}

data <- Kronos.scRT::load_multiple_df(opt$RTs,delim = '\t',col_types = readr::cols(chr='c'))


if('group_filter' %in% names(opt)){
    opt$group_filter=stringr::str_split(opt$group_filter, ',', simplify = F)[[1]]
    data=data%>%
        dplyr::filter(!group %in% opt$group_filter)
}

if('colors' %in% names(opt)){

  opt$colors=stringr::str_split(opt$colors, ',', simplify = F)[[1]]

  if(length(opt$colors)<3){
    warning('Not enough colors were provided. Standard settings will be used.')
    opt$colors=NULL

  }else if(!all(sapply(opt$colors[1:3],function(x) tryCatch(is.matrix(col2rgb(x)),
           error=function(e) F )))){
    warning('Provided colors are not valid. Standard settings will be used.')
    opt$colors=NULL
  }

}

#run RT_changes_plot
plot=Kronos.scRT::RT_changes_plot(data,colors = opt$colors,deltaRT = opt$deltaRT_threshold)

if('deltaRT_threshold' %in% names(opt)){
  suppressMessages(
    ggplot2::ggsave(
      paste0(file.path(opt$out, 'ChangingRT_deltaRT_'),opt$deltaRT_threshold,
             '.pdf'),
      plot = plot,
      limitsize = FALSE,
      device = grDevices::cairo_pdf
    )
  )
}else{
  suppressMessages(
    ggplot2::ggsave(
      file.path(opt$out, 'ChangingRT.pdf'),
      plot = plot,
      limitsize = FALSE,
      device = grDevices::cairo_pdf
    )
  )
}

print('done')


