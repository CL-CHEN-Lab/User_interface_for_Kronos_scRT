#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

#options
option_list = list(
    optparse::make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Per cell stat file , if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    optparse::make_option(
        c("-T", "--tracks"),
        type = "character",
        default = NULL,
        help = "Tracks file,  if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    optparse::make_option(
      c("-b", "--base_name"),
      type = "character",
      default = NULL,
      help = "Basename of an experiment, if multiple basenames are provided they have to be separated by a comma",
      metavar = "character"
    ),
    optparse::make_option(
      c("-g", "--groups"),
      type = "character",
      default = NULL,
      help = "Group name of an experiment, if multiple groups are provided they have to be separated by a comma",
      metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--out"),
        type = "character",
        default = "output",
        help = "Output directory [default= %default]",
        metavar = "character"
    )
)

#recover inputs
opt = optparse::parse_args(object = optparse::OptionParser(option_list = option_list))

#create directory
if(!dir.exists(opt$out)){
    dir.create(opt$out,recursive = T)
}

#recover multiple paths
opt$file = stringr::str_split(opt$file, ',')[[1]]

opt$tracks = stringr::str_split(opt$tracks, ',')[[1]]

opt$base_name = stringr::str_split(opt$base_name, ',')[[1]]

opt$groups = stringr::str_split(opt$groups, ',')[[1]]

#load operators
`%do%`=foreach::`%do%`
`%>%`=tidyr::`%>%`


#check length inputs and proceed with the conversion
if(length(opt$file)==length(opt$tracks) &
length(opt$file)==length(opt$base_name) &
length(opt$file)==length(opt$groups)){
.tmp=foreach::foreach(file=opt$file,tracks=opt$tracks,basename=opt$base_name,group=opt$groups)%do%{

  results=tryCatch(Kronos.scRT::TenXtoKronos(PerCell = file,CNV = tracks,basename = basename,group = group),
                   error=function(e) stop(paste('One of the following input lead to an error:',file,tracks,basename,group,sep = '\n')))

  results$PerCell%>%
    readr::write_csv(paste0(file.path(opt$out, 'Kronos_format_'), basename(file)))

  results$CNV%>%
    readr::write_tsv(paste0(file.path(opt$out, 'Kronos_format_'), stringr::str_replace(string = basename(tracks),pattern = '.bed',replacement = '.tsv')))

  basename
}
}else{
  stop('inputs have different length')
}

print('done')



