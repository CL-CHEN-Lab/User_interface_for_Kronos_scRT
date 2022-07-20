#options
options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

option_list = list(
  optparse::make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Dataset file name",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--out"),
        type = "character",
        default = "./output",
        help = "Output directory [default= %default]",
        metavar = "character"
    ),
    optparse::make_option(
        c("-C", "--Auto_correct"),
        type = "logical",
        action = "store_true",
        default = F,
        help = "If True diagnostic tries to automatically find the best parameters to reconstitute the S phase [default= %default]",
        metavar = "logical"
    ),
    optparse::make_option(
        c("-S", "--threshold_Sphase"),
        type = "double",
        help = "Threshold to identify S-phase cells",
        metavar = "double",
        default = NULL
    ),
    optparse::make_option(
        c("-G", "--threshold_G1G2phase"),
        type = "double",
        help = "Threshold to identify G1-phase cells. -S has to be selected and has to be bigger than -G",
        metavar = "double",
        default = NULL
    ),
    optparse::make_option(
        c("-f", "--Sphase_first_part"),
        type = "double",
        help = "Correction parameter for the first part of the S-phase [0.95,1]",
        metavar = "double",
        default = NULL
    ),
    optparse::make_option(
        c("-s", "--Sphase_second_part"),
        type = "double",
        help = "Correction parameter for the second part of the S-phase [0.5,0.55]",
        metavar = "double",
        default = NULL
    ),
    optparse::make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        help = "Number of parallel jobs to run [default= %default] ",
        metavar = "integer"
    ),
    optparse::make_option(
        c("-m", "--min_n_reads"),
        type = "double",
        default = 160,
        action = 'store',
        help = "Min n of reads per million per haploid genome to keep a cell in the analysis [default= %default]",
        metavar = "double"
    )
)

opt = optparse::parse_args( optparse::OptionParser(option_list=option_list))

#check inputs
if (!'file' %in% names(opt)) {
    stop("PerCell stat file must be provided. See script usage (--help)")
} else{
    # check file
    if (!file.exists(opt$file)) {
        stop("The provided PerCell stat file does not exist.")
        #check format
    } else if (!Kronos.scRT::right_format(
        file_path = opt$file,
        columns_to_check = c(
            "Cell",
            "normalized_dimapd",
            "mean_ploidy",
            "ploidy_confidence",
            "is_high_dimapd",
            "is_noisy",
            "coverage_per_1Mbp",
            "basename",
            "group"
        ),
        delim = ',',
        logical = T
    )) {
        stop("The wrong file has been provided as PerCell stat file.")
    }
}

#create output directory
if(!dir.exists(opt$out)){
    dir.create(opt$out,recursive = T)
}

#load data
Results=Kronos.scRT::diagnostic(PerCell = readr::read_csv(opt$file,col_types = readr::cols()),
                        interactive = F,
                        min_RPMPH = opt$min_n_reads,
                        G1G2_th = opt$threshold_G1G2phase,
                        S_th = opt$threshold_Sphase,
                        First_half_factor = ifelse(!is.null(opt$Sphase_first_part),opt$Sphase_first_part,1),
                        Second_half_factor = ifelse(!is.null(opt$Sphase_second_part),opt$Sphase_second_part,1),
                        Automatic_correction =opt$Auto_correct,
                        cores = opt$cores)

#save outputs
readr::write_tsv(x = Results$Settings, file = file.path(opt$out, paste0(
  Results$Settings$basename, '_settings.txt'
)))
suppressMessages(ggplot2::ggsave(
  plot = Results$all_cells_plot,
  filename = file.path(opt$out, paste0(
    Results$Settings$basename, '_all_cells.pdf'
  )),
  device = grDevices::cairo_pdf
))
suppressMessages(ggplot2::ggsave(
  plot = Results$first_filtering_plot,
  filename = file.path(
    opt$out,
    paste0(Results$Settings$basename, '_filtered_cells.pdf')
  ),
  device = grDevices::cairo_pdf
))
suppressMessages(ggplot2::ggsave(
  plot = Results$first_filtering_plot,
  filename = file.path(
    opt$out,
    paste0(Results$Settings$basename, '_filtered_cells.pdf')
  ),
  device = grDevices::cairo_pdf
))
if(!is.na(Results$Settings$threshold_Sphase)){
  suppressMessages(ggplot2::ggsave(
  plot = Results$selected_G1_S_cells_plot,
  filename = file.path(
    opt$out,
    paste0(Results$Settings$basename, '_Final_G1G2_S_phase_staging.pdf')
  ),
  device = grDevices::cairo_pdf
))
}
suppressMessages(ggplot2::ggsave(
  plot = Results$fixed_S_phase_plot,
  filename = file.path(
    opt$out,
    paste0(Results$Settings$basename, '_final_staging',
          ifelse(is.na(Results$Settings$threshold_G1G2phase),'',paste0('_G1G2_phase_th_',Results$Settings$threshold_G1G2phase)),
          ifelse(is.na(Results$Settings$threshold_Sphase),'',paste0('_S_phase_th_',Results$Settings$threshold_Sphase)),
          '_Sphase_first_part_parameter_',Results$Settings$Sphase_first_part,
          '_Sphase_second_part_parameter_',Results$Settings$Sphase_second_part,
          '.pdf')
  ),
  device = grDevices::cairo_pdf
))

suppressMessages(ggplot2::ggsave(
  plot = Results$S_phase_cell_distribution_plot,
  filename = file.path(
    opt$out,
    paste0(Results$Settings$basename, '_final_Sphase_distribution',
           ifelse(is.na(Results$Settings$threshold_G1G2phase),'',paste0('_G1G2_phase_th_',Results$Settings$threshold_G1G2phase)),
           ifelse(is.na(Results$Settings$threshold_Sphase),'',paste0('_S_phase_th_',Results$Settings$threshold_Sphase)),
           '_Sphase_first_part_param_',Results$Settings$Sphase_first_part,
           '_Sphase_second_part_par_',Results$Settings$Sphase_second_part,
           '.pdf')
  ),
  device = grDevices::cairo_pdf
))

print('done')
