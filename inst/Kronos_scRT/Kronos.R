#!/usr/local/bin/Rscript
#recover arguments
args <- commandArgs(trailingOnly = FALSE)

#recover path to folder
Path_Kronos=args[stringr::str_detect(string = args,pattern = '--file=')]
Path_Kronos=dirname(stringr::str_remove_all(string = Path_Kronos,pattern = '--file='))

#recover arguments to pass to other scripts
args_position=which(stringr::str_detect(string = args,pattern = '--args'))


#dictionary check call and return file name
dictionary=function(call){

    result=dplyr::case_when(
                     call[1]=='binning' ~ list('Kronos_binning.R',2),
                     call[1]=='RT' ~ list('Kronos_RT.R',2),
                     call[1]=='WhoIsWho' ~ list('Kronos_whoswho.R',2),
                     call[1]=='CNV' ~ list('Kronos_cnv.R',2),
                     call[1]=='10XtoKronos' ~ list('10XtoKronos.R',2),
                     call[1]=='diagnostic' ~ list('Kronos_diagnostic.R',2),
                     call[1]=='fastqtoBAM' ~ list('Kronos_fastqtoBAM.R',2),
                     call[1]=='scPlots' ~ list('Kronos_scPlots.R',2),
                     call[1]=='Corr' ~ list('Kronos_Corr.R',2),
                     call[1]=='DRed' ~ list('Kronos_DRed.R',2),
                     call[1]=='BinRepProb' ~ list('Kronos_BinRepProb.R',2),
                     call[1]=='RefRT' ~ list('Kronos_RefRT.R',2),
                     call[1]=='cluster'& call[2] == 'RT' ~ list('Kronos_cluster_RT.R',3),
                     call[1]=='compare'& call[2] == 'RT' ~ list('Kronos_compare_RT.R',3),
                     call[1]=='compare'& call[2] == 'TW' ~ list('Kronos_compare_TW.R',3),
                     T ~ list(NA,NA))

    names(result)=c('module','arg_start')

    return(result)

}

if(length(length(args) >= (args_position + 1))) {
    args = args[(args_position + 1):length(args)]
    called_module = dictionary(args)
} else{
    called_module=list(module=NA)
}
#if module does not exist, send help message
if(is.na(called_module$module)){
message("\nProgram: Kronos (Tools for the single cell analysis of the replication timing program)\n
Usage: Kronos <command> [options]\n
commands:\n
\t 10XtoKronos \t\t Converts 10XGenomics output to Kronos scRT format
\t fastqtoBAM \t\t Trims and maps reads
\t binning\t\t Calculates mappability and gc content for bins to be used with Kronos CNV
\t CNV\t\t\t Calculates copy number variation
\t WhoIsWho \t\t Manually assign cell cycle stage
\t diagnostic \t\t Plotting tools to identify appropriate thresholds for Kronos RT
\t RT \t\t\t Calculates scReplication profiles and scRT
\t RefRT \t\t\t Rebins Reference RT
\t Corr \t\t\t Calculates pairwise spearman correlation between multipe pseudo-bulk replication timing/ rescaled bulk RT files
\t cluster RT \t\t Identify RT signatures between different experiments
\t compare RT \t\t Compares RT results from multiple experiments
\t compare TW \t\t Compares variability from multiple experiments and/or over multiple regions
\t scPlots \t\t scRT plots
\t DRed \t\t\t Performs Dimension Reduction using TSNE and UMAP
\t BinRepProb \t\t Calculates bin replication probablitity
")

}else{
    # if no arguments have been provided show help
    if (length(args) < called_module$arg_start) {
        out = system(intern = T,
                     command = paste(
                         'Rscript',
                         file.path(Path_Kronos, 'modules', called_module$module),
                         '-h'
                     ))
        out[1] = paste('Usage: Kronos',paste(args[1:(called_module$arg_start-1)],collapse=' '),'[options]')
        message(paste(out, collapse = '\n'))
    } else{
        args = paste(args[called_module$arg_start:length(args)], collapse = ' ')
        system(command = paste(
            'Rscript',
            file.path(Path_Kronos, 'modules', called_module$module),
            args
        ))
    }

}
