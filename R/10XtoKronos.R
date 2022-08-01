#' It converts 10X Genomics data into Kronos format
#'
#' @return list
#'
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate select
#' @importFrom readr read_csv read_tsv cols
#'
#' @param PerCell, PerCell 10XGenomics file path
#' @param CNV, CNV 10XGenomics file path
#' @param basename, sample basename
#' @param group, experimental group
#'
#' @export
#'
TenXtoKronos = function(PerCell,CNV,basename,group) {

  #load operators
  `%>%`=tidyr::`%>%`

  #check percell file format,load and convert
  if (Kronos.scRT::right_format(
    file_path = PerCell,
    delim = ',',
    columns_to_check = c(
      "barcode",
      "cell_id",
      "total_num_reads",
      "num_unmapped_reads" ,
      "num_lowmapq_reads",
      "num_duplicate_reads" ,
      "num_mapped_dedup_reads",
      "frac_mapped_duplicates",
      "effective_depth_of_coverage",
      "effective_reads_per_1Mbp"  ,
      "raw_mapd",
      "normalized_mapd" ,
      "raw_dimapd",
      "normalized_dimapd" ,
      "mean_ploidy",
      "ploidy_confidence" ,
      "is_high_dimapd",
      "is_noisy"
    ),
    logical = T
  )){

    PerCell=readr::read_csv(PerCell,col_types = readr::cols())%>%
      dplyr::select(
        'Cell'=cell_id,
        normalized_dimapd,
        mean_ploidy,
        ploidy_confidence,
        is_high_dimapd,
        is_noisy,
        'coverage_per_1Mbp'=effective_reads_per_1Mbp
      )%>%
      dplyr::mutate(is_high_dimapd=as.logical(is_high_dimapd),
                    is_noisy=as.logical(is_noisy),
                    basename=basename,
                    group=group)

  }else{
    stop('Provided percell file does not have the expected format')
  }


  #check percell file format,load and convert
  if (Kronos.scRT::right_format(
    file_path = CNV,
    delim = '\t',skip = 2,
    columns_to_check = c(
      "#chrom",
      "start",
      "end",
      "id",
      "copy_number",
      "event_confidence"
    ),
    logical = T
  )){

    CNV=readr::read_tsv(CNV,col_types = readr::cols(),skip = 2)%>%
      dplyr::select(
        'Cell'=id,
        'chr'=`#chrom`,
        start,
        end,
        copy_number
      )%>%
      dplyr::mutate(reads = '10X',
                    basename=basename,
                    group=group)

  }else{
    stop('Provided CNV file does not have the expected format')
  }
  #return
  return(list(PerCell = PerCell, CNV = CNV))
  }

