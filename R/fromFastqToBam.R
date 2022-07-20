#' Identifies adapters present in a fastq
#'
#' @importFrom Biostrings vcountPattern DNAStringSet
#' @importFrom stringr str_sub
#' @importFrom readr read_lines
#' @param File, path to a fastq file
#' @param Adapter_sequences, a vector containing the adapters to test
#'
#' @export
#'

IdAdapters = function(File, Adapter_sequences) {
  #read first 10^6 reads
  reads = readr::read_lines(File, n_max = 4 * 10 ^ 6)
  #select sequences and convert into DNAstringset
  reads = Biostrings::DNAStringSet(reads[base::seq(2, length(reads), 4)])

  #recover short sequence
  short = sapply(Adapter_sequences, function(x)
    stringr::str_sub(x, start = 1, end = 13), simplify = T)

  #count occurrences
  counts = sapply(short, function(x)
    mean(Biostrings::vcountPattern(
      subject = reads, pattern = x
    )), simplify = T)

  #return adapter sequence if at least 10% of reads presents such adapter
  counts = ceiling(round(counts, 2)) == 1
  return(Adapter_sequences[counts])
}

#' Sequentially trim reads
#'
#' @importFrom tidyr %>% drop_na
#' @importFrom dplyr tibble lead mutate
#' @importFrom stringr str_split
#' @importFrom Rbowtie2 remove_adapters
#'
#' @param File, path to a fastq file
#' @param Adapter_sequences, a vector containing the adapters to test
#' @param outputdir, output directory
#' @export
#'
sequential_trimming = function(File, Read_adaters, outputdir,min_size=25,cores=1) {
  #load operator
  `%>%` = tidyr::`%>%`

  #create folder if it does not exits
  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = T)
  }

  #count adapters
  n_adapters = length(Read_adaters)

  #create df with directory and final file output
  temporary_files_definition = dplyr::tibble(IN = File,
                                             OUT = file.path(outputdir, paste0(
                                               stringr::str_split(
                                                 string = basename(File),
                                                 pattern = '.fastq|.fq',
                                                 simplify = T
                                               )[1, 1],
                                               '_trimmed.fq'
                                             )),
                                             Delete = NA)
  if (n_adapters == 0) {
    #if no adapters is found copy file into exit folder
    file.copy(temporary_files_definition$IN,
              temporary_files_definition$OUT)

  } else{
    #if you have multiple adapters add tmp files
    if (n_adapters != 1) {
      temporary_files_definition = dplyr::tibble(
        IN = c(
          temporary_files_definition$IN,
          paste0(
            temporary_files_definition$OUT,
            1:(n_adapters - 1),
            'tmp.fq'
          ),
          temporary_files_definition$OUT
        )
      )

      temporary_files_definition = temporary_files_definition %>%
        dplyr::mutate(OUT = dplyr::lead(IN, 1)) %>%
        tidyr::drop_na() %>%
        dplyr::mutate(Delete = IN)

      temporary_files_definition$Delete[1] = NA

    }

    #run trimming
    for (i in 1:nrow(temporary_files_definition)) {
      Rbowtie2::remove_adapters(
        file1 = temporary_files_definition$IN[i],
        adapter1 = Read_adaters[i],
        output1 = temporary_files_definition$OUT[i],... = paste('--minlength',min_size,'--threads',cores)
      )

      #delete tmp files
      if (!is.na(temporary_files_definition$Delete[i])) {
        file.remove(temporary_files_definition$Delete[i])
      }
    }
  }
}


#' Trims and maps reads based on Rbowtie2
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom stringr str_split str_replace
#' @importFrom doSNOW registerDoSNOW
#' @importFrom snow makeCluster stopCluster
#' @importFrom Rsamtools ScanBamParam scanBam asBam
#' @importFrom Rbowtie2 bowtie2
#' @importFrom dplyr tibble
#'
#' @param bowtie2_index, path to bowtie2 index
#' @param File1, list of paths to fastq files (R1 files for PE)
#' @param File2, list of paths to R2 fastq files for PE
#' @param adapters_list, list of adapters to trim, if not provided general illumina adaters will be used
#' @param outputdir, output path
#' @param cores, number of parallel jobs
#' @param trim, logical value, if False reads won't be trimmed
#' @param phred33, logical value, if False phred64 will be used
#' @param Read_min_size_after_trimming, read minimum size to be kept
#'
#' @export
#'
FastqToBam = function(bowtie2_index,
                      File1,
                      File2 = NULL,
                      adapters_list = NULL,
                      outputdir = file.path(getwd(), 'FastqToBam'),
                      cores = 1,
                      trim = T,
                      phred33 = T,
                      Read_min_size_after_trimming=25) {
  #load operator
  `%dopar%` = foreach::`%dopar%`

  #set cluster and its stop at the exit
  cl = snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))

  info = foreach::foreach(Cell = 1:length(File1),.combine = 'rbind') %dopar% {
    #if the user wants to trim
    if (trim) {
      #define the list of adapters to look for
      if (is.null(adapters_list)) {
        adapters = c('CTGTCTCTTATACACATCT',
                     'ATGTGTATAAGAGACA',
                     'AGATGTGTATAAGAGACAG')
      } else{
        adapters = adapters_list
      }

      #count occurrences
      Read1_adaters = Kronos.scRT::IdAdapters(File = File1[Cell], Adapter_sequences = adapters)

      #check that if both 'ATGTGTATAAGAGACA','AGATGTGTATAAGAGACAG' are present, only 'AGATGTGTATAAGAGACAG' will go on
      if (is.null(adapters_list) &
          all(c('ATGTGTATAAGAGACA', 'AGATGTGTATAAGAGACAG') %in% Read1_adaters)) {
        Read1_adaters = Read1_adaters[Read1_adaters != 'ATGTGTATAAGAGACA']
      }

      #trim
      Kronos.scRT::sequential_trimming(
        File = File1[Cell],
        Read_adaters = Read1_adaters,
        outputdir = file.path(outputdir, 'trimmed'),
        min_size = Read_min_size_after_trimming
      )
      #reconstruct same name given to the output of sequential_trimming
      timmedF1 = paste0(
        stringr::str_split(
          string = basename(File1[Cell]),
          pattern = '.fastq|.fq',
          simplify = T
        )[1, 1],
        '_trimmed.fq'
      )

      #do the same if file2 exits
      if (!is.null(File2[Cell])) {
        #count occurrences
        Read2_adaters = Kronos.scRT::IdAdapters(File = File2[Cell], Adapter_sequences = adapters)

        #check that if both 'ATGTGTATAAGAGACA','AGATGTGTATAAGAGACAG' are present, only 'AGATGTGTATAAGAGACAG' will go on
        if (is.null(adapters_list) &
            all(c('ATGTGTATAAGAGACA', 'AGATGTGTATAAGAGACAG') %in% Read2_adaters)) {
          Read2_adaters = Read2_adaters[Read2_adaters != 'ATGTGTATAAGAGACA']
        }

        #trim
        Kronos.scRT::sequential_trimming(
          File = File2[Cell] ,
          Read_adaters = Read2_adaters,
          outputdir = file.path(outputdir, 'trimmed'),
          min_size = Read_min_size_after_trimming
        )
        #reconstruct same name given to the output of sequential_trimming
        timmedF2 = paste0(
          stringr::str_split(
            string = basename(File2[Cell]),
            pattern = '.fastq|.fq',
            simplify = T
          )[1, 1],
          '_trimmed.fq'
        )
      }
    } else{
      #if trim is not required
      timmedF1 = File1[Cell]
      timmedF2 = File2[Cell]
    }
    #create BAM folder if it doesn't exist
    if (!dir.exists(file.path(outputdir, 'BAM'))) {
      dir.create(file.path(outputdir, 'BAM'), recursive = T)
    }

    if (is.null(File2[Cell])) {
      #defene output name
      SamFile = file.path(
        outputdir,
        'BAM',
        stringr::str_replace(
          string = timmedF1,
          pattern = ".fq|.fastq|.fastq.gz|fastq.zip",
          replacement = '.sam'
        )
      )
      #align with bowtie2
      Rbowtie2::bowtie2(
        bt2Index = bowtie2_index,
        samOutput = SamFile,
        seq1 = file.path(outputdir, 'trimmed', timmedF1),
        ... = ifelse(phred33, '--phred33', '--phred64'),
        overwrite = TRUE
      )
    } else{
      #define output name
      SamFile = file.path(
        outputdir,
        'BAM',
        stringr::str_replace(
          string = stringr::str_replace(
            string = timmedF1,
            pattern = ".fq|.fastq|.fastq.gz|fastq.zip",
            replacement = '.sam'
          ),
          pattern = '_R1_',
          replacement = '_'
        )
      )
      #align with bowtie2
      Rbowtie2::bowtie2(
        bt2Index = bowtie2_index,
        samOutput = SamFile,
        seq1 = file.path(outputdir, 'trimmed', timmedF1),
        seq2 = file.path(outputdir, 'trimmed', timmedF2),
        ... = ifelse(phred33, '--phred33', '--phred64'),
        overwrite = TRUE
      )
    }

    #convert to BAM, sort and index
    Bamfile<-Rsamtools::asBam(SamFile)
    #remove sam
    file.remove(SamFile)

    Kronos.scRT::BamMetrics(Bamfile = Bamfile ,isPE = !is.null(File2[Cell]))

  }

  return(info)
}

#' merges a list of fastq files into one
#'
#' @importFrom stringr str_sub
#'
#' @param FileList, a list of files to merge
#' @param NewName, output file name if not provided this function looks for a common base between files
#' @param output, output directory
#'
#' @export
#'

mergefastq = function(FileList,
                      NewName = NULL,
                      output = file.path(getwd(), 'mergedFastq')) {
  #find common base name
  if (is.null(NewName)) {
    f = sapply(basename(FileList), utf8ToInt, simplify = T)
    f = which.min(apply(f, 1, function(x)
      length(unique(x)) == 1)) - 1
    NewName = paste0(stringr::str_sub(
      string = basename(FileList)[1],
      start = 1,
      end = f
    ),
    '_combined.fastq')
  }
  #if directory does not exists create id
  if (!dir.exists(output)) {
    dir.create(output, recursive = T)
  }
  #create destination file
  NewName = file.path(output, NewName)
  file.create(NewName)

  #merge
  for (File in FileList) {
    dFile = Kronos.scRT::Decompress_if_needed(File = File, outdir = file.path(output, paste0('tmp_', basename(NewName))))

    file.append(file1 = NewName, file2 = dFile)

  }

  unlink(x = file.path(output, paste0('tmp_', basename(NewName))), recursive = T)

  return(NewName)
}


#' Checks if files are compressed and decompressed them in a destination folder returning the full path of the decompressed file or original file if no decompression was needed
#'
#' @importFrom R.utils isBzipped bunzip2 isGzipped gunzip
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom utils unzip
#'
#' @param File, file to test
#' @param outdir, output directory
#'
#' @export
#'

Decompress_if_needed = function(File, outdir = NULL) {
  #if no output dir save in the same folder
  if (is.null(outdir)) {
    outdir = dirname(File)
  }

  #this should be the final decompressed name
  exit = File
  destination = file.path(outdir, tools::file_path_sans_ext(basename(File)))

  #check for different types of compression
  if (R.utils::isBzipped(File)) {
    R.utils::bunzip2(File, destname = destination, remove = F)
    exit = destination

  } else if (R.utils::isGzipped(File, destname = destination)) {
    R.utils::gunzip(File, destname = destination, remove = F)
    exit = destination

  } else if (tools::file_ext(File) == 'zip') {
    exit = destination
    exit = tryCatch(
      expr = unzip(zipfile = File, exdir = outdir),
      warning = function(x) {
        return(File)
      }
    )
  }

  return(exit)
}

#' recover metrics from bam files
#'
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom dplyr tibble
#' @importFrom tidyr %>%
#'
#' @param Bamfile, a bam file directory
#' @param isPE, logical value indicating whether reads are PE (TRUE) or SE (FALSE)
#'
#' @export
#'
BamMetrics=function(Bamfile,isPE=T){
  #return info
  if(isPE){
    #First in a pair proper paired reads
    param1 <-
      Rsamtools::ScanBamParam(
        what = c('qname','mapq', 'qwidth', 'isize','rname'),
        flag = Rsamtools::scanBamFlag(
          isPaired = T,
          isUnmappedQuery = F,
          isFirstMateRead = T,
          isDuplicate = F
        )
      )

    #Second in a pair proper paired reads
    param2 <-
      Rsamtools::ScanBamParam(
        what = c('qname','mapq', 'qwidth', 'isize','rname'),
        flag = Rsamtools::scanBamFlag(
          isPaired = T,
          isUnmappedQuery = F,
          isSecondMateRead = T,
          isDuplicate = F
        )
      )

    info=full_join(dplyr::as_tibble(Rsamtools::scanBam(file = Bamfile,param =param1 )[[1]]),
                   dplyr::as_tibble(Rsamtools::scanBam(file = Bamfile,param =param2 )[[1]]),by='qname')


    info=dplyr::tibble(
      File = basename(Bamfile),
      Parameter =c(rep('MappingQuality', 11),
                   rep('ReadLength',11),
                   rep('Fragment size',11),
                   'TotalReads (log10)',
                   'ProperPair (%)',
                   'Unmapped Fraction (%)'),
      Values = c(
        quantile(
          x = c(info$mapq.x,info$mapq.y),
          probs = c(0.05,seq(0.1, 0.9, 0.1),0.95),
          na.rm = T
        ),
        quantile(
          x = c(info$qwidth.x,info$qwidth.y),
          probs = c(0.05,seq(0.1, 0.9, 0.1),0.95),
          na.rm = T
        ),
        quantile(
          x = abs(info$isize.x),
          probs = c(0.05,seq(0.1, 0.9, 0.1),0.95),
          na.rm = T
        ),
        log10(length(info$mapq.x)+length(info$mapq.y)),
        round(100*sum(!is.na(info$rname.x) & !is.na(info$rname.y))/length(info$mapq.x),2),
        round(100*(sum(is.na(info$rname.x))+sum(is.na(info$rname.y)))/ (length(info$mapq.x)+length(info$mapq.y)),2)),
      Category = c(rep(paste0(c(0.05,seq(0.1, 0.9, 0.1),0.95)*100,'th percentile'),3), 'Total Reads','Proper Pair','Unmapped Fraction')
    )
  }else{
  param <- Rsamtools::ScanBamParam(what=c('qwidth','mapq','rname'))
  info=Rsamtools::scanBam(file = Bamfile,param =param )[[1]]
  info=dplyr::tibble(
    File = basename(Bamfile),
    Parameter =c(rep('MappingQuality', 11),
                 rep('ReadLength',11),
                 'TotalReads (log10)',
                 'Unmapped Fraction (%)'),
    Values = c(
      quantile(
        x = info$mapq,
        probs = c(0.05,seq(0.1, 0.9, 0.1),0.95),
        na.rm = T
      ),
      quantile(
        x = info$qwidth,
        probs = c(0.05,seq(0.1, 0.9, 0.1),0.95),
        na.rm = T
      ),
      log10(length(info$mapq)),
      round(100*sum(is.na(info$rname)) / length(info$mapq),1)),
    Category = c(rep(paste0(c(0.05,seq(0.1, 0.9, 0.1),0.95)*100,'th percentile'),2), 'Total Reads', 'Unmapped Fraction')
  )
  }

  return(info)
}
