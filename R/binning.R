#'  Given a genomic sequence it returns the relative coordinates where the sequence is known
#'
#' @return tibble
#'
#' @importFrom stringr str_locate_all
#' @importFrom dplyr tibble
#' @param x, a genomic sequence
#'
#' @export
#'

find_known_sequences = function(x) {
  y = stringr::str_locate_all(x, '.[TACG]+')
  return(dplyr::tibble(start = y[[1]][, 1], end = y[[1]][, 2]))
}

#'  Randomly mutates the 4 nucleotides
#'
#' @return vector
#'
#' @param Base, either a list of characters or integers corresponding to the 4 nucleotides (A,T,C,G)
#'
#' @export
#'

mutate_sequence = Vectorize(function(Base) {
  #if the function got a character as input convert it to int
  if (is.character(Base)) {
    B = utf8ToInt(Base)
  } else{
    B = Base
  }

  #apply mutation
  bases = c(65, 67, 84, 71)
  bases = bases[bases != B]

  mut = sample(bases, 1)

  #return the mutation with the same data format
  return(if (is.character(Base)) {
    intToUtf8(mut)
  } else{
    mut
  })

}, vectorize.args = 'Base')

#'  Given a genomic start and end recover a genomic sequence and simulate mutations due to sequencing
#'
#' @return tibble
#'
#' @importFrom stringr str_length str_remove str_sub
#' @importFrom dplyr tibble mutate
#' @param Reads, tibble with start and end coordinates
#' @param Sequence, genomic sequence from which recover reads
#' @param errorRate, error rate to simulate mutation
#'
#' @export
#'
#'
recover_and_mutate = function(Reads, Sequence, errorRate) {
  #load operator
  `%>%` = tidyr::`%>%`

  #if no order is provided assign the input order
  if (!'order' %in% names(Reads)) {
    Reads$order = 1:nrow(Reads)
  }

  len_Seq = stringr::str_length(Sequence)
  mutate = sample(1:len_Seq,
                  errorRate * len_Seq)

  #convert characters into int
  REF_mutated = utf8ToInt(as.character(Sequence[[1]]))

  #mutation
  REF_mutated[mutate] = Kronos.scRT::mutate_sequence(REF_mutated[mutate])

  #convert back to string
  REF_mutated = intToUtf8(REF_mutated)

  #recover reads sequences
  Reads = dplyr::tibble(
    reads = stringr::str_sub(
      string = REF_mutated,
      start = Reads$start,
      end = Reads$end - 1
    ),
    order = Reads$order
  ) %>%
    dplyr::mutate(
      reads = stringr::str_remove(reads, 'N{1,1000}$'),
      reads = stringr::str_remove(reads, '^N{1,1000}')
    )

  return(Reads)
}


#'  Reshapes data into fastq format and save them into a file
#'
#' @return NULL
#'
#' @importFrom readr write_delim
#' @importFrom dplyr select mutate arrange n
#' @importFrom tidyr gather
#' @importFrom stringr str_count
#'
#' @param Reads, output of recover_and_mutate
#' @param path, path where to save the file
#' @param Chr, Chromosome from which the reads are originated
#' @param ID, read identifier
#'
#' @export
#'
#'
reshape_and_save = function(Reads, Chrom, ID, path = './') {
  #load operator
  `%>%` = tidyr::`%>%`

  Reads = Reads %>%
    dplyr::arrange(order)

  colnames(Reads) = c('2', 'order')


  Reads %>%
    dplyr::mutate(
      n = stringr::str_count(`2`),
      `1` = paste0('@read', 1:dplyr::n(), Chrom, 'round', ID),
      `3` = '+',
      `4` = unlist(lapply(n, function(x)
        paste0(
          rep('D', x), collapse = ''
        )))
    ) %>%
    dplyr::select(-n) %>%
    tidyr::gather('pos', 'towrite', -order) %>%
    dplyr::arrange(order, pos) %>%
    dplyr::select(towrite) %>%
    readr::write_delim(file = path,
                       col_names = F,
                       append = T)
  return(NULL)
}

#'  Given a genomic sequence it calculates its reverse complementary
#'
#' @return string
#'
#' @importFrom stringr str_extract_all
#' @importFrom XVector rev
#'
#' @param sequence, a genomic sequence
#'
#' @export
#'

rev_com = function(sequence) {
  dict = list(
    A = 'T',
    `T` = 'A',
    C = 'G',
    G = 'C',
    N = 'N'
  )
  sequence = stringr::str_extract_all(sequence, '.')
  sequence = unlist(lapply(sequence , function(x)
    paste0(
      sapply(XVector::rev(x), function(y)
        dict[[y]], simplify = T),
      collapse = ''
    )))
  return(sequence)
}


#'  Calculates genomic bins with GC content and mappability
#'
#' @return dataframe
#'
#' @importFrom  dplyr filter as_tibble tibble arrange lead mutate right_join row_number select summarise ungroup
#' @importFrom  tidyr drop_na %>%
#' @importFrom  foreach %do% %dopar% foreach
#' @importFrom stringr str_count str_length str_split str_sub
#' @importFrom Biostrings readDNAStringSet width
#' @importFrom snow makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom Rsamtools scanBam ScanBamParam asBam scanBamFlag
#' @importFrom Rbowtie2 bowtie2
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges findOverlaps pintersect
#' @importFrom GenomicRanges makeGRangesListFromDataFrame width
#'
#' @param RefGenome, path to .fa file
#' @param bowtie2_index, bowtie2 index
#' @param bin_size, size of the final bin in bp: default = 20 000
#' @param read_size, size of a read in bp: default= 40bp
#' @param fragment_size, size of a fragment in bp if paired_ends=T: default= 200bp
#' @param paired_ends, logical, whether reads have to be PE or SE: default= False
#' @param directory_to_bamfiles, path to bam file folder from which estimate read_size, fragment_size and paired_ends.
#' @param tmp_dir, path where to create a temporary folder: default= bins in the wd
#' @param upper_mappability_th, max mappability above which a bin is considered unsuitable for the analysis
#' @param lower_mappability_th, min mappability below which a bin is considered unsuitable for the analysis
#' @param black_list, path to a blacklist file
#' @param coverage, an integer corresponding to the simulated coverage
#' @param errorRate, simulated sequencer percentage error rate
#' @param chr_prefix, prefix identifying a chromosome in the .fa file: default='chr'
#' @param chr_range, chromosomes to consider in the analysis. Use ':' to identify a range and ',' for isolated chromosomes. default='1:22'
#' @param cores, number of cores to use
#' @param return_estimated_param, logical value. If True and directory_to_bamfiles has been passed it returns a list wit the estimated parameters and the calculated bins
#' @export
#'
#'
binning = function(RefGenome,
                   bowtie2_index,
                   bin_size = 20 * 10 ^ 3,
                   read_size = 40,
                   fragment_size = 200,
                   paired_ends = F,
                   directory_to_bamfiles = NULL,
                   tmp_dir = file.path(getwd(), 'bins'),
                   upper_mappability_th = 1.5,
                   lower_mappability_th = 0.8,
                   black_list = NULL,
                   coverage = 1,
                   errorRate = 10 ^ -3,
                   chr_prefix = 'chr',
                   chr_range = NULL,
                   cores = 1,
                   return_estimated_param = F) {
  #load operators
  `%>%` = tidyr::`%>%`
  `%dopar%` = foreach::`%dopar%`
  `%do%` = foreach::`%do%`

  options(scipen = 9999)

  # create temp directory
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }

  #loading reference fa
  reference = Biostrings::readDNAStringSet(RefGenome)

  #declare clusters
  cl = snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))


  #sample 30 files (if available) to estimate parameters
  if (!is.null(directory_to_bamfiles)) {
    list = list.files(path = directory_to_bamfiles,
                      pattern = 'bam$',
                      full.names = T)

    if (length(list) > 30) {
      list = sample(list, 30)
    }

    parameters = foreach::foreach(file = list,
                                  .combine = 'rbind') %dopar% {
                                    sapply(Rsamtools::scanBam(
                                      file = file,
                                      param = Rsamtools::ScanBamParam(what =
                                                                        c('isize', 'qwidth'))
                                    )[[1]],
                                    function(x)
                                      median(abs(x), na.rm = T))
                                  }

    parameters = dplyr::as_tibble(parameters) %>%
      dplyr::summarise(qwidth = round(median(qwidth)),
                       isize = round(median(isize)))
    if (parameters$isize != 0) {
      paired_ends = T
      fragment_size = parameters$isize
      read_size = parameters$qwidth
    } else{
      paired_ends = F
      read_size = parameters$qwidth
    }

  }

  if (return_estimated_param) {
    if (paired_ends) {
      param = list(
        paired_ends = T,
        fragment_size = fragment_size,
        read_size = read_size
      )
    } else{
      param = list(paired_ends = F,
                   read_size = read_size)
    }
  }

  # select chrs of interest
  # convert string into range
  if (is.null(chr_range)) {
    chr_list = unique(names(reference))
  } else{
    chr_list = paste0(ifelse(is.null(chr_prefix), '', chr_prefix),
                      unlist(Kronos.scRT::String_to_Range(stringr::str_split(chr_range, ',',simplify = T))))
    chr_list = chr_list[chr_list %in% unique(names(reference))]
  }
  genome.Chromsizes = foreach::foreach(Chr = chr_list,
                                       .combine = 'rbind') %dopar% {
                                         #load operators
                                         `%>%` = tidyr::`%>%`
                                         `%do%` = foreach::`%do%`

                                         #genome size
                                         genome.Chromsizes = dplyr::tibble(chr = Chr,
                                                                           size = Biostrings::width(reference[Chr]))

                                         position = Kronos.scRT::find_known_sequences(reference[Chr])

                                         #look for seeds
                                         if (paired_ends) {
                                           #initialize simulated reads
                                           size =  fragment_size

                                         } else{
                                           size = read_size

                                         }

                                         Variability = round(seq(-size, size, by = 2 * size / (coverage + 1)))[(1:coverage) +
                                                                                                                 1]

                                         tmp = foreach::foreach(variability = Variability) %do% {
                                           #look for seeds
                                           if (paired_ends) {
                                             #initialize simulated reads
                                             position = position %>%
                                               dplyr::filter(end - start > size)

                                             simulated_reads_1 = foreach::foreach(i = 1:length(position$start),
                                                                                  .combine = 'c') %do% {
                                                                                    seq(position$start[i] + variability, position$end[i], by = size)
                                                                                  }
                                             simulated_reads_1 = dplyr::tibble(start = unlist(simulated_reads_1),
                                                                               end = start + read_size) %>%
                                               dplyr::mutate(order =  dplyr::row_number())
                                             simulated_reads_2 = simulated_reads_1 %>%
                                               dplyr::mutate(end = start + fragment_size,
                                                             start = end - read_size)

                                           } else{
                                             position = position %>%
                                               dplyr::filter(end - start > size)
                                             #initialize simulated reads
                                             simulated_reads = foreach::foreach(i = 1:length(position$start)) %do% {
                                               seq(position$start[i] + variability, position$end[i], by = read_size)
                                             }
                                             simulated_reads = dplyr::tibble(start = unlist(simulated_reads),
                                                                             end = start + read_size) %>%
                                               dplyr::mutate(order =  dplyr::row_number())
                                           }

                                           if (paired_ends) {
                                             #recover strings reads and mutate some of them
                                             simulated_reads_1 = Kronos.scRT::recover_and_mutate(Reads = simulated_reads_1,
                                                                                                 Sequence = reference[Chr],
                                                                                                 errorRate = errorRate)
                                             simulated_reads_2 = Kronos.scRT::recover_and_mutate(Reads = simulated_reads_2,
                                                                                                 Sequence = reference[Chr],
                                                                                                 errorRate = errorRate)


                                             simulated_reads_2 = simulated_reads_2 %>%
                                               dplyr::mutate(reads = Kronos.scRT::rev_com(reads))
                                           } else{
                                             #recover strings reads and mutate some of them
                                             simulated_reads = Kronos.scRT::recover_and_mutate(Reads = simulated_reads,
                                                                                               Sequence = reference[Chr],
                                                                                               errorRate = errorRate)
                                           }

                                           if (paired_ends) {
                                             #save fastq files
                                             Kronos.scRT::reshape_and_save(
                                               Reads = simulated_reads_1,
                                               path = file.path(
                                                 tmp_dir,
                                                 paste0(
                                                   basename(bowtie2_index),
                                                   '_',
                                                   Chr,
                                                   '_simulated_reads_1.fq'
                                                 )
                                               ),
                                               Chrom = Chr,
                                               ID = variability
                                             )
                                             Kronos.scRT::reshape_and_save(
                                               Reads = simulated_reads_2,
                                               path = file.path(
                                                 tmp_dir,
                                                 paste0(
                                                   basename(bowtie2_index),
                                                   '_',
                                                   Chr,
                                                   '_simulated_reads_2.fq'
                                                 )
                                               ),
                                               Chrom = Chr,
                                               ID = variability
                                             )
                                             #remove from memory
                                             rm('simulated_reads_1')
                                             rm('simulated_reads_2')

                                           } else{
                                             #save fastq file
                                             Kronos.scRT::reshape_and_save(
                                               Reads = simulated_reads,
                                               path =  file.path(
                                                 tmp_dir,
                                                 paste0(
                                                   basename(bowtie2_index),
                                                   '_',
                                                   Chr,
                                                   '_simulated_reads.fq'
                                                 )
                                               ),
                                               Chrom = Chr,
                                               ID = variability
                                             )
                                             #remove from memory
                                             rm('simulated_reads')
                                           }
                                           variability
                                         }
                                         genome.Chromsizes
                                       }



  if (paired_ends) {
    #create new file

    Reads1 = file.path(tmp_dir, paste0(basename(bowtie2_index),
                                      '_simulated_reads_1.fq'))
    file.create(Reads1)

    Reads2 = file.path(tmp_dir, paste0(basename(bowtie2_index),
                                       '_simulated_reads_2.fq'))
    file.create(Reads2)

    #combine files
    trash = foreach::foreach(Chr = chr_list) %do% {
      ReadsTMP1 = file.path(tmp_dir,
                            paste0(
                              basename(bowtie2_index),
                              '_',
                              Chr,
                              '_simulated_reads_1.fq'
                            ))
      ReadsTMP2 = file.path(tmp_dir,
                            paste0(
                              basename(bowtie2_index),
                              '_',
                              Chr,
                              '_simulated_reads_2.fq'
                            ))

      file.append(Reads1, ReadsTMP1)
      file.append(Reads2, ReadsTMP2)

      # remove from hd simulated_reads.fq
      file.remove(ReadsTMP1)
      file.remove(ReadsTMP2)
    }
    rm('trash')
    #align with bowtie2
    Rbowtie2::bowtie2(
      bt2Index = bowtie2_index,
      samOutput = file.path(tmp_dir, paste0(
        basename(bowtie2_index), '_simulated_reads.sam'
      )),
      seq1 = Reads1,
      seq2 = Reads2,
      ... = paste0('--phred33 --ignore-quals -p ', cores),
      overwrite = TRUE
    )

  } else{
    #create new file

    Reads = file.path(tmp_dir, paste0(basename(bowtie2_index),
                                      '_simulated_reads.fq'))
    file.create(Reads)


    #combine files
    trash = foreach::foreach(Chr = chr_list) %do% {
      ReadsTMP = file.path(tmp_dir,
                           paste0(
                             basename(bowtie2_index),
                             '_',
                             Chr,
                             '_simulated_reads.fq'
                           ))

      file.append(Reads, ReadsTMP)

      # remove from hd simulated_reads.fq
      file.remove(ReadsTMP)
    }
    rm('trash')
    #align with bowtie2
    Rbowtie2::bowtie2(
      bt2Index = bowtie2_index,
      samOutput = file.path(tmp_dir, paste0(
        basename(bowtie2_index), '_simulated_reads.sam'
      )),
      seq1 = Reads,
      ... = paste0('--phred33 --ignore-quals -p ', cores),
      overwrite = TRUE
    )

  }

  # calculate bins
  dir.bam <-
    Rsamtools::asBam(file.path(tmp_dir, paste0(
      basename(bowtie2_index), '_simulated_reads.sam'
    )))

  # remove from hd simulated_reads.fq
  file.remove(file.path(tmp_dir, paste0(
    basename(bowtie2_index), '_simulated_reads.sam'
  )))

  if (paired_ends) {
    param1 <-
      Rsamtools::ScanBamParam(
        what = c('rname', 'pos', 'isize', 'mapq'),
        flag = Rsamtools::scanBamFlag(hasUnmappedMate = T, isUnmappedQuery = F)
      )
    param2 <-
      Rsamtools::ScanBamParam(
        what = c('rname', 'pos', 'isize', 'mapq', 'mrnm'),
        flag = Rsamtools::scanBamFlag(isPaired = T, isUnmappedQuery = F)
      )
    bins = rbind(
      dplyr::as_tibble(Rsamtools::scanBam(file = dir.bam, param = param1)[[1]]) %>%
        dplyr::filter(mapq >= 30) %>%
        dplyr::select('chr' = rname, pos) %>%
        dplyr::mutate(read = 1),
      dplyr::as_tibble(Rsamtools::scanBam(file = dir.bam, param = param2)[[1]]) %>%
        dplyr::filter(mapq >= 30) %>%
        dplyr::mutate(read = ifelse(
          rname == mrnm &
            abs(isize) < bin_size, 0.5, 1
        )) %>%
        dplyr::select('chr' = rname, pos, read)
    ) %>%
      tidyr::drop_na()
    #parameter used to estimate mappability th
    theoretical_reads = bin_size / fragment_size

  } else{
    param <- Rsamtools::ScanBamParam(
      what = c('rname', 'pos', 'mapq'),
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = F)
    )

    bins = dplyr::as_tibble(Rsamtools::scanBam(file = dir.bam, param = param)[[1]]) %>%
      dplyr::filter(mapq >= 30) %>%
      dplyr::mutate(read = 1) %>%
      dplyr::select('chr' = rname, pos, read)

    #parameter used to estimate mappability th
    theoretical_reads = bin_size / read_size

  }

  bins = foreach::foreach (Chr = genome.Chromsizes$chr,
                           .combine = 'rbind') %do% {
                             size = genome.Chromsizes$size[genome.Chromsizes$chr == Chr]

                             bins_chr = dplyr::tibble(chr = Chr,
                                                      start = seq(0, size, by = bin_size)) %>%
                               dplyr::mutate(end = dplyr::lead(start, n = 1, default =  size))
                             ## calculate reads per bin
                             reads_proper <-
                               bins$pos[bins$chr == Chr &
                                          bins$read == 0.5]
                             reads_notproper <-
                               bins$pos[bins$chr == Chr &
                                          bins$read == 1]
                             if (length(reads_proper) != 0) {
                               reads_proper[reads_proper <= 0] <- 1
                               reads_proper <-
                                 hist(reads_proper,
                                      breaks =  c(1, bins_chr$end),
                                      plot = F)
                               reads_proper <-
                                 reads_proper$counts / 2
                             }
                             if (length(reads_notproper) != 0) {
                               reads_notproper[reads_notproper <= 0] <- 1
                               reads_notproper <-
                                 hist(reads_notproper,
                                      breaks =  c(1, bins_chr$end),
                                      plot = F)
                               reads_notproper <-
                                 reads_notproper$counts
                             }
                             if (length(reads_proper) != 0 &
                                 length(reads_notproper) != 0) {
                               Reads = reads_notproper + reads_proper
                             } else if (length(reads_proper) != 0 &
                                        length(reads_notproper) == 0) {
                               Reads = reads_proper
                             } else if (length(reads_proper) == 0 &
                                        length(reads_notproper) != 0) {
                               Reads = reads_notproper
                             } else{
                               Reads = 0
                             }
                             ## Concatenate
                             bins_chr %>%
                               dplyr::mutate(reads = Reads)
                           }

  bins = bins %>%
    dplyr::mutate(
      mappability = reads / (coverage * theoretical_reads),
      mappability_th = ifelse(
        mappability >= lower_mappability_th &
          mappability <= upper_mappability_th ,
        T,
        F
      )
    ) %>%
    dplyr::select(chr, start, end, mappability, mappability_th)

  #calculate gc % per bin
  bins = foreach::foreach(i = unique(bins$chr), .combine = 'rbind') %dopar% {
    #load operator
    `%>%` = tidyr::`%>%`

    bins %>%
      dplyr::filter(chr == i) %>%
      dplyr::mutate(
        seq = stringr::str_sub(
          string =  reference[names(reference) == i],
          start = start + 1,
          end = end
        ),
        gc_frequency = stringr::str_count(seq, 'G|C') / stringr::str_length(seq)
      ) %>%
      dplyr::select(-seq)
  }



  bins = bins %>%
    dplyr::mutate(type = ifelse(paired_ends, 'PE', 'SE')) %>%
    dplyr::ungroup()

  if (!is.null(black_list)) {
    black_list = readr::read_tsv(black_list,col_names = c('chr','start','end')) %>%
      GenomicRanges::makeGRangesFromDataFrame()

    tbins = bins %>% GenomicRanges::makeGRangesFromDataFrame()

    hits = IRanges::findOverlaps(query = tbins, subject = black_list)
    overlaps = IRanges::pintersect(black_list[S4Vectors::subjectHits(hits)], tbins[S4Vectors::queryHits(hits)])
    overlaps = overlaps[GenomicRanges::width(overlaps) > read_size]
    bins = overlaps %>%
      dplyr::as_tibble() %>%
      dplyr::select('chr' = seqnames, start, end, hit) %>%
      dplyr::right_join(bins, by = c("chr", "start", "end"))

    bins = bins %>%
      dplyr::mutate(mappability_th = ifelse(!is.na(hit), F, mappability_th)) %>%
      dplyr::select(-hit) %>%
      dplyr::arrange(chr, start)

  }

  #remove folder
  unlink(x = tmp_dir, recursive = T)

  if (return_estimated_param) {
    return(list(bins = bins, param = param))
  } else{
    return(bins)
  }

}
