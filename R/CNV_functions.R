#' Calls CNV over genomic bins defined by binning
#'
#' @return  list of dataframes
#'
#' @importFrom  tidyr %>% drop_na
#' @importFrom  dplyr n filter lag rename as_tibble tibble arrange bind_cols full_join group_by inner_join left_join mutate nth pull rowwise summarise ungroup
#' @importFrom  foreach %do% %:% %dopar% foreach
#' @importFrom  doSNOW registerDoSNOW
#' @importFrom  snow stopCluster makeCluster
#' @importFrom  stringr str_remove str_split
#' @importFrom  Rsamtools ScanBamParam scanBamFlag scanBam
#' @importFrom  DNAcopy CNA smooth.CNA segment
#' @importFrom  MASS fitdistr
#' @importFrom readr cols read_tsv write_tsv
#' @importFrom matrixStats weightedSd
#'
#' @export
#'
#' @param directory, Single cell Bamfiles directory
#' @param bins, bins produced by binning
#' @param chrom_size, a dataframe (chr=chrom name, size=chrom size)
#' @param basename, sample name
#' @param group, experimental group name
#' @param tmp_dir, temporary folder.
#' @param min_n_reads, Min n of reads to keep a cell in the analysis
#' @param mim_mean_CN_accepted, Min mean CN accepted as result
#' @param max_mean_CN_accepted, Max mean CN accepted as result
#' @param ploidy, user extimated ploidy
#' @param chr_prefix, Chromosome prefix, if there is no prefix use NULL
#' @param chr_range, Chromosomes to consider in the analysis (example 1:5,8,15:18,X)
#' @param cores, Number of cores to use
#'
#'

CallCNV = function(directory,
                   bins,
                   chrom_size,
                   basename = 'Exp',
                   group = NULL,
                   tmp_dir = file.path(getwd(), 'CNV'),
                   min_n_reads = 200000,
                   mim_mean_CN_accepted = 2,
                   max_mean_CN_accepted = 8,
                   ploidy = NULL,
                   chr_prefix = 'chr',
                   chr_range = NULL,
                   cores = 1) {
  #load operators
  `%>%` = tidyr::`%>%`
  `%do%` = foreach::`%do%`
  `%:%` = foreach::`%:%`
  `%dopar%` = foreach::`%dopar%`

  #if group is null group=basename
  if (is.null(group)) {
    group = basename
  }

  #stop if binning file has not the right format
  if (!all(
    c(
      'chr',
      'start',
      'end',
      'mappability',
      'mappability_th',
      'gc_frequency',
      'type'
    ) %in% colnames(bins)
  )) {
    stop('a non appropriate bins dataframe was provided')
  }

  #check directory and load files paths
  if (dir.exists(directory)) {
    files = list.files(directory, pattern = '.bam$', full.names = T)
    if (length(files) == 0) {
      stop(paste0(directory, " does not contain bam files."))
    }
  } else{
    stop(paste0(directory, " does not exist."))
  }

  #create temp folder
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }

  #calcualte genome size
  genome_size = sum(chrom_size$size)

  #select chrs
  if (is.null(chr_range)) {
    chromosome = unique(bins$chr)
  } else{
    chromosome = paste0(ifelse(is.null(chr_prefix), '', chr_prefix),
                        unlist(Kronos.scRT::String_to_Range(stringr::str_split(chr_range, ',')[[1]])))
    bins = bins[bins$chr %in% chromosome,]
  }

  #select chromosomes of interest
  chrom_size = chrom_size %>%
    dplyr::mutate(chr = factor(chr, levels = chromosome)) %>%
    tidyr::drop_na()

  # SE or PE ?
  type = bins %>%
    dplyr::select(type) %>%
    unique() %>%
    dplyr::pull()

  #declair cluster
  cl = snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))

  bins = bins %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(bin = 1:dplyr::n())

  bins_median_size = median(bins$end - bins$start)

  # calculating profile
  files = foreach::foreach (file = files,
                            .combine = 'rbind') %dopar% {
                              #load operators
                              `%>%` = tidyr::`%>%`
                              `%do%` = foreach::`%do%`

                              if (type == 'PE') {
                                # Single reads
                                param1 <- Rsamtools::ScanBamParam(
                                  what = c('qname','rname', 'pos', 'mapq', 'qwidth', 'strand'),
                                  flag = Rsamtools::scanBamFlag(
                                    isPaired = F,
                                    isUnmappedQuery = F,
                                    isDuplicate = F
                                  ),
                                  mapqFilter = 30
                                )

                                #First in a pair proper paired reads
                                param2 <-
                                  Rsamtools::ScanBamParam(
                                    what = c('qname', 'rname', 'pos', 'mapq', 'qwidth', 'strand'),
                                    flag = Rsamtools::scanBamFlag(
                                      isPaired = T,
                                      isUnmappedQuery = F,
                                      isFirstMateRead = T,
                                      isDuplicate = F
                                    ),
                                    mapqFilter = 30
                                  )

                                #Second in a pair proper paired reads
                                param3 <-
                                  Rsamtools::ScanBamParam(
                                    what = c('qname', 'rname', 'pos', 'mapq', 'qwidth', 'strand'),
                                    flag = Rsamtools::scanBamFlag(
                                      isPaired = T,
                                      isUnmappedQuery = F,
                                      isSecondMateRead = T,
                                      isDuplicate = F
                                    ),
                                    mapqFilter = 30
                                  )
                                # load first in a pair
                                FP = dplyr::as_tibble(Rsamtools::scanBam(file, param = param2)[[1]]) %>%
                                  tidyr::drop_na() %>%
                                  dplyr::filter(rname %in% chromosome) %>%
                                  dplyr::group_by(rname) %>%
                                  dplyr::mutate(bin = ceiling(pos / bins_median_size))
                                # load second in a pair
                                SP = as.data.frame(Rsamtools::scanBam(file, param = param3)[[1]]) %>%
                                  tidyr::drop_na() %>%
                                  dplyr::filter(rname %in% chromosome) %>%
                                  dplyr::group_by(rname) %>%
                                  dplyr::mutate(bin = ceiling(pos / bins_median_size))

                                # remube duplicates
                                #calculate the cumulative size based on chr order. merge with data and add it to the position of the read
                                if (nrow(FP) > 0 | nrow(SP) > 0) {
                                  Deduplicated_reads_list = dplyr::full_join(
                                    FP %>% dplyr::inner_join(
                                      chrom_size %>% dplyr::mutate(size = cumsum(size),
                                                                   size =
                                                                     size - min(size)),
                                      by = c('rname' = 'chr')
                                    ),
                                    SP %>% dplyr::inner_join(
                                      chrom_size %>% dplyr::mutate(size = cumsum(size),
                                                                   size =
                                                                     size - min(size)),
                                      by = c('rname' = 'chr')
                                    ),
                                    by = c('qname')
                                  ) %>%
                                    #if position is sin the minus strand it is actully the last base of the read
                                    dplyr::mutate(
                                      pos.x = ifelse(strand.x == '+', pos.x + size.x, pos.x + size.x + qwidth.x),
                                      pos.y = ifelse(strand.y == '+', pos.y + size.y, pos.y +
                                                       size.y + qwidth.y)
                                    ) %>%
                                    dplyr::rowwise() %>%
                                    dplyr::mutate(start = min(pos.x, pos.y, na.rm = T),
                                                  #if the mate was not mapped, use qwidth to ide the end of the read
                                                  end = max(pos.x, pos.y, na.rm = T)) %>%
                                    dplyr::group_by(start, end) %>%
                                    #for duplicated reads keep either of the two
                                    dplyr::summarise(read_to_keep = sample(qname, 1)) %>% dplyr::pull(read_to_keep)

                                  # select identifiers reads
                                  qname_fp = FP %>%
                                    dplyr::ungroup() %>%
                                    dplyr::select(qname, rname,  'mate_bin' = bin)
                                  qname_sp = SP %>%
                                    dplyr::ungroup() %>%
                                    dplyr::select(qname, rname, 'mate_bin' = bin)

                                  # if a read in a pair has it's paired in the same bin it counts as 1/2 a read eles one read.
                                  FP = FP %>%
                                    dplyr::filter(qname %in% Deduplicated_reads_list) %>%
                                    dplyr::left_join(qname_sp, by = c("qname", 'rname')) %>%
                                    dplyr::mutate(read = ifelse(bin != mate_bin |
                                                                  is.na(mate_bin), 1, 0.5)) %>%
                                    dplyr::select('chr' = rname, pos, read, bin) %>%
                                    dplyr::ungroup()

                                  SP = SP %>%
                                    dplyr::filter(qname %in% Deduplicated_reads_list) %>%
                                    dplyr::left_join(qname_fp, by = c("qname", 'rname')) %>%
                                    dplyr::mutate(read = ifelse(bin != mate_bin |
                                                                  is.na(mate_bin), 1, 0.5)) %>%
                                    dplyr::select('chr' = rname, pos, read, bin) %>%
                                    dplyr::ungroup()

                                }
                                #load non paired reads
                                SR = dplyr::as_tibble(Rsamtools::scanBam(file, param = param1)[[1]]) %>%
                                  tidyr::drop_na() %>%
                                  dplyr::filter(rname %in% chromosome)

                                Deduplicated_SR_list = SR %>%
                                  dplyr::mutate(pos = ifelse(strand == '+', pos, pos + qwidth)) %>%
                                  dplyr::group_by(rname, pos) %>%
                                  #for duplicated reads keep either of the two
                                  dplyr::summarise(read_to_keep = sample(qname, 1)) %>% dplyr::pull(read_to_keep)

                                SR = SR %>%
                                  dplyr::filter(qname %in% Deduplicated_SR_list) %>%
                                  dplyr::mutate(read = 1) %>%
                                  dplyr::rename('chr' = rname) %>%
                                  dplyr::ungroup()


                                #total reads
                                count_reads = nrow(FP) + nrow(SP) + nrow(SR)

                                # include non paired reads if existing
                                # sum all the reads in a bin
                                if (length(SR$chr) != 0) {
                                  sam = rbind(SR = SR  %>%
                                                dplyr::group_by(chr) %>%
                                                dplyr::mutate(bin = ceiling(pos / bins_median_size)),
                                              FP ,
                                              SP) %>%
                                    dplyr::group_by(chr, bin) %>%
                                    dplyr::summarise(reads = sum(read)) %>%
                                    dplyr::ungroup()

                                } else{
                                  sam = rbind(FP ,
                                              SP) %>%
                                    dplyr::group_by(chr, bin) %>%
                                    dplyr::summarise(reads = sum(read)) %>%
                                    dplyr::ungroup()
                                }
                              } else if (type == 'SE') {
                                param <- Rsamtools::ScanBamParam(
                                  what = c('qname','rname', 'pos', 'mapq', 'strand','qwidth'),
                                  flag = Rsamtools::scanBamFlag(isUnmappedQuery = F,
                                                                isDuplicate = F),
                                  mapqFilter = 30
                                )

                                # calculate reads in each bin
                                sam = dplyr::as_tibble(Rsamtools::scanBam(file, param =
                                                                            param)[[1]]) %>%
                                  dplyr::filter(rname %in% chromosome)

                                Deduplicated_sam_list = sam %>%
                                  dplyr::mutate(pos = ifelse(strand == '+', pos, pos + qwidth)) %>%
                                  dplyr::group_by(rname, pos) %>%
                                  #for duplicated reads keep either of the two
                                  dplyr::summarise(read_to_keep = sample(qname, 1)) %>% dplyr::pull(read_to_keep)

                                sam = sam %>%
                                  dplyr::filter(qname %in% Deduplicated_sam_list) %>%
                                  dplyr::mutate(read = 1) %>%
                                  dplyr::rename('chr' = rname) %>%
                                  dplyr::group_by(chr) %>%
                                  dplyr::mutate(bin = ceiling(pos / bins_median_size)) %>%
                                  dplyr::group_by(chr, bin) %>%
                                  dplyr::summarise(reads = sum(read)) %>%
                                  dplyr::ungroup()

                                count_reads = sum(sam$reads)

                              } else{
                                stop('Bins file does not contain a correct sequencing type (SE/PE)')
                              }
                              if (count_reads >= min_n_reads) {
                                # save one file per cell with a coverage track
                                bins %>%
                                  dplyr::left_join(sam %>%
                                                     dplyr::mutate(chr = as.character(chr)),
                                                   by = c('chr', 'bin')) %>%
                                  dplyr::mutate(reads = ifelse(is.na(reads), 0, reads),
                                                Cell = file) %>%
                                  dplyr::select(-bin) %>%
                                  readr::write_tsv(file = file.path(tmp_dir, paste0(basename(file), '.tmp')))
                                dplyr::tibble(file = paste0(basename(file), '.tmp'),
                                              count_reads = count_reads)
                              }else{
                                dplyr::tibble()
                              }

                            }

  # correct for mappability and normalize for GC content based on all the cells
  #(norm reads= reads * median reads  / median reads per interval of GC (rounded at 2 digits))
  data = bins

  data$reads = foreach::foreach(file = files$file,
                                .combine = '+') %dopar% {
                                  #load operators
                                  `%>%` = tidyr::`%>%`

                                  readr::read_tsv(file.path(tmp_dir, file), col_types = readr::cols(chr = 'c')) %>%
                                    dplyr::pull(reads)
                                }


  gc_correction_value = data %>%
    dplyr::filter(mappability_th) %>%
    dplyr::group_by(chr, start, end, mappability)  %>%
    dplyr::mutate(
      gc_frequency = round(gc_frequency, 1),
      reads_mappability = reads / mappability
    ) %>%
    dplyr::group_by(gc_frequency) %>%
    dplyr::mutate(mGC = median(reads_mappability)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      m = median(reads_mappability),
      gc_corretion_values = ifelse(mGC == 0, NA, m / mGC)
    ) %>%
    dplyr::select(chr, start, end, gc_corretion_values)

  mapd = foreach::foreach (file = files$file,
                           .combine = 'rbind',
                           .errorhandling =  "remove") %dopar% {
                             #load operators
                             `%>%` = tidyr::`%>%`
                             `%do%`=foreach::`%do%`

                             data = readr::read_tsv(file.path(tmp_dir, file), col_types = readr::cols(chr = 'c'))
                             data = dplyr::left_join(data, gc_correction_value, by = c('chr', 'start', 'end')) %>%
                               dplyr::mutate(
                                 reads_mappability = reads / mappability,
                                 gc_corrected_reads = reads_mappability * gc_corretion_values
                               )

                             data_500Kb = data %>%
                               dplyr::group_by(chr) %>%
                               dplyr::mutate(range = ceiling(end / 500000)) %>%
                               dplyr::group_by(chr, range) %>%
                               dplyr::summarise(
                                 start = min(start),
                                 end = max(end),
                                 gc_corrected_reads = sum(gc_corrected_reads, na.rm = T)
                               )

                             CovReadsMega = files[files$file == file,] %>%
                               dplyr::summarise(coverage = 1000000 * count_reads / genome_size) %>%
                               dplyr::pull(coverage)

                             #calculate normalized_MAPD and mapd
                             mapd = data_500Kb %>%
                               dplyr::group_by(chr) %>%
                               dplyr::mutate(read_1n = dplyr::lag(gc_corrected_reads, 1)) %>%
                               tidyr::drop_na() %>%
                               dplyr::ungroup() %>%
                               dplyr::mutate(
                                 mean_n = mean(gc_corrected_reads),
                                 normalized_mapd = (gc_corrected_reads - read_1n) / mean_n
                               ) %>%
                               dplyr::summarise(
                                 normalized_mapd = median(abs(
                                   normalized_mapd - median(normalized_mapd)
                                 )),
                                 coverage = 1000000 * sum(gc_corrected_reads) / genome_size,
                                 normalized_dimapd = normalized_mapd * sqrt(coverage)
                               ) %>%
                               dplyr::mutate(CovReadsMega =  CovReadsMega)
                             #spread data for segmentation
                             data = data %>%
                               dplyr::filter(mappability_th) %>%
                               dplyr::select(chr, start, end, gc_corrected_reads, Cell) %>%
                               dplyr::ungroup()

                             # create object
                             CNA.object <-
                               DNAcopy::CNA(as.matrix(data$gc_corrected_reads),
                                            data$chr,
                                            data$start,
                                            sampleid = file)

                             # smooth data
                             smoothed.CNA.object <-
                               DNAcopy::smooth.CNA(CNA.object)

                             # free memory
                             rm('CNA.object')

                             # segment
                             segment.smoothed.CNA.object <-
                               DNAcopy::segment(smoothed.CNA.object)

                             #free memory
                             rm('smoothed.CNA.object')

                             segment.smoothed.CNA.object = dplyr::as_tibble(segment.smoothed.CNA.object$output)

                             #free memory
                             rm('data')

                             bin_size = bins[1,] %>%
                               dplyr::mutate(bs = end - start) %>%
                               dplyr::pull(bs)

                             # identify CN based on minimum of the target function

                             weitghts = (segment.smoothed.CNA.object$loc.end - segment.smoothed.CNA.object$loc.start) / bin_size

                             possible_factors = foreach::foreach(i = seq(0.1, 1000 , 0.1),
                                                                 .combine = 'rbind') %do% {
                                                                   TargetF = sqrt(sum((
                                                                     weitghts * sinpi(segment.smoothed.CNA.object$seg.mean / i) ^ 2
                                                                   )))

                                                                   mean_cn = weighted.mean(round(segment.smoothed.CNA.object$seg.mean / i),
                                                                                           weitghts)

                                                                   Variability = 100 * matrixStats::weightedSd(x=segment.smoothed.CNA.object$seg.mean,w=weitghts) / weighted.mean(segment.smoothed.CNA.object$seg.mean, weitghts)
                                                                   TargetF = dplyr::tibble(
                                                                     possible_factors = TargetF,
                                                                     X = i,
                                                                     mean_cn = mean_cn,
                                                                     Variability = Variability
                                                                   )

                                                                 }

                             Var = unique(possible_factors$Variability)
                             min = possible_factors$possible_factors[which(diff(sign(
                               diff(possible_factors$possible_factors)
                             )) == 2) + 1]


                             if (!is.null(ploidy)) {
                               possible_factors = possible_factors %>%
                                 dplyr::filter(possible_factors %in% min)


                               if (sum((
                                 possible_factors$mean_cn >= ploidy / 1.5 &
                                 possible_factors$mean_cn <= ploidy * 2
                               )
                               ) == 0) {
                                 selected = possible_factors$X[which(abs(possible_factors$mean_cn -
                                                                           ploidy) == min(abs(
                                                                             possible_factors$mean_cn - ploidy
                                                                           )))]

                                 mean_cn = possible_factors$mean_cn[which(abs(possible_factors$mean_cn -
                                                                                ploidy) == min(abs(
                                                                                  possible_factors$mean_cn - ploidy
                                                                                )))]
                                 PloConf = -200
                               } else{
                                 possible_factors = possible_factors %>%
                                   dplyr::filter(possible_factors %in% min,
                                                 mean_cn <= ploidy * 2,
                                                 mean_cn >= ploidy / 1.5)

                                 selected = min(possible_factors$possible_factors)
                                 mean_cn = possible_factors$mean_cn[possible_factors$possible_factors ==
                                                                      selected]
                                 selected = possible_factors$X[possible_factors$possible_factors ==
                                                                 selected]
                                 PloConf = -100
                               }

                             } else{
                               possible_factors = possible_factors %>%
                                 dplyr::filter(
                                   possible_factors %in% min,
                                   mean_cn <= max_mean_CN_accepted,
                                   mean_cn >= mim_mean_CN_accepted
                                 )

                               if (Var < 5) {
                                 selected = possible_factors$X[which(abs(possible_factors$mean_cn -
                                                                           2) == min(abs(possible_factors$mean_cn - 2)))]
                                 mean_cn = possible_factors$mean_cn[which(abs(possible_factors$mean_cn -
                                                                                2) == min(abs(possible_factors$mean_cn - 2)))]
                                 PloConf = -2
                               } else{
                                 selected = min(possible_factors$possible_factors)
                                 if(nrow(possible_factors)>1){
                                 PloConf = dplyr::nth(possible_factors$possible_factors[base::order(possible_factors$possible_factors, decreasing = T)], n = -2) -
                                   selected
                                 }else{
                                   PloConf = 0
                                 }
                                 mean_cn = possible_factors$mean_cn[possible_factors$possible_factors ==
                                                                      selected]
                                 selected = possible_factors$X[possible_factors$possible_factors ==
                                                                 selected]

                               }
                             }
                             CNV_correction = dplyr::tibble(
                               ID = unique(segment.smoothed.CNA.object$ID),
                               Cell = stringr::str_remove(file, '.tmp$'),
                               X = selected,
                               ploidy_confidence = PloConf,
                               mean_ploidy = mean_cn
                             )

                             # called CNV
                             CNV = segment.smoothed.CNA.object %>%
                               dplyr::inner_join(CNV_correction, by = 'ID') %>%
                               dplyr::mutate(CNV = round(seg.mean / X, 0)) %>%
                               dplyr::select(
                                 Cell,
                                 'chr' = chrom,
                                 'start' = loc.start,
                                 'end' = loc.end,
                                 'copy_number' = CNV,
                                 'reads' = seg.mean
                               )


                             CNV_correction = CNV_correction %>%
                               dplyr::select(-X)


                             CNV %>%
                               readr::write_tsv(file.path(tmp_dir, paste0(file, '_cnv_calls.bed')),
                                                col_names = T)

                             file.remove(file.path(tmp_dir, file))

                             # calculate mean CNV
                             mapd = mapd %>%
                               dplyr::bind_cols(CNV_correction) %>%
                               dplyr::select(-ID)

                             mapd

                           }

  #calculate DiMApd
  mapd = mapd %>%
    dplyr::mutate(d = coverage - median(coverage))

  LM_mapd_coverage = lm(formula = normalized_dimapd ~ d , data = mapd)

  mapd = mapd %>%
    dplyr::mutate(
      normalized_dimapd = normalized_dimapd - d * LM_mapd_coverage$coefficients[[2]],
      is_high_dimapd = F
    ) %>%
    dplyr::select(-d)

  #fit dimapd to gaussian dist
  mem = 0
  while (T) {
    fit <-
      MASS::fitdistr(mapd$normalized_dimapd[!mapd$is_high_dimapd], 'normal')

    mapd = mapd %>%
      dplyr::mutate(
        is_high_dimapd = pnorm(
          normalized_dimapd,
          mean = fit$estimate[[1]],
          sd = fit$estimate[[2]],
          lower.tail = F
        ) <= 0.01
      )

    if (mem == sum(mapd$is_high_dimapd)) {
      break
    } else{
      mem = sum(mapd$is_high_dimapd)
    }
  }

  #load CNV files
  CNV = foreach::foreach(file = files$file,
                         .combine = 'rbind',
                         .inorder = T) %dopar% {
                           readr::read_tsv(file.path(tmp_dir, paste0(file, '_cnv_calls.bed'))) %>%
                             dplyr::arrange(chr, start, end)

                         }

  CNV = CNV %>%
    dplyr::mutate(basename = basename,
                  group = group)


  #delete tmp files
  sapply(file.path(tmp_dir, paste0(files$file, '_cnv_calls.bed')), file.remove)

  #reshpe percell DF
  mapd = mapd %>%
    dplyr::mutate(
      is_noisy = ifelse(
        is_high_dimapd |
          (ploidy_confidence < 2 & ploidy_confidence != -100),
        T,
        F
      ),
      coverage_per_1Mbp = CovReadsMega,
      Cell = stringr::str_remove(Cell, '.tmp$')
    ) %>%
    dplyr::select(
      Cell,
      normalized_dimapd,
      mean_ploidy,
      ploidy_confidence,
      is_high_dimapd,
      is_noisy,
      coverage_per_1Mbp
    ) %>%
    dplyr::mutate(basename = basename,
                  group = group)

  #delete temp folder
  unlink(x = tmp_dir, recursive = T)

  #return
  return(list(PerCell = mapd, CNV = CNV))

}
