#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for MCF7-subpop1 cells. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 1192 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 1366469 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"MCF7_subpop1"

#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for MCF7-subpop2 cells. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 503 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 562621 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"MCF7_subpop2"

#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for MCF7 cells before identifying its two subpopulations. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 2321 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 2501060 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"MCF7"

#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for Jeff cell replicate 1. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 952 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 579410 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"Jeff_rep1"

#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for Jeff cell replicate 2. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 503 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 562621 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"Jeff_rep2"

#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for HeLa S3 cell replicate 1. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 459 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 291039 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"HeLa_rep1"

#' @title SingneCell RT Test object
#'
#' @description A \code{list} with two \code{tibble}: PerCell and CNV for HeLa S3 cell replicate 2. Assembly Hg38.
#'
#' @format A list containing two tibbles variables:
#' \describe{
#' \itemize{
#'   \item{PerCell}{Per Cell metrics: a tibble containing 293 rows and 6 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{normalized_dimapd}{normalized Depth Independent Median Absolute deviation of Pairwise Differences}
#'   \item{mean_ploidy}{Single Cell mean ploidy}
#'   \item{ploidy_confidence}{Confidence in the CN calling}
#'   \item{is_high_dimapd}{High DIMAPD cells have high internal variability and could be S-phase cells}
#'   \item{is_noisy}{noisy cells have either poor ploidy confidence or high DIMAPD}
#'   \item{coverage_per_1Mbp}{Reads coverage per megabase}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#'   }
#'   \itemize{
#'   \item{CNV}{Single cell copy number 310431 rows and 8 columns}
#'   \itemize{
#'   \item{Cell}{Single Cell identifier}
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{copy_number}{Copy number Called for a bin}
#'   \item{reads}{Normalised read count inside a bin}
#'   \item{basename}{Sample identifier}
#'   \item{group}{Experimental group identifier}
#'   }
#' }
#' }
#' @source \url{https://www.nature.com/articles/s41467-022-30043-x}
"HeLa_rep2"


#' @title MCF7 Reference RT
#'
#' @description A \code{tibble} containing bulk replication timing data derived from a Hg19 -> Hg38 liftover.
#'
#' @format tibble with 2831095 rows and 4 columns.
#' \describe{
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{RT}{Replication timing data for a certain bin}
#'   }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM923442}
"MCF7_Reference"

#' @title GM12878 Reference RT
#'
#' @description A \code{tibble} containing bulk replication timing data derived from a Hg19 -> Hg38 liftover.
#'
#' @format tibble with 2831101 rows and 4 columns.
#' \describe{
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{RT}{Replication timing data for a certain bin}
#'   }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM923451}
"GM12878_Reference"

#' @title HeLa S3 Reference RT
#'
#' @description A \code{tibble} containing bulk replication timing data derived from a Hg19 -> Hg38 liftover.
#'
#' @format tibble with 2831129 rows and 4 columns.
#' \describe{
#'   \item{chr}{Chromosome position of a bin}
#'   \item{start}{Start position of a bin}
#'   \item{end}{End position of a bin}
#'   \item{RT}{Replication timing data for a certain bin}
#'   }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM923449}
"HeLa_Reference"
