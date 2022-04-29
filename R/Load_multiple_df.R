#'  Loads and binds a dataframe divided in multiple files
#'
#' @return a string informing the user file format is the required one
#'
#' @importFrom  readr read_delim
#' @importFrom  dplyr tibble
#'
#' @param dirs, list of directories
#' @param delim, delimiter for read_delim
#' @param ..., other parameters to pass to read_delim
#'
#' @export

load_multiple_df = function(dirs, delim = '\t', ...) {
  df = lapply(dirs,
              function(x)
                tryCatch(
                  readr::read_delim(file = x, delim = delim, ...),
                  error = function(x)
                    dplyr::tibble()
                ))


  df = do.call(rbind, df)
  return(df)
}
