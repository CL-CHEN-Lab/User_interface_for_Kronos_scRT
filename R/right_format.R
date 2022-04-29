#'  checks whether a file contains the right number of columns and/or the right colums names
#'
#' @return a string informing the user whether a file format is the required one
#'
#' @importFrom  readr read_delim cols
#' @importFrom dplyr tibble
#'
#' @param file_path, path to a file of interest
#' @param columns_to_check, either a list of columns name or the number of colums that are supposed to be in the file
#' @param delim, delimiter for read_delim
#' @param skip, number of rows to skip for read_delim
#' @param wrong_message, message returned if the format is not the one requested
#' @param rigth_message, message returned if the format is requested one
#'
#' @export




right_format = Vectorize(function(file_path,
                                  columns_to_check,
                                  delim = '\t',
                                  skip = 0,
                                  wrong_message = 'does not have the proper format \n',
                                  rigth_message = '',
                                  logical = F) {
  # if columns_to_check is numeric check the number of colums
  if (is.numeric(columns_to_check)) {
    if (ncol(tryCatch(
      expr =  readr::read_delim(
        file_path,
        col_types = readr::cols(),
        n_max = 0,
        delim = delim,
        skip = skip
      ),
      error = function(x)
        dplyr::tibble()
    )) != columns_to_check) {
      return(ifelse(logical, F, wrong_message))
    } else{
      return(ifelse(logical, T, rigth_message))
    }

    # if columns_to_check is not numeric check columns names
  } else{
    #checks if it has the right format
    if (!all(columns_to_check %in% colnames(tryCatch(
      expr =  readr::read_delim(
        file_path,
        col_types = readr::cols(),
        n_max = 0,
        delim = delim,
        skip = skip
      ),
      error = function(x)
        dplyr::tibble()
    )))) {
      return(ifelse(logical, F, wrong_message))
    } else{
      return(ifelse(logical, T, rigth_message))
    }
  }

}, vectorize.args = 'file_path')
