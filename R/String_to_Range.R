#' converts a strin representing a range of values into an array of integers
#'
#' @return an array of integers ranging from the first to the last number in the string
#'
#' @importFrom  stringr str_detect str_split
#'
#' @param x, a string
#' @param div, symbol used to indicate a ragnge
#'
#'@export
#'
String_to_Range = Vectorize(function(x, div = ':') {
  if (stringr::str_detect(x, div)) {
    x = stringr::str_split(x, div)[[1]]
    return(as.numeric(x[1]):as.numeric(x[2]))
  } else{
    return(x)
  }
}, vectorize.args = 'x')
