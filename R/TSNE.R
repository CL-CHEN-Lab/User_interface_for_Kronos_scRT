#' adds an option to set a seed in Rtsnee
#'
#' @return Rtsne object
#' @importFrom Rtsne Rtsne
#'
#' @export
#'
#' @param X, matrix; Data matrix (each row is an observation, each column is a variable)
#' @param seed, an integer
#' @param ..., other parameters that can be passed to Rtsne
#'
#' @examples
#' \dontrun{
#' TSNE(X,seed,...)
#' }
#'

TSNE = function(X, seed, ...) {
  #reset seed on exit
  old_seed = .Random.seed
  on.exit(.Random.seed <<- old_seed)
  # set seed locally
  set.seed(as.integer(seed))
  #run tSNE
  tsne <-
    Rtsne::Rtsne(X = X, ...)
  return(tsne)
}
