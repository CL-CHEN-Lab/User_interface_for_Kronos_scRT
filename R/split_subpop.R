#'  Splits PerCell data frame based on identified sub-populations
#'
#' @return list
#'
#' @importFrom  tidyr %>%
#' @importFrom  dplyr select group_split inner_join
#'
#' @param PerCell, Kornos PerCell data frame produced by Kronos
#' @param Supbop, Sub-population grouping created by Kronos
#'
#' @export
#'
split_subpop = function(PerCell,
                        Supbop,
                        group = NULL,
                        basename = NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`

  if (!is.null(group)) {
    Supbop = Supbop[Supbop$group == group, ]
  }

  if (!is.null(basename)) {
    Supbop = Supbop[Supbop$basename == basename, ]
  }

  result = dplyr::inner_join(PerCell, Supbop, by = 'Cell') %>%
    dplyr::select(-group, -basename) %>%
    dplyr::group_split(subpopulation)

  return(result)

}

#'  Splits single cell CN dataframe  based on identified sub-populations
#'
#' @return list
#'
#' @importFrom  tidyr %>%
#' @importFrom  dplyr mutate select inner_join
#'
#' @param scCN, single cell CN data frame produced by Kronos
#' @param Supbop, Sub-population grouping created by Kronos
#'
#' @export
#'
Assign_subpop = function(scCN,
                         Supbop,
                         group = NULL,
                         basename = NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`

  if (!is.null(group)) {
    Supbop = Supbop[Supbop$group == group, ]
  }

  if (!is.null(basename)) {
    Supbop = Supbop[Supbop$basename == basename, ]
  }

  result = dplyr::inner_join(PerCell, Supbop, by = c('Cell', 'group', 'basename')) %>%
    dplyr::mutate(group = paste(group, subpopulation, sep = '_')) %>%
    dplyr::select(-subpopulation)

  return(result)

}
