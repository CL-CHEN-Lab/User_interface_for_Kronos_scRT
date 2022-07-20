#'  Hclustering to identify RT signatures
#'
#' @return list
#'
#' @importFrom stringr str_locate_all
#' @importFrom dplyr select mutate group_by n all_of summarise
#' @importFrom tidyr spread drop_na %>% gather unite
#' @importFrom foreach foreach %:% %do%
#'
#' @param ..., bulk/pseudo-bulk RTs dataframes (chr,start,end,RT,group)
#' @param deltaRT_th, min RT change for a region to be considered
#' @param CrossingRT, IF true, only changes that are crossing the mid RT value will be considered
#' @param n_clusters, number of clusters
#' @param colors, min, midpoint and max colors
#'
#' @export
#'

RT_clustering=function(...,deltaRT_th=0.1,CrossingRT=T,n_clusters=NULL,colors=NULL){

  #define operator
  `%>%`=tidyr::`%>%`
  `%:%`=foreach::`%:%`
  `%do%`=foreach::`%do%`

  #set colors
  if(is.null(colors)){
    colors=c('#005095', 'white','#a7001b')
  }
  #recover inputs
  data=list(...)
  data <- do.call('rbind',data)

  data = tryCatch(
  data %>%
    dplyr::select(chr, start, end, RT, basename) %>%
    tidyr::spread(basename, RT) %>%
    tidyr::drop_na(),
  error= function(e) data %>%
    tidyr::unite(name,basename,group,sep = ' - ')%>%
    dplyr::select(chr, start, end, RT, name) %>%
    tidyr::spread(name, RT) %>%
    tidyr::drop_na())

  Groups = names(data)[!names(data) %in% c('chr', 'start', 'end')]

  #calculate deltaRT
  distance = foreach::foreach(names1 = 1:(length(Groups) - 1), .combine = cbind) %:%
    foreach::foreach(names2 = (names1 + 1):length(Groups),
            .combine = cbind) %do% {
              data[Groups[names1]] - data[Groups[names2]]
            }


keep=abs(distance)>=deltaRT_th
keep=rowSums(keep)>0

if (CrossingRT){

  crossing = foreach::foreach(names1 = 1:(length(Groups) - 1), .combine = cbind) %:%
    foreach::foreach(names2 = (names1 + 1):length(Groups),
                     .combine = cbind) %do% {
                       sign((0.5 - data[Groups[names1]]) / (0.5 - data[Groups[names2]]))
                     }

  crossing=crossing==-1
  crossing=rowSums(crossing)>0
  keep=  crossing & keep

}

#if the number of clusters is not defined
if (is.null(n_clusters)){
  n_clusters=(2^length(Groups))-2
}

data=data[keep,]
clusters <- hclust(dist(data[Groups]))
data$clusters <-  paste('Cluster', cutree(clusters, n_clusters))


data=data%>%
  dplyr::mutate(clusters=factor(clusters, levels =paste('Cluster', 1:n_clusters) ))%>%
  dplyr::group_by(clusters)%>%
  dplyr::mutate(n=1:dplyr::n())%>%
  tidyr::gather(name,RT,dplyr::all_of(Groups))

p=data%>%
  ggplot2::ggplot(ggplot2::aes(x=name,y=n,hight=1,width=1,fill=RT))+
  ggplot2::geom_raster()+
  ggplot2::facet_grid(clusters ~ name,scales = 'free',space = 'free')+
  ggplot2::scale_fill_gradient2(low = colors[1],midpoint = 0.5,high = colors[3] ,mid = colors[2])+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 axis.title = ggplot2::element_blank(),
                 panel.spacing = ggplot2::unit(0,'line'),
                 strip.text.y = ggplot2::element_text(angle = 0))+
  ggplot2::scale_x_discrete( expand = c(0, 0)) +
  ggplot2::scale_y_continuous( expand = c(0, 0))

clusters=data%>%
  dplyr::group_by(clusters,chr,start,end)%>%
  dplyr::summarise(.groups = 'drop')

return(list(clusters = clusters, plot = p))
}


#'  Gene ontology on RT clusters using clusterprofiler
#'
#' @return ggplot
#'
#' @importFrom stringr str_locate_all str_split
#' @importFrom dplyr select mutate group_by n all_of summarise tibble arrange filter ungroup
#' @importFrom tidyr spread drop_na %>% gather
#' @importFrom foreach foreach %:% %do%
#' @importFrom clusterProfiler enrichGO
#' @importFrom GenomicRanges makeGRangesFromDataFrame split
#' @importFrom IRanges findOverlaps pintersect width
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom ggplot2 ggplot aes geom_col facet_grid theme_bw
#'
#' @param clusters, clusters from RT_clustering
#' @param genes, either a dataframe containing chr,start,end,gene column, where gene is an ENSEMBL identifier or the path to a gtf file
#' @param orgDB, OrgDb.
#' @param ontology, One of "BP", "MF", and "CC" sub-ontologies, or "ALL" for all three.
#' @param showCategory, number of categories to show per cluster
#' @param qvalueCutoff, cut-off for qval
#' @param ..., other parameters to pass to enrichGO
#'
#' @export
#'
Cluster_enrichment=function(clusters,genes,orgDB,ontology=c("BP", "MF", "CC", "ALL"),showCategory=10,qvalueCutoff=0.05,...){
  #operator
  `%>%`=tidyr::`%>%`
  #if a path to a gtf is provided load it and select info of interest
  if (class(genes) == 'character') {
    genes = readr::read_tsv(genes, col_names = F) %>%
      dplyr::mutate(gene = stringr::str_extract(string = X9, 'ENSG[0-9]{11}')) %>%
      dplyr::select('chr' = X1,
                    'start' = X4,
                    'end' = X5,
                    gene) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(chr=unique(chr),start = min(start), end = max(end))
  } else{
    #remove suffix from gene names
    genes$gene = sapply(strsplit(x = genes$gene, split = '\\.'), function(x)
      x[1])
  }

  convert_to_granges_if_needed=function(x) {
    if (class(x)[1] != 'GRanges') {
      y = tryCatch(
        expr = GenomicRanges::makeGRangesFromDataFrame(df = x, keep.extra.columns = T),
        error = function(e)
          stop("Provaded objects are not GRanges nor DataFrames")
      )
      return(y)
    } else{
      return(x)
    }
  }

  clusters=convert_to_granges_if_needed(clusters)
  genes=convert_to_granges_if_needed(genes)

  #look for overlap
  hits = IRanges::findOverlaps(query = genes, subject = clusters)

  overlaps = IRanges::pintersect(genes[S4Vectors::queryHits(hits)], clusters[S4Vectors::subjectHits(hits)])
  # check percentage of overlap of the first group region
  overlaps$hit <-
    IRanges::width(overlaps) / IRanges::width(genes[S4Vectors::queryHits(hits)]) >= 0.5

  overlaps$clusters = clusters$clusters[S4Vectors::subjectHits(hits)]

  overlaps=overlaps[overlaps$hit]

  #split into lists
  overlaps=GenomicRanges::split(overlaps,overlaps$clusters)

  overlaps = lapply(overlaps, function(x) {
    clusterProfiler::enrichGO(
      gene = x$gene,
      OrgDb= orgDB,
      keyType= 'ENSEMBL',
      ont= ontology[1],
      qvalueCutoff =qvalueCutoff,
      ...
    )
  })

  # recover inf
  overlaps=lapply(1:length(overlaps), function(x){
    results=overlaps[[x]]@result[overlaps[[x]]@result$p.adjust < overlaps[[x]]@qvalueCutoff,]
    if(nrow(results)>0){
    results$cluster=names(overlaps)[x]
    }else{
      results=dplyr::tibble()
    }
    return(results)
  })

  #merge dataframes
  overlaps=do.call(what = 'rbind',overlaps)

  #plot
  plot=overlaps%>%
    dplyr::arrange(-p.adjust)%>%
    dplyr::group_by(cluster)%>%
    dplyr::mutate(n=1:dplyr::n())%>%
    dplyr::filter(n<=showCategory)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(GeneRatio=sapply(GeneRatio,function(GeneRatio){
      x=stringr::str_split(GeneRatio,'/',simplify = T)
      return(as.numeric(x[1])/as.numeric(x[2]))
    }))%>%
    ggplot2::ggplot(ggplot2::aes(GeneRatio,Description,fill=p.adjust))+
    ggplot2::geom_col()+
    ggplot2::facet_grid(cluster~.,space = 'free',scales = 'free')+
    ggplot2::theme_bw()
  return(plot)
}


