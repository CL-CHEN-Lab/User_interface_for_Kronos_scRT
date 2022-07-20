#' calculates twidth
#'
#' @return tibble
#'
#' @importFrom tidyr %>% gather spread
#' @importFrom dplyr filter select group_by mutate summarise ungroup
#'
#' @export
#'
#' @param df, output of Twidth_fit_data
#'
#'

Twidth = function(df) {
  #load required functions
  `%>%` = tidyr::`%>%`

  t = df %>%
    dplyr::filter(t75 | t25) %>%
    tidyr::gather('t', 'value', t25, t75) %>%
    dplyr::filter(value) %>%
    dplyr::select(-percentage, -value) %>%
    dplyr::group_by(group, category, t) %>%
    dplyr::summarise(time = min(time)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(t, time) %>%
    dplyr::mutate(Twidth = abs(t75 - t25))

  return(t)
}

#' calculates twidth
#'
#' @return tibble
#' @importFrom  foreach foreach %dopar%
#' @importFrom  doSNOW registerDoSNOW
#' @importFrom snow makeCluster stopCluster
#'
#' @export
#'
#' @param df, dataframe containing the columns group and category
#' @param cores, number of cores for parallelization
#'
#'

Twidth_fit_data = function(df, ncores = 1) {
  df = Kronos.scRT::AverageVariability(df)

  #load required operators
  `%dopar%` = foreach::`%dopar%`
  `%:%` = foreach::`%:%`
  #declare clusters
  cl <- snow::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))

  fitted_data = foreach::foreach(
    cat = unique(df$category),
    .combine = 'rbind',
    .errorhandling = 'remove'
  ) %:%foreach::foreach(
    Group = unique(df$group),
    .combine = 'rbind',
    .errorhandling = 'remove'
  )%dopar% {
    Kronos.scRT::T25_75(df = df[df$category == cat & df$group == Group, ], Group, cat)
  }

  return(fitted_data)

}

#' fit RT data into sigmoid function and return T25 and T75
#'
#' @return tibble
#' @importFrom dplyr mutate select add_row
#' @importFrom tidyr %>%
#'
#' @export
#'
#' @param df, dataframe of group type
#' @param sample, the name of a sample
#' @param category, Genomic Category
#'
#'

T25_75 = function(df, sample, category) {
  #load required functions
  `%>%` = tidyr::`%>%`

  model = tryCatch(
    stats::nls(
      percentage ~ stats::SSlogis(time, Asym, xmid, scal),
      data = df[, c('percentage', 'time')] %>%
        dplyr::add_row(percentage = 1, time = -10) %>%
        dplyr::add_row(percentage = 0, time = 10),
      control = nls.control(maxiter = 100),
      algorithm = 'port',
      start = c(
        Asym = 1,
        xmid = 0,
        scal = -0.5
      )
    ),
    #If the data cannot be fitted with a Gauss-Newton algorithm, try the
    #Golub and Pereyra algorithm for the solution of a nonlinear least squares
    #problem which assumes a number of the parameters are linear.
    #Also, add a higher tolerance (1e-04 Vs 1e-05).
    error = tryCatch(
      function(e)
        stats::nls(
          percentage ~ stats::SSlogis(time, Asym, xmid, scal),
          data = df[, c('percentage', 'time')] %>%
            dplyr::add_row(percentage = 1, time = -10) %>%
            dplyr::add_row(percentage = 0, time = 10),
          algorithm = 'plinear',
          control = nls.control(
            maxiter = 100,
            tol = 1e-04,
            warnOnly = T
          )
        ),
      error = function(e)
        print('Try to reduce the number of RT groups')
    )
  )
  min = min(df$time)
  max = max(df$time)
  data = stats::predict(model,
                        newdata = data.frame(time = seq(min, max, 0.01)),
                        type = "l")
  result = data.frame(
    time = seq(min, max, 0.01),
    percentage = data,
    group = sample
  )

  t = result %>%
    dplyr::mutate(
      distance75 = abs(percentage - 0.75),
      distance25 = abs(percentage - 0.25)
    ) %>%
    dplyr::mutate(
      min75 = min(distance75),
      min25 = min(distance25),
      t75 = distance75 == min75,
      t25 = distance25 == min25
    ) %>%
    dplyr::select(group, time, percentage, t75, t25) %>%
    dplyr::mutate(category = category)

  return(t)
}

#' given an RT value (Early == 1 and Late == 0), this function returns its RT category base on the number of categories to return (2,3 or 5)
#'
#' @return a string or a vector of strings
#' @importFrom dplyr case_when
#'
#' @export
#'
#' @param RT, a number of a vector of numbers
#' @param number, number of categories to create
#' @param ..., other parameters that can be passed to Rtsne
#'
#'

split_into_categoreis = Vectorize(function(RT, number = 2) {
  if (number == 3) {
    return(
      dplyr::case_when(
        RT > 0.6666667 ~ 'Early',
        RT > 0.3333333 & RT <= 0.6666667 ~ 'Mid',
        RT <= 0.3333333 ~ 'Late'
      )
    )
  } else if (number == 5) {
    return(
      dplyr::case_when(
        RT > 0.8 ~ 'Very Early',
        RT <= 0.8 & RT > 0.6 ~ 'Early',
        RT <= 0.6 & RT > 0.4 ~ 'Mid',
        RT <= 0.4 & RT > 0.2 ~ 'Late',
        RT <= 0.2 ~ 'Very Late'
      )
    )
  } else {
    return(dplyr::case_when(RT > 0.5 ~ 'Early',
                            RT <= 0.5 ~ 'Late'))
  }
}, vectorize.args = 'RT')

#' Order in which RT categories should be plot
#'
#' @return character vector
#'
#' @export
#'
#' @param number, number of categories (2,3 or 5)
#'
#'

cat_levels = function(number) {
  if (number == 3) {
    return(c('Early',
             'Mid',
             'Late'))
  } else if (number == 5) {
    return(c('Very Early',
             'Early',
             'Mid',
             'Late',
             'Very Late'))
  } else if (number == 2) {
    return(c('Early',
             'Late'))
  }else if (number == 1) {
    return(c('All'))
  }else{
    stop('Allowed number of categories: 1,2,3 and 5.')
  }
}

#' Extended Twidth plot
#'
#' @return ggplot element
#' @importFrom ggplot2 aes facet_grid geom_line geom_point geom_text geom_vline ggplot scale_x_reverse scale_y_continuous theme
#' @export
#'
#' @param Variability, dataframe created by either TW_RTAnnotation or TW_GenomeAnnotation
#' @param Fitted_data, dataframe creted by Twidth_fit_data
#' @param Twidth, dataframe created by Twidth
#' @param Color, plot color (s)
#'
#'

Twidth_extended_plot = function(Variability, Fitted_data, Twidth, Color='red') {
  #calculate average variability
  Variability = Kronos.scRT::AverageVariability(Variability)

  #plot
  plot = ggplot2::ggplot(Variability) +
    ggplot2::geom_point(ggplot2::aes(time, percentage), color = Color) +
    ggplot2::geom_line(data = Fitted_data,
                       ggplot2::aes(time, percentage),
                       color = 'blue') +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::geom_vline(data = Twidth,
                        ggplot2::aes(xintercept = t25),
                        color = 'red') +
    ggplot2::geom_vline(data = Twidth,
                        ggplot2::aes(xintercept = t75),
                        color = 'red') +
    ggplot2::geom_text(
      data = Twidth,
      ggplot2::aes(label = paste('TW\n', Twidth)),
      x = Inf,
      y = 0.5,
      hjust = 1.25
    ) +
    ggplot2::facet_grid( ~ category)

  return(plot)
}


#' Twidth barplotplot
#'
#' @return ggplot element
#' @importFrom ggplot2 aes element_text geom_col geom_text ggplot theme xlab ylab
#' @importFrom dplyr n group_by summarise select inner_join
#' @importFrom tidyr %>%
#' @export
#'
#' @param Variability, dataframe created by either TW_RTAnnotation or TW_GenomeAnnotation
#' @param Twidth, dataframe created by Twidth
#' @param Color, barplot color(s)
#' @param pval, dataframe produced by twidth_pval
#'
#'

Twidth_barplot = function(Variability, Twidth,Color='lightblue',pval=NULL) {
  #load operator
  `%>%` = tidyr::`%>%`

  #calculate bins info
  Twidth = dplyr::inner_join(
    Twidth,
    Variability %>%
      dplyr::select(group, chr, start, end, category) %>%
      unique() %>%
      dplyr::group_by(group, category) %>%
      dplyr::summarise(`N of bins` = n()),
    by = c('group', 'category')
  )
  #plot
  p = ggplot2::ggplot(Twidth) +
    ggplot2::geom_col(ggplot2::aes(category, Twidth),
                      position = 'dodge',
                      fill = Color) +
    ggplot2::ylab('Twidth') + ggplot2::xlab('') +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = 'none') +
    ggplot2::geom_text(
      ggplot2::aes(category, Twidth / 2, label = paste0('n bins:\n', `N of bins`)),
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ) +
    ggplot2::geom_text(aes(category, Twidth, label = paste0('Twidth: ', Twidth)),
                       vjust = -0.25)

  if(!is.null(pval)){
    pval=pval%>%
      dplyr::mutate(a=as.numeric(category1),
                    b=as.numeric(category2),
                    position=(a+b)/2,
                    x=ifelse(a>b,b,a),
                    xend=ifelse(a<b,b,a),
                    pval_label = tryCatch(
                      format(
                        adj_pval,
                        digits = 2,
                        scientific = T
                      ),
                      error=function(x)
                        format(
                          pval,
                          digits = 2,
                          scientific = T
                        )
                    ),
                    Statistically=ifelse(
                      tryCatch(
                        adj_pval < 0.05 & iterations > 2000,
                        error=function(x)
                          pval < 0.05 & iterations > 2000
                      ),
                      'Significant',
                      'Non-significant')
      )%>%
      dplyr::group_by(group)%>%
      dplyr::mutate( y=seq(1.1*max(Twidth$Twidth),(1+dplyr::n()/10)*max(Twidth$Twidth),length.out = dplyr::n()))%>%
      dplyr::ungroup()

    p=p+
      ggplot2::geom_text(data=pval, ggplot2::aes(x=position,y=y,label=pval_label,color=Statistically),vjust=-0.25)+
      ggplot2::geom_segment(data=pval,ggplot2::aes(
        x=x,
        xend=xend,
        y=y,
        yend=y,color=Statistically))+
      ggplot2::scale_color_manual(values = c('Significant'='red','Non-significant'='black'))+
      ggplot2::theme(legend.position = 'top')+
      ggplot2::labs(color='')
    }

  return(p)
}

#' Calculates the proportion of replicated bins in function of time
#'
#' @return ggplot element
#'
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom tidyr %>%
#' @export
#'
#' @param Variability, dataframe created by either TW_RTAnnotation or TW_GenomeAnnotation
#'
#'

AverageVariability = function(Variability) {
  #load operator
  `%>%` = tidyr::`%>%`

  Variability %>%
    dplyr::group_by(group, time, category) %>%
    dplyr::summarise(percentage = mean(percentage)) %>%
    dplyr::ungroup() %>%
    return()
}


#' Annotate Variability dataframe with a customized annotation
#'
#' @importFrom dplyr group_by mutate summarise ungroup select as_tibble bind_rows inner_join
#' @importFrom tidyr %>%
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges findOverlaps pintersect width
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @export
#'
#' @param Variability, dataframe created by Kronos
#' @param GenomeAnnotation, a dataFrame containing chr, start, end and annotation columns.
#'


TW_GenomeAnnotation = function(Variability,
                               GenomeAnnotation) {
  #load operator
  `%>%` = tibble::`%>%`

  #convert into Grange
  GenomeAnnotation = GenomeAnnotation %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = T,
      seqnames.field = 'chr',
      end.field = 'end',
      start.field = 'start'
    )

  # create GRange Bins
  Bin = Variability %>%
    dplyr::select(chr, start, end) %>%
    unique() %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = T,
      seqnames.field = 'chr',
      end.field = 'end',
      start.field = 'start'
    )

  #find overlaps
  hits = IRanges::findOverlaps(Bin, GenomeAnnotation)

  #info about overlapping regions
  overlaps <-
    IRanges::pintersect(GenomeAnnotation[S4Vectors::subjectHits(hits)], Bin[S4Vectors::queryHits(hits)])

  #Add annotation
  add = dplyr::as_tibble(Bin[S4Vectors::queryHits(hits)])
  not_add = dplyr::as_tibble(Bin[-S4Vectors::queryHits(hits)])

  #Based on the overlap define the predominant notation of each bin.
  #A category to be chosen has to be predominant (at least 60% of the total overlaps in the bin)
  Annotation = dplyr::bind_rows(
    add %>%
      dplyr::mutate(
        size = IRanges::width(overlaps),
        category = overlaps$annotation
      ) %>%
      dplyr::group_by(seqnames, start, end, category) %>%
      dplyr::summarise(n = sum(size)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(seqnames, start, end) %>%
      dplyr::mutate(
        select = (n / sum(n) >= 0.6),
        category = ifelse(any(select), category, '_Unknown_')
      ) %>%
      dplyr::select(seqnames, start, end, category),
    not_add %>%
      dplyr::mutate(category = '_Unknown_') %>%
      dplyr::select(seqnames, start, end, category)
  )

  Annotation = Annotation %>%
    dplyr::ungroup() %>%
    dplyr::mutate(seqnames = as.character(seqnames))

  Variability = Variability %>%
    dplyr::inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))

  return(Variability)
}


#' Annotate Variability dataframe based on RT Categories
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr %>%
#' @export
#'
#' @param Variability, dataframe created by Kronos
#' @param RT_Groups, number of RT groups (1,2,3,5)
#'

TW_RTAnnotation = function(Variability, RT_Groups = 2) {
  #load operator
  `%>%` = tibble::`%>%`

  #assign bins to RT groups
  if (RT_Groups == 1) {
    Variability = Variability %>%
      dplyr::mutate(category = factor('All'))
  } else if (RT_Groups %in% c(2, 3, 5)) {
    Variability = Variability %>%
      dplyr::mutate(
        category = Kronos.scRT::split_into_categoreis(RT, number = RT_Groups),
        category = factor(category,
                          levels = Kronos.scRT::cat_levels(RT_Groups))
      )
  }
  return(Variability)
}

#' Hypothesis test between different RT categories inside a sample
#'
#'
#' @importFrom dplyr mutate tibble filter pull select inner_join
#' @importFrom tidyr %>%
#' @importFrom foreach %do% %:% %dopar% foreach
#' @importFrom  doSNOW registerDoSNOW
#' @importFrom snow makeCluster stopCluster
#' @export
#'
#' @param Variability, dataframe created by Kronos
#' @param twidth,  dataframe created by Twidth
#' @param pairs_to_test, a dataframe containing the columns Category1, Category2 to test
#' @param adjust.methods, correction method. Can be abbreviated. by default it is set as none
#' @param alternative, Hypothesis test. If pairs_to_test is null it is always two.sided
#' @param nIterations, number of iterations
#' @param ncores, number of cores for parallelization
#'


Twidth_pval = function(variability,
                       twidth,
                       pairs_to_test = NULL,
                       adjust.methods = c("none","holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                        "fdr"),
                       alternative = c('two.sided', 'greater', 'lower'),
                       nIterations = 10 ^ 4,
                       ncores = 3) {

  #define operators
  `%>%` = tidyr::`%>%`
  `%do%` = foreach::`%do%`
  `%dopar%` = foreach::`%dopar%`
  `%:%` = foreach::`%:%`
  #identify pairs to test
  if (is.null(pairs_to_test)) {
    #if pairs to test is not set alternative must be two.sided
    if (alternative != 'two.sided') {
      alternative = 'two.sided'
      warning('Since pairs_to_test is null, alternative has been set to two.sided')
    }
    # if not provided create all possible combinations
    pairs_to_test = sort(unique(twidth$category))

    pairs_to_test = foreach::foreach(i = 1:(length(pairs_to_test) - 1), .combine = 'rbind') %:%
      foreach::foreach(h = (i + 1):length(pairs_to_test),
                       .combine = 'rbind') %do%
      {
        dplyr::tibble(Category1 = pairs_to_test[i],
                      Category2 = pairs_to_test[h])

      }

  }

  # base of the input invert the order of element or calculate absolute
  stat_type = function(what, x, y) {
    return(switch (
      what,
      'greater' = x - y,
      'lower' = y - x,
      'two.sided' = abs(y - x),
      stop('wrong alternative hypotesys provided')
    ))
  }

  pval = foreach::foreach(g = unique(variability$group),
                          .combine = 'rbind') %:%
    foreach::foreach (
      i = 1:nrow(pairs_to_test),
      .combine = 'rbind' ,
      .errorhandling = "remove"
    ) %do% {
      #calculate real difference between two categories inside a group
      Real_difference = stat_type(
        alternative,
        twidth %>%
          dplyr::filter(category == pairs_to_test$Category1[i],
                        group == g) %>%
          dplyr::pull(Twidth),
        twidth %>%
          dplyr::filter(category == pairs_to_test$Category2[i],
                        group == g) %>%
          dplyr::pull(Twidth)
      )

      Variabilit_to_test = variability %>%
        dplyr::filter(
          category == pairs_to_test$Category1[i] |
            category == pairs_to_test$Category2[i],
          group == g
        )

      Bins = Variabilit_to_test %>%
        dplyr::select(chr, start, end, category) %>%
        unique()

      Variabilit_to_test$category = NULL

      #n of bins in each category
      nCat1= Bins%>%dplyr::filter(category==pairs_to_test$Category1[i])%>%nrow()
      nCat2= Bins%>%dplyr::filter(category==pairs_to_test$Category2[i])%>%nrow()

      #start cluster and permutate
      cl <- snow::makeCluster(ncores)
      doSNOW::registerDoSNOW(cl)
      on.exit(snow::stopCluster(cl))

      Boot_Strapped = foreach::foreach(
        index = 1:(nIterations-1),
        .combine = '+',
        .errorhandling = "remove"
      ) %dopar% {
        #define operators
        `%>%`=tidyr::`%>%`

        #shuffle categories
        Bins_perm = rbind(Bins[sample(1:nrow(Bins),size = nCat1,replace = T),]%>%
            dplyr::mutate(category=pairs_to_test$Category1[i]),
            Bins[sample(1:nrow(Bins),size = nCat2,replace = T),]%>%
              dplyr::mutate(category=pairs_to_test$Category2[i])
          )

        Variabilit_to_test_perm = Variabilit_to_test %>%
          dplyr::inner_join(Bins_perm, by = c("chr", "start", "end"))


        #fit data
        twidth_fitted_data_perm = Kronos.scRT::Twidth_fit_data(df = Variabilit_to_test_perm,
                                                          ncores = 1)

        #calculate TW
        twidth_perm = Kronos.scRT::Twidth(twidth_fitted_data_perm)

        c(
          stat_type(
            alternative,
            twidth_perm %>%
              dplyr::filter(category == pairs_to_test$Category1[i],
                            group == g) %>%
              dplyr::pull(Twidth),
            twidth_perm %>%
              dplyr::filter(category == pairs_to_test$Category2[i],
                            group == g) %>%
              dplyr::pull(Twidth)
          ) >=
            Real_difference,
          1
        )


      }

      #return pvalues
      dplyr::tibble(
        group = g,
        category1 = pairs_to_test$Category1[i],
        category2 = pairs_to_test$Category2[i],
        #calculate pvalue over effective iterations
        pval = (Boot_Strapped[1] + 1) / (Boot_Strapped[2] + 1),
        # actual iterations
        iterations = Boot_Strapped[2] + 1
      )
    }

  if(adjust.methods=='none'){
      return(pval)
  }else{

    pval%>%
      dplyr::mutate(adj.pval=p.adjust(pval,method = adjust.methods))%>%
      dplyr::select(group,category1,category2,pval,adj.pval,iterations)%>%
      return()
    }

}
