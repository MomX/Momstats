#' Plot confusion matrices
#'
#' Plots as heatmap and with useful cosmetics
#'
#' @param x confusion matrix or `CV_tbl` returned by lda
#' @param prop `logical` if `TRUE` (default), use accuracies (within an actual class); if `FALSE` use counts
#' @param percent `logical` whether to use percentages (default to `TRUE`)
#' @param drop_zeros `logical` whether to drop zeros in labels (default to `TRUE`)
#' @param text `logical` whether display cell information (default to `TRUE`)
#' @param text_colour `logical` whether display cell information (default to `white`)
#' @param text_size `numeric` size for cell information, in points (default to `8`)
#' @param text_signif `integer` to feed [signif] (default to `2``)
#' @param axis_size `numeric` size for axes labels
#' @param axis_x_angle `numeric` angle for x-axis labels
#'
#' @details when a `matrix` is passed, it assumes that it is:
#'
#' * squared
#' * rows are for actual classes
#' * columns for predicted classes
#' * class accuracies are calculated per actual classes, ie rowwise
#'
#' @examples
#' # works fine on matrices too
#' # set.seed(2329) # for the sake of replicability
#' m <- matrix(sample(1:100), 10)
#' gg_CV(m)
#'
#' @export
gg_CV <- function(x,
                  prop         = TRUE,
                  percent      = TRUE,
                  drop_zeros   = TRUE,
                  text         = TRUE,
                  text_colour  = "white",
                  text_size    = 8,
                  text_signif  = 2,
                  axis_size    = 10,
                  axis_x_angle = 0){
  UseMethod("gg_CV")
}


#' @rdname gg_CV
#' @export
gg_CV.default <- function(x,
                          prop         = TRUE,
                          percent      = TRUE,
                          drop_zeros   = TRUE,
                          text         = TRUE,
                          text_colour  = "white",
                          text_size    = 8,
                          text_signif  = 2,
                          axis_size    = 10,
                          axis_x_angle = 0){

  df <- x

  # early steps for gg
  fills <- range(df$prop)
  accuracy_legend <- "accuracy"

  # handles cell info
  # make a copy from appropriate column to label
  # percentage if required
  if (prop){
    if (percent){
      df <- df %>% dplyr::mutate(prop=100*prop)
      accuracy_legend <- paste0(accuracy_legend, " (%)")
    }
    df <- dplyr::mutate(df, label=signif(.data$prop, text_signif))
  } else {
    df <- dplyr::mutate(df, label=.data$n)
  }

  # empty str if drop_zeros
  if (drop_zeros)
    df <- dplyr::mutate(df, label=ifelse(.data$label==0, "", .data$label))

  # build the plot
  gg <- df %>%
    ggplot() +
    aes(x     = .data$predicted,
        y     = forcats::fct_rev(.data$actual), # for sensible plotting
        fill  = .data$prop) +
    # square tiles
    geom_tile() +
    coord_equal() +
    # scales is forced to [0; 1] and put on top
    scale_fill_viridis_c(begin = fills[1], end=fills[2]) +
    guides(fill=guide_legend(accuracy_legend, reverse = TRUE)) +
    # cosmetics
    labs(x="predicted as", y="actual class") +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = axis_size),
      axis.text.x = element_text(angle = axis_x_angle, hjust = 0)
    )

  # add values if required
  if (text)
    gg <- gg + geom_text(aes(label=.data$label),
                         colour=text_colour, size=text_size/.point)

  print(df)
  # return this beauty
  gg
}

#' @rdname gg_CV
#' @export
gg_CV.stat_lda_full <- function(x,
                                prop         = TRUE,
                                percent      = TRUE,
                                drop_zeros   = TRUE,
                                text         = TRUE,
                                text_colour  = "white",
                                text_size    = 8,
                                text_signif  = 2,
                                axis_size    = 10,
                                axis_x_angle = 0){
  gg_CV(x            = x$CV_tbl,
        percent      = percent,
        drop_zeros   = drop_zeros,
        text         = text,
        text_colour  = text_colour,
        text_size    = text_size,
        text_signif  = text_signif,
        axis_size    = axis_size,
        axis_x_angle = axis_x_angle)
}

#' @rdname gg_CV
#' @export
gg_CV.matrix <- function(x,
                         prop         = TRUE,
                         percent      = TRUE,
                         drop_zeros   = TRUE,
                         text         = TRUE,
                         text_colour  = "white",
                         text_size    = 8,
                         text_signif  = 2,
                         axis_size    = 10,
                         axis_x_angle = 0){
  gg_CV(x = .CV_tbl(x),
        percent      = percent,
        drop_zeros   = drop_zeros,
        text         = text,
        text_colour  = text_colour,
        text_size    = text_size,
        text_signif  = text_signif,
        axis_size    = axis_size,
        axis_x_angle = axis_x_angle)
}

#' Plot stat_lda_bootstrap digest
#'
#' Plots as heatmap and with useful cosmetics
#'
#' @param x object returned by [stat_lda_bootstrap]
#'
#' @examples
#' dummy_df %>%
#'     stat_lda_bootstrap(foo2_NA, a:e, k=10) %>%
#'     gg_stat_lda_bootstrap()
#' @export
gg_stat_lda_bootstrap <- function(x){

  # extract what we need
  df       <- x$CV_accuracy
  observed <- x$observed
  k        <- x$k

  # let's go
  df %>%
    ggplot2::ggplot() +
    ggplot2::aes(x=.data$what, y=.data$accuracy) +
    ggplot2::geom_hline(yintercept = observed, linetype="dashed", colour="grey50") +
    ggplot2::geom_violin(alpha=0.75) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_discrete(labels=c("balanced", "random", "random and balanced")) +
    ggplot2::annotate("text", x = 0, y = observed, hjust=0, vjust=-0.5,
                      label = paste0("observed=", signif(observed, 3)),
                      size=8/.point, colour="grey20") +
    ggplot2::labs(x="resampled datasets", y="accuracy", subtitle = paste0(k, " permutations"))
}

