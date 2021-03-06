% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_lda_gr.R
\name{gg_stat_lda_bootstrap}
\alias{gg_stat_lda_bootstrap}
\title{Plot stat_lda_bootstrap digest}
\usage{
gg_stat_lda_bootstrap(x)
}
\arguments{
\item{x}{object returned by \link{stat_lda_bootstrap}}
}
\description{
Plots as heatmap and with useful cosmetics
}
\examples{
dummy_df \%>\%
    stat_lda_bootstrap(foo2_NA, a:e, k=10) \%>\%
    gg_stat_lda_bootstrap()
}
