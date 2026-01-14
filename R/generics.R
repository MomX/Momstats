# Generics ----

#' Collect results from statistical analyses
#'
#' Generic function to extract and add results from statistical analyses
#' to tibbles. Each method adds appropriate columns based on the analysis type.
#'
#' @param x An object from a statistical analysis (e.g., `stat_pca`, `stat_lda`)
#' @param ... Additional arguments passed to methods
#'
#' @return A tibble with results added
#'
#' @details
#' Different statistical methods add different types of results:
#' * `stat_pca`: PC scores
#' * `stat_lda`: Discriminant scores, predictions, posterior probabilities
#' * `stat_kmeans`: Cluster assignments
#'
#' @seealso [collect.stat_pca()]
#'
#' @export
collect <- function(x, ...) {
  UseMethod("collect")
}

