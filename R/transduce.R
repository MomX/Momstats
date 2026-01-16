# transduce -----

#' Transduce model positions to coefficients and shapes
#'
#' Predict coefficients and reconstruct shapes at specified positions in model space.
#'
#' @param object A stat_* object (stat_pca, stat_lda, stat_lm, ...)
#' @param positions A tibble with position coordinates. Column names must match
#'   the model's coordinate system (e.g., PC1, PC2 for PCA). Each row represents
#'   one position where coefficients will be reconstructed.
#'
#' @return A tibble with columns:
#'   * Position column(s): from input positions tibble
#'   * Coefficient column(s): coe (or coe_A, coe_B, ...) with proper classes
#'   * Inverse shape column(s): coe_i (or coe_A_i, coe_B_i, ...) as matrices
#'
#' @details
#' `transduce()` reconstructs coefficient vectors at specified positions in model
#' space and applies inverse transformations to obtain shapes.
#'
#' The `positions` tibble defines where to reconstruct shapes. For PCA, columns
#' should be named PC1, PC2, etc. Unspecified PCs default to 0 (mean position).
#'
#' Shapes are reconstructed using the appropriate inverse method based on
#' coefficient class (eft_i, rft_i, opoly_i, etc.).
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#'
#' # Single axis
#' pos <- tibble(PC1 = c(-2, 0, 2))
#' trans <- transduce(pca, pos)
#'
#' # Two axes (manual grid)
#' pos <- tibble(
#'   PC1 = c(-2, 0, 2, -2, 0, 2),
#'   PC2 = c(-1, -1, -1, 1, 1, 1)
#' )
#' trans <- transduce(pca, pos)
#'
#' # Or use expand_grid for convenience
#' pos <- expand_grid(PC1 = c(-2, 0, 2), PC2 = c(-1, 1))
#' trans <- transduce(pca, pos)
#'
#' # Access reconstructed shapes
#' plot(trans$coe_i[[1]], type = "l")
#' }
#'
#' @seealso [tidyr::expand_grid()] for creating position grids
#'
#' @export
transduce <- function(object, positions) {
  UseMethod("transduce")
}
