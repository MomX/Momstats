# stat_pca.R

#' Principal Component Analysis for morphometric data
#'
#' Perform PCA on coefficient data from morphometric analyses.
#'
#' @param data A tibble with coefficient columns
#' @param formula A formula specifying predictors. Can be:
#'   * Missing: auto-detects single coe column
#'   * Bare column name: `coe`
#'   * Formula: `~ coe`, `~ coe + size`, `~ coe1 + coe2`
#'   * Use `coe` in formula to auto-detect coefficient columns
#' @param center Logical. Should data be centered? Default `TRUE`.
#' @param scale Logical or NULL. Should data be scaled to unit variance?
#'   If `NULL` (default), automatically determined based on predictor types.
#' @param ... Additional arguments passed to [stats::prcomp()]
#'
#' @return An object of class `c("stat_pca", "momstats")` containing:
#' * `data`: Original tibble (unchanged)
#' * `model`: The [stats::prcomp()] object
#' * `method`: "pca"
#' * `call`: The function call
#' * `formula`: Formula used (if any)
#' * `coe_cols`: Names of coefficient columns used
#' * `covariate_cols`: Names of covariate columns used (if any)
#' * `predictor_cols`: All predictor column names
#' * `group_col`: NULL (no grouping for PCA)
#' * `variance_explained`: Proportion of variance explained by each PC
#' * `cumvar_explained`: Cumulative proportion of variance explained
#' * `center`: Logical, was centering applied
#' * `scale`: Logical, was scaling applied
#'
#' @details
#' `stat_pca()` provides a unified interface for PCA on morphometric coefficient
#' data. It handles:
#'
#' * **Automatic coefficient detection**: Finds `coe` columns via class
#' * **Coefficient unfolding**: Expands list-columns to individual harmonic columns
#' * **Multiple coefficient types**: Can combine different coefficient columns
#' * **Covariates**: Add non-morphometric predictors (size, environmental variables)
#' * **Smart defaults**: Automatically determines centering and scaling
#'
#' ## Formula syntax
#'
#' The formula specifies which predictors to use:
#' * `~ coe`: Use auto-detected coefficient column(s)
#' * `~ coe1 + coe2`: Use specific coefficient columns
#' * `~ coe + size`: Coefficient column plus a covariate
#' * Bare name or missing: auto-detect single coe column
#'
#' ## Centering and scaling
#'
#' * `center = TRUE` (default): Variables are centered to mean = 0
#' * `scale = NULL` (default): Automatic determination:
#'   - Pure EFT coefficients: `scale = FALSE`
#'   - Mixed predictors or non-EFT: `scale = TRUE`
#' * `scale = TRUE/FALSE`: Override automatic behavior
#'
#' ## Getting results
#'
#' Use [collect()] to add PC scores to your data:
#' ```r
#' pca <- boteft %>% stat_pca()
#' boteft_with_pcs <- collect(pca, retain = 5)  # Add first 5 PCs
#' ```
#'
#' @examples
#' \dontrun{
#' # Basic PCA on coefficients
#' pca1 <- boteft %>% stat_pca()
#'
#' # Explicit column
#' pca2 <- boteft %>% stat_pca(coe)
#'
#' # Formula syntax
#' pca3 <- boteft %>% stat_pca(~ coe)
#'
#' # With covariate
#' pca4 <- boteft %>% stat_pca(~ coe + length)
#'
#' # Force scaling
#' pca5 <- boteft %>% stat_pca(~ coe, scale = TRUE)
#'
#' # Add PC scores to data
#' boteft_pca <- collect(pca1, retain = 5)
#'
#' # Or fold into single column
#' boteft_pca <- collect(pca1, retain = 5, fold = "pca_scores")
#'
#' # Plot results
#' plot(pca1)  # Score plot (default)
#' plot(pca1, type = "scree")
#' plot(pca1, type = "loadings")
#' plot(pca1, color = type)  # Color by grouping variable
#' plot(pca1, labels = type)  # Text labels
#' plot(pca1, color = type, chull = TRUE)  # With convex hulls
#' }
#'
#' @seealso [stats::prcomp()], [collect.stat_pca()], [plot.stat_pca()]
#'
#' @export
stat_pca <- function(data, formula = NULL, center = TRUE, scale = NULL, ...) {

  # Store call
  call <- match.call()

  # Check input
  if (!is.data.frame(data)) {
    stop("data must be a tibble or data frame")
  }

  # Parse formula to get predictor columns (captures messages for later)
  formula_parsed <- parse_unsupervised_formula(data, rlang::enquo(formula))

  coe_cols <- formula_parsed$coe_cols
  covariate_cols <- formula_parsed$covariate_cols
  predictor_cols <- formula_parsed$predictor_cols

  # Determine scale default if NULL
  if (is.null(scale)) {
    scale <- determine_scale_default(data, coe_cols, covariate_cols)
  }

  # Build predictor matrix
  predictor_matrix <- build_predictor_matrix(data, coe_cols, covariate_cols)

  # Run PCA
  pca_model <- stats::prcomp(predictor_matrix, center = center, scale. = scale, ...)

  # Calculate variance explained
  variance_explained <- pca_model$sdev^2 / sum(pca_model$sdev^2)
  cumvar_explained <- cumsum(variance_explained)

  # Build result object (minimal structure)
  result <- list(
    data = data,  # Original, not augmented
    model = pca_model,
    method = "pca",
    call = call,
    formula = formula,
    coe_cols = coe_cols,
    covariate_cols = covariate_cols,
    predictor_cols = predictor_cols,
    group_col = NULL,
    variance_explained = variance_explained,
    cumvar_explained = cumvar_explained,
    center = center,
    scale = scale
  )

  class(result) <- c("stat_pca", "momstats")
  result
}


# Smart defaults ----

#' Determine default for scale parameter
#'
#' @param data A tibble
#' @param coe_cols Coefficient column names
#' @param covariate_cols Covariate column names
#'
#' @return Logical
#'
#' @keywords internal
#' @noRd
determine_scale_default <- function(data, coe_cols, covariate_cols) {

  # If covariates present, scale
  if (length(covariate_cols) > 0) {
    return(TRUE)
  }

  # If multiple coe columns, scale
  if (length(coe_cols) > 1) {
    return(TRUE)
  }

  # If single coe column, check if it's EFT
  if (length(coe_cols) == 1) {
    col <- data[[coe_cols[1]]]
    # Check if eft class
    if ("eft" %in% class(col)) {
      return(FALSE)  # Don't scale pure EFT
    } else {
      return(TRUE)   # Scale other coefficient types
    }
  }

  # Default to TRUE
  TRUE
}


# Collect method ----

#' Collect PC scores from PCA results
#'
#' Extract PC scores from a PCA result and add them to a tibble.
#'
#' @param x A `stat_pca` object
#' @param data A tibble. If NULL, uses the original data from the PCA.
#' @param retain How many PCs to retain:
#'   * `NULL` (default): All PCs
#'   * Integer (e.g., `5`): First N PCs
#'   * Numeric 0-1 (e.g., `0.95`): PCs explaining this proportion of variance
#' @param fold How to add PC scores:
#'   * `FALSE` (default): Add as separate columns (`PC1`, `PC2`, ...)
#'   * `TRUE`: Fold into single list-column named `"pca"`
#'   * Character: Fold into single list-column with this name
#' @param ... Additional arguments (reserved)
#'
#' @return A tibble with PC scores added
#'
#' @details
#' The `retain` parameter allows flexible selection of PCs:
#' * `retain = NULL`: Keep all PCs
#' * `retain = 5`: Keep first 5 PCs
#' * `retain = 0.95`: Keep PCs explaining 95% of variance
#' * `retain = 1`: Keep only PC1
#'
#' When `fold = TRUE` or a character name, PC scores are stored as a list-column
#' with class `c("pca", "coe")`, making them usable in downstream analyses.
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#'
#' # Add all PCs
#' collect(pca)
#'
#' # Add first 5 PCs
#' collect(pca, retain = 5)
#'
#' # Add PCs explaining 95% variance
#' collect(pca, retain = 0.95)
#'
#' # Fold into list-column
#' collect(pca, retain = 5, fold = TRUE)
#' collect(pca, retain = 5, fold = "pca_scores")
#' }
#'
#' @export
collect.stat_pca <- function(x, data = NULL, retain = NULL, fold = FALSE, ...) {

  # Use original data if not provided
  if (is.null(data)) {
    data <- x$data
  }

  # Get PC scores
  pc_scores <- x$model$x

  # Determine which PCs to retain
  if (is.null(retain)) {
    # Keep all
    pcs_to_keep <- seq_len(ncol(pc_scores))
  } else if (retain >= 1) {
    # Integer: first N PCs
    pcs_to_keep <- seq_len(min(retain, ncol(pc_scores)))
  } else if (retain > 0 && retain < 1) {
    # Proportion: PCs explaining this much variance
    pcs_to_keep <- seq_len(which(x$cumvar_explained >= retain)[1])
  } else {
    stop("retain must be NULL, an integer >= 1, or a proportion between 0 and 1")
  }

  # Extract selected PCs
  selected_scores <- pc_scores[, pcs_to_keep, drop = FALSE]

  # Add to data
  if (isFALSE(fold)) {
    # Unfold: add as separate columns
    for (i in seq_len(ncol(selected_scores))) {
      col_name <- paste0("PC", pcs_to_keep[i])
      data[[col_name]] <- selected_scores[, i]
    }
  } else {
    # Fold: add as list-column
    if (isTRUE(fold)) {
      col_name <- "pca"
    } else {
      col_name <- as.character(fold)
    }

    # Check if column exists
    if (col_name %in% names(data)) {
      stop(sprintf("Column '%s' already exists", col_name))
    }

    # Create list-column
    pc_list <- lapply(seq_len(nrow(selected_scores)), function(i) {
      vec <- selected_scores[i, ]
      names(vec) <- paste0("PC", pcs_to_keep)
      class(vec) <- c("pca", "coe", "numeric")
      vec
    })

    class(pc_list) <- c("pca", "coe", "list")
    data[[col_name]] <- pc_list
  }

  data
}


# Print method ----

#' @export
print.stat_pca <- function(x, ...) {
  cat("Principal Component Analysis\n")
  cat("---------------------------------------\n")

  # Build predictor description
  pred_parts <- character()
  if (length(x$coe_cols) > 0) {
    pred_parts <- c(pred_parts, paste(x$coe_cols, collapse = " + "))
  }
  if (length(x$covariate_cols) > 0) {
    pred_parts <- c(pred_parts, paste(x$covariate_cols, collapse = " + "))
  }
  pred_desc <- paste(pred_parts, collapse = " + ")

  cat(sprintf("* %d observations, %d predictors (%s)\n",
              nrow(x$data), ncol(x$model$x), pred_desc))

  # Centering/scaling info
  center_str <- if (x$center) "TRUE" else "FALSE"
  scale_str <- if (x$scale) "TRUE" else "FALSE"

  # Build reason message
  if (!x$scale) {
    reason <- "(pure EFT coefficients)"
  } else if (length(x$covariate_cols) > 0) {
    reason <- "(mixed predictors)"
  } else if (length(x$coe_cols) > 1) {
    reason <- "(multiple coefficient types)"
  } else {
    reason <- "(non-EFT coefficients)"
  }

  cat(sprintf("* Centering: %s, Scaling: %s %s\n", center_str, scale_str, reason))

  cat(sprintf("* %d principal components\n\n", length(x$model$sdev)))

  cat("Variance explained (first 5 PCs):\n")
  n_show <- min(5, length(x$variance_explained))
  var_df <- data.frame(
    PC = paste0("PC", seq_len(n_show)),
    Variance = round(x$variance_explained[seq_len(n_show)] * 100, 2),
    Cumulative = round(x$cumvar_explained[seq_len(n_show)] * 100, 2)
  )
  colnames(var_df) <- c("", "Var (%)", "Cum (%)")
  print(var_df, row.names = FALSE)

  cat("\n* calculated using prcomp ($model)\n")
  cat("* collect() to add PC scores to data\n")
  cat("* plot() to visualize results (?plot.stat_pca)\n")

  invisible(x)
}


# Plot method ----

#' Plot PCA results
#'
#' Visualize PCA results with score plots, scree plots, or loading plots.
#'
#' @param x A `stat_pca` object
#' @param type Character. Type of plot:
#'   * `"scores"`: PC score plot (default)
#'   * `"scree"`: Variance explained by each PC
#'   * `"loadings"`: Loading plot for first two PCs
#' @param pcs Integer vector of length 2. Which PCs to plot for score/loading plots.
#'   Default is `c(1, 2)`.
#' @param color Column name (bare or quoted) for coloring points in score plot.
#' @param labels Column name (bare or quoted) for text labels in score plot.
#'   For loadings, logical: should variable names be shown? Default `TRUE`.
#' @param chull Logical. Draw convex hulls around groups? Only works when `color`
#'   is a factor. Default `TRUE`.
#' @param n_pcs Integer. For scree plot, how many PCs to show? Default is all.
#' @param cex Numeric. Character expansion for labels. Default is 0.7.
#' @param legend Logical. Show legend for colors? Default `TRUE`.
#' @param ... Additional arguments (reserved for future use)
#'
#' @return NULL (invisibly). Draws plot as side effect.
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#'
#' # Score plot (default)
#' plot(pca)
#' plot(pca, color = type)
#' plot(pca, labels = type, color = type)
#' plot(pca, color = type, chull = TRUE)
#'
#' # Different PCs
#' plot(pca, pcs = c(2, 3))
#'
#' # Scree plot
#' plot(pca, type = "scree")
#' plot(pca, type = "scree", n_pcs = 10)  # First 10 only
#'
#' # Loading plot
#' plot(pca, type = "loadings")
#' plot(pca, type = "loadings", labels = FALSE)
#' }
#'
#' @export
plot.stat_pca <- function(x, type = c("scores", "scree", "loadings"),
                          pcs = c(1, 2), color = NULL, labels = NULL,
                          chull = TRUE, n_pcs = NULL, cex = 0.7,
                          legend = TRUE, ...) {

  type <- match.arg(type)

  if (type == "scree") {
    plot_pca_scree(x, n_pcs)
  } else if (type == "scores") {
    plot_pca_scores(x, pcs, rlang::enquo(color), rlang::enquo(labels),
                    chull, cex, legend)
  } else if (type == "loadings") {
    # Default labels to TRUE for loadings
    if (is.null(labels)) labels <- TRUE
    plot_pca_loadings(x, pcs, labels, cex)
  }

  invisible(NULL)
}


#' @keywords internal
#' @noRd
plot_pca_scree <- function(x, n_pcs = NULL) {

  var_pct <- x$variance_explained * 100
  cumvar_pct <- x$cumvar_explained * 100

  # Determine how many to show
  if (is.null(n_pcs)) {
    n_show <- length(var_pct)
  } else {
    n_show <- min(n_pcs, length(var_pct))
  }

  var_pct <- var_pct[seq_len(n_show)]
  cumvar_pct <- cumvar_pct[seq_len(n_show)]

  # Main plot
  plot(seq_along(var_pct), var_pct,
       type = "b", pch = 19, col = rgb(0, 0, 0, 0.5),
       xlab = "Principal Component",
       ylab = "Variance Explained (%)",
       main = "Scree Plot",
       las = 1,
       ylim = c(0, max(var_pct) * 1.15))

  abline(h = 0, col = "gray", lty = 2)

  # Add cumulative variance text ABOVE each point
  text(seq_along(var_pct), var_pct,
       labels = sprintf("%.1f%%", cumvar_pct),
       pos = 3,  # above
       cex = 0.7,
       col = "gray40")
}


#' @keywords internal
#' @noRd
plot_pca_scores <- function(x, pcs, color_quo, labels_quo, chull, cex, legend) {

  pc_scores <- x$model$x
  pc1 <- pc_scores[, pcs[1]]
  pc2 <- pc_scores[, pcs[2]]

  var1 <- round(x$variance_explained[pcs[1]] * 100, 1)
  var2 <- round(x$variance_explained[pcs[2]] * 100, 1)

  # Get color values
  color_vals <- get_plot_variable(x$data, color_quo)

  # Get label values
  label_vals <- get_plot_variable(x$data, labels_quo)

  # Determine colors
  if (!is.null(color_vals)) {
    if (is.factor(color_vals) || is.character(color_vals)) {
      # Categorical
      color_vals <- as.factor(color_vals)
      cols <- rainbow(nlevels(color_vals), alpha = 0.5)[as.integer(color_vals)]
    } else if (is.numeric(color_vals)) {
      # Continuous - use color gradient
      col_ramp <- colorRampPalette(c("blue", "red"))
      color_indices <- cut(color_vals, breaks = 100, labels = FALSE)
      cols <- col_ramp(100)[color_indices]
      cols <- paste0(cols, "80")  # Add transparency
    } else {
      cols <- rgb(0, 0, 0, 0.5)
    }
  } else {
    cols <- rgb(0, 0, 0, 0.5)
  }

  # Create plot
  plot(pc1, pc2,
       xlab = sprintf("PC%d (%s%%)", pcs[1], var1),
       ylab = sprintf("PC%d (%s%%)", pcs[2], var2),
       main = "PCA Score Plot",
       type = "n",  # Don't plot points yet
       asp = 1,
       las = 1)

  abline(h = 0, v = 0, col = "gray", lty = 2)

  # Draw convex hulls if requested (behind points)
  if (chull && !is.null(color_vals) &&
      (is.factor(color_vals) || is.character(color_vals))) {
    draw_chulls(pc1, pc2, color_vals, rainbow(nlevels(as.factor(color_vals)), alpha = 1))
  }

  # Now add points
  points(pc1, pc2, pch = 19, col = cols)

  # Add labels if requested
  if (!is.null(label_vals)) {
    text(pc1, pc2, labels = label_vals, cex = cex, pos = 3)
  }

  # Add legend if requested and we have colors
  if (legend && !is.null(color_vals)) {
    add_plot_legend(color_vals, cols, position = "topright")
  }
}


#' @keywords internal
#' @noRd
plot_pca_loadings <- function(x, pcs, labels, cex) {

  loadings <- x$model$rotation[, pcs]

  # Create empty plot (no points)
  plot(loadings[, 1], loadings[, 2],
       xlab = sprintf("PC%d loading", pcs[1]),
       ylab = sprintf("PC%d loading", pcs[2]),
       main = "PCA Loading Plot",
       type = "n",  # No points
       asp = 1,
       las = 1)

  abline(h = 0, v = 0, col = "gray", lty = 2)

  # Add arrows from origin (no points at tips)
  arrows(0, 0, loadings[, 1], loadings[, 2],
         length = 0.1, col = rgb(0, 0, 0, 0.5), lwd = 1.5)

  # Add labels if requested (default TRUE)
  if (isTRUE(labels)) {
    text(loadings[, 1], loadings[, 2],
         labels = rownames(loadings),
         cex = cex, pos = 3)
  }
}


# Plot helpers (shared across methods) ----

#' Add legend to plot
#'
#' @param color_vals Values used for coloring
#' @param cols Colors used
#' @param position Legend position
#'
#' @keywords internal
#' @noRd
add_plot_legend <- function(color_vals, cols, position = "topright") {

  if (is.null(color_vals)) return(invisible())

  if (is.factor(color_vals) || is.character(color_vals)) {
    # Categorical legend
    color_vals <- as.factor(color_vals)
    lvls <- levels(color_vals)
    # Get unique colors for each level
    unique_cols <- rainbow(nlevels(color_vals), alpha = 0.5)

    legend(position,
           legend = lvls,
           pch = 19,
           col = unique_cols,
           bty = "n",
           cex = 0.8)

  } else if (is.numeric(color_vals)) {
    # Continuous legend - show gradient
    range_vals <- range(color_vals, na.rm = TRUE)

    # Create simple gradient legend
    legend_text <- c(
      sprintf("%.2f", range_vals[2]),  # max
      "",
      sprintf("%.2f", range_vals[1])   # min
    )

    col_ramp <- colorRampPalette(c("blue", "red"))
    legend_cols <- col_ramp(3)

    legend(position,
           legend = legend_text,
           pch = 15,
           col = legend_cols,
           bty = "n",
           cex = 0.8,
           title = "Value")
  }
}


#' Draw convex hulls for groups
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param groups Factor or character vector of group assignments
#' @param cols Colors for each group
#'
#' @keywords internal
#' @noRd
draw_chulls <- function(x, y, groups, cols = NULL) {

  if (is.null(groups)) return(invisible())

  groups <- as.factor(groups)
  lvls <- levels(groups)

  if (is.null(cols)) {
    cols <- rainbow(length(lvls), alpha = 1)
  }

  for (i in seq_along(lvls)) {
    grp <- lvls[i]
    idx <- which(groups == grp)

    if (length(idx) < 3) next  # Need at least 3 points for hull

    x_grp <- x[idx]
    y_grp <- y[idx]

    # Get convex hull
    hull_idx <- grDevices::chull(x_grp, y_grp)
    hull_idx <- c(hull_idx, hull_idx[1])  # Close the polygon

    # Draw outline only (no fill)
    lines(x_grp[hull_idx], y_grp[hull_idx],
          col = cols[i],
          lwd = 2)
  }
}


# Summary method ----

#' @export
summary.stat_pca <- function(object, ...) {
  cat("Principal Component Analysis\n")
  cat(sprintf("\n%d observations, %d variables, %d PCs\n",
              nrow(object$data), ncol(object$model$rotation),
              length(object$model$sdev)))

  cat("\nImportance of components:\n")
  importance <- rbind(
    "Standard deviation" = object$model$sdev,
    "Proportion of Variance" = object$variance_explained,
    "Cumulative Proportion" = object$cumvar_explained
  )
  colnames(importance) <- paste0("PC", seq_len(ncol(importance)))

  # Show first 10 PCs
  n_show <- min(10, ncol(importance))
  print(round(importance[, seq_len(n_show)], 4))

  if (ncol(importance) > n_show) {
    cat(sprintf("\n... and %d more PCs\n", ncol(importance) - n_show))
  }

  invisible(object)
}
