# stat_kmeans ------

#' K-means clustering for morphometric data
#'
#' Perform k-means clustering on coefficient data.
#'
#' @param data A tibble with coefficient columns
#' @param formula A formula specifying predictors. Can be:
#'   * Missing: auto-detects single coe column
#'   * Bare column name: `coe`
#'   * Formula: `~ coe`, `~ coe + size`, `~ coe1 + coe2`
#'   * Use `coe` in formula to auto-detect coefficient columns
#' @param k Integer. Number of clusters. Default is 3.
#' @param nstart Integer. Number of random starts for k-means. Default is 25.
#' @param ... Additional arguments passed to [stats::kmeans()]
#'
#' @return An object of class `c("stat_kmeans", "momstats")` containing:
#' * `data`: Original tibble (unchanged)
#' * `model`: The [stats::kmeans()] object
#' * `method`: "kmeans"
#' * `call`: The function call
#' * `formula`: Formula used (if any)
#' * `predictor_cols`: All predictor column names
#' * `k`: Number of clusters
#' * `centers`: Cluster centers in predictor space (matrix)
#' * `cluster_sizes`: Size of each cluster
#' * `withinss`: Within-cluster sum of squares
#' * `tot_withinss`: Total within-cluster sum of squares
#' * `betweenss`: Between-cluster sum of squares
#' * `variance_explained`: Proportion of variance explained by clustering
#'
#' @details
#' `stat_kmeans()` provides a unified interface for k-means clustering on
#' morphometric coefficient data.
#'
#' ## Formula syntax
#'
#' The formula specifies which predictors to use:
#' * `~ coe`: Use auto-detected coefficient column(s)
#' * `~ coe1 + coe2`: Use specific coefficient columns
#' * `~ coe + size`: Coefficient column plus a covariate
#' * Bare name or missing: auto-detect single coe column
#'
#' ## Getting results
#'
#' Use [collect()] to add cluster assignments to your data:
#' ```r
#' km <- boteft %>% stat_kmeans(k = 3)
#' boteft_clustered <- collect(km)  # Adds 'cluster' column
#' ```
#'
#' Use [transduce()] to get shapes at cluster centers:
#' ```r
#' center_shapes <- transduce(km, tibble(cluster = 1:3))
#' ```
#'
#' @examples
#' \dontrun{
#' # Basic k-means with 3 clusters
#' km1 <- boteft %>% stat_kmeans(k = 3)
#'
#' # More clusters
#' km2 <- boteft %>% stat_kmeans(k = 5)
#'
#' # With covariate
#' km3 <- boteft %>% stat_kmeans(~ coe + length, k = 4)
#'
#' # Add cluster assignments
#' boteft_clustered <- collect(km1)
#'
#' # Get cluster center shapes
#' centers <- transduce(km1, tibble(cluster = 1:3))
#'
#' # Plot results
#' plot(km1)  # Cluster visualization
#' plot(km1, color = type)  # Color by original grouping
#' }
#'
#' @seealso [stats::kmeans()], [collect.stat_kmeans()], [plot.stat_kmeans()], [transduce.stat_kmeans()]
#'
#' @export
stat_kmeans <- function(data, formula = NULL, k = 3, nstart = 25, ...) {

  # Store call
  call <- match.call()

  # Check input
  if (!is.data.frame(data)) {
    stop("data must be a tibble or data frame")
  }

  if (k < 2) {
    stop("k must be at least 2")
  }

  # Parse formula to get predictor columns
  formula_parsed <- parse_unsupervised_formula(data, rlang::enquo(formula))

  coe_cols <- formula_parsed$coe_cols
  covariate_cols <- formula_parsed$covariate_cols
  predictor_cols <- formula_parsed$predictor_cols

  # Build predictor matrix
  predictor_matrix <- build_predictor_matrix(data, coe_cols, covariate_cols)

  # Run k-means
  km_model <- stats::kmeans(predictor_matrix, centers = k, nstart = nstart, ...)

  # Calculate variance explained
  total_ss <- km_model$totss
  between_ss <- km_model$betweenss
  within_ss <- km_model$tot.withinss
  variance_explained <- between_ss / total_ss

  # Build result object
  result <- list(
    data = data,
    model = km_model,
    method = "kmeans",
    call = call,
    formula = formula,
    predictor_cols = predictor_cols,
    k = k,
    centers = km_model$centers,
    cluster_sizes = km_model$size,
    withinss = km_model$withinss,
    tot_withinss = within_ss,
    betweenss = between_ss,
    variance_explained = variance_explained
  )

  class(result) <- c("stat_kmeans", "momstats")
  result
}


# Print method ----

#' @export
print.stat_kmeans <- function(x, ...) {
  cat("K-means Clustering\n")
  cat("---------------------------------------\n")

  # Identify coe columns from predictor_cols
  coe_cols <- x$predictor_cols[
    sapply(x$data[x$predictor_cols], function(x) "coe" %in% class(x))
  ]
  other_cols <- setdiff(x$predictor_cols, coe_cols)

  # Build predictor description
  pred_parts <- character()
  if (length(coe_cols) > 0) {
    pred_parts <- c(pred_parts, paste(coe_cols, collapse = " + "))
  }
  if (length(other_cols) > 0) {
    pred_parts <- c(pred_parts, paste(other_cols, collapse = " + "))
  }
  pred_desc <- paste(pred_parts, collapse = " + ")

  cat(sprintf("* %d observations, %d predictors (%s)\n",
              nrow(x$data), ncol(x$centers), pred_desc))

  cat(sprintf("* %d clusters\n", x$k))
  cat(sprintf("* Variance explained: %.2f%%\n\n", x$variance_explained * 100))

  cat("Cluster sizes:\n")
  cluster_df <- data.frame(
    Cluster = seq_len(x$k),
    Size = x$cluster_sizes,
    Within_SS = round(x$withinss, 2)
  )
  print(cluster_df, row.names = FALSE)

  cat(sprintf("\n* Total within-cluster SS: %.2f\n", x$tot_withinss))
  cat(sprintf("* Between-cluster SS: %.2f\n", x$betweenss))

  cat("\n* collect() to add cluster assignments to data\n")
  cat("* transduce() to get shapes at cluster centers\n")
  cat("* plot() to visualize clusters\n")

  invisible(x)
}


# Summary method ----

#' @export
summary.stat_kmeans <- function(object, ...) {
  cat("K-means Clustering Summary\n")
  cat("==================================================\n\n")

  cat(sprintf("Number of clusters: %d\n", object$k))
  cat(sprintf("Number of observations: %d\n", nrow(object$data)))
  cat(sprintf("Number of predictors: %d\n\n", ncol(object$centers)))

  cat("Cluster characteristics:\n")
  for (i in seq_len(object$k)) {
    cat(sprintf("\nCluster %d:\n", i))
    cat(sprintf("  Size: %d (%.1f%%)\n",
                object$cluster_sizes[i],
                object$cluster_sizes[i] / sum(object$cluster_sizes) * 100))
    cat(sprintf("  Within SS: %.4f\n", object$withinss[i]))
  }

  cat(sprintf("\nTotal within-cluster SS: %.4f\n", object$tot_withinss))
  cat(sprintf("Between-cluster SS: %.4f\n", object$betweenss))
  cat(sprintf("Variance explained: %.2f%%\n", object$variance_explained * 100))

  invisible(object)
}


# Collect method ----

#' Collect cluster assignments from k-means
#'
#' Extract cluster assignments from k-means results and add them to a tibble.
#'
#' @param x A `stat_kmeans` object
#' @param data A tibble. If NULL, uses the original data from the clustering.
#' @param name Character. Name for the cluster assignment column. Default is "cluster".
#' @param ... Additional arguments (reserved)
#'
#' @return A tibble with cluster assignments added
#'
#' @examples
#' \dontrun{
#' km <- boteft %>% stat_kmeans(k = 3)
#'
#' # Add cluster assignments
#' collect(km)
#'
#' # Custom column name
#' collect(km, name = "k3_cluster")
#' }
#'
#' @export
collect.stat_kmeans <- function(x, data = NULL, name = "cluster", ...) {

  # Use original data if not provided
  if (is.null(data)) {
    data <- x$data
  }

  # Check if column exists
  if (name %in% names(data)) {
    stop(sprintf("Column '%s' already exists", name))
  }

  # Add cluster assignments as factor
  data[[name]] <- factor(x$model$cluster)

  data
}


# Transduce method ----

#' Transduce k-means cluster centers to shapes
#'
#' Reconstruct shapes at cluster centers or at positions relative to centers.
#'
#' @param object A `stat_kmeans` object
#' @param positions A tibble with cluster specifications. Must contain a `cluster`
#'   column with cluster numbers (1 to k).
#'
#' @return A tibble with reconstructed coefficients and shapes
#'
#' @details
#' For k-means, transduce simply retrieves the cluster centers from the model
#' and reconstructs the corresponding shapes. Since cluster centers are already
#' in the original predictor space, no inverse transformation is needed.
#'
#' @examples
#' \dontrun{
#' km <- boteft %>% stat_kmeans(k = 3)
#'
#' # Get shapes at all cluster centers
#' centers <- transduce(km, tibble(cluster = 1:3))
#'
#' # Get shape at specific cluster
#' center1 <- transduce(km, tibble(cluster = 1))
#' }
#'
#' @export
transduce.stat_kmeans <- function(object, positions) {

  if (!is.data.frame(positions)) {
    stop("positions must be a tibble or data frame")
  }

  if (nrow(positions) == 0) {
    stop("positions tibble is empty")
  }

  # Check for cluster column
  if (!"cluster" %in% names(positions)) {
    stop("positions must contain a 'cluster' column")
  }

  # Validate cluster numbers
  clusters <- positions$cluster
  if (any(clusters < 1 | clusters > object$k)) {
    stop(sprintf("cluster numbers must be between 1 and %d", object$k))
  }

  # Identify coe columns from predictor_cols
  coe_cols <- object$predictor_cols[
    sapply(object$data[object$predictor_cols], function(x) "coe" %in% class(x))
  ]

  if (length(coe_cols) == 0) {
    stop("No coe columns found in model predictors")
  }

  # Get coe classes
  coe_classes <- sapply(object$data[coe_cols], function(col) class(col)[1])

  # Reconstruct coefficients for each position
  n_positions <- nrow(positions)

  # Pre-allocate lists for coe columns
  coe_lists <- vector("list", length(coe_cols))
  names(coe_lists) <- coe_cols

  for (i in seq_along(coe_cols)) {
    coe_lists[[i]] <- vector("list", n_positions)
  }

  # For each position (cluster)
  for (i in seq_len(n_positions)) {

    cluster_num <- clusters[i]

    # Get center for this cluster (already in predictor space!)
    center <- object$centers[cluster_num, ]

    # Split center back into coe columns
    col_start <- 1
    for (j in seq_along(coe_cols)) {
      # Get number of coefficients for this coe column
      example_coe <- object$data[[coe_cols[j]]][[1]]
      n_coe <- length(example_coe)
      col_end <- col_start + n_coe - 1

      # Extract coefficients
      coe_vec <- center[col_start:col_end]
      names(coe_vec) <- names(example_coe)
      class(coe_vec) <- c(coe_classes[j], "coe", "numeric")

      # Store
      coe_lists[[j]][[i]] <- coe_vec

      col_start <- col_end + 1
    }
  }

  # Build result tibble: positions + coe columns
  result <- positions

  # Add coe columns
  for (j in seq_along(coe_cols)) {
    coe_col <- coe_cols[j]
    class(coe_lists[[j]]) <- c(coe_classes[j], "coe", "list")
    result[[coe_col]] <- coe_lists[[j]]
  }

  # Add inverse columns by calling the proper inverse function
  for (j in seq_along(coe_cols)) {
    coe_col <- coe_cols[j]
    shape_col <- paste0(coe_col, "_i")
    coe_class <- coe_classes[j]

    # Get the inverse function name
    inverse_fn_name <- paste0(coe_class, "_i")

    # Get the function
    inverse_fn <- tryCatch(
      get(inverse_fn_name, envir = asNamespace("Momocs2")),
      error = function(e) {
        stop(sprintf("Inverse function '%s' not found for class '%s'",
                     inverse_fn_name, coe_class))
      }
    )

    # Call it on the coe column
    result[[shape_col]] <- inverse_fn(result[[coe_col]])
  }

  tibble::as_tibble(result)
}


# Plot method ----

#' Plot k-means clustering results
#'
#' Visualize k-means clustering with various plot types.
#'
#' @param x A `stat_kmeans` object
#' @param type Character. Type of plot:
#'   * `"clusters"`: Cluster visualization (requires PCA projection)
#'   * `"withinss"`: Within-cluster sum of squares by cluster
#' @param color Column name (bare or quoted) for coloring points (for cluster plot).
#'   If NULL, colors by cluster assignment.
#' @param labels Column name (bare or quoted) for text labels.
#' @param ... Additional arguments (reserved)
#'
#' @return NULL (invisibly). Draws plot as side effect.
#'
#' @examples
#' \dontrun{
#' km <- boteft %>% stat_kmeans(k = 3)
#'
#' # Cluster plot (colored by cluster)
#' plot(km)
#'
#' # Color by original grouping
#' plot(km, color = type)
#'
#' # Within-cluster SS plot
#' plot(km, type = "withinss")
#' }
#'
#' @export
plot.stat_kmeans <- function(x, type = c("clusters", "withinss"),
                             color = NULL, labels = NULL, ...) {

  type <- match.arg(type)

  if (type == "withinss") {
    # Bar plot of within-cluster SS
    barplot(x$withinss,
            names.arg = paste0("C", seq_len(x$k)),
            xlab = "Cluster",
            ylab = "Within-cluster SS",
            main = "Within-cluster Sum of Squares",
            col = rainbow(x$k, alpha = 0.5),
            las = 1)

  } else if (type == "clusters") {
    # For cluster visualization, we need to project to 2D
    # Use PCA on the predictor matrix

    coe_cols <- x$predictor_cols[
      sapply(x$data[x$predictor_cols], function(x) "coe" %in% class(x))
    ]
    other_cols <- setdiff(x$predictor_cols, coe_cols)

    predictor_matrix <- build_predictor_matrix(x$data, coe_cols, other_cols)

    # Quick PCA for visualization
    pca_temp <- stats::prcomp(predictor_matrix, center = TRUE, scale. = FALSE)
    pc_scores <- pca_temp$x[, 1:2]

    # Project cluster centers
    center_pcs <- predict(pca_temp, newdata = x$centers)[, 1:2]

    var1 <- round(summary(pca_temp)$importance[2, 1] * 100, 1)
    var2 <- round(summary(pca_temp)$importance[2, 2] * 100, 1)

    # Get color values
    color_quo <- rlang::enquo(color)
    if (rlang::quo_is_null(color_quo)) {
      # Color by cluster
      color_vals <- factor(x$model$cluster)
    } else {
      color_vals <- get_plot_variable(x$data, color_quo)
    }

    # Get label values
    label_quo <- rlang::enquo(labels)
    label_vals <- get_plot_variable(x$data, label_quo)

    # Determine colors
    if (!is.null(color_vals)) {
      if (is.factor(color_vals) || is.character(color_vals)) {
        color_vals <- as.factor(color_vals)
        cols <- rainbow(nlevels(color_vals), alpha = 0.5)[as.integer(color_vals)]
      } else if (is.numeric(color_vals)) {
        col_ramp <- colorRampPalette(c("blue", "red"))
        color_indices <- cut(color_vals, breaks = 100, labels = FALSE)
        cols <- col_ramp(100)[color_indices]
        cols <- paste0(cols, "80")
      } else {
        cols <- rgb(0, 0, 0, 0.5)
      }
    } else {
      cols <- rgb(0, 0, 0, 0.5)
    }

    # Create plot
    plot(pc_scores[, 1], pc_scores[, 2],
         xlab = sprintf("PC1 (%s%%)", var1),
         ylab = sprintf("PC2 (%s%%)", var2),
         main = sprintf("K-means Clusters (k=%d)", x$k),
         type = "n",
         asp = 1,
         las = 1)

    abline(h = 0, v = 0, col = "gray", lty = 2)

    # Add points
    points(pc_scores[, 1], pc_scores[, 2], pch = 19, col = cols)

    # Add cluster centers
    points(center_pcs[, 1], center_pcs[, 2],
           pch = 3, cex = 2, lwd = 3, col = "black")
    text(center_pcs[, 1], center_pcs[, 2],
         labels = seq_len(x$k), pos = 3, cex = 1.2, font = 2)

    # Add labels if requested
    if (!is.null(label_vals)) {
      text(pc_scores[, 1], pc_scores[, 2], labels = label_vals,
           cex = 0.7, pos = 3)
    }
  }

  invisible(NULL)
}
