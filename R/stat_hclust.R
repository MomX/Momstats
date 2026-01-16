# stat_hclust --------

#' Hierarchical clustering for morphometric data
#'
#' Perform hierarchical clustering on coefficient data.
#'
#' @param data A tibble with coefficient columns
#' @param formula A formula specifying predictors. Can be:
#'   * Missing: auto-detects single coe column
#'   * Bare column name: `coe`
#'   * Formula: `~ coe`, `~ coe + size`, `~ coe1 + coe2`
#' @param method Character. Agglomeration method for hierarchical clustering.
#'   One of "ward.D2" (default), "single", "complete", "average", "mcquitty",
#'   "median", "centroid". See [stats::hclust()].
#' @param dist_method Character. Distance metric. Default is "euclidean".
#'   Other options: "manhattan", "maximum", "canberra", "binary", "minkowski".
#'   See [stats::dist()].
#' @param center Logical. Should data be centered? Default `TRUE`.
#' @param scale Logical or NULL. Should data be scaled to unit variance?
#'   If `NULL` (default), automatically determined based on predictor types.
#' @param k Integer. Optional. If provided, cuts tree to k clusters.
#' @param h Numeric. Optional. If provided, cuts tree at height h.
#' @param ... Additional arguments passed to [stats::hclust()]
#'
#' @return An object of class `c("stat_hclust", "momstats")` containing:
#' * `data`: Original tibble (unchanged)
#' * `model`: The [stats::hclust()] object
#' * `dist_matrix`: Distance matrix used for clustering
#' * `method`: Agglomeration method
#' * `dist_method`: Distance metric used
#' * `call`: The function call
#' * `formula`: Formula used (if any)
#' * `predictor_cols`: All predictor column names
#' * `center`: Logical, was centering applied
#' * `scale`: Logical, was scaling applied
#' * `k`: Number of clusters (if tree was cut)
#' * `h`: Height threshold (if tree was cut by height)
#' * `clusters`: Cluster assignments (if tree was cut)
#'
#' @details
#' `stat_hclust()` provides hierarchical clustering for morphometric data with
#' proper handling of coefficient columns and optional covariates.
#'
#' ## Agglomeration methods
#'
#' * `"ward.D2"` (default): Ward's minimum variance method - typically best for
#'   morphometric data as it minimizes within-cluster variance
#' * `"complete"`: Maximum distance between clusters
#' * `"average"`: UPGMA - average distance between clusters
#' * `"single"`: Minimum distance (tends to chain)
#'
#' ## Distance metrics
#'
#' * `"euclidean"` (default): Standard L2 distance
#' * `"manhattan"`: L1 distance, more robust to outliers
#' * `"maximum"`: Chebyshev distance
#'
#' ## Cutting the tree
#'
#' The tree can be cut during creation or later:
#' ```r
#' # Cut during creation
#' hc <- stat_hclust(data, k = 4)
#'
#' # Cut later via collect
#' hc <- stat_hclust(data)
#' data_clustered <- collect(hc, k = 4)
#' ```
#'
#' ## Getting results
#'
#' Use [collect()] to add cluster assignments:
#' ```r
#' hc <- boteft %>% stat_hclust()
#' boteft_clustered <- collect(hc, k = 4)
#' ```
#'
#' Use [transduce()] to reconstruct shapes:
#' ```r
#' # Cluster centers (if tree was cut)
#' centers <- transduce(hc, tibble(cluster = 1:4))
#'
#' # Any node in the tree
#' node_shapes <- transduce(hc, tibble(node = c(45, 50, 55)))
#' ```
#'
#' @examples
#' \dontrun{
#' # Basic hierarchical clustering
#' hc1 <- boteft %>% stat_hclust()
#'
#' # With specific method
#' hc2 <- boteft %>% stat_hclust(method = "average")
#'
#' # Cut tree during creation
#' hc3 <- boteft %>% stat_hclust(k = 4)
#'
#' # Different distance
#' hc4 <- boteft %>% stat_hclust(dist_method = "manhattan")
#'
#' # Add cluster assignments
#' boteft_clustered <- collect(hc1, k = 4)
#'
#' # Get cluster center shapes
#' centers <- transduce(hc3, tibble(cluster = 1:4))
#'
#' # Get shapes at internal nodes
#' nodes <- transduce(hc1, tibble(node = c(45, 50)))
#'
#' # Plot (requires ape package)
#' plot(hc1)  # Unrooted phylogram (default)
#' plot(hc1, color = type)  # Color by grouping
#' plot(hc1, type = "dendrogram")  # Classic dendrogram
#' }
#'
#' @seealso [stats::hclust()], [stats::dist()], [collect.stat_hclust()],
#'   [plot.stat_hclust()], [transduce.stat_hclust()]
#'
#' @export
stat_hclust <- function(data, formula = NULL, method = "ward.D2",
                        dist_method = "euclidean",
                        center = TRUE, scale = NULL,
                        k = NULL, h = NULL, ...) {

  # Store call
  call <- match.call()

  # Check input
  if (!is.data.frame(data)) {
    stop("data must be a tibble or data frame")
  }

  # Check k and h not both specified
  if (!is.null(k) && !is.null(h)) {
    stop("Specify either k or h, not both")
  }

  # Parse formula to get predictor columns
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

  # Center and scale
  if (center) {
    predictor_matrix <- scale(predictor_matrix, center = TRUE, scale = FALSE)
  }
  if (scale) {
    predictor_matrix <- scale(predictor_matrix, center = FALSE, scale = TRUE)
  }

  # Compute distance matrix
  dist_matrix <- stats::dist(predictor_matrix, method = dist_method)

  # Run hierarchical clustering
  hc_model <- stats::hclust(dist_matrix, method = method, ...)

  # Cut tree if requested
  clusters <- NULL
  if (!is.null(k)) {
    clusters <- stats::cutree(hc_model, k = k)
  } else if (!is.null(h)) {
    clusters <- stats::cutree(hc_model, h = h)
    k <- length(unique(clusters))
  }

  # Build result object
  result <- list(
    data = data,
    model = hc_model,
    dist_matrix = dist_matrix,
    method = method,
    dist_method = dist_method,
    call = call,
    formula = formula,
    predictor_cols = predictor_cols,
    center = center,
    scale = scale,
    k = k,
    h = h,
    clusters = clusters
  )

  class(result) <- c("stat_hclust", "momstats")
  result
}


# Print method ----

#' @export
print.stat_hclust <- function(x, ...) {
  cat("Hierarchical Clustering\n")
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
              nrow(x$data), length(x$predictor_cols), pred_desc))

  cat(sprintf("* Method: %s\n", x$method))
  cat(sprintf("* Distance: %s\n", x$dist_method))

  # Centering/scaling info
  center_str <- if (x$center) "TRUE" else "FALSE"
  scale_str <- if (x$scale) "TRUE" else "FALSE"

  # Build reason message
  if (!x$scale) {
    reason <- "(pure coefficient data)"
  } else if (length(other_cols) > 0) {
    reason <- "(mixed predictors)"
  } else if (length(coe_cols) > 1) {
    reason <- "(multiple coefficient types)"
  } else {
    reason <- "(non-EFT coefficients)"
  }

  cat(sprintf("* Centering: %s, Scaling: %s %s\n", center_str, scale_str, reason))

  # Clustering info
  if (!is.null(x$clusters)) {
    cat(sprintf("\n* Tree cut into %d clusters\n", x$k))
    cluster_sizes <- table(x$clusters)
    cat("  Cluster sizes:", paste(cluster_sizes, collapse = ", "), "\n")
  } else {
    cat("\n* Tree not cut (use collect(hc, k = ...) to cut)\n")
  }

  cat("\n* collect() to add cluster assignments to data\n")
  cat("* transduce() to reconstruct shapes at nodes/clusters\n")
  cat("* plot() to visualize tree (requires ape package)\n")

  invisible(x)
}


# Summary method ----

#' @export
summary.stat_hclust <- function(object, ...) {
  cat("Hierarchical Clustering Summary\n")
  cat("==================================================\n\n")

  cat(sprintf("Number of observations: %d\n", nrow(object$data)))
  cat(sprintf("Number of predictors: %d\n", length(object$predictor_cols)))
  cat(sprintf("Agglomeration method: %s\n", object$method))
  cat(sprintf("Distance metric: %s\n\n", object$dist_method))

  if (!is.null(object$clusters)) {
    cat(sprintf("Tree cut into %d clusters:\n", object$k))
    cluster_table <- table(object$clusters)
    for (i in seq_along(cluster_table)) {
      cat(sprintf("  Cluster %d: %d observations (%.1f%%)\n",
                  i, cluster_table[i],
                  cluster_table[i] / sum(cluster_table) * 100))
    }
  } else {
    cat("Tree not yet cut into clusters\n")
  }

  cat("\nTree structure:\n")
  cat(sprintf("  Height range: [%.4f, %.4f]\n",
              min(object$model$height), max(object$model$height)))
  cat(sprintf("  Number of merges: %d\n", nrow(object$model$merge)))

  invisible(object)
}


# Collect method ----

#' Collect cluster assignments from hierarchical clustering
#'
#' Extract cluster assignments from hierarchical clustering and add them to a tibble.
#'
#' @param x A `stat_hclust` object
#' @param data A tibble. If NULL, uses the original data.
#' @param k Integer. Number of clusters to cut tree into.
#' @param h Numeric. Height at which to cut tree.
#' @param name Character. Name for cluster assignment column. Default is "cluster".
#' @param ... Additional arguments (reserved)
#'
#' @return A tibble with cluster assignments added
#'
#' @details
#' If the tree was already cut during `stat_hclust()`, those clusters are used
#' unless `k` or `h` is specified here.
#'
#' @examples
#' \dontrun{
#' hc <- boteft %>% stat_hclust()
#'
#' # Cut and collect
#' collect(hc, k = 4)
#'
#' # Different k
#' collect(hc, k = 5)
#'
#' # Cut by height
#' collect(hc, h = 10)
#'
#' # Custom column name
#' collect(hc, k = 3, name = "hc_cluster")
#' }
#'
#' @export
collect.stat_hclust <- function(x, data = NULL, k = NULL, h = NULL,
                                name = "cluster", ...) {

  # Use original data if not provided
  if (is.null(data)) {
    data <- x$data
  }

  # Check if column exists
  if (name %in% names(data)) {
    stop(sprintf("Column '%s' already exists", name))
  }

  # Determine clusters
  if (!is.null(k)) {
    # Cut at k
    clusters <- stats::cutree(x$model, k = k)
  } else if (!is.null(h)) {
    # Cut at height h
    clusters <- stats::cutree(x$model, h = h)
  } else if (!is.null(x$clusters)) {
    # Use existing cut
    clusters <- x$clusters
  } else {
    stop("Specify k or h to cut tree, or cut tree during stat_hclust()")
  }

  # Add cluster assignments as factor
  data[[name]] <- factor(clusters)

  data
}


# Transduce method ----

#' Transduce hierarchical clustering nodes to shapes
#'
#' Reconstruct shapes at cluster centers or at any node in the tree.
#'
#' @param object A `stat_hclust` object
#' @param positions A tibble specifying positions. Can contain:
#'   * `cluster` column: cluster numbers (if tree was cut)
#'   * `node` column: node numbers in the tree (tips are 1:n, internal nodes are n+1:...)
#'
#' @return A tibble with reconstructed coefficients and shapes
#'
#' @details
#' For hierarchical clustering, transduce can reconstruct shapes at:
#'
#' **Cluster centers** (if tree was cut with k or h):
#' ```r
#' transduce(hc, tibble(cluster = 1:4))
#' ```
#' Returns the mean shape of each cluster.
#'
#' **Any node in the tree**:
#' ```r
#' transduce(hc, tibble(node = c(45, 50, 55)))
#' ```
#' Returns the mean shape of all descendants of that node. This allows
#' visualization of "ancestral" shapes at internal nodes.
#'
#' Node numbering follows [stats::hclust()] convention:
#' * Nodes 1 to n: original observations (tips)
#' * Nodes n+1 to 2n-1: internal nodes from merges
#'
#' @examples
#' \dontrun{
#' hc <- boteft %>% stat_hclust(k = 4)
#'
#' # Cluster centers
#' centers <- transduce(hc, tibble(cluster = 1:4))
#'
#' # Shapes at internal nodes
#' hc2 <- boteft %>% stat_hclust()
#' internal <- transduce(hc2, tibble(node = c(45, 50, 55, 60)))
#'
#' # Root node (all data)
#' root <- transduce(hc2, tibble(node = 2 * nrow(boteft) - 1))
#' }
#'
#' @export
transduce.stat_hclust <- function(object, positions) {

  if (!is.data.frame(positions)) {
    stop("positions must be a tibble or data frame")
  }

  if (nrow(positions) == 0) {
    stop("positions tibble is empty")
  }

  # Check what kind of positions we have
  has_cluster <- "cluster" %in% names(positions)
  has_node <- "node" %in% names(positions)

  if (!has_cluster && !has_node) {
    stop("positions must contain either 'cluster' or 'node' column")
  }

  if (has_cluster && has_node) {
    stop("positions should contain either 'cluster' OR 'node', not both")
  }

  # Identify coe columns
  coe_cols <- object$predictor_cols[
    sapply(object$data[object$predictor_cols], function(x) "coe" %in% class(x))
  ]

  if (length(coe_cols) == 0) {
    stop("No coe columns found in model predictors")
  }

  # Get coe classes
  coe_classes <- sapply(object$data[coe_cols], function(col) class(col)[1])

  # Build predictor matrix (original, unscaled)
  predictor_matrix <- build_predictor_matrix(object$data, coe_cols,
                                             setdiff(object$predictor_cols, coe_cols))

  n_obs <- nrow(predictor_matrix)

  # Get descendants for each requested position
  if (has_cluster) {
    # Cluster mode
    if (is.null(object$clusters)) {
      stop("Tree was not cut. Use node-based transduce or cut tree first.")
    }

    clusters_requested <- positions$cluster
    if (any(clusters_requested < 1 | clusters_requested > object$k)) {
      stop(sprintf("Cluster numbers must be between 1 and %d", object$k))
    }

    # For each cluster, get observations
    descendant_lists <- lapply(clusters_requested, function(cl) {
      which(object$clusters == cl)
    })

  } else {
    # Node mode
    nodes_requested <- positions$node

    # Validate nodes
    max_node <- 2 * n_obs - 1
    if (any(nodes_requested < 1 | nodes_requested > max_node)) {
      stop(sprintf("Node numbers must be between 1 and %d", max_node))
    }

    # For each node, get all descendant tips
    descendant_lists <- lapply(nodes_requested, function(node) {
      get_node_descendants(object$model, node, n_obs)
    })
  }

  # Reconstruct coefficients for each position
  n_positions <- nrow(positions)

  # Pre-allocate lists for coe columns
  coe_lists <- vector("list", length(coe_cols))
  names(coe_lists) <- coe_cols

  for (i in seq_along(coe_cols)) {
    coe_lists[[i]] <- vector("list", n_positions)
  }

  # For each position
  for (i in seq_len(n_positions)) {

    # Get descendants for this position
    descendants <- descendant_lists[[i]]

    # Calculate mean of descendants
    if (length(descendants) == 0) {
      stop(sprintf("No descendants found for position %d", i))
    }

    mean_predictors <- colMeans(predictor_matrix[descendants, , drop = FALSE])

    # Split mean back into coe columns
    col_start <- 1
    for (j in seq_along(coe_cols)) {
      # Get number of coefficients for this coe column
      example_coe <- object$data[[coe_cols[j]]][[1]]
      n_coe <- length(example_coe)
      col_end <- col_start + n_coe - 1

      # Extract coefficients
      coe_vec <- mean_predictors[col_start:col_end]
      names(coe_vec) <- names(example_coe)
      class(coe_vec) <- c(coe_classes[j], "coe", "numeric")

      # Store
      coe_lists[[j]][[i]] <- coe_vec

      col_start <- col_end + 1
    }
  }

  # Build result tibble
  result <- positions

  # Add coe columns
  for (j in seq_along(coe_cols)) {
    coe_col <- coe_cols[j]
    class(coe_lists[[j]]) <- c(coe_classes[j], "coe", "list")
    result[[coe_col]] <- coe_lists[[j]]
  }

  # Add inverse columns
  for (j in seq_along(coe_cols)) {
    coe_col <- coe_cols[j]
    shape_col <- paste0(coe_col, "_i")
    coe_class <- coe_classes[j]

    # Get inverse function
    inverse_fn_name <- paste0(coe_class, "_i")
    inverse_fn <- tryCatch(
      get(inverse_fn_name, envir = asNamespace("Momocs2")),
      error = function(e) {
        stop(sprintf("Inverse function '%s' not found for class '%s'",
                     inverse_fn_name, coe_class))
      }
    )

    # Call it
    result[[shape_col]] <- inverse_fn(result[[coe_col]])
  }

  tibble::as_tibble(result)
}


# Helper: Get all descendant tips of a node ----

#' Get descendant tips for a node in hclust tree
#'
#' @param hc An hclust object
#' @param node Node number
#' @param n_obs Number of observations (tips)
#'
#' @return Vector of tip (observation) indices
#'
#' @keywords internal
#' @noRd
get_node_descendants <- function(hc, node, n_obs) {

  # Tips (1 to n_obs) return themselves
  if (node <= n_obs) {
    return(node)
  }

  # Internal nodes: recursively get descendants
  merge_idx <- node - n_obs
  left <- hc$merge[merge_idx, 1]
  right <- hc$merge[merge_idx, 2]

  # Negative values in merge are tips, positive are internal nodes
  left_descendants <- if (left < 0) {
    abs(left)
  } else {
    get_node_descendants(hc, left + n_obs, n_obs)
  }

  right_descendants <- if (right < 0) {
    abs(right)
  } else {
    get_node_descendants(hc, right + n_obs, n_obs)
  }

  c(left_descendants, right_descendants)
}


# Plot method ----

#' Plot hierarchical clustering tree
#'
#' Visualize hierarchical clustering with various tree layouts.
#'
#' @param x A `stat_hclust` object
#' @param type Character. Type of plot:
#'   * `"unrooted"`: Unrooted phylogram (default, requires ape)
#'   * `"dendrogram"`: Classic dendrogram
#'   * `"fan"`: Fan/circular tree (requires ape)
#'   * `"phylogram"`: Rectangular phylogram (requires ape)
#' @param color Column name (bare or quoted) for coloring tips. If NULL and tree
#'   was cut, colors by cluster assignment.
#' @param tip_labels Logical. Show tip labels? Default TRUE for small trees (<50 tips).
#' @param ... Additional arguments passed to plotting functions
#'
#' @return NULL (invisibly). Draws plot as side effect.
#'
#' @details
#' Most plot types require the `ape` package. If not installed, falls back to
#' base dendrogram plot.
#'
#' @examples
#' \dontrun{
#' hc <- boteft %>% stat_hclust(k = 4)
#'
#' # Unrooted tree (default)
#' plot(hc)
#'
#' # Color by original grouping
#' plot(hc, color = type)
#'
#' # Color by clusters
#' plot(hc)  # Automatic if tree was cut
#'
#' # Classic dendrogram
#' plot(hc, type = "dendrogram")
#'
#' # Fan tree
#' plot(hc, type = "fan")
#' }
#'
#' @export
plot.stat_hclust <- function(x, type = c("unrooted", "dendrogram", "fan", "phylogram"),
                             color = NULL, tip_labels = NULL, ...) {

  type <- match.arg(type)

  # Determine tip labels
  if (is.null(tip_labels)) {
    tip_labels <- nrow(x$data) < 50
  }

  # Get color values
  color_quo <- rlang::enquo(color)
  if (rlang::quo_is_null(color_quo)) {
    # Default: color by cluster if tree was cut
    if (!is.null(x$clusters)) {
      color_vals <- factor(x$clusters)
    } else {
      color_vals <- NULL
    }
  } else {
    color_vals <- get_plot_variable(x$data, color_quo)
  }

  # Determine colors
  if (!is.null(color_vals)) {
    if (is.factor(color_vals) || is.character(color_vals)) {
      color_vals <- as.factor(color_vals)
      tip_colors <- rainbow(nlevels(color_vals), alpha = 0.8)[as.integer(color_vals)]
    } else if (is.numeric(color_vals)) {
      col_ramp <- colorRampPalette(c("blue", "red"))
      color_indices <- cut(color_vals, breaks = 100, labels = FALSE)
      tip_colors <- col_ramp(100)[color_indices]
    } else {
      tip_colors <- "black"
    }
  } else {
    tip_colors <- "black"
  }

  # Plot based on type
  if (type == "dendrogram") {
    # Base graphics dendrogram
    plot(x$model, labels = if (tip_labels) x$model$labels else FALSE, ...)

    # Color tips if requested
    if (!is.null(color_vals) && !is.null(x$clusters)) {
      # Add colored rectangles around clusters
      rect.hclust(x$model, k = x$k, border = rainbow(x$k, alpha = 0.8))
    }

  } else {
    # Phylo-style plots require ape
    if (!requireNamespace("ape", quietly = TRUE)) {
      message("Package 'ape' required for phylo-style plots. Falling back to dendrogram.")
      message("Install with: install.packages('ape')")
      plot(x$model, labels = if (tip_labels) x$model$labels else FALSE, ...)
      return(invisible(NULL))
    }

    # Convert to phylo
    phylo_tree <- ape::as.phylo(x$model)

    # Plot with ape
    if (type == "unrooted") {
      ape::plot.phylo(phylo_tree, type = "unrooted",
                      show.tip.label = tip_labels,
                      tip.color = tip_colors,
                      ...)
    } else if (type == "fan") {
      ape::plot.phylo(phylo_tree, type = "fan",
                      show.tip.label = tip_labels,
                      tip.color = tip_colors,
                      ...)
    } else if (type == "phylogram") {
      ape::plot.phylo(phylo_tree, type = "phylogram",
                      show.tip.label = tip_labels,
                      tip.color = tip_colors,
                      ...)
    }
  }

  invisible(NULL)
}
