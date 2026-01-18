# stat_meanshape.R

#' Calculate mean (or aggregated) shapes
#'
#' Compute mean or other aggregated shapes, either globally or by groups.
#'
#' @param data A tibble with coefficient columns
#' @param formula A formula or bare column name specifying grouping:
#'   * Missing/NULL: Global mean across all observations
#'   * Bare name: `status` → group by this column, use all coe columns
#'   * Formula: `status ~ .` → explicit all coe columns
#'   * Formula: `status ~ VD` → specific coe column(s)
#'   * Formula: `~ status` → one-sided formula (group by status)
#' @param .fn Aggregation function. Default `mean`. Can be `median` or any
#'   function that takes a matrix column and returns a scalar.
#' @param ... Additional arguments (reserved for future use)
#'
#' @return An object of class `c("stat_meanshape", "momstats")` containing:
#' * `data`: Original tibble (unchanged)
#' * `summary`: Tibble with group(s), n, and mean coe column(s)
#' * `method`: "meanshape"
#' * `call`: The function call
#' * `formula`: Formula used (if any)
#' * `group_cols`: Grouping column name(s) (NULL if global)
#' * `coe_cols`: Coefficient column names
#' * `.fn`: Aggregation function name
#'
#' @details
#' `stat_meanshape()` provides a unified interface for aggregating shapes.
#' It handles:
#'
#' * **Global aggregation**: No grouping, single mean shape
#' * **Grouped aggregation**: By one or more grouping variables
#' * **Flexible syntax**: Bare names or formula syntax
#' * **Multiple coe columns**: Aggregates all coe columns simultaneously
#' * **Custom aggregation**: mean (default), median, or custom function
#' * **Automatic inverse**: Computes reconstructed shapes via `coe_class_i()`
#'
#' ## Formula syntax
#'
#' The formula specifies grouping and optionally which coe columns:
#' * Missing: Global mean of all coe columns
#' * `status`: Group by status, all coe columns
#' * `~ status`: Same as above
#' * `status ~ .`: Explicit all coe columns
#' * `status ~ VD`: Only VD coe column
#' * `status ~ VD + VL`: Multiple specific coe columns
#'
#' ## Aggregation functions
#'
#' * `.fn = mean` (default): Arithmetic mean
#' * `.fn = median`: Median shape
#' * Custom: Any function taking a numeric vector and returning a scalar
#'
#' ## Getting results
#'
#' Use [collect()] to add mean shapes back to original data (matched by group):
#' ```r
#' ms <- vitis %>% stat_meanshape(status)
#' vitis_with_means <- collect(ms)  # Adds mean_VD, mean_VL, etc.
#' ```
#'
#' @examples
#' \dontrun{
#' # Global mean
#' ms1 <- vitis %>% stat_meanshape()
#'
#' # By group
#' ms2 <- vitis %>% stat_meanshape(status)
#'
#' # Formula syntax
#' ms3 <- vitis %>% stat_meanshape(~ status)
#' ms4 <- vitis %>% stat_meanshape(status ~ .)
#'
#' # Specific coe column
#' ms5 <- vitis %>% stat_meanshape(status ~ VD)
#'
#' # Median instead of mean
#' ms6 <- vitis %>% stat_meanshape(status, .fn = median)
#'
#' # Add mean shapes to data
#' vitis_means <- collect(ms2)
#'
#' # Plot mean shapes
#' plot(ms2)  # Mean shapes only
#' plot(ms2, with_data = TRUE)  # Mean shapes + contributing data
#' }
#'
#' @seealso [collect.stat_meanshape()], [plot.stat_meanshape()]
#'
#' @export
stat_meanshape <- function(data, formula = NULL, .fn = mean, ...) {

  # Store call
  call <- match.call()

  # Check input
  if (!is.data.frame(data)) {
    stop("data must be a tibble or data frame")
  }

  # Get all coe columns in data
  all_coe_cols <- Momocs2::get_all_coe_cols(data)

  if (length(all_coe_cols) == 0) {
    stop("No coe columns found in data")
  }

  # Parse formula (handle missing, bare name, formula)
  formula_expr <- substitute(formula)

  # Case 1: Missing/NULL → global mean of all coe
  if (is.null(formula_expr) || identical(formula_expr, quote(expr = ))) {
    group_cols <- NULL
    coe_cols <- all_coe_cols

  } else if (!is.call(formula_expr) || !identical(formula_expr[[1]], as.name("~"))) {
    # Case 2: Bare name → group by this, use all coe
    group_cols <- as.character(formula_expr)

    # Check if it exists
    if (!group_cols %in% names(data)) {
      stop(sprintf("Column '%s' not found in data", group_cols))
    }

    coe_cols <- all_coe_cols

  } else {
    # Case 3: Formula
    # Check if one-sided or two-sided
    if (length(formula_expr) == 2) {
      # One-sided: ~ status
      lhs <- NULL
      rhs <- formula_expr[[2]]
    } else if (length(formula_expr) == 3) {
      # Two-sided: status ~ .
      lhs <- formula_expr[[2]]
      rhs <- formula_expr[[3]]
    } else {
      stop("Invalid formula")
    }

    # Parse grouping (LHS or RHS if one-sided)
    if (is.null(lhs)) {
      # One-sided: grouping from RHS
      group_cols <- as.character(rhs)
      if (!group_cols %in% names(data)) {
        stop(sprintf("Column '%s' not found in data", group_cols))
      }
      coe_cols <- all_coe_cols
    } else {
      # Two-sided: grouping from LHS
      group_cols <- as.character(lhs)
      if (!group_cols %in% names(data)) {
        stop(sprintf("Column '%s' not found in data", group_cols))
      }

      # Parse RHS for coe columns
      if (identical(rhs, as.name("."))) {
        # "." → all coe
        coe_cols <- all_coe_cols
      } else {
        # Specific coe columns
        coe_cols <- parse_coe_from_rhs(rhs, data)
      }
    }
  }

  # Get function name for printing
  fn_name <- deparse(substitute(.fn))

  # Compute aggregated shapes
  if (is.null(group_cols)) {
    # Global aggregation
    summary_tbl <- compute_global_meanshape(data, coe_cols, .fn)
  } else {
    # Grouped aggregation
    summary_tbl <- compute_grouped_meanshape(data, group_cols, coe_cols, .fn)
  }

  # Add inverse shapes
  summary_tbl <- add_inverse_shapes(summary_tbl, coe_cols)

  # Build result object
  result <- list(
    data = data,
    summary = summary_tbl,
    method = "meanshape",
    call = call,
    formula = formula_expr,
    group_cols = group_cols,
    coe_cols = coe_cols,
    .fn = fn_name
  )

  class(result) <- c("stat_meanshape", "momstats")
  result
}


# Helper functions ----

#' Parse coe columns from formula RHS
#'
#' @param rhs Right-hand side expression
#' @param data A tibble
#'
#' @return Character vector of coe column names
#'
#' @keywords internal
#' @noRd
parse_coe_from_rhs <- function(rhs, data) {

  # Handle addition: coe1 + coe2
  if (is.call(rhs) && identical(rhs[[1]], as.name("+"))) {
    left <- parse_coe_from_rhs(rhs[[2]], data)
    right <- parse_coe_from_rhs(rhs[[3]], data)
    return(c(left, right))
  }

  # Simple name
  col_name <- as.character(rhs)

  if (!col_name %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col_name))
  }

  if (!"coe" %in% class(data[[col_name]])) {
    stop(sprintf("Column '%s' is not a coe column", col_name))
  }

  col_name
}


#' Compute global mean shape
#'
#' @param data A tibble
#' @param coe_cols Coefficient column names
#' @param .fn Aggregation function
#'
#' @return A tibble with one row and mean coe columns
#'
#' @keywords internal
#' @noRd
compute_global_meanshape <- function(data, coe_cols, .fn) {

  n <- nrow(data)

  # Build result tibble
  result <- tibble::tibble(n = n)

  for (col in coe_cols) {
    coe_list <- data[[col]]
    coe_class <- class(coe_list)[1]

    # Unfold to matrix
    coe_matrix <- Momocs2::unfold(coe_list)

    # Apply aggregation function column-wise
    mean_coe <- apply(coe_matrix, 2, .fn)

    # Restore class and names
    class(mean_coe) <- c(coe_class, "coe", "numeric")

    # Wrap in list-column
    result[[col]] <- list(mean_coe)
  }

  # Fix list-column classes
  for (col in coe_cols) {
    coe_class <- class(data[[col]])[1]
    class(result[[col]]) <- c(coe_class, "coe", "list")
  }

  result
}


#' Compute grouped mean shapes
#'
#' @param data A tibble
#' @param group_cols Grouping column name(s)
#' @param coe_cols Coefficient column names
#' @param .fn Aggregation function
#'
#' @return A tibble with one row per group and mean coe columns
#'
#' @keywords internal
#' @noRd
compute_grouped_meanshape <- function(data, group_cols, coe_cols, .fn) {

  # Get unique groups
  groups <- unique(data[[group_cols]])

  # Initialize result lists
  result_list <- vector("list", length(groups))

  for (i in seq_along(groups)) {
    grp <- groups[i]

    # Subset data for this group
    idx <- which(data[[group_cols]] == grp)
    n_grp <- length(idx)

    # Build row for this group
    row <- tibble::tibble(
      !!group_cols := grp,
      n = n_grp
    )

    for (col in coe_cols) {
      coe_list <- data[[col]][idx]
      coe_class <- class(coe_list)[1]

      # Unfold to matrix
      coe_matrix <- Momocs2::unfold(coe_list)

      # Apply aggregation function column-wise
      mean_coe <- apply(coe_matrix, 2, .fn)

      # Restore class and names
      class(mean_coe) <- c(coe_class, "coe", "numeric")

      # Add to row
      row[[col]] <- list(mean_coe)
    }

    result_list[[i]] <- row
  }

  # Combine rows
  result <- dplyr::bind_rows(result_list)

  # Fix list-column classes
  for (col in coe_cols) {
    coe_class <- class(data[[col]])[1]
    class(result[[col]]) <- c(coe_class, "coe", "list")
  }

  result
}


#' Add inverse shapes to summary tibble
#'
#' @param summary_tbl A tibble with coe columns
#' @param coe_cols Coefficient column names
#'
#' @return The tibble with added `_i` columns
#'
#' @keywords internal
#' @noRd
add_inverse_shapes <- function(summary_tbl, coe_cols) {

  for (col in coe_cols) {
    coe_class <- class(summary_tbl[[col]])[1]
    shape_col <- paste0(col, "_i")

    # Get inverse function
    inverse_fn_name <- paste0(coe_class, "_i")

    inverse_fn <- tryCatch(
      get(inverse_fn_name, envir = asNamespace("Momocs2")),
      error = function(e) {
        stop(sprintf("Inverse function '%s' not found for class '%s'",
                     inverse_fn_name, coe_class))
      }
    )

    # Apply inverse
    summary_tbl[[shape_col]] <- inverse_fn(summary_tbl[[col]])
  }

  summary_tbl
}


# Methods ----

#' @export
print.stat_meanshape <- function(x, ...) {
  cat("Mean Shape Calculation\n")
  cat("---------------------------------------\n")

  # Aggregation function
  cat(sprintf("* Aggregation: %s\n", x$.fn))

  # Grouping info
  if (is.null(x$group_cols)) {
    cat(sprintf("* Global aggregation (%d observations)\n", nrow(x$data)))
  } else {
    n_groups <- nrow(x$summary)
    cat(sprintf("* Grouped by: '%s' (%d groups)\n", x$group_cols, n_groups))

    # Show group sizes
    group_info <- sprintf("%s (n=%d)",
                          x$summary[[x$group_cols]],
                          x$summary$n)
    cat(sprintf("* Groups: %s\n", paste(group_info, collapse = ", ")))
  }

  # Coe columns
  cat(sprintf("* Coefficient columns: %s\n\n", paste(x$coe_cols, collapse = ", ")))

  # Show summary tibble (first few rows)
  cat("Mean shapes:\n")
  print(x$summary, n = min(10, nrow(x$summary)))

  cat("\n* collect() to add mean shapes to data\n")
  cat("* plot() to visualize mean shapes (?plot.stat_meanshape)\n")

  invisible(x)
}


#' @export
summary.stat_meanshape <- function(object, ...) {
  cat("Mean Shape Calculation\n")
  cat(sprintf("\nAggregation function: %s\n", object$.fn))

  if (is.null(object$group_cols)) {
    cat(sprintf("Global aggregation: %d observations\n", nrow(object$data)))
  } else {
    cat(sprintf("Grouped by: '%s'\n", object$group_cols))
    cat(sprintf("Number of groups: %d\n", nrow(object$summary)))
  }

  cat(sprintf("Coefficient columns: %s\n\n", paste(object$coe_cols, collapse = ", ")))

  cat("Mean shapes:\n")
  print(object$summary)

  invisible(object)
}


#' Collect mean shapes and add to original data
#'
#' Add mean shapes back to the original data, matched by group.
#'
#' @param x A `stat_meanshape` object
#' @param data A tibble. If NULL, uses the original data from the analysis.
#' @param prefix Character. Prefix for mean shape columns. Default "mean_".
#' @param ... Additional arguments (reserved)
#'
#' @return A tibble with mean shape columns added (and their `_i` inverses)
#'
#' @details
#' For grouped mean shapes, mean shapes are matched back to original rows
#' by group membership. For global mean shapes, the same mean is added to
#' all rows.
#'
#' New columns added:
#' * `{prefix}{coe_col}`: Mean coefficient column (e.g., `mean_VD`)
#' * `{prefix}{coe_col}_i`: Mean shape (e.g., `mean_VD_i`)
#'
#' @examples
#' \dontrun{
#' ms <- vitis %>% stat_meanshape(status)
#'
#' # Add mean shapes to data
#' vitis_with_means <- collect(ms)
#'
#' # Custom prefix
#' vitis_with_means <- collect(ms, prefix = "grp_mean_")
#' }
#'
#' @export
collect.stat_meanshape <- function(x, data = NULL, prefix = "mean_", ...) {

  # Use original data if not provided
  if (is.null(data)) {
    data <- x$data
  }

  # Global or grouped?
  if (is.null(x$group_cols)) {
    # Global: add same mean to all rows
    for (col in x$coe_cols) {
      mean_col_name <- paste0(prefix, col)
      shape_col_name <- paste0(prefix, col, "_i")

      # Replicate mean for all rows
      n <- nrow(data)
      coe_class <- class(x$summary[[col]])[1]

      # Extract the single mean coe
      mean_coe <- x$summary[[col]][[1]]

      # Replicate
      coe_list <- rep(list(mean_coe), n)
      class(coe_list) <- c(coe_class, "coe", "list")
      data[[mean_col_name]] <- coe_list

      # Replicate inverse
      mean_shape <- x$summary[[paste0(col, "_i")]][[1]]
      shape_list <- rep(list(mean_shape), n)
      class(shape_list) <- c("coo", "list")
      data[[shape_col_name]] <- shape_list
    }

  } else {
    # Grouped: match by group
    for (col in x$coe_cols) {
      mean_col_name <- paste0(prefix, col)
      shape_col_name <- paste0(prefix, col, "_i")

      coe_class <- class(x$summary[[col]])[1]

      # Initialize list-columns
      coe_list <- vector("list", nrow(data))
      shape_list <- vector("list", nrow(data))

      # Match by group
      for (i in seq_len(nrow(x$summary))) {
        grp <- x$summary[[x$group_cols]][i]
        idx <- which(data[[x$group_cols]] == grp)

        mean_coe <- x$summary[[col]][[i]]
        mean_shape <- x$summary[[paste0(col, "_i")]][[i]]

        # Assign to matching rows
        for (j in idx) {
          coe_list[[j]] <- mean_coe
          shape_list[[j]] <- mean_shape
        }
      }

      class(coe_list) <- c(coe_class, "coe", "list")
      class(shape_list) <- c("coo", "list")

      data[[mean_col_name]] <- coe_list
      data[[shape_col_name]] <- shape_list
    }
  }

  data
}


#' Plot mean shapes
#'
#' Visualize mean shapes, optionally with contributing data.
#'
#' @param x A `stat_meanshape` object
#' @param all Logical. Should all contributing shapes be plotted behind the mean?
#'   Default FALSE.
#' @param facet Logical. For grouped data, should groups be offset vertically?
#'   Default TRUE. If FALSE, groups are overlaid with colors + legend.
#' @param col Character. Color(s) for mean shapes. Default "red". For facet=FALSE,
#'   vector of colors per group can be provided.
#' @param lwd Numeric. Line width for mean shapes (or cex for ldk). Default 2.
#' @param all_col Character. Color for contributing shapes when facet=TRUE.
#'   Default "gray". When facet=FALSE, ignored (uses parent colors with transparency).
#' @param all_lwd Numeric. Line width for contributing shapes (or cex for ldk). Default 0.5.
#' @param labels Logical. Show group labels/legend? Default TRUE.
#' @param label_cex Numeric. Label size. Default 0.7.
#' @param ... Additional graphical parameters (reserved)
#'
#' @export
plot.stat_meanshape <- function(x, all = FALSE, facet = TRUE,
                                col = "red", lwd = 2,
                                all_col = "#000000E6", all_lwd = 0.1,
                                labels = TRUE, label_cex = 0.7,
                                ...) {

  n_groups <- if (is.null(x$group_cols)) 1 else nrow(x$summary)
  shape_cols <- paste0(x$coe_cols, "_i")
  n_coe <- length(x$coe_cols)

  # Prepare colors for mean shapes
  if (length(col) == 1 && !facet && n_groups > 1) {
    group_cols <- rainbow(n_groups)
  } else {
    group_cols <- rep(col, length.out = n_groups)
  }

  # Prepare colors for contributing shapes
  if (!facet && n_groups > 1) {
    # Use parent colors with E6 transparency (90% transparent)
    all_group_cols <- paste0(group_cols, "E6")
  } else {
    # Faceted: use uniform color
    all_group_cols <- rep(all_col, n_groups)
  }

  # Collect all shapes
  all_shapes <- list()

  if (is.null(x$group_cols)) {
    # Mean shapes (use _i columns from summary)
    for (col_name in shape_cols) {
      all_shapes <- c(all_shapes, x$summary[[col_name]])
    }

    # Contributing shapes (use original coe columns from data, then inverse)
    if (all) {
      for (coe_col in x$coe_cols) {
        # Get inverse function
        coe_class <- class(x$data[[coe_col]])[1]
        inverse_fn_name <- paste0(coe_class, "_i")
        inverse_fn <- get(inverse_fn_name, envir = asNamespace("Momocs2"))

        # Apply inverse to get shapes
        shapes_list <- inverse_fn(x$data[[coe_col]])
        all_shapes <- c(all_shapes, shapes_list)
      }
    }
  } else {
    for (i in seq_len(n_groups)) {
      grp <- x$summary[[x$group_cols]][i]

      # Mean shapes (use _i columns from summary)
      for (col_name in shape_cols) {
        all_shapes <- c(all_shapes, list(x$summary[[col_name]][[i]]))
      }

      # Contributing shapes (use original coe columns from data)
      if (all) {
        idx <- which(x$data[[x$group_cols]] == grp)
        for (j in idx) {
          for (coe_col in x$coe_cols) {
            # Get inverse function
            coe_class <- class(x$data[[coe_col]])[1]
            inverse_fn_name <- paste0(coe_class, "_i")
            inverse_fn <- get(inverse_fn_name, envir = asNamespace("Momocs2"))

            # Apply inverse to single coe
            shape <- inverse_fn(x$data[[coe_col]][j])[[1]]
            all_shapes <- c(all_shapes, list(shape))
          }
        }
      }
    }
  }

  # Template
  templated <- Momocs2::coo_template_relatively(all_shapes, size = 0.9)

  # Plot extent
  xlim <- c(-0.5, n_coe - 0.5)
  ylim <- if (facet) c(-0.6, n_groups - 0.4) else c(-0.6, 0.6)

  plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", asp = 1, axes = FALSE)

  # Draw
  shape_idx <- 1

  for (i in seq_len(n_groups)) {
    group_y <- if (facet) n_groups - i else 0

    # Draw contributing shapes FIRST (background)
    if (all) {
      n_contrib <- if (is.null(x$group_cols)) nrow(x$data) else
        sum(x$data[[x$group_cols]] == x$summary[[x$group_cols]][i])

      mean_end <- shape_idx + n_coe - 1
      contrib_start <- mean_end + 1

      for (k in seq_len(n_contrib)) {
        for (j in seq_len(n_coe)) {
          idx <- contrib_start + (k - 1) * n_coe + (j - 1)
          shape <- templated[[idx]]

          if (!is.matrix(shape) || ncol(shape) < 2) next

          x_offset <- j - 1
          shape_trans <- shape
          shape_trans[, 1] <- shape[, 1] + x_offset
          shape_trans[, 2] <- shape[, 2] + group_y

          if ("ldk" %in% class(x$data[[x$coe_cols[j]]])) {
            points(shape_trans[, 1], shape_trans[, 2], col = all_group_cols[i], pch = 20, cex = all_lwd)
          } else {
            lines(shape_trans[, 1], shape_trans[, 2], col = all_group_cols[i], lwd = all_lwd)
          }
        }
      }
    }

    # Draw mean shapes (foreground)
    for (j in seq_len(n_coe)) {
      shape <- templated[[shape_idx]]
      shape_idx <- shape_idx + 1

      if (!is.matrix(shape) || ncol(shape) < 2) next

      x_offset <- j - 1
      shape_trans <- shape
      shape_trans[, 1] <- shape[, 1] + x_offset
      shape_trans[, 2] <- shape[, 2] + group_y

      if ("ldk" %in% class(x$summary[[shape_cols[j]]])) {
        points(shape_trans[, 1], shape_trans[, 2], col = group_cols[i], pch = 20, cex = lwd)
      } else {
        lines(shape_trans[, 1], shape_trans[, 2], col = group_cols[i], lwd = lwd)
      }
    }

    # Skip contributing shapes indices
    if (all) {
      n_contrib <- if (is.null(x$group_cols)) nrow(x$data) else
        sum(x$data[[x$group_cols]] == x$summary[[x$group_cols]][i])
      shape_idx <- shape_idx + n_contrib * n_coe
    }

    # Labels (only if faceted)
    if (labels && facet && !is.null(x$group_cols)) {
      grp_label <- as.character(x$summary[[x$group_cols]][i])
      text((n_coe - 1) / 2, group_y - 0.5, labels = grp_label, cex = label_cex)
    }
  }

  # Legend (only if overlaid)
  if (!facet && labels && !is.null(x$group_cols)) {
    grp_labels <- as.character(x$summary[[x$group_cols]])
    legend("topright", legend = grp_labels, col = group_cols, lwd = lwd, bty = "n", cex = label_cex)
  }

  # Title
  if (is.null(x$group_cols)) {
    title(main = sprintf("Mean shape (%s)", x$.fn))
  } else {
    title(main = sprintf("Mean shapes by %s (%s)", x$group_cols, x$.fn))
  }

  invisible(NULL)
}
