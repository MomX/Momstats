
# Morphospace helpers ----

#' Morphospace visualization: shapes along axes
#'
#' Create shapes along specified axes.
#'
#' @param object A stat_* object (stat_pca, stat_lda, etc.)
#' @param pcs Integer vector of length 2. Which axes to use. Default c(1, 2).
#' @param n Integer. Number of shapes per axis direction. Default 1.
#' @param extent Numeric. Extent factor (< 1 = zoom in, > 1 = extrapolate). Default 1.0.
#' @param transduce Logical. Should shapes be transduced? Default TRUE.
#' @param template Logical or numeric. Templating option. Default TRUE.
#' @param draw Logical. Should shapes be drawn? Default TRUE.
#' @param ... Additional graphical parameters (lwd, col, etc.)
#'
#' @return A tibble with positions, coefficients, and shapes
#'
#' @details
#' Creates shapes along each axis. For n=1, creates 4 shapes at extremes.
#' For n>1, creates more intermediate shapes along each axis.
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#' plot(pca, morphospace = "axes")
#'
#' # More shapes
#' plot(pca, morphospace = FALSE)
#' morphospace_axes(pca, n = 3)
#' }
#'
#' @export
morphospace_axes <- function(object, pcs = c(1, 2), n = 1, extent = 1.0,
                             transduce = TRUE, template = TRUE,
                             draw = TRUE, ...) {

  axis_names <- get_axis_names(object, pcs)
  pc_scores <- get_pc_scores(object, pcs)
  pc1_range <- range(pc_scores[, 1])
  pc2_range <- range(pc_scores[, 2])

  pc1_center <- mean(pc1_range)
  pc2_center <- mean(pc2_range)

  pc1_half_width <- diff(pc1_range) / 2 * extent
  pc2_half_width <- diff(pc2_range) / 2 * extent

  pc1_limits <- c(pc1_center - pc1_half_width, pc1_center + pc1_half_width)
  pc2_limits <- c(pc2_center - pc2_half_width, pc2_center + pc2_half_width)

  if (n == 1) {
    positions <- tibble::tibble(
      !!axis_names[1] := c(pc1_limits[1], pc1_limits[2], 0, 0),
      !!axis_names[2] := c(0, 0, pc2_limits[1], pc2_limits[2])
    )
  } else {
    pc1_seq <- seq(pc1_limits[1], pc1_limits[2], length.out = 2*n + 1)
    pc2_seq <- seq(pc2_limits[1], pc2_limits[2], length.out = 2*n + 1)

    pc1_axis <- tibble::tibble(
      !!axis_names[1] := pc1_seq,
      !!axis_names[2] := 0
    )

    pc2_seq_no_origin <- pc2_seq[pc2_seq != 0]
    pc2_axis <- tibble::tibble(
      !!axis_names[1] := 0,
      !!axis_names[2] := pc2_seq_no_origin
    )

    positions <- rbind(pc1_axis, pc2_axis)
  }

  if (transduce) {
    positions <- transduce(object, positions)
  }

  # Template and translate even when draw=FALSE (only skip if template=FALSE)
  if (!isFALSE(template) && transduce) {
    positions <- template_morphospace(object, positions, pcs,
                                      template_size = template,
                                      grid_info = NULL)
  }

  if (draw && transduce) {
    draw_morphospace(positions,  ...)
  }

  if (draw) invisible(positions) else positions
}


#' Morphospace visualization: uniform grid
#'
#' Create shapes on a uniform grid across the full PC plane.
#'
#' @param object A stat_* object (stat_pca, stat_lda, etc.)
#' @param pcs Integer vector of length 2. Which axes to use. Default c(1, 2).
#' @param nrow Integer. Number of rows in grid. Default 5.
#' @param ncol Integer. Number of columns in grid. Default 5.
#' @param extent Numeric. Extent factor. Default 1.0.
#' @param transduce Logical. Should shapes be transduced? Default TRUE.
#' @param template Logical or numeric. Templating option. Default TRUE.
#' @param draw Logical. Should shapes be drawn? Default TRUE.
#' @param ... Additional graphical parameters (lwd, col, etc.)
#'
#' @return A tibble with positions, coefficients, and shapes
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#' plot(pca, morphospace = "grid", nrow = 3, ncol = 3)
#' }
#'
#' @export
morphospace_grid <- function(object, pcs = c(1, 2), nrow = 5, ncol = 5,
                             extent = 1.0, transduce = TRUE, template = TRUE,
                             draw = TRUE, ...) {

  axis_names <- get_axis_names(object, pcs)
  pc_scores <- get_pc_scores(object, pcs)
  pc1_range <- range(pc_scores[, 1])
  pc2_range <- range(pc_scores[, 2])

  pc1_center <- mean(pc1_range)
  pc2_center <- mean(pc2_range)

  pc1_half_width <- diff(pc1_range) / 2 * extent
  pc2_half_width <- diff(pc2_range) / 2 * extent

  pc1_limits <- c(pc1_center - pc1_half_width, pc1_center + pc1_half_width)
  pc2_limits <- c(pc2_center - pc2_half_width, pc2_center + pc2_half_width)

  pc1_seq <- seq(pc1_limits[1], pc1_limits[2], length.out = ncol)
  pc2_seq <- seq(pc2_limits[1], pc2_limits[2], length.out = nrow)

  # Create grid with base R
  grid_df <- expand.grid(
    pc1 = pc1_seq,
    pc2 = pc2_seq
  )

  positions <- tibble::tibble(
    !!axis_names[1] := grid_df$pc1,
    !!axis_names[2] := grid_df$pc2
  )

  if (transduce) {
    positions <- transduce(object, positions)
  }

  if (!isFALSE(template) && transduce) {
    grid_info <- list(
      pc1_limits = pc1_limits,
      pc2_limits = pc2_limits,
      nrow = nrow,
      ncol = ncol
    )

    positions <- template_morphospace(object, positions, pcs,
                                      template_size = template,
                                      grid_info = grid_info)
  }

  if (draw && transduce) {
    draw_morphospace(positions,  ...)
  }

  if (draw) invisible(positions) else positions
}


#' Morphospace visualization: data range grid
#'
#' Create shapes on a grid within the actual data bounding box.
#'
#' @param object A stat_* object (stat_pca, stat_lda, etc.)
#' @param pcs Integer vector of length 2. Which axes to use. Default c(1, 2).
#' @param nrow Integer. Number of rows in grid. Default 5.
#' @param ncol Integer. Number of columns in grid. Default 5.
#' @param extent Numeric. Extent factor (not used for grid placement). Default 1.0.
#' @param transduce Logical. Should shapes be transduced? Default TRUE.
#' @param template Logical or numeric. Templating option. Default TRUE.
#' @param draw Logical. Should shapes be drawn? Default TRUE.
#' @param ... Additional graphical parameters (lwd, col, etc.)
#'
#' @return A tibble with positions, coefficients, and shapes
#'
#' @details
#' Creates a grid within the actual data bounding box.
#' Grid always stays within data range regardless of extent.
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#' plot(pca, morphospace = "range")
#' }
#'
#' @export
morphospace_range <- function(object, pcs = c(1, 2), nrow = 5, ncol = 5,
                              extent = 1.0, transduce = TRUE, template = TRUE,
                              draw = TRUE, ...) {

  axis_names <- get_axis_names(object, pcs)
  pc_scores <- get_pc_scores(object, pcs)
  pc1_bbox <- range(pc_scores[, 1])
  pc2_bbox <- range(pc_scores[, 2])

  pc1_seq <- seq(pc1_bbox[1], pc1_bbox[2], length.out = ncol)
  pc2_seq <- seq(pc2_bbox[1], pc2_bbox[2], length.out = nrow)

  # Create grid with base R
  grid_df <- expand.grid(
    pc1 = pc1_seq,
    pc2 = pc2_seq
  )

  positions <- tibble::tibble(
    !!axis_names[1] := grid_df$pc1,
    !!axis_names[2] := grid_df$pc2
  )

  if (transduce) {
    positions <- transduce(object, positions)
  }

  if (!isFALSE(template) && transduce) {
    grid_info <- list(
      pc1_limits = pc1_bbox,
      pc2_limits = pc2_bbox,
      nrow = nrow,
      ncol = ncol
    )

    positions <- template_morphospace(object, positions, pcs,
                                      template_size = template,
                                      grid_info = grid_info)
  }

  if (draw && transduce) {
    draw_morphospace(positions,  ...)
  }

  if (draw) invisible(positions) else positions
}


#' Morphospace visualization: group centroids
#'
#' Create shapes at group centroid positions.
#'
#' @param object A stat_* object (stat_pca, stat_lda, etc.)
#' @param group Column name for grouping (bare or quoted). For PCA, required.
#'   For LDA, optional (uses discriminant groups if NULL).
#' @param group_vals Pre-evaluated group values (used internally by plot).
#' @param transduce Logical. Should shapes be transduced? Default TRUE.
#' @param template Logical or numeric. Templating option. Default TRUE.
#' @param draw Logical. Should shapes be drawn? Default TRUE.
#' @param ... Additional graphical parameters (lwd, col, etc.)
#'
#' @return A tibble with positions, coefficients, and shapes
#'
#' @examples
#' \dontrun{
#' pca <- boteft %>% stat_pca()
#' plot(pca, color = type, morphospace = "centroids")
#'
#' # Or call directly
#' plot(pca, morphospace = FALSE)
#' morphospace_centroids(pca, group = type)
#' }
#'
#' @export
morphospace_centroids <- function(object, group = NULL, group_vals = NULL,
                                  transduce = TRUE, template = TRUE,
                                  draw = TRUE, ...) {

  if (inherits(object, "stat_lda")) {
    positions <- get_lda_centroids(object)
    pcs <- c(1, 2)
  } else if (inherits(object, "stat_pca")) {
    # Accept either group (quosure) or group_vals (pre-evaluated)
    if (!is.null(group_vals)) {
      # Use pre-evaluated values (from plot call)
      group_data <- group_vals
    } else {
      # Use quosure (from direct call)
      group_quo <- rlang::enquo(group)
      if (rlang::quo_is_null(group_quo)) {
        stop("For PCA, group must be specified")
      }
      group_data <- get_plot_variable(object$data, group_quo)
      if (is.null(group_data)) {
        stop("Group variable not found")
      }
    }

    positions <- get_pca_centroids(object, group_data, pcs = c(1, 2))
    pcs <- c(1, 2)
  } else {
    stop("morphospace_centroids not implemented for this model type")
  }

  if (transduce) {
    positions <- transduce(object, positions)
  }

  if (!isFALSE(template) && transduce) {
    positions <- template_morphospace(object, positions, pcs,
                                      template_size = template,
                                      grid_info = NULL)
  }

  if (draw && transduce) {
    draw_morphospace(positions,  ...)
  }

  if (draw) invisible(positions) else positions
}


# Helpers ----

#' @keywords internal
#' @noRd
get_axis_names <- function(object, pcs) {
  if (inherits(object, "stat_pca")) {
    paste0("PC", pcs)
  } else if (inherits(object, "stat_lda")) {
    paste0("LD", pcs)
  } else {
    paste0("Axis", pcs)
  }
}

#' @keywords internal
#' @noRd
get_pc_scores <- function(object, pcs) {
  if (inherits(object, "stat_pca")) {
    object$model$x[, pcs, drop = FALSE]
  } else if (inherits(object, "stat_lda")) {
    predictor_matrix <- build_predictor_matrix(
      object$data,
      object$predictor_cols[sapply(object$data[object$predictor_cols],
                                   function(x) "coe" %in% class(x))],
      character(0))
    centered <- scale(predictor_matrix,
                      center = colMeans(object$model$means),
                      scale = FALSE)
    ld_scores <- centered %*% object$model$scaling
    ld_scores[, pcs, drop = FALSE]
  } else {
    stop("get_pc_scores not implemented")
  }
}

#' @keywords internal
#' @noRd
get_lda_centroids <- function(object) {
  group_means <- object$model$means
  grand_mean <- colMeans(group_means)
  centered_means <- scale(group_means, center = grand_mean, scale = FALSE)
  ld_positions <- centered_means %*% object$model$scaling

  tibble::tibble(
    group = rownames(group_means),
    LD1 = ld_positions[, 1],
    LD2 = if (ncol(ld_positions) > 1) ld_positions[, 2] else 0
  )
}

#' @keywords internal
#' @noRd
get_pca_centroids <- function(object, group_vals, pcs) {
  pc_scores <- object$model$x[, pcs, drop = FALSE]
  group_factor <- as.factor(group_vals)

  centroids <- lapply(levels(group_factor), function(lvl) {
    idx <- which(group_factor == lvl)
    colMeans(pc_scores[idx, , drop = FALSE])
  })

  centroids_matrix <- do.call(rbind, centroids)

  # Create tibble with ONLY PC columns (no group column)
  result <- tibble::tibble(
    !!paste0("PC", pcs[1]) := centroids_matrix[, 1],
    !!paste0("PC", pcs[2]) := centroids_matrix[, 2]
  )

  # Store group names as attribute for later use if needed
  attr(result, "groups") <- levels(group_factor)

  result
}

#' @keywords internal
#' @noRd
template_morphospace <- function(object, positions, pcs, template_size, grid_info) {
  axis_names <- get_axis_names(object, pcs)
  coe_cols <- names(positions)[sapply(positions, function(x) "coe" %in% class(x))]

  if (length(coe_cols) == 0) {
    warning("No coe columns for templating")
    return(positions)
  }

  n_coe <- length(coe_cols)

  if (isTRUE(template_size)) {
    if (!is.null(grid_info)) {
      grid_width <- diff(grid_info$pc1_limits)
      grid_height <- diff(grid_info$pc2_limits)
      cell_width <- grid_width / (grid_info$ncol - 1)
      cell_height <- grid_height / (grid_info$nrow - 1)
      base_size <- min(cell_width, cell_height) * 0.4
    } else {
      pc_scores <- get_pc_scores(object, pcs)
      data_width <- diff(range(pc_scores[, 1]))
      data_height <- diff(range(pc_scores[, 2]))
      base_size <- min(data_width, data_height) * 0.08
    }
    size_factor <- base_size
  } else if (is.numeric(template_size)) {
    if (!is.null(grid_info)) {
      grid_width <- diff(grid_info$pc1_limits)
      grid_height <- diff(grid_info$pc2_limits)
      cell_width <- grid_width / (grid_info$ncol - 1)
      cell_height <- grid_height / (grid_info$nrow - 1)
      base_size <- min(cell_width, cell_height) * 0.4
    } else {
      pc_scores <- get_pc_scores(object, pcs)
      data_width <- diff(range(pc_scores[, 1]))
      data_height <- diff(range(pc_scores[, 2]))
      base_size <- min(data_width, data_height) * 0.08
    }
    size_factor <- base_size * template_size
  } else {
    return(positions)
  }

  x_pos <- positions[[axis_names[1]]]
  y_pos <- positions[[axis_names[2]]]

  view_offsets <- calculate_view_offsets(n_coe, size_factor)

  for (j in seq_along(coe_cols)) {
    coe_col <- coe_cols[j]
    shape_col <- paste0(coe_col, "_i")

    if (!shape_col %in% names(positions)) next

    offset_x <- view_offsets$x[j]
    offset_y <- view_offsets$y[j]

    for (i in seq_len(nrow(positions))) {
      shape <- positions[[shape_col]][[i]]
      if (!is.matrix(shape) || ncol(shape) < 2) next

      shape_centered <- scale(shape, center = TRUE, scale = FALSE)
      shape_scaled <- shape_centered * size_factor
      shape_scaled[, 1] <- shape_scaled[, 1] + x_pos[i] + offset_x
      shape_scaled[, 2] <- shape_scaled[, 2] + y_pos[i] + offset_y

      positions[[shape_col]][[i]] <- shape_scaled
    }
  }

  positions
}

#' @keywords internal
#' @noRd
calculate_view_offsets <- function(n_views, size_factor) {
  if (n_views == 1) {
    return(list(x = 0, y = 0))
  }
  if (n_views == 2) {
    spacing <- size_factor * 1.2
    return(list(x = c(-spacing/2, spacing/2), y = c(0, 0)))
  }
  if (n_views == 3) {
    spacing <- size_factor * 1.2
    return(list(
      x = c(0, -spacing/2, spacing/2),
      y = c(spacing * 0.6, -spacing * 0.3, -spacing * 0.3)
    ))
  }
  if (n_views == 4) {
    spacing <- size_factor * 1.2
    return(list(
      x = c(-spacing/2, spacing/2, -spacing/2, spacing/2),
      y = c(spacing/2, spacing/2, -spacing/2, -spacing/2)
    ))
  }
  angles <- seq(0, 2 * pi, length.out = n_views + 1)[1:n_views]
  radius <- size_factor * 1.2
  list(x = cos(angles) * radius, y = sin(angles) * radius)
}


#' @keywords internal
#' @noRd
draw_morphospace <- function(positions, lwd = 1, col = "black", ...) {
  shape_cols <- names(positions)[grepl("_i$", names(positions))]
  if (length(shape_cols) == 0) {
    warning("No shape columns for drawing")
    return(invisible(NULL))
  }

  for (shape_col in shape_cols) {
    shapes <- positions[[shape_col]]
    for (i in seq_along(shapes)) {
      shape <- shapes[[i]]
      if (!is.matrix(shape) || ncol(shape) < 2) next
      lines(shape[, 1], shape[, 2], lwd = lwd, col = col, ...)
    }
  }

  invisible(NULL)
}

#' @keywords internal
#' @noRd
add_plot_legend <- function(color_vals, cols, position = "topright") {
  if (is.null(color_vals)) return(invisible())

  if (is.factor(color_vals) || is.character(color_vals)) {
    color_vals <- as.factor(color_vals)
    lvls <- levels(color_vals)
    unique_cols <- rainbow(nlevels(color_vals), alpha = 0.5)
    legend(position, legend = lvls, pch = 19, col = unique_cols,
           bty = "n", cex = 0.8)
  } else if (is.numeric(color_vals)) {
    range_vals <- range(color_vals, na.rm = TRUE)
    legend_text <- c(sprintf("%.2f", range_vals[2]), "",
                     sprintf("%.2f", range_vals[1]))
    col_ramp <- colorRampPalette(c("blue", "red"))
    legend_cols <- col_ramp(3)
    legend(position, legend = legend_text, pch = 15, col = legend_cols,
           bty = "n", cex = 0.8, title = "Value")
  }
}

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
    if (length(idx) < 3) next

    x_grp <- x[idx]
    y_grp <- y[idx]
    hull_idx <- grDevices::chull(x_grp, y_grp)
    hull_idx <- c(hull_idx, hull_idx[1])
    lines(x_grp[hull_idx], y_grp[hull_idx], col = cols[i], lwd = 2)
  }
}
