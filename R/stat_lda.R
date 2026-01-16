# stat_lda.R

#' Linear Discriminant Analysis for morphometric data
#'
#' Perform LDA on coefficient data from morphometric analyses. Supports both
#' bare response names (auto-detects all coe columns) and full formula syntax.
#'
#' @param data A tibble with coefficient columns
#' @param formula A formula or bare column name specifying response and predictors:
#'   * Bare name: `species` → auto-detects all coe columns
#'   * Formula: `species ~ .` → all coe columns
#'   * Formula: `species ~ . + length` → all coe + covariate
#'   * Formula: `species ~ coe1` → specific coe column
#'   * Formula: `species ~ coe1 + coe2 + length` → multiple specific
#' @param ... Additional arguments passed to [MASS::lda()]
#'
#' @return An object of class `c("stat_lda", "momstats")` containing:
#' * `data`: Original tibble (unchanged)
#' * `model`: List containing both the LDA model and CV predictions
#' * `method`: "lda"
#' * `call`: The function call
#' * `formula`: Formula used
#' * `response_col`: Name of response column
#' * `predictor_cols`: All predictor column names (coe + others)
#' * `kept_predictor_cols`: Predictor columns after removing collinear
#' * `n_groups`: Number of groups
#' * `n_ld`: Number of discriminant functions
#'
#' @details
#' `stat_lda()` provides a unified interface for LDA on morphometric coefficient
#' data. It handles:
#'
#' * **Flexible formula syntax**: Bare response names or full formulas
#' * **Automatic coefficient detection**: Finds all `coe` columns via class
#' * **Coefficient unfolding**: Expands list-columns to individual harmonic columns
#' * **Multiple coefficient types**: Can combine different coefficient columns
#' * **Covariates**: Add non-morphometric predictors (size, environmental variables)
#' * **Cross-validation**: Always runs with CV=TRUE for reliable predictions
#'
#' ## Formula syntax
#'
#' The formula specifies response and predictors:
#' * `species` → Response only, auto-detect all coe columns
#' * `species ~ .` → All coe columns (explicit)
#' * `species ~ coe1` → Specific coe column
#' * `species ~ coe1 + coe2` → Multiple specific coe columns
#' * `species ~ . + length` → All coe plus a covariate
#'
#' ## Centering and scaling
#'
#' LDA internally centers and scales the data as part of the algorithm.
#' No manual preprocessing is needed.
#'
#' ## Cross-validation
#'
#' The function runs LDA twice:
#' 1. `CV = FALSE`: Builds the full model (discriminant functions, loadings)
#' 2. `CV = TRUE`: Generates cross-validated predictions
#'
#' Both are stored in `$model` for easy access.
#'
#' ## Getting results
#'
#' Use [collect()] to add predictions to your data:
#' ```r
#' lda <- boteft %>% stat_lda(type)
#' boteft_pred <- collect(lda)  # Adds pred, prob columns
#' ```
#'
#' @examples
#' \dontrun{
#' # Basic LDA - auto-detects all coe columns
#' lda1 <- boteft %>% stat_lda(type)
#'
#' # Explicit formula
#' lda2 <- boteft %>% stat_lda(type ~ .)
#'
#' # Specific coe column
#' lda3 <- boteft %>% stat_lda(type ~ coe)
#'
#' # With covariate
#' lda4 <- boteft %>% stat_lda(type ~ . + length)
#'
#' # Add predictions to data
#' boteft_pred <- collect(lda1)  # Adds pred, prob
#'
#' # Add LD scores too
#' boteft_ld <- collect(lda1, retain = 2)  # Adds pred, prob, LD1, LD2
#'
#' # Plot results
#' plot(lda1)  # Discrimination plot (default)
#' plot(lda1, type = "loadings")
#' plot(lda1, color = type)  # Color by true groups
#' plot(lda1, labels = pred)  # Label by predictions
#' }
#'
#' @seealso [MASS::lda()], [collect.stat_lda()], [plot.stat_lda()]
#'
#' @export
stat_lda <- function(data, formula, ...) {

  # Store call
  call <- match.call()

  # Check input
  if (!is.data.frame(data)) {
    stop("data must be a tibble or data frame")
  }

  # CRITICAL: Use substitute() for bare name support in pipes
  # enquo() doesn't work with bare names in pipes
  formula_expr <- substitute(formula)

  # Parse formula to get response and predictor columns
  formula_parsed <- parse_supervised_formula(data, formula_expr)

  response_col <- formula_parsed$response_col
  coe_cols <- formula_parsed$coe_cols
  covariate_cols <- formula_parsed$covariate_cols
  predictor_cols <- formula_parsed$predictor_cols

  # Get response variable
  response <- data[[response_col]]

  # Check response is factor or can be coerced
  if (!is.factor(response)) {
    response <- as.factor(response)
    message(sprintf("Converting response '%s' to factor", response_col))
  }

  # Build predictor matrix
  predictor_matrix <- build_predictor_matrix(data, coe_cols, covariate_cols)

  # Check for and remove constant/collinear columns
  # LDA doesn't like constant columns or perfect collinearity
  col_vars <- apply(predictor_matrix, 2, stats::var, na.rm = TRUE)
  constant_cols <- which(col_vars == 0 | is.na(col_vars))

  if (length(constant_cols) > 0) {
    message(sprintf("Removing %d constant column(s): %s",
                    length(constant_cols),
                    paste(colnames(predictor_matrix)[constant_cols], collapse = ", ")))
    predictor_matrix <- predictor_matrix[, -constant_cols, drop = FALSE]
  }

  # Check for perfect collinearity using QR decomposition
  qr_decomp <- qr(predictor_matrix)
  if (qr_decomp$rank < ncol(predictor_matrix)) {
    # Some columns are collinear
    n_collinear <- ncol(predictor_matrix) - qr_decomp$rank
    # Keep only independent columns
    keep_cols <- qr_decomp$pivot[seq_len(qr_decomp$rank)]
    removed_cols <- colnames(predictor_matrix)[-keep_cols]

    message(sprintf("Removing %d collinear column(s): %s",
                    length(removed_cols),
                    paste(removed_cols, collapse = ", ")))

    predictor_matrix <- predictor_matrix[, keep_cols, drop = FALSE]
  }

  # Store which columns were kept for prediction later
  kept_predictor_cols <- colnames(predictor_matrix)

  # Run LDA twice: once for model, once for CV predictions
  # Model (CV = FALSE) - gives discriminant functions, loadings, etc.
  lda_model <- MASS::lda(predictor_matrix, grouping = response, ...)

  # CV predictions (CV = TRUE) - gives cross-validated predictions
  lda_cv <- MASS::lda(predictor_matrix, grouping = response, CV = TRUE, ...)

  # Combine into single model object (Option B: merge but preserve structure)
  lda_model$cv_class <- lda_cv$class
  lda_model$cv_posterior <- lda_cv$posterior

  # Calculate number of discriminant functions
  n_groups <- length(lda_model$prior)
  n_ld <- min(n_groups - 1, ncol(predictor_matrix))

  # Build result object
  result <- list(
    data = data,  # Original, not augmented
    model = lda_model,
    method = "lda",
    call = call,
    formula = formula_expr,
    response_col = response_col,
    predictor_cols = predictor_cols,
    kept_predictor_cols = kept_predictor_cols,  # After removing constant/collinear
    n_groups = n_groups,
    n_ld = n_ld
  )

  class(result) <- c("stat_lda", "momstats")
  result
}


# Helper functions ----

#' Get coefficient columns from predictor columns
#'
#' @param data A tibble
#' @param predictor_cols Character vector of predictor column names
#'
#' @return Character vector of coe column names
#'
#' @keywords internal
#' @noRd
get_coe_from_predictors <- function(data, predictor_cols) {
  predictor_cols[
    sapply(data[predictor_cols], function(x) "coe" %in% class(x))
  ]
}


# Methods ----

#' Collect predictions and LD scores from LDA results
#'
#' Extract predictions, posterior probabilities, and optionally LD scores from
#' an LDA result and add them to a tibble.
#'
#' @param x A `stat_lda` object
#' @param data A tibble. If NULL, uses the original data from the LDA.
#' @param retain How many LDs to retain:
#'   * `FALSE` (default): No LD scores, only predictions
#'   * `TRUE`: All LDs
#'   * Integer (e.g., `2`): First N LDs
#'   * Numeric 0-1 (e.g., `0.95`): LDs explaining this proportion of variance
#' @param fold How to add LD scores (only if retain != FALSE):
#'   * `FALSE` (default): Add as separate columns (`LD1`, `LD2`, ...)
#'   * `TRUE`: Fold into single list-column named `"lda"`
#'   * Character: Fold into single list-column with this name
#' @param ... Additional arguments (reserved)
#'
#' @return A tibble with predictions added:
#' * `pred`: Predicted class (factor)
#' * `prob`: Posterior probability of predicted class (numeric)
#' * `LD1`, `LD2`, ... : LD scores (if `retain` is not FALSE)
#'
#' @details
#' The function always adds cross-validated predictions (`pred`) and the
#' posterior probability of the predicted class (`prob`). These come from
#' leave-one-out cross-validation performed during model fitting.
#'
#' Optionally, discriminant scores can be added via the `retain` parameter:
#' * `retain = FALSE`: No LD scores (default, fastest)
#' * `retain = TRUE`: All LD scores
#' * `retain = 2`: First 2 LDs
#' * `retain = 0.95`: LDs explaining 95% of between-group variance
#'
#' When `fold = TRUE` or a character name, LD scores are stored as a list-column
#' with class `c("lda", "coe")`, making them usable in downstream analyses.
#'
#' @examples
#' \dontrun{
#' lda <- boteft %>% stat_lda(type)
#'
#' # Add predictions only (default)
#' collect(lda)
#'
#' # Add predictions + all LD scores
#' collect(lda, retain = TRUE)
#'
#' # Add predictions + first 2 LDs
#' collect(lda, retain = 2)
#'
#' # Add predictions + LDs explaining 95% variance
#' collect(lda, retain = 0.95)
#'
#' # Fold LD scores into list-column
#' collect(lda, retain = 2, fold = TRUE)
#' collect(lda, retain = 2, fold = "lda_scores")
#' }
#'
#' @export
collect.stat_lda <- function(x, data = NULL, retain = FALSE, fold = FALSE, ...) {

  # Use original data if not provided
  if (is.null(data)) {
    data <- x$data
  }

  # Add CV predictions and probabilities
  data$pred <- x$model$cv_class
  data$prob <- apply(x$model$cv_posterior, 1, max)  # Prob of predicted class

  # If retain is FALSE, we're done
  if (isFALSE(retain)) {
    return(data)
  }

  # Otherwise, compute LD scores and add them
  # Predict on the same data to get LD scores
  coe_cols <- get_coe_from_predictors(x$data, x$predictor_cols)
  other_cols <- setdiff(x$predictor_cols, coe_cols)

  pred_matrix <- build_predictor_matrix(data, coe_cols, other_cols)

  # Keep only columns that were used in training (after removing collinear)
  pred_matrix <- pred_matrix[, x$kept_predictor_cols, drop = FALSE]

  # Get LD scores using predict.lda
  ld_scores <- predict(x$model, newdata = pred_matrix)$x

  # Determine which LDs to retain
  if (isTRUE(retain)) {
    # Keep all
    lds_to_keep <- seq_len(ncol(ld_scores))
  } else if (retain >= 1) {
    # Integer: first N LDs
    lds_to_keep <- seq_len(min(retain, ncol(ld_scores)))
  } else if (retain > 0 && retain < 1) {
    # Proportion: LDs explaining this much variance
    svd_vals <- x$model$svd
    prop_var <- svd_vals^2 / sum(svd_vals^2)
    cumvar <- cumsum(prop_var)
    lds_to_keep <- seq_len(which(cumvar >= retain)[1])
  } else {
    stop("retain must be FALSE, TRUE, an integer >= 1, or a proportion between 0 and 1")
  }

  # Extract selected LDs
  selected_scores <- ld_scores[, lds_to_keep, drop = FALSE]

  # Add to data
  if (isFALSE(fold)) {
    # Unfold: add as separate columns
    for (i in seq_len(ncol(selected_scores))) {
      col_name <- paste0("LD", lds_to_keep[i])
      data[[col_name]] <- selected_scores[, i]
    }
  } else {
    # Fold: add as list-column
    if (isTRUE(fold)) {
      col_name <- "lda"
    } else {
      col_name <- as.character(fold)
    }

    # Check if column exists
    if (col_name %in% names(data)) {
      stop(sprintf("Column '%s' already exists", col_name))
    }

    # Create list-column
    ld_list <- lapply(seq_len(nrow(selected_scores)), function(i) {
      vec <- selected_scores[i, ]
      names(vec) <- paste0("LD", lds_to_keep)
      class(vec) <- c("lda", "coe", "numeric")
      vec
    })

    class(ld_list) <- c("lda", "coe", "list")
    data[[col_name]] <- ld_list
  }

  data
}


#' @export
print.stat_lda <- function(x, ...) {
  cat("Linear Discriminant Analysis\n")
  cat("---------------------------------------\n")

  # Identify coe columns from predictor_cols
  coe_cols <- get_coe_from_predictors(x$data, x$predictor_cols)
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

  # Response info
  cat(sprintf("* Response: '%s' (%d groups)\n", x$response_col, x$n_groups))

  # Group names and sizes
  group_info <- sprintf("%s (n=%d)", names(x$model$prior), x$model$counts)
  cat(sprintf("* Groups: %s\n", paste(group_info, collapse = ", ")))
  cat(sprintf("* %d discriminant function%s\n\n", x$n_ld, if (x$n_ld > 1) "s" else ""))

  # Proportion of trace (variance explained by each LD)
  if (x$n_ld > 0) {
    svd_vals <- x$model$svd
    prop_var <- svd_vals^2 / sum(svd_vals^2)
    cumvar <- cumsum(prop_var)

    cat("Proportion of trace:\n")
    n_show <- min(5, x$n_ld)
    var_df <- data.frame(
      LD = paste0("LD", seq_len(n_show)),
      Proportion = round(prop_var[seq_len(n_show)] * 100, 2),
      Cumulative = round(cumvar[seq_len(n_show)] * 100, 2)
    )
    colnames(var_df) <- c("", "Prop (%)", "Cum (%)")
    print(var_df, row.names = FALSE)
  }

  # CV accuracy
  cv_correct <- sum(x$model$cv_class == x$data[[x$response_col]])
  cv_accuracy <- cv_correct / nrow(x$data) * 100

  cat(sprintf("\n* CV accuracy: %.1f%% (%d/%d correct)\n",
              cv_accuracy, cv_correct, nrow(x$data)))

  cat("* calculated using MASS::lda with CV=TRUE ($model)\n")
  cat("* collect() to add predictions to data\n")
  cat("* plot() to visualize results (?plot.stat_lda)\n")

  invisible(x)
}


#' @export
summary.stat_lda <- function(object, ...) {
  cat("Linear Discriminant Analysis\n")
  cat(sprintf("\n%d observations, %d variables, %d groups\n",
              nrow(object$data), length(object$predictor_cols), object$n_groups))

  cat("\nGroup prior probabilities:\n")
  print(round(object$model$prior, 4))

  cat("\nGroup means (first 5 predictors):\n")
  n_show <- min(5, ncol(object$model$means))
  print(round(object$model$means[, seq_len(n_show)], 4))

  if (object$n_ld > 0) {
    cat("\nCoefficients of linear discriminants (first 5 predictors):\n")
    n_pred_show <- min(5, nrow(object$model$scaling))
    n_ld_show <- min(3, ncol(object$model$scaling))
    print(round(object$model$scaling[seq_len(n_pred_show), seq_len(n_ld_show), drop = FALSE], 4))

    cat("\nProportion of trace:\n")
    svd_vals <- object$model$svd
    prop_var <- svd_vals^2 / sum(svd_vals^2)
    names(prop_var) <- paste0("LD", seq_along(prop_var))
    print(round(prop_var, 4))
  }

  invisible(object)
}


# Plot method ----

#' Plot LDA results
#'
#' Visualize LDA results with discrimination plots, loading plots, or histograms
#' (for binary classification).
#'
#' @param x A `stat_lda` object
#' @param type Character. Type of plot:
#'   * `"scores"`: LD score plot (default) - scatter for 2+ LDs, boxplot for 1 LD
#'   * `"loadings"`: Loading plot for first two LDs
#' @param lds Integer vector of length 2. Which LDs to plot for score/loading plots.
#'   Default is `c(1, 2)`. For binary classification (1 LD), only first element used.
#' @param color Column name (bare or quoted) for coloring points. Default uses
#'   the response variable.
#' @param labels Column name (bare or quoted) for text labels in score plot.
#'   For loadings, logical: should variable names be shown? Default `TRUE`.
#' @param chull Logical. Draw convex hulls around groups? Only works when `color`
#'   is a factor. Default `TRUE`.
#' @param cex Numeric. Character expansion for labels. Default is 0.7.
#' @param legend Logical. Show legend for colors? Default `TRUE`.
#' @param ... Additional arguments (reserved for future use)
#'
#' @return NULL (invisibly). Draws plot as side effect.
#'
#' @examples
#' \dontrun{
#' lda <- boteft %>% stat_lda(type)
#'
#' # Score plot (default)
#' plot(lda)
#' plot(lda, color = type)
#' plot(lda, labels = pred)
#' plot(lda, color = type, chull = TRUE)
#'
#' # Different LDs
#' plot(lda, lds = c(2, 3))
#'
#' # Loading plot
#' plot(lda, type = "loadings")
#' plot(lda, type = "loadings", labels = FALSE)
#' }
#'
#' @export
plot.stat_lda <- function(x, type = c("scores", "loadings"),
                          lds = c(1, 2), color = NULL, labels = NULL,
                          chull = TRUE, cex = 0.7, legend = TRUE, ...) {

  type <- match.arg(type)

  # Capture color argument properly
  color_expr <- substitute(color)

  # Default color to response variable if NULL
  if (is.null(color_expr) || identical(color_expr, quote(NULL))) {
    color_quo <- rlang::sym(x$response_col)
  } else {
    color_quo <- rlang::enquo(color)
  }

  if (type == "scores") {
    plot_lda_scores(x, lds, color_quo, rlang::enquo(labels), chull, cex, legend)
  } else if (type == "loadings") {
    # Default labels to TRUE for loadings
    if (is.null(labels)) labels <- TRUE
    plot_lda_loadings(x, lds, labels, cex)
  }

  # return the object to continue piping
  invisible(x)
}


#' @keywords internal
#' @noRd
plot_lda_scores <- function(x, lds, color_quo, labels_quo, chull, cex, legend) {

  # Get LD scores
  # Predict on same data to get LD scores
  coe_cols <- get_coe_from_predictors(x$data, x$predictor_cols)
  other_cols <- setdiff(x$predictor_cols, coe_cols)

  pred_matrix <- build_predictor_matrix(x$data, coe_cols, other_cols)

  # Keep only columns that were used in training (after removing collinear)
  pred_matrix <- pred_matrix[, x$kept_predictor_cols, drop = FALSE]

  ld_scores <- predict(x$model, newdata = pred_matrix)$x

  # Get color values
  color_vals <- get_plot_variable(x$data, color_quo)

  # Get label values
  label_vals <- get_plot_variable(x$data, labels_quo)

  # Binary classification (1 LD) - use boxplot
  if (x$n_ld == 1) {
    plot_lda_scores_1d(ld_scores[, 1], color_vals, label_vals, x, cex, legend)
    return(invisible(x))
  }

  # Multi-class (2+ LDs) - use scatter plot
  ld1 <- ld_scores[, lds[1]]
  ld2 <- ld_scores[, lds[2]]

  # Calculate proportion of trace for axes
  svd_vals <- x$model$svd
  prop_var <- svd_vals^2 / sum(svd_vals^2)
  var1 <- round(prop_var[lds[1]] * 100, 1)
  var2 <- round(prop_var[lds[2]] * 100, 1)

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
  plot(ld1, ld2,
       xlab = sprintf("LD%d (%s%%)", lds[1], var1),
       ylab = sprintf("LD%d (%s%%)", lds[2], var2),
       main = "LDA Score Plot",
       type = "n",  # Don't plot points yet
       asp = 1,
       las = 1)

  abline(h = 0, v = 0, col = "gray", lty = 2)

  # Draw convex hulls if requested (behind points)
  if (chull && !is.null(color_vals) &&
      (is.factor(color_vals) || is.character(color_vals))) {
    draw_chulls(ld1, ld2, color_vals, rainbow(nlevels(as.factor(color_vals)), alpha = 1))
  }

  # Now add points
  points(ld1, ld2, pch = 19, col = cols)

  # Add labels if requested
  if (!is.null(label_vals)) {
    text(ld1, ld2, labels = label_vals, cex = cex, pos = 3)
  }

  # Add legend if requested and we have colors
  if (legend && !is.null(color_vals)) {
    add_plot_legend(color_vals, cols, position = "topright")
  }

  # pass the object so that we can continue the pipe
  invisible(x)
}


#' @keywords internal
#' @noRd
plot_lda_scores_1d <- function(ld_scores, color_vals, label_vals, x, cex, legend) {

  # Binary classification - boxplot with jittered points
  if (is.null(color_vals)) {
    color_vals <- x$data[[x$response_col]]
  }

  color_vals <- as.factor(color_vals)
  lvls <- levels(color_vals)

  # Create boxplot
  boxplot(ld_scores ~ color_vals,
          xlab = "Group",
          ylab = "LD1",
          main = "LDA Discrimination (Binary)",
          col = rainbow(length(lvls), alpha = 0.3),
          las = 1)

  abline(h = 0, col = "gray", lty = 2)

  # Add jittered points
  for (i in seq_along(lvls)) {
    grp <- lvls[i]
    idx <- which(color_vals == grp)
    y_vals <- ld_scores[idx]
    x_vals <- jitter(rep(i, length(y_vals)), amount = 0.1)

    points(x_vals, y_vals, pch = 19,
           col = rainbow(length(lvls), alpha = 0.5)[i])

    # Add labels if requested
    if (!is.null(label_vals)) {
      text(x_vals, y_vals, labels = label_vals[idx], cex = cex, pos = 3)
    }
  }
}


#' @keywords internal
#' @noRd
plot_lda_loadings <- function(x, lds, labels, cex) {

  # Binary classification - only 1 LD
  if (x$n_ld == 1) {
    # Show loadings as barplot
    loadings <- x$model$scaling[, 1]

    # Show top 20 loadings by magnitude
    top_n <- min(20, length(loadings))
    top_idx <- order(abs(loadings), decreasing = TRUE)[seq_len(top_n)]
    top_loadings <- loadings[top_idx]

    barplot(top_loadings,
            main = "LD1 Loadings (Top 20)",
            ylab = "Loading",
            las = 2,
            cex.names = 0.6,
            col = ifelse(top_loadings > 0, "steelblue", "coral"))
    abline(h = 0)

    return(invisible(x))
  }

  # Multi-class - loading plot for 2 LDs
  loadings <- x$model$scaling[, lds]

  # Create empty plot (no points)
  plot(loadings[, 1], loadings[, 2],
       xlab = sprintf("LD%d loading", lds[1]),
       ylab = sprintf("LD%d loading", lds[2]),
       main = "LDA Loading Plot",
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
  # pass the object so that we can continue the pipe
  invisible(x)
}


# Plot helpers (shared with stat_pca) ----

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


# Predict method ----

#' Predict method for LDA
#'
#' Predict class membership and discriminant scores for new data using a fitted LDA model.
#'
#' @param object A `stat_lda` object
#' @param newdata A tibble with the same predictor columns as training data
#' @param retain How many LDs to return:
#'   * `FALSE` (default): No LD scores, only predictions
#'   * `TRUE`: All LDs
#'   * Integer (e.g., `2`): First N LDs
#'   * Numeric 0-1 (e.g., `0.95`): LDs explaining this proportion of variance
#' @param fold How to return LD scores (only if retain != FALSE):
#'   * `FALSE` (default): Add as separate columns (`LD1`, `LD2`, ...)
#'   * `TRUE`: Fold into single list-column named `"lda"`
#'   * Character: Fold into single list-column with this name
#' @param .collect Logical. Should predictions be added to `newdata` (TRUE, default)
#'   or returned as a standalone tibble (FALSE)?
#' @param ... Additional arguments (reserved)
#'
#' @return If `.collect = TRUE`, returns `newdata` with predictions added. If
#'   `.collect = FALSE`, returns a tibble with identifier columns (non-predictors) and
#'   predictions only. Always includes:
#' * `pred`: Predicted class (factor)
#' * `prob`: Posterior probability of predicted class (numeric)
#' * `LD1`, `LD2`, ... : LD scores (if `retain` is not FALSE)
#'
#' @examples
#' \dontrun{
#' # Train LDA
#' lda <- boteft %>% stat_lda(type)
#'
#' # Predict on new data
#' new_preds <- predict(lda, new_data)
#'
#' # Pipe the model
#' new_preds <- training %>%
#'   stat_lda(type) %>%
#'   predict(testing)
#' }
#'
#' @export
predict.stat_lda <- function(object, newdata, retain = FALSE, fold = FALSE,
                             .collect = TRUE, ...) {

  if (!is.data.frame(newdata)) {
    stop("newdata must be a tibble or data frame")
  }

  # Identify coe columns from predictor_cols
  coe_cols <- get_coe_from_predictors(object$data, object$predictor_cols)
  other_cols <- setdiff(object$predictor_cols, coe_cols)

  # Check all required columns present
  missing_cols <- setdiff(object$predictor_cols, names(newdata))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in newdata: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Build predictor matrix
  pred_matrix <- build_predictor_matrix(newdata, coe_cols, other_cols)

  # Keep only columns used in training
  if (!all(object$kept_predictor_cols %in% colnames(pred_matrix))) {
    stop("Predictor structure doesn't match training data")
  }

  pred_matrix <- pred_matrix[, object$kept_predictor_cols, drop = FALSE]

  # Predict using predict.lda
  lda_pred <- predict(object$model, newdata = pred_matrix)

  # Prepare output
  if (.collect) {
    result <- newdata
  } else {
    non_pred_cols <- setdiff(names(newdata), object$predictor_cols)
    if (length(non_pred_cols) > 0) {
      result <- newdata[, non_pred_cols, drop = FALSE]
    } else {
      result <- tibble::tibble(.rows = nrow(newdata))
    }
  }

  # Add predictions
  result$pred <- lda_pred$class
  result$prob <- apply(lda_pred$posterior, 1, max)

  # If retain is FALSE, we're done
  if (isFALSE(retain)) {
    return(result)
  }

  # Otherwise, add LD scores
  ld_scores <- lda_pred$x

  # Determine which LDs to retain
  if (isTRUE(retain)) {
    lds_to_keep <- seq_len(ncol(ld_scores))
  } else if (retain >= 1) {
    lds_to_keep <- seq_len(min(retain, ncol(ld_scores)))
  } else if (retain > 0 && retain < 1) {
    svd_vals <- object$model$svd
    prop_var <- svd_vals^2 / sum(svd_vals^2)
    cumvar <- cumsum(prop_var)
    lds_to_keep <- seq_len(which(cumvar >= retain)[1])
  } else {
    stop("retain must be FALSE, TRUE, an integer >= 1, or a proportion between 0 and 1")
  }

  selected_scores <- ld_scores[, lds_to_keep, drop = FALSE]

  # Add LD scores
  if (isFALSE(fold)) {
    for (i in seq_len(ncol(selected_scores))) {
      col_name <- paste0("LD", lds_to_keep[i])
      result[[col_name]] <- selected_scores[, i]
    }
  } else {
    col_name <- if (isTRUE(fold)) "lda" else as.character(fold)

    if (col_name %in% names(result)) {
      stop(sprintf("Column '%s' already exists", col_name))
    }

    ld_list <- lapply(seq_len(nrow(selected_scores)), function(i) {
      vec <- selected_scores[i, ]
      names(vec) <- paste0("LD", lds_to_keep)
      class(vec) <- c("lda", "coe", "numeric")
      vec
    })

    class(ld_list) <- c("lda", "coe", "list")
    result[[col_name]] <- ld_list
  }

  result
}


#' @rdname transduce
#' @export
transduce.stat_lda <- function(object, positions) {

  if (!is.data.frame(positions)) {
    stop("positions must be a tibble or data frame")
  }

  if (nrow(positions) == 0) {
    stop("positions tibble is empty")
  }

  # Validate column names are LDs
  pos_names <- names(positions)

  if (!all(grepl("^LD[0-9]+$", pos_names))) {
    stop("positions columns must be LD axes (e.g., LD1, LD2, ...)")
  }

  # Extract LD numbers and validate they exist
  ld_nums <- as.integer(sub("LD", "", pos_names))
  max_ld <- ncol(object$model$scaling)

  if (any(ld_nums < 1 | ld_nums > max_ld)) {
    stop(sprintf("LD numbers must be between 1 and %d", max_ld))
  }

  # Identify coe columns from predictor_cols by checking classes in data
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

  # Get LDA components
  scaling <- object$model$scaling  # p × k matrix
  grand_mean <- colMeans(object$model$means)  # Grand mean in predictor space

  # Get the scaling matrix for the specified LDs only
  scaling_subset <- scaling[, ld_nums, drop = FALSE]

  # Compute pseudoinverse of the (subset) scaling matrix
  scaling_pinv <- MASS::ginv(scaling_subset)  # k × p

  # For each position (row in positions tibble)
  for (i in seq_len(n_positions)) {

    # Get LD scores for this position as a row vector
    ld_scores <- as.numeric(positions[i, ])

    # Reconstruct: X = Z %*% W_pinv + grand_mean
    # ld_scores is 1 × k, scaling_pinv is k × p
    reconstructed <- as.vector(ld_scores %*% scaling_pinv) + grand_mean

    # Split reconstructed predictors back into coe columns
    col_start <- 1
    for (j in seq_along(coe_cols)) {
      # Get number of coefficients for this coe column
      n_coe <- nrow(object$model$scaling)
      col_end <- n_coe

      # Extract coefficients
      coe_vec <- reconstructed[col_start:col_end]

      # Get harmonic names from kept_predictor_cols
      names(coe_vec) <- object$kept_predictor_cols

      # Add classes
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
