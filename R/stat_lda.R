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
#' * `model`: List containing both the LDA model and CV predictions:
#'   - All components from [MASS::lda()] (prior, counts, means, scaling, svd, etc.)
#'   - `cv_class`: Cross-validated predicted classes
#'   - `cv_posterior`: Cross-validated posterior probabilities
#' * `method`: "lda"
#' * `call`: The function call
#' * `formula`: Formula used
#' * `response_col`: Name of response column
#' * `coe_cols`: Names of coefficient columns used
#' * `covariate_cols`: Names of covariate columns used (if any)
#' * `predictor_cols`: All predictor column names
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
    coe_cols = coe_cols,
    covariate_cols = covariate_cols,
    predictor_cols = predictor_cols,
    kept_predictor_cols = kept_predictor_cols,  # After removing constant/collinear
    n_groups = n_groups,
    n_ld = n_ld
  )

  class(result) <- c("stat_lda", "momstats")
  result
}


# Collect method ----

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
  pred_matrix <- build_predictor_matrix(data, x$coe_cols, x$covariate_cols)

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


# Print method ----

#' @export
print.stat_lda <- function(x, ...) {
  cat("Linear Discriminant Analysis\n")
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


# Summary method ----

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

  invisible(NULL)
}


#' @keywords internal
#' @noRd
plot_lda_scores <- function(x, lds, color_quo, labels_quo, chull, cex, legend) {

  # Get LD scores
  # Predict on same data to get LD scores
  pred_matrix <- build_predictor_matrix(x$data, x$coe_cols, x$covariate_cols)

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
    return(invisible(NULL))
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

    return(invisible(NULL))
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
}


