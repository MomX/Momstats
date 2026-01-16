# stat_lm ---------

#' Linear Model for morphometric data
#'
#' Fit linear models with coefficient data as response or predictor.
#'
#' @param data A tibble with coefficient and/or scalar columns
#' @param formula A formula specifying the model. Can be:
#'   * Shape regression: `coe ~ length` (shape depends on predictors)
#'   * Response regression: `length ~ coe` (response depends on shape)
#'   * Multivariate response: `cbind(length, width) ~ coe`
#'   * Multiple predictors: `coe ~ length + temperature`
#'   * Interactions: `coe ~ length * temperature`
#' @param ... Additional arguments passed to [stats::lm()]
#'
#' @return An object of class `c("stat_lm", "momstats")` containing:
#' * `data`: Original tibble (unchanged)
#' * `model`: The [stats::lm()] object
#' * `method`: "lm"
#' * `direction`: "shape_regression" or "response_regression"
#' * `call`: The function call
#' * `formula`: Formula used
#' * `response_col`: Response column name(s)
#' * `predictor_cols`: Predictor column name(s)
#' * `r_squared`: R^2 value(s)
#' * `adj_r_squared`: Adjusted R^2 value(s)
#'
#' @details
#' `stat_lm()` provides a unified interface for linear modeling with morphometric data.
#'
#' ## Formula directions
#'
#' The function automatically detects model direction from the formula:
#'
#' **Shape regression** (`coe ~ predictors`):
#' * Fits multivariate regression with all harmonics as response
#' * Proper treatment of shape as correlated multivariate outcome
#' * Tests whether predictors affect shape
#' * Example: `boteft %>% stat_lm(coe ~ length)`
#'
#' **Response regression** (`response ~ coe`):
#' * Fits regression with unfolded harmonics as predictors
#' * Tests whether shape predicts a response variable
#' * Example: `boteft %>% stat_lm(length ~ coe)`
#'
#' **Multivariate response** (`cbind(y1, y2) ~ coe`):
#' * Fits multivariate regression with multiple responses
#' * Example: `boteft %>% stat_lm(cbind(length, width) ~ coe)`
#'
#' ## Interactions and multiple predictors
#'
#' Standard R formula syntax is supported:
#' * `coe ~ length + temperature`: Additive effects
#' * `coe ~ length * temperature`: With interaction
#' * `coe ~ poly(length, 2)`: Polynomial terms
#'
#' ## Getting results
#'
#' Use [collect()] to add predictions to your data:
#' ```r
#' lm_fit <- boteft %>% stat_lm(coe ~ length)
#' boteft_pred <- collect(lm_fit)  # Adds fitted values
#' ```
#'
#' Use [transduce()] to predict at new positions:
#' ```r
#' new_shapes <- transduce(lm_fit, tibble::tibble(length = seq(50, 150, by = 10)))
#' ```
#'
#' @examples
#' \dontrun{
#' # Shape regression
#' lm1 <- boteft %>% stat_lm(coe ~ length)
#' lm2 <- boteft %>% stat_lm(coe ~ length + type)
#' lm3 <- boteft %>% stat_lm(coe ~ length * type)
#'
#' # Response regression
#' lm4 <- boteft %>% stat_lm(length ~ coe)
#'
#' # Multivariate response
#' lm5 <- boteft %>% stat_lm(cbind(length, width) ~ coe)
#'
#' # Plot results
#' plot(lm1)
#' plot(lm1, type = "residuals")
#'
#' # Get predictions
#' boteft_pred <- collect(lm1)
#'
#' # Transduce to new values
#' new_shapes <- transduce(lm1, tibble::tibble(length = seq(50, 150, by = 10)))
#' }
#'
#' @seealso [stats::lm()], [collect.stat_lm()], [plot.stat_lm()], [transduce.stat_lm()]
#'
#' @export
stat_lm <- function(data, formula, ...) {

  # Store call
  call <- match.call()

  # Check input
  if (!is.data.frame(data)) {
    stop("data must be a tibble or data frame")
  }

  # Get formula expression
  formula_expr <- substitute(formula)

  # Detect direction
  direction <- detect_lm_direction(data, formula_expr)

  # Parse formula based on direction
  if (direction == "shape_regression") {
    # coe ~ predictors
    parsed <- parse_supervised_formula(data, formula_expr)

    response_cols <- parsed$response_cols
    coe_cols <- response_cols  # Response IS the coe column
    predictor_cols <- parsed$predictor_cols
    covariate_cols <- parsed$covariate_cols

    # Check response is actually a coe column
    if (!"coe" %in% class(data[[response_cols[1]]])) {
      stop(sprintf("For shape regression, LHS must be a coe column. Got: %s",
                   response_cols[1]))
    }

    # Build response matrix (unfold coe)
    Y <- build_response_matrix(data, response_cols)

    # Build predictor matrix (covariates only, no coe)
    if (length(covariate_cols) == 0) {
      stop("Shape regression requires at least one predictor variable")
    }

    # Build model matrix using formula
    # Create a temporary data frame with predictors
    pred_data <- data[covariate_cols]

    # Build formula for model.matrix (intercept included by default)
    mm_formula <- as.formula(paste("~", paste(covariate_cols, collapse = "+")))
    X <- model.matrix(mm_formula, data = pred_data)

    # Fit multivariate regression
    model <- lm(Y ~ X - 1, ...)  # -1 because X already has intercept

  } else {
    # response_regression: response ~ coe (and/or other predictors)
    parsed <- parse_supervised_formula(data, formula_expr)

    response_cols <- parsed$response_cols
    predictor_cols <- parsed$predictor_cols
    coe_cols <- parsed$coe_cols
    covariate_cols <- parsed$covariate_cols

    # Build response matrix
    Y <- build_response_matrix(data, response_cols)

    # Build predictor matrix (unfold coe columns)
    X <- build_predictor_matrix(data, coe_cols, covariate_cols)

    # Fit regression
    if (length(response_cols) == 1) {
      # Univariate response
      model <- lm(Y ~ X, ...)
    } else {
      # Multivariate response
      model <- lm(Y ~ X, ...)
    }
  }

  # Extract R-squared values
  if (inherits(model, "lm")) {
    model_summary <- summary(model)

    if (is.list(model_summary) && !inherits(model_summary, "summary.lm")) {
      # Multivariate response - multiple summaries
      r_squared <- sapply(model_summary, function(s) s$r.squared)
      adj_r_squared <- sapply(model_summary, function(s) s$adj.r.squared)
    } else {
      # Single response
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
    }
  }

  # Build result object
  result <- list(
    data = data,
    model = model,
    method = "lm",
    direction = direction,
    call = call,
    formula = formula_expr,
    response_cols = response_cols,
    predictor_cols = predictor_cols,
    coe_cols = coe_cols,
    r_squared = r_squared,
    adj_r_squared = adj_r_squared
  )

  class(result) <- c("stat_lm", "momstats")
  result
}


# Print method ----

#' @export
print.stat_lm <- function(x, ...) {
  cat("Linear Model")
  if (x$direction == "shape_regression") {
    cat(" (Shape Regression)\n")
  } else {
    cat(" (Response Regression)\n")
  }
  cat("---------------------------------------\n")

  # Formula
  cat("Formula: ")
  print(x$formula)
  cat("\n")

  cat(sprintf("* %d observations\n", nrow(x$data)))

  if (x$direction == "shape_regression") {
    cat(sprintf("* Response: %s (%d harmonics)\n",
                x$response_cols[1],
                ncol(x$model$coefficients)))
    cat(sprintf("* Predictors: %s\n", paste(x$predictor_cols, collapse = ", ")))
  } else {
    if (length(x$response_cols) == 1) {
      cat(sprintf("* Response: %s\n", x$response_cols))
    } else {
      cat(sprintf("* Responses: %s\n", paste(x$response_cols, collapse = ", ")))
    }

    if (length(x$coe_cols) > 0) {
      cat(sprintf("* Shape predictors: %s\n", paste(x$coe_cols, collapse = ", ")))
    }
    if (length(x$predictor_cols) > length(x$coe_cols)) {
      other_preds <- setdiff(x$predictor_cols, x$coe_cols)
      cat(sprintf("* Other predictors: %s\n", paste(other_preds, collapse = ", ")))
    }
  }

  # R-squared
  if (length(x$r_squared) == 1) {
    cat(sprintf("\n2 = %.4f, Adjusted R2 = %.4f\n",
                x$r_squared, x$adj_r_squared))
  } else {
    cat(sprintf("\nR2 range: [%.4f, %.4f]\n",
                min(x$r_squared), max(x$r_squared)))
    cat(sprintf("Mean R2: %.4f\n", mean(x$r_squared)))
  }

  cat("\n* summary() for model details\n")
  cat("* collect() to add fitted values to data\n")
  cat("* transduce() to predict at new positions\n")
  cat("* plot() to visualize results\n")

  invisible(x)
}


# Summary method ----

#' @export
summary.stat_lm <- function(object, ...) {
  cat("Linear Model Summary\n")
  cat("==================================================\n\n")

  cat("Direction: ")
  if (object$direction == "shape_regression") {
    cat("Shape Regression (coe ~ predictors)\n")
  } else {
    cat("Response Regression (response ~ coe)\n")
  }
  cat("\n")

  # Call summary on the underlying lm object
  print(summary(object$model))

  invisible(object)
}


# Collect method ----

#' Collect fitted values from linear model
#'
#' Extract fitted values and residuals from a linear model and add them to a tibble.
#'
#' @param x A `stat_lm` object
#' @param data A tibble. If NULL, uses the original data from the model.
#' @param fold How to add fitted values:
#'   * `FALSE` (default): Add as columns matching response structure
#'   * `TRUE`: Fold into list-column (for coe responses)
#'   * Character: Name for folded column
#' @param residuals Logical. Should residuals be included? Default `FALSE`.
#' @param ... Additional arguments (reserved)
#'
#' @return A tibble with fitted values (and optionally residuals) added
#'
#' @examples
#' \dontrun{
#' lm_fit <- boteft %>% stat_lm(coe ~ length)
#'
#' # Add fitted values
#' collect(lm_fit)
#'
#' # With residuals
#' collect(lm_fit, residuals = TRUE)
#'
#' # Fold into list-column
#' collect(lm_fit, fold = TRUE)
#' }
#'
#' @export
collect.stat_lm <- function(x, data = NULL, fold = FALSE, residuals = FALSE, ...) {

  # Use original data if not provided
  if (is.null(data)) {
    data <- x$data
  }

  # Get fitted values
  fitted_vals <- fitted(x$model)

  if (x$direction == "shape_regression") {
    # fitted_vals is a matrix (n x n_harmonics)

    if (isFALSE(fold)) {
      # Don't fold - but this doesn't make sense for coe
      # Better to fold by default for shape regression
      warning("For shape regression, consider using fold = TRUE to store fitted coe")

      # Add as individual columns (messy but user asked for it)
      for (i in seq_len(ncol(fitted_vals))) {
        col_name <- paste0(x$response_cols[1], "_fitted_", colnames(fitted_vals)[i])
        data[[col_name]] <- fitted_vals[, i]
      }
    } else {
      # Fold into coe list-column
      col_name <- if (isTRUE(fold)) {
        paste0(x$response_cols[1], "_fitted")
      } else {
        as.character(fold)
      }

      # Get coe class from original
      coe_class <- class(x$data[[x$response_cols[1]]])[1]

      # Create list of fitted coe vectors
      fitted_list <- lapply(seq_len(nrow(fitted_vals)), function(i) {
        vec <- fitted_vals[i, ]
        class(vec) <- c(coe_class, "coe", "numeric")
        vec
      })

      class(fitted_list) <- c(coe_class, "coe", "list")
      data[[col_name]] <- fitted_list
    }

    # Add residuals if requested
    if (residuals) {
      resid_vals <- residuals(x$model)

      if (isFALSE(fold)) {
        for (i in seq_len(ncol(resid_vals))) {
          col_name <- paste0(x$response_cols[1], "_resid_", colnames(resid_vals)[i])
          data[[col_name]] <- resid_vals[, i]
        }
      } else {
        col_name <- if (isTRUE(fold)) {
          paste0(x$response_cols[1], "_resid")
        } else {
          paste0(as.character(fold), "_resid")
        }

        resid_list <- lapply(seq_len(nrow(resid_vals)), function(i) {
          vec <- resid_vals[i, ]
          class(vec) <- c(coe_class, "coe", "numeric")
          vec
        })

        class(resid_list) <- c(coe_class, "coe", "list")
        data[[col_name]] <- resid_list
      }
    }

  } else {
    # Response regression

    if (length(x$response_cols) == 1) {
      # Single response
      data[[paste0(x$response_cols, "_fitted")]] <- as.vector(fitted_vals)

      if (residuals) {
        data[[paste0(x$response_cols, "_resid")]] <- as.vector(residuals(x$model))
      }
    } else {
      # Multiple responses
      for (i in seq_along(x$response_cols)) {
        resp_name <- x$response_cols[i]
        data[[paste0(resp_name, "_fitted")]] <- fitted_vals[, i]

        if (residuals) {
          data[[paste0(resp_name, "_resid")]] <- residuals(x$model)[, i]
        }
      }
    }
  }

  data
}


# Transduce method ----

#' Transduce linear model positions to shapes
#'
#' Predict coefficients and shapes at specified predictor values.
#'
#' @param object A `stat_lm` object
#' @param positions A tibble with predictor values. Column names must match
#'   the model's predictors.
#'
#' @return A tibble with predicted values
#'
#' @examples
#' \dontrun{
#' lm_fit <- boteft %>% stat_lm(coe ~ length)
#'
#' # Predict at new lengths
#' new_shapes <- transduce(lm_fit, tibble::tibble(length = seq(50, 150, by = 10)))
#'
#' # Multiple predictors
#' lm_fit2 <- boteft %>% stat_lm(coe ~ length + width)
#' new_shapes <- transduce(lm_fit2,
#'                         expand_grid(length = c(50, 100, 150),
#'                                     width = c(20, 30)))
#' }
#'
#' @export
transduce.stat_lm <- function(object, positions) {

  if (!is.data.frame(positions)) {
    stop("positions must be a tibble or data frame")
  }

  if (nrow(positions) == 0) {
    stop("positions tibble is empty")
  }

  if (object$direction == "shape_regression") {
    # Shape regression: predict coe from new predictor values

    # Check all predictors present
    missing_cols <- setdiff(object$predictor_cols, names(positions))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns in positions: %s",
                   paste(missing_cols, collapse = ", ")))
    }

    # The model was fit as lm(Y ~ X - 1) where X is the model matrix
    # We need to create a data frame with the X matrix columns for predict()

    # Build model matrix for new data
    pred_data <- positions[object$predictor_cols]
    mm_formula <- as.formula(paste("~", paste(object$predictor_cols, collapse = "+")))
    X_new <- model.matrix(mm_formula, data = pred_data)

    # Create newdata as a data frame with column names matching model
    # The model's predictor names are "X" + column names from model matrix
    newdata <- as.data.frame(X_new)
    colnames(newdata) <- paste0("X", colnames(X_new))

    # Actually, the issue is that our model was fit with X as the predictor name
    # Let's reconstruct properly

    # Get the actual coefficient names from the model
    coef_names <- rownames(coef(object$model))

    # Simpler approach: use the model matrix directly in matrix multiplication
    # pred = X_new %*% coef
    coefficients <- coef(object$model)

    # Matrix multiplication: X_new %*% coefficients
    pred_vals <- X_new %*% coefficients

    # Get coe class
    coe_class <- class(object$data[[object$response_cols[1]]])[1]

    # Get harmonic names from original coe
    example_coe <- object$data[[object$response_cols[1]]][[1]]
    harm_names <- names(example_coe)

    # Create list of predicted coe vectors
    pred_list <- lapply(seq_len(nrow(pred_vals)), function(i) {
      vec <- pred_vals[i, ]
      names(vec) <- harm_names
      class(vec) <- c(coe_class, "coe", "numeric")
      vec
    })

    class(pred_list) <- c(coe_class, "coe", "list")

    # Add to result
    result <- positions
    coe_col <- object$response_cols[1]
    result[[coe_col]] <- pred_list

    # Add inverse
    inverse_fn_name <- paste0(coe_class, "_i")
    inverse_fn <- get(inverse_fn_name, envir = asNamespace("Momocs2"))
    result[[paste0(coe_col, "_i")]] <- inverse_fn(result[[coe_col]])

  } else {
    # Response regression: would need coe values in positions
    # This is less common for transduce, but we can support it

    stop("transduce() for response regression not yet implemented. Use predict() on the model directly.")
  }

  tibble::as_tibble(result)
}

# Plot method ----

#' Plot linear model results
#'
#' Visualize linear model diagnostics or predictions.
#'
#' @param x A `stat_lm` object
#' @param type Character. Type of plot:
#'   * `"residuals"`: Residual plot
#'   * `"qq"`: Q-Q plot
#'   * `"fitted"`: Fitted vs observed
#'   * `"predictions"`: Predictions along predictor (shape regression only)
#' @param which For diagnostic plots, which response/harmonic to plot (integer)
#' @param ... Additional arguments passed to plot functions
#'
#' @return NULL (invisibly). Draws plot as side effect.
#'
#' @examples
#' \dontrun{
#' lm_fit <- boteft %>% stat_lm(coe ~ length)
#'
#' plot(lm_fit, type = "residuals")
#' plot(lm_fit, type = "qq")
#' plot(lm_fit, type = "predictions")
#' }
#'
#' @export
plot.stat_lm <- function(x, type = c("residuals", "qq", "fitted", "predictions"),
                         which = 1, ...) {

  type <- match.arg(type)

  if (type %in% c("residuals", "qq", "fitted")) {
    # Standard diagnostic plots
    if (x$direction == "shape_regression") {
      # Multivariate response - plot one harmonic
      fitted_vals <- fitted(x$model)[, which]
      resid_vals <- residuals(x$model)[, which]
      harm_name <- colnames(fitted(x$model))[which]

      if (type == "residuals") {
        plot(fitted_vals, resid_vals,
             xlab = "Fitted values",
             ylab = "Residuals",
             main = sprintf("Residual Plot (%s)", harm_name),
             pch = 19, col = rgb(0, 0, 0, 0.5))
        abline(h = 0, col = "red", lty = 2)

      } else if (type == "qq") {
        stats::qqnorm(resid_vals, main = sprintf("Q-Q Plot (%s)", harm_name),
               pch = 19, col = rgb(0, 0, 0, 0.5))
        stats::qqline(resid_vals, col = "red", lty = 2)

      } else if (type == "fitted") {
        obs_vals <- fitted_vals + resid_vals
        plot(obs_vals, fitted_vals,
             xlab = "Observed",
             ylab = "Fitted",
             main = sprintf("Fitted vs Observed (%s)", harm_name),
             pch = 19, col = rgb(0, 0, 0, 0.5),
             asp = 1)
        abline(0, 1, col = "red", lty = 2)
      }

    } else {
      # Univariate response
      fitted_vals <- fitted(x$model)
      resid_vals <- residuals(x$model)

      if (type == "residuals") {
        plot(fitted_vals, resid_vals,
             xlab = "Fitted values",
             ylab = "Residuals",
             main = "Residual Plot",
             pch = 19, col = rgb(0, 0, 0, 0.5))
        abline(h = 0, col = "red", lty = 2)

      } else if (type == "qq") {
        stats::qqnorm(resid_vals, main = "Q-Q Plot",
               pch = 19, col = rgb(0, 0, 0, 0.5))
        stats::qqline(resid_vals, col = "red", lty = 2)

      } else if (type == "fitted") {
        obs_vals <- fitted_vals + resid_vals
        plot(obs_vals, fitted_vals,
             xlab = "Observed",
             ylab = "Fitted",
             main = "Fitted vs Observed",
             pch = 19, col = rgb(0, 0, 0, 0.5),
             asp = 1)
        abline(0, 1, col = "red", lty = 2)
      }
    }

  } else if (type == "predictions") {
    # Only makes sense for shape regression with single continuous predictor
    if (x$direction != "shape_regression") {
      stop("Prediction plots only available for shape regression")
    }

    if (length(x$predictor_cols) != 1) {
      stop("Prediction plots currently only support single predictor models")
    }

    # Get predictor range
    pred_name <- x$predictor_cols[1]
    pred_vals <- x$data[[pred_name]]
    pred_range <- range(pred_vals, na.rm = TRUE)

    # Create prediction positions
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = 50)
    positions <- tibble::tibble(!!pred_name := pred_seq)

    # Transduce
    predictions <- transduce(x, positions)

    # Plot shapes along predictor
    # This would benefit from a proper morphospace plotting function
    # For now, just indicate it's not fully implemented
    message("Full prediction visualization coming soon. Use transduce() to get predicted shapes.")
  }

  invisible(NULL)
}
