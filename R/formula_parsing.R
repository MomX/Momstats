# formula_parsing.R
# Formula parsing and predictor matrix building for statistical methods

# Formula parsing ----

#' Parse unsupervised formula
#'
#' Extract predictor columns from formula for unsupervised methods (PCA, clustering, etc.)
#'
#' @param data A tibble
#' @param formula_quo A quosure or expression containing the formula or bare column name
#'
#' @return A list with:
#' * `coe_cols`: Coefficient column names
#' * `covariate_cols`: Covariate column names
#' * `predictor_cols`: All predictor column names
#' * `messages`: Messages about column detection
#'
#' @details
#' Supports multiple syntax styles:
#' * Bare column: `coe` or `length`
#' * Formula: `~ coe`, `~ coe + length`
#' * Auto-detection: `~ .` (all coe columns)
#' * Tidyselect: `~ starts_with("Petal")`, `~ contains("Width")`
#' * Combinations: `~ starts_with("Petal") + length`
#'
#' @keywords internal
#' @noRd
parse_unsupervised_formula <- function(data, formula_quo) {

  # Extract expression (works with both enquo and raw expressions)
  if (inherits(formula_quo, "quosure")) {
    formula_expr <- rlang::quo_get_expr(formula_quo)
  } else {
    formula_expr <- formula_quo
  }

  # If missing/NULL, auto-detect single or multiple coe
  if (is.null(formula_expr) || identical(formula_expr, quote(expr = ))) {
    coe_cols <- get_all_coe_cols(data)
    return(list(
      coe_cols = coe_cols,
      covariate_cols = character(0),
      predictor_cols = coe_cols,
      messages = sprintf("Using coe column: '%s'", coe_cols)
    ))
  }

  # Case 1: Bare column name (not a formula)
  if (!is.call(formula_expr) || !identical(formula_expr[[1]], as.name("~"))) {
    # It's a bare name like `coe`
    col_name <- as.character(formula_expr)

    # Check if it exists
    if (!col_name %in% names(data)) {
      stop(sprintf("Column '%s' not found in data", col_name))
    }

    # Check if it's a coe column
    if ("coe" %in% class(data[[col_name]])) {
      return(list(
        coe_cols = col_name,
        covariate_cols = character(0),
        predictor_cols = col_name,
        messages = sprintf("Using coe column: '%s'", col_name)
      ))
    } else {
      # It's a covariate
      return(list(
        coe_cols = character(0),
        covariate_cols = col_name,
        predictor_cols = col_name,
        messages = sprintf("Using covariate column: '%s'", col_name)
      ))
    }
  }

  # Case 2: Formula
  # Check if one-sided or two-sided
  if (length(formula_expr) == 2) {
    # One-sided: ~ rhs
    rhs <- formula_expr[[2]]
  } else if (length(formula_expr) == 3) {
    # Two-sided: lhs ~ rhs (for unsupervised, just use rhs)
    rhs <- formula_expr[[3]]
  } else {
    stop("Invalid formula")
  }

  # Parse RHS to get predictor terms
  predictor_terms <- parse_formula_terms(rhs, data)

  list(
    coe_cols = predictor_terms$coe_cols,
    covariate_cols = predictor_terms$covariate_cols,
    predictor_cols = predictor_terms$predictor_cols,
    messages = predictor_terms$messages
  )
}


#' Parse supervised formula
#'
#' Extract response and predictor columns from formula for supervised methods (LDA, MANOVA, lm, etc.)
#'
#' @param data A tibble
#' @param formula_expr An expression from substitute() or enquo()
#'
#' @return A list with:
#' * `response_col`: Response column name (first response, for backwards compatibility)
#' * `response_cols`: Response column names (for multivariate response, or wraps single response)
#' * `coe_cols`: Coefficient column names
#' * `covariate_cols`: Covariate column names
#' * `predictor_cols`: All predictor column names
#' * `messages`: Messages about column detection
#'
#' @details
#' Supports multiple syntax styles:
#' * Bare response: `species` (auto-detects all coe as predictors)
#' * Simple formula: `species ~ coe`, `species ~ coe + length`
#' * Multivariate response: `cbind(length, width) ~ coe`
#' * Auto-detection: `species ~ .` (all coe columns)
#' * Tidyselect: `species ~ starts_with("harm")`
#'
#' @keywords internal
#' @noRd
parse_supervised_formula <- function(data, formula_expr) {

  # Case 1: Bare column name → auto-grab all coe columns as predictors
  if (!is.call(formula_expr) || !identical(formula_expr[[1]], as.name("~"))) {
    # It's a bare name like `species`
    response_col <- as.character(formula_expr)

    # Check if it exists
    if (!response_col %in% names(data)) {
      stop(sprintf("Column '%s' not found in data", response_col))
    }

    # Auto-grab all coe columns
    coe_cols <- Momocs2::get_all_coe_cols(data)

    return(list(
      response_col = response_col,
      response_cols = response_col,
      coe_cols = coe_cols,
      covariate_cols = character(0),
      predictor_cols = coe_cols,
      messages = sprintf("Response: '%s', using all coe columns: %s",
                         response_col, paste(coe_cols, collapse = ", "))
    ))
  }

  # Case 2: Formula with response ~ predictors
  # Extract LHS (response)
  lhs <- formula_expr[[2]]
  if (is.null(lhs)) {
    stop("Formula must have a response variable (e.g., species ~ coe or coe ~ length)")
  }

  # Handle multivariate response (cbind(y1, y2) ~ x)
  response_cols <- if (is.call(lhs) && identical(lhs[[1]], as.name("cbind"))) {
    # Multivariate response
    as.character(lhs[-1])
  } else {
    # Single response
    as.character(lhs)
  }

  # For backwards compatibility, keep response_col as first/only response
  response_col <- response_cols[1]

  # Check all responses exist
  missing <- setdiff(response_cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("Response column(s) not found: %s", paste(missing, collapse = ", ")))
  }

  # Extract RHS (predictors)
  rhs <- formula_expr[[3]]

  # Parse RHS to get predictor terms
  predictor_terms <- parse_formula_terms(rhs, data)

  # Build messages
  resp_msg <- if (length(response_cols) == 1) {
    sprintf("Response: '%s'", response_col)
  } else {
    sprintf("Responses: %s", paste(response_cols, collapse = ", "))
  }

  list(
    response_col = response_col,
    response_cols = response_cols,
    coe_cols = predictor_terms$coe_cols,
    covariate_cols = predictor_terms$covariate_cols,
    predictor_cols = predictor_terms$predictor_cols,
    messages = c(resp_msg, predictor_terms$messages)
  )
}


#' Detect formula direction for stat_lm
#'
#' Determine if formula is shape regression (coe ~ x) or response regression (y ~ coe)
#'
#' @param data A tibble
#' @param formula_expr Formula expression
#'
#' @return Character: "shape_regression" or "response_regression"
#'
#' @details
#' * Shape regression: `coe ~ length` - predicts coefficients from predictors
#' * Response regression: `length ~ coe` - predicts response from coefficients
#'
#' @keywords internal
#' @noRd
detect_lm_direction <- function(data, formula_expr) {

  # Must be a two-sided formula
  if (!is.call(formula_expr) || !identical(formula_expr[[1]], as.name("~"))) {
    stop("stat_lm requires a two-sided formula (e.g., coe ~ length or length ~ coe)")
  }

  if (length(formula_expr) != 3) {
    stop("stat_lm requires a two-sided formula (e.g., coe ~ length or length ~ coe)")
  }

  # Extract LHS
  lhs <- formula_expr[[2]]
  lhs_name <- as.character(lhs)

  # Check if LHS is a coe column
  if (lhs_name %in% names(data) && "coe" %in% class(data[[lhs_name]])) {
    return("shape_regression")  # coe ~ predictors
  } else {
    return("response_regression")  # response ~ coe (and/or other predictors)
  }
}


#' Parse formula terms
#'
#' Extract coe and covariate columns from formula RHS
#' Handles "." expansion to all coe columns and tidyselect helpers
#'
#' @param rhs Right-hand side of formula (expression)
#' @param data A tibble
#'
#' @return List with coe_cols, covariate_cols, predictor_cols, messages
#'
#' @details
#' Supports:
#' * Simple names: `coe`, `length`
#' * Dot expansion: `.` (all coe columns)
#' * Tidyselect helpers: `starts_with("Petal")`, `contains("Width")`, etc.
#' * Addition: `coe + length`, `starts_with("Petal") + length`
#'
#' Note: Exclusion syntax (`-Species`) is not currently supported.
#'
#' @keywords internal
#' @noRd
parse_formula_terms <- function(rhs, data) {

  coe_cols <- character(0)
  covariate_cols <- character(0)
  messages <- character(0)

  # Special case: "." expands to all coe columns
  if (identical(rhs, as.name("."))) {
    detected <- Momocs2::get_all_coe_cols(data)
    coe_cols <- c(coe_cols, detected)
    messages <- c(messages, sprintf("Expanded '.' to all coe columns: %s",
                                    paste(detected, collapse = ", ")))

    return(list(
      coe_cols = unique(coe_cols),
      covariate_cols = unique(covariate_cols),
      predictor_cols = unique(c(coe_cols, covariate_cols)),
      messages = messages
    ))
  }

  # Check if it's a tidyselect expression
  if (is.call(rhs)) {
    fn_name <- as.character(rhs[[1]])

    # Tidyselect helpers
    if (fn_name %in% c("all_of", "any_of", "starts_with", "ends_with",
                       "contains", "matches", "num_range", "everything",
                       "last_col", "where")) {

      # Evaluate tidyselect expression
      selected_cols <- tryCatch(
        names(tidyselect::eval_select(rhs, data)),
        error = function(e) {
          stop(sprintf("Error evaluating tidyselect expression: %s", e$message))
        }
      )

      # Classify each selected column
      for (col in selected_cols) {
        if ("coe" %in% class(data[[col]])) {
          coe_cols <- c(coe_cols, col)
        } else {
          covariate_cols <- c(covariate_cols, col)
        }
      }

      messages <- c(messages, sprintf("Selected columns: %s",
                                      paste(selected_cols, collapse = ", ")))

      return(list(
        coe_cols = unique(coe_cols),
        covariate_cols = unique(covariate_cols),
        predictor_cols = unique(c(coe_cols, covariate_cols)),
        messages = messages
      ))
    }

    # Addition operator: handle recursively
    if (fn_name == "+") {
      # Parse left and right sides
      left_terms <- parse_formula_terms(rhs[[2]], data)
      right_terms <- parse_formula_terms(rhs[[3]], data)

      # Combine results
      coe_cols <- c(left_terms$coe_cols, right_terms$coe_cols)
      covariate_cols <- c(left_terms$covariate_cols, right_terms$covariate_cols)
      messages <- c(left_terms$messages, right_terms$messages)

      return(list(
        coe_cols = unique(coe_cols),
        covariate_cols = unique(covariate_cols),
        predictor_cols = unique(c(coe_cols, covariate_cols)),
        messages = messages
      ))
    }
  }

  # Simple variable names
  terms_list <- all.vars(rhs)

  for (term in terms_list) {
    if (term %in% names(data)) {
      # Check if it's a coe column by class
      if ("coe" %in% class(data[[term]])) {
        coe_cols <- c(coe_cols, term)
        messages <- c(messages, sprintf("Using coe column: '%s'", term))
      } else {
        # It's a covariate
        covariate_cols <- c(covariate_cols, term)
        messages <- c(messages, sprintf("Using covariate: '%s'", term))
      }
    } else {
      stop(sprintf("Column '%s' not found in data", term))
    }
  }

  # Check we have at least something
  if (length(coe_cols) == 0 && length(covariate_cols) == 0) {
    stop("No valid predictors found in formula")
  }

  list(
    coe_cols = unique(coe_cols),
    covariate_cols = unique(covariate_cols),
    predictor_cols = unique(c(coe_cols, covariate_cols)),
    messages = messages
  )
}


# Matrix building ----

#' Build predictor matrix
#'
#' Combine coefficient and covariate columns into a single matrix for analysis
#'
#' @param data A tibble
#' @param coe_cols Coefficient column names
#' @param covariate_cols Covariate column names
#'
#' @return A numeric matrix with one row per observation
#'
#' @details
#' Coefficient columns are unfolded using [Momocs2::unfold()] to expand list-columns
#' into individual harmonic columns. Covariate columns are added as-is and must be numeric.
#'
#' @keywords internal
#' @noRd
build_predictor_matrix <- function(data, coe_cols, covariate_cols) {

  matrices <- list()

  # Unfold coefficient columns using Momocs2::unfold
  for (col in coe_cols) {
    unfolded <- Momocs2::unfold(data[[col]])
    matrices[[col]] <- unfolded
  }

  # Add covariate columns
  for (col in covariate_cols) {
    # Ensure it's numeric
    if (!is.numeric(data[[col]])) {
      stop(sprintf("Covariate '%s' must be numeric", col))
    }
    matrices[[col]] <- matrix(data[[col]], ncol = 1,
                              dimnames = list(NULL, col))
  }

  # Combine all matrices
  do.call(cbind, matrices)
}


#' Build response matrix
#'
#' Prepare response variable(s) for modeling
#'
#' @param data A tibble
#' @param response_cols Response column name(s)
#'
#' @return A numeric matrix or vector
#'
#' @details
#' For single responses, returns a vector. For multiple responses, returns a matrix.
#' If response is a coe column, it is unfolded first.
#'
#' @keywords internal
#' @noRd
build_response_matrix <- function(data, response_cols) {

  if (length(response_cols) == 1) {
    # Single response - return as vector
    y <- data[[response_cols]]

    # Check if it's a coe column (needs unfolding)
    if ("coe" %in% class(y)) {
      return(Momocs2::unfold(y))
    }

    # Otherwise ensure numeric
    if (!is.numeric(y)) {
      stop(sprintf("Response '%s' must be numeric", response_cols))
    }

    return(y)
  } else {
    # Multiple responses - return as matrix
    matrices <- list()

    for (col in response_cols) {
      y <- data[[col]]

      # Check if coe column
      if ("coe" %in% class(y)) {
        matrices[[col]] <- Momocs2::unfold(y)
      } else if (is.numeric(y)) {
        matrices[[col]] <- matrix(y, ncol = 1, dimnames = list(NULL, col))
      } else {
        stop(sprintf("Response '%s' must be numeric", col))
      }
    }

    return(do.call(cbind, matrices))
  }
}


# Plot helpers ----

#' Get plot variable (color or labels)
#'
#' Extract column for plotting, handling both bare and quoted names
#'
#' @param data A data frame
#' @param var_quo A quosure containing column name or NULL
#'
#' @return Vector of values or NULL
#'
#' @keywords internal
#' @noRd
get_plot_variable <- function(data, var_quo) {

  # If NULL or missing, return NULL
  if (missing(var_quo) || is.null(var_quo)) {
    return(NULL)
  }

  # Handle quosures
  if (inherits(var_quo, "quosure")) {
    if (rlang::quo_is_null(var_quo)) {
      return(NULL)
    }

    # Try as bare name first
    var_name <- tryCatch(
      rlang::as_name(var_quo),
      error = function(e) NULL
    )

    # If that failed, try evaluating (for quoted names)
    if (is.null(var_name)) {
      var_name <- rlang::eval_tidy(var_quo)
    }
  } else {
    # Not a quosure, just use as is
    var_name <- as.character(var_quo)
  }

  # Check if column exists
  if (!var_name %in% names(data)) {
    warning(sprintf("Column '%s' not found in data", var_name))
    return(NULL)
  }

  data[[var_name]]
}
