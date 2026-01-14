# formula_parsing.R
# Formula parsing and predictor matrix building for statistical methods

# Formula parsing ----

#' Parse unsupervised formula
#'
#' Extract predictor columns from formula for unsupervised methods (PCA, clustering, etc.)
#'
#' @param data A tibble
#' @param formula_quo A quosure containing the formula or bare column name
#'
#' @return A list with:
#' * `coe_cols`: Coefficient column names
#' * `covariate_cols`: Covariate column names
#' * `predictor_cols`: All predictor column names
#'
#' @keywords internal
#' @noRd
parse_unsupervised_formula <- function(data, formula_quo) {

  # If missing/NULL, auto-detect
  if (rlang::quo_is_null(formula_quo)) {
    coe_cols <- get_coe_cols(data)
    return(list(
      coe_cols = coe_cols,
      covariate_cols = character(0),
      predictor_cols = coe_cols,
      messages = sprintf("Using coe column: '%s'", coe_cols)
    ))
  }

  # Get the actual expression
  formula_expr <- rlang::quo_get_expr(formula_quo)

  # Case 1: Bare column name (not a formula)
  if (!rlang::is_call(formula_expr, "~")) {
    # It's a bare name like `coe`
    col_name <- rlang::as_name(formula_quo)

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
  # Extract RHS (there's no LHS for unsupervised)
  rhs <- rlang::f_rhs(rlang::eval_tidy(formula_quo))

  # Parse RHS to get predictor terms
  predictor_terms <- parse_formula_terms(rhs, data)

  list(
    coe_cols = predictor_terms$coe_cols,
    covariate_cols = predictor_terms$covariate_cols,
    predictor_cols = predictor_terms$predictor_cols,
    messages = predictor_terms$messages
  )
}


#' Parse formula terms
#'
#' Extract coe and covariate columns from formula RHS
#'
#' @param rhs Right-hand side of formula
#' @param data A tibble
#'
#' @return List with coe_cols, covariate_cols, predictor_cols, messages
#'
#' @keywords internal
#' @noRd
parse_formula_terms <- function(rhs, data) {

  # Get all term names from RHS
  terms_list <- all.vars(rhs)

  coe_cols <- character(0)
  covariate_cols <- character(0)
  messages <- character(0)

  for (term in terms_list) {
    # Special case: "coe" triggers auto-detection
    if (term == "coe") {
      detected <- get_coe_cols(data)
      coe_cols <- c(coe_cols, detected)
      messages <- c(messages, sprintf("Auto-detected coe column: '%s'", detected))
    } else if (term %in% names(data)) {
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


# Plot helpers ----

#' Get plot variable (color or labels)
#'
#' Extract column for plotting, handling both bare and quoted names
#'
#' @param data A data frame
#' @param var_quo A quosure containing column name
#'
#' @return Vector of values or NULL
#'
#' @keywords internal
#' @noRd
get_plot_variable <- function(data, var_quo) {

  # If NULL, return NULL
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

  # Check if column exists
  if (!var_name %in% names(data)) {
    warning(sprintf("Column '%s' not found in data", var_name))
    return(NULL)
  }

  data[[var_name]]
}
