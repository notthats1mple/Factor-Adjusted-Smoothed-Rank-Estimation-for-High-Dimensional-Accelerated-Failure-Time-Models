# Phase 2a protocol Section 7: experiment-layer feature schema guard.
# No core API changes, fitting, outcome access, column reordering, or RNG.

.phase2a_schema_stop <- function(code, message, block = NULL,
                                 expected = NULL, observed = NULL) {
  stop(structure(list(message = message, call = NULL, code = code,
                      block = block, expected = expected, observed = observed),
                 class = c("phase2a_schema_error", "error", "condition")))
}

.phase2a_names <- function(size, block) {
  if (block == "X") sprintf("gene_%04d", seq_len(size)) else
    sprintf("clin_%02d", seq_len(size))
}

.phase2a_matrix_columns <- function(x, block) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x))
    .phase2a_schema_stop("invalid_matrix", paste(block, "must be a numeric matrix"), block)
  if (any(!is.finite(x)))
    .phase2a_schema_stop("nonfinite_input", paste(block, "contains a non-finite value"), block)
  if (ncol(x) == 0L) return(character())
  columns <- colnames(x)
  if (is.null(columns) || length(columns) != ncol(x) ||
      anyNA(columns) || any(!nzchar(trimws(columns))))
    .phase2a_schema_stop("missing_column_names",
      paste(block, "requires a nonempty name for every column"), block,
      observed = columns)
  if (anyDuplicated(columns))
    .phase2a_schema_stop("duplicate_column_names",
      paste(block, "column names must be unique"), block, observed = columns)
  columns
}

.phase2a_check_saved_schema <- function(schema) {
  if (!is.list(schema) || !identical(schema$format, "phase2a-training-schema") ||
      !identical(schema$version, 1L) || length(schema$training_rows) != 1L ||
      !is.numeric(schema$training_rows) || !is.finite(schema$training_rows) ||
      schema$training_rows < 2 || schema$training_rows != floor(schema$training_rows))
    .phase2a_schema_stop("invalid_saved_schema", "Invalid or unsupported saved training schema")
  for (block in c("X", "Z")) {
    entry <- schema[[block]]
    if (!is.list(entry) || length(entry$ncol) != 1L || !is.numeric(entry$ncol) ||
        !is.finite(entry$ncol) || entry$ncol < 0 || entry$ncol != floor(entry$ncol) ||
        !is.character(entry$columns) || length(entry$columns) != entry$ncol ||
        anyNA(entry$columns) || anyDuplicated(entry$columns) ||
        !identical(entry$columns, .phase2a_names(entry$ncol, block)))
      .phase2a_schema_stop("invalid_saved_schema",
        paste("Invalid saved", block, "column schema"), block)
  }
  invisible(TRUE)
}

# The returned plain list can be saved with saveRDS()/readRDS(), or serialized
# as JSON. The runner persists it once with the training-data record. Protocol
# dimensions (p=300,d=2) are enforced by the locked experiment configuration;
# this helper supports smaller named matrices for deterministic guard tests.
phase2a_schema <- function(X, Z = NULL) {
  x_columns <- .phase2a_matrix_columns(X, "X")
  if (nrow(X) < 2L)
    .phase2a_schema_stop("invalid_training_rows", "Training schema requires at least two rows")
  if (is.null(Z)) Z <- matrix(numeric(), nrow(X), 0L)
  z_columns <- .phase2a_matrix_columns(Z, "Z")
  if (nrow(Z) != nrow(X))
    .phase2a_schema_stop("row_count_mismatch", "Training X and Z row counts differ", "Z",
                        expected = nrow(X), observed = nrow(Z))
  for (block in c("X", "Z")) {
    columns <- if (block == "X") x_columns else z_columns
    expected <- .phase2a_names(length(columns), block)
    if (!identical(columns, expected))
      .phase2a_schema_stop("noncanonical_training_names",
        paste("Training", block, "names/order do not match the Phase 2a convention"),
        block, expected = expected, observed = columns)
  }
  list(format = "phase2a-training-schema", version = 1L, training_rows = nrow(X),
       X = list(ncol = ncol(X), columns = x_columns),
       Z = list(ncol = ncol(Z), columns = z_columns),
       matching = "exact column names, order, and count; duplicates forbidden",
       automatic_reordering = FALSE)
}

# Validation never changes the supplied matrices or the saved schema.
phase2a_validate_schema <- function(schema, X, Z = NULL) {
  .phase2a_check_saved_schema(schema)
  x_columns <- .phase2a_matrix_columns(X, "X")
  if (is.null(Z)) {
    if (schema$Z$ncol != 0L)
      .phase2a_schema_stop("missing_clinical_matrix", "Z is required by the saved training schema", "Z")
    Z <- matrix(numeric(), nrow(X), 0L)
  }
  z_columns <- .phase2a_matrix_columns(Z, "Z")
  if (nrow(Z) != nrow(X))
    .phase2a_schema_stop("row_count_mismatch", "Prediction X and Z row counts differ", "Z",
                        expected = nrow(X), observed = nrow(Z))
  for (block in c("X", "Z")) {
    actual <- if (block == "X") x_columns else z_columns
    expected <- schema[[block]]$columns
    if (length(actual) != schema[[block]]$ncol)
      .phase2a_schema_stop("column_count_mismatch",
        paste(block, "column count differs from the saved training schema"), block,
        expected = schema[[block]]$ncol, observed = length(actual))
    if (!identical(actual, expected))
      .phase2a_schema_stop("column_names_or_order_mismatch",
        paste(block, "column names/order differ from the saved training schema; no reordering is performed"),
        block, expected = expected, observed = actual)
  }
  invisible(TRUE)
}

phase2a_predict <- function(fit, schema, Xnew, Znew = NULL) {
  if (!inherits(fit, "aft_fit") || !inherits(fit$design, "aft_design") ||
      !identical(fit$design$kind, "pca"))
    .phase2a_schema_stop("invalid_prediction_fit", "Prediction requires a core aft_fit with training preprocessing")
  .phase2a_check_saved_schema(schema)
  if (!identical(as.integer(fit$design$n), as.integer(schema$training_rows)))
    .phase2a_schema_stop("training_schema_mismatch", "Fit training row count differs from the saved schema")
  # Also prevent attaching a different schema to an already fitted design.
  phase2a_validate_schema(schema, fit$design$Xc, fit$design$Zc)
  phase2a_validate_schema(schema, Xnew, Znew)
  # Preserve the core's success check and frozen training mapping exactly.
  predict.aft_fit(fit, Xnew = Xnew, Znew = Znew)
}
