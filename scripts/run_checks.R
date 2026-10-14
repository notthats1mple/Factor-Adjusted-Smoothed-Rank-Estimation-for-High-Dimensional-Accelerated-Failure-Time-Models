#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1L])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "load_project.R"))
source(file.path(aft_project_root(), "tests", "public", "checks.R"))
options(warn = 1)
results <- list()
check <- function(name, passed) {
  success <- isTRUE(passed)
  results[[length(results) + 1L]] <<- list(name = name, pass = success)
  cat(if (success) "PASS" else "FAIL", name, "\n")
  invisible(success)
}
metadata <- tryCatch(public_checks(check), error = function(e) {
  results[[length(results) + 1L]] <<- list(name = "unexpected_exception", pass = FALSE,
                                          message = conditionMessage(e))
  cat("FAIL unexpected_exception:", conditionMessage(e), "\n")
  NULL
})
success <- length(results) > 0 && all(vapply(results, `[[`, logical(1), "pass"))
report <- list(status = if (success) "PASS" else "FAIL", checks = results, metadata = metadata,
               R_version = as.character(getRversion()), platform = R.version$platform)
trailing <- commandArgs(trailingOnly = TRUE)
if (length(trailing)) {
  if (length(trailing) != 2L || trailing[1] != "--output") stop("Usage: Rscript scripts/run_checks.R [--output result.json]")
  jsonlite::write_json(report, trailing[2], auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
cat(sprintf("%s: %d checks; fixed-array tests only.\n", report$status, length(results)))
if (!success) quit(status = 1L)
