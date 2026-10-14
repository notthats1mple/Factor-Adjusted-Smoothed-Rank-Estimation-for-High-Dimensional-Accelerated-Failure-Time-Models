# Source this file by its full path from any working directory.
# All algorithm definitions are loaded into the source() target environment.
.aft_loader_env <- environment()
.aft_loader_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
if (!length(.aft_loader_files))
  stop("Load this entry point with source('/path/to/scripts/load_project.R')")
.aft_loader_root <- normalizePath(file.path(dirname(tail(.aft_loader_files, 1L)[[1L]]), ".."),
                                  winslash = "/", mustWork = TRUE)
# restore_dependencies.R installs locked packages here in a fresh checkout.
if (dir.exists(file.path(.aft_loader_root, ".library")))
  .libPaths(c(normalizePath(file.path(.aft_loader_root, ".library"), winslash = "/"), .libPaths()))
aft_project_root <- local({ root <- .aft_loader_root; function() root })
aft_public_config <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Install the jsonlite version recorded in renv.lock before loading configuration")
  jsonlite::fromJSON(file.path(aft_project_root(), "config", "solver_adapter.json"),
                     simplifyVector = FALSE)
}
.aft_loader_sources <- c(
  "R/preprocess.R", "R/loss.R", "R/diagnostics.R", "R/solver.R", "R/path.R",
  "diagnostics/phase2b/endpoint.R", "optim/phase2c/arithmetic.R", "optim/phase2c/fista.R",
  "optim/active_face_polish/arithmetic.R", "optim/active_face_polish/polish.R",
  "pilot/schema.R", "tuning/phase3a/data_validation.R", "tuning/phase3a/adapter.R")
for (.aft_loader_file in .aft_loader_sources)
  sys.source(file.path(.aft_loader_root, .aft_loader_file), envir = .aft_loader_env)
phase2c_arithmetic_load(.aft_loader_root)
rm(list = c(".aft_loader_env", ".aft_loader_files", ".aft_loader_root",
            ".aft_loader_sources", ".aft_loader_file"))
