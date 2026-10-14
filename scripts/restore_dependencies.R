# Explicit, project-local restoration of the archived dependency versions.
# Run once with: Rscript --vanilla scripts/restore_dependencies.R
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script) != 1L) stop("Run this script with Rscript")
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
if (getRversion() != package_version("4.5.2"))
  stop("The lockfile requires R 4.5.2. Select that R version before restoration.")
library_path <- file.path(root, ".library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(library_path), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600,
        renv.config.auto.snapshot = FALSE)
bootstrap_version <- "1.2.4"
installed_bootstrap <- tryCatch(packageVersion("renv"), error = function(e) NULL)
if (is.null(installed_bootstrap) ||
    installed_bootstrap != package_version(bootstrap_version)) {
  if ("renv" %in% loadedNamespaces())
    stop("A different renv version is already loaded; restart with Rscript --vanilla.")
  tarball <- tempfile(fileext = ".tar.gz")
  bases <- c("https://cloud.r-project.org/src/contrib/Archive/renv/",
             "https://cloud.r-project.org/src/contrib/")
  downloaded <- FALSE
  for (base in bases) {
    url <- paste0(base, "renv_", bootstrap_version, ".tar.gz")
    downloaded <- tryCatch(isTRUE(download.file(url, tarball, mode = "wb") == 0),
                           error = function(e) FALSE)
    if (downloaded) break
  }
  if (!downloaded) stop("Could not retrieve renv 1.2.4; no alternate version was installed.")
  install.packages(tarball, repos = NULL, type = "source", lib = library_path)
  unlink(tarball)
}
if (packageVersion("renv") != package_version(bootstrap_version))
  stop("renv bootstrap version mismatch")
renv::restore(project = root, library = normalizePath(library_path),
              lockfile = file.path(root, "renv.lock"), prompt = FALSE,
              retry = FALSE, clean = FALSE, transactional = TRUE)
lock <- jsonlite::read_json(file.path(root, "renv.lock"), simplifyVector = FALSE)
stopifnot(getRversion() == package_version(lock$R$Version))
for (name in names(lock$Packages)) {
  expected <- package_version(lock$Packages[[name]]$Version)
  actual <- packageVersion(name, lib.loc = library_path)
  if (actual != expected) stop(paste("Locked version mismatch:", name))
}
cat("Restored and checked", length(lock$Packages), "locked packages.\n")
