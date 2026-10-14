#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1L])
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
source(file.path(dirname(normalizePath(script)), "..", "scripts", "load_project.R"))

# A balanced four-row array has an analytic partial-penalty solution.
X <- matrix(c(-1, 1, -1, 1), 4L, 1L, dimnames = list(NULL, "x1"))
Z <- matrix(c(-1, -1, 1, 1), 4L, 1L, dimnames = list(NULL, "z1"))
Y <- as.vector(Z + .4 * X)
delta <- rep(1, 4)
design <- aft_preprocess(X, Z, K = 0L)
expected <- c(1, .4 - .3 * qnorm(.575))
config <- aft_public_config()
config$fista$max_iter <- 2000L # Small example budget; saved configuration stays unchanged.
fit <- aft_solve_cv_v012(Y, delta, design, initial = expected + .01,
                          h = .6, lambda = .1, config = config)
stopifnot(isTRUE(fit$success), fit$standard$kkt$residual_with_underflow_budget <= 1e-7)
print(data.frame(term = c("unpenalized z1", "penalized x1"),
                 analytic = expected, fitted = fit$standard$coefficients))
cat("Status:", fit$status, "\nAccepted FISTA updates:", fit$fista$iterations,
    "\nFull KKT + budget:", fit$standard$kkt$residual_with_underflow_budget, "\n")

# Factor-adjusted prediction uses means and PCA learned from six fixed training rows.
# The final two rows illustrate new observations; they never enter preprocessing.
Xfa <- rbind(c(1, 3, -1), c(2, -2, 1), c(4, 0, 3), c(-1, 1, -2),
             c(0, 4, 2), c(3, -3, 0), c(5, -4, 2), c(-3, 2, -1))
colnames(Xfa) <- sprintf("gene_%04d", 1:3)
Zfa <- matrix(c(2, 0, -1, 3, 1, -2, 4, -3), 8L, 1L, dimnames = list(NULL, "clin_01"))
prepared <- p3a_prepare_split(Xfa, Zfa, c(-2, -1, 0, 1, 2, 3, 4, -4),
                              c(1, 1, 0, 1, 1, 0, 1, 1), 1:8, 1:6, 7:8, K = 1L)
illustrative_coefficients <- c(.2, .3, .4, 0, -.2)
scores <- p3a_predict_frozen(prepared, illustrative_coefficients,
                            Xfa[7:8, , drop = FALSE], Zfa[7:8, , drop = FALSE])
cat("Frozen FA prediction with explicitly supplied illustrative coefficients:\n")
print(data.frame(row = 7:8, score = scores))
