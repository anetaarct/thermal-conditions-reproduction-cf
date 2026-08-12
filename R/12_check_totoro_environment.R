#!/usr/bin/env Rscript

packages <- c("brms", "cmdstanr", "rstan", "posterior", "readr", "dplyr", "ggplot2")
for (pkg in packages) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(packageVersion(pkg)) else NA_character_
  cat(pkg, "available=", available, "version=", version, "\n")
}
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("cmdstan_path=", tryCatch(cmdstanr::cmdstan_path(),
                                error = function(e) conditionMessage(e)), "\n")
}
cat("libPaths=", paste(.libPaths(), collapse = ";"), "\n")
