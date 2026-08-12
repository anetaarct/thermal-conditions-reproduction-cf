#!/usr/bin/env Rscript

options(repos = c("https://stan-dev.r-universe.dev",
                  CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr")
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr installation failed")
}

path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
if (!nzchar(path) || !dir.exists(path)) {
  cmdstanr::install_cmdstan(cores = min(8L, parallel::detectCores()),
                           overwrite = FALSE)
}

cat("cmdstanr=", as.character(packageVersion("cmdstanr")), "\n")
cat("cmdstan_path=", cmdstanr::cmdstan_path(), "\n")
cat("cmdstan_version=", as.character(cmdstanr::cmdstan_version()), "\n")
