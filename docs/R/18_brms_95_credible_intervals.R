#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(brms); library(readr); library(tibble)})

summarise_model <- function(path, output) {
  model <- readRDS(path)
  result <- posterior_summary(model, probs = c(0.025, 0.975)) |>
    as.data.frame() |>
    rownames_to_column("parameter") |>
    as_tibble()
  names(result) <- c("parameter", "estimate", "est_error",
                     "CrI_2.5", "CrI_97.5")
  write_csv(result, output)
  result
}

cs <- summarise_model(
  "models/brms/clutch_size_full.rds",
  "tables/brms_cs_full_posterior_95.csv"
)
fs <- summarise_model(
  "models/brms/fledging_positive_full.rds",
  "tables/brms_fs_full_posterior_95.csv"
)

keep <- function(x) {
  x[grepl("^(b_|sd_).*", x$parameter), ]
}
print(keep(cs), n = Inf)
print(keep(fs), n = Inf)
