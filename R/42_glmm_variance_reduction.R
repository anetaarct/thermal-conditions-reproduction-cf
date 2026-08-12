#!/usr/bin/env Rscript
.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(glmmTMB); library(readr); library(tibble)})
args <- commandArgs(trailingOnly = TRUE)

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))
ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 4)
sd_year <- function(m) attr(VarCorr(m)$cond$year_f, "stddev")[[1]]
sd_female <- function(m) attr(VarCorr(m)$cond$female_id, "stddev")[[1]]
reduction <- function(full, reduced) 100 * (1 - sd_year(full)^2 / sd_year(reduced)^2)

fit_cached <- function(path, formula, disp, zi, family, start = NULL) {
  if (file.exists(path)) return(readRDS(path))
  m <- glmmTMB(formula, dispformula = disp, ziformula = zi, family = family,
               data = dat, control = ctrl, start = start)
  stopifnot(m$fit$convergence == 0L, isTRUE(m$sdr$pdHess))
  saveRDS(m, path); m
}

fs_full <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")
fs_disp <- ~ rel_LD + temp_mean_c + year_c + age_class
fs_start <- function(formula) {
  fixed_formula <- lme4::nobars(formula)
  wanted <- colnames(model.matrix(fixed_formula, data = dat))
  beta <- fixef(fs_full)$cond[wanted]
  stopifnot(!anyNA(beta))
  list(
    beta = unname(beta),
    betadisp = unname(fixef(fs_full)$disp),
    theta = unname(getME(fs_full, "theta"))
  )
}
fs_null_formula <- cbind(fledged, failed) ~ rel_LD * age_class +
  (1 | year_f) + (1 | female_id)
fs_temp_formula <- cbind(fledged, failed) ~ rel_LD * age_class + temp_mean_c +
  (1 | year_f) + (1 | female_id)
fs_null <- fit_cached(
  "models/glmmTMB/fledging_variance_no_seasonal.rds",
  fs_null_formula, fs_disp, ~0, betabinomial(), fs_start(fs_null_formula))
fs_temp <- fit_cached(
  "models/glmmTMB/fledging_variance_temperature_only.rds",
  fs_temp_formula, fs_disp, ~0, betabinomial(), fs_start(fs_temp_formula))
fs_share <- 100 * c(year = sd_year(fs_full)^2, female = sd_female(fs_full)^2) /
  (sd_year(fs_full)^2 + sd_female(fs_full)^2)
fs_out <- tribble(
  ~component, ~quantity, ~estimate, ~unit,
  "Location", "Year random-intercept SD", sd_year(fs_full), "logit scale",
  "Location", "Female random-intercept SD", sd_female(fs_full), "logit scale",
  "Location", "Share of modelled random-intercept variance: year", fs_share["year"], "%",
  "Location", "Share of modelled random-intercept variance: female", fs_share["female"], "%",
  "Location", "Reduction in year random-intercept variance: temperature and year",
  reduction(fs_full, fs_null), "%",
  "Location", "Reduction in year random-intercept variance: temperature alone",
  reduction(fs_temp, fs_null), "%",
  "Precision", "Random effects", NA_real_, "none specified")
write_csv(fs_out, "tables/fledging_variance_decomposition.csv")

if ("fledging" %in% args) {
  print(fs_out, n = Inf)
  quit(save = "no", status = 0)
}

hu_full <- readRDS("models/glmmTMB/fledged_count_hurdle_poisson.rds")
hu_zi <- ~ rel_LD + temp_mean_c + year_c + age_class
hu_null <- fit_cached(
  "models/glmmTMB/realised_variance_no_seasonal_nofemale.rds",
  fledged ~ rel_LD * age_class + (1 | year_f), ~1, hu_zi, truncated_poisson())
hu_temp <- fit_cached(
  "models/glmmTMB/realised_variance_temperature_only_nofemale.rds",
  fledged ~ rel_LD * age_class + temp_mean_c + (1 | year_f),
  ~1, hu_zi, truncated_poisson())
hu_out <- tribble(
  ~component, ~quantity, ~estimate, ~unit,
  "Positive count", "Year random-intercept SD", sd_year(hu_full), "log scale",
  "Positive count", "Reduction in year random-intercept variance: temperature and year",
  reduction(hu_full, hu_null), "%",
  "Positive count", "Reduction in year random-intercept variance: temperature alone",
  reduction(hu_temp, hu_null), "%",
  "Complete failure", "Random effects", NA_real_, "none specified")
write_csv(hu_out, "tables/realised_fledgling_variance_decomposition.csv")

print(fs_out, n = Inf); print(hu_out, n = Inf)
