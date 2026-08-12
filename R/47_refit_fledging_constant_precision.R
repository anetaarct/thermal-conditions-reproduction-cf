#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(glmmTMB); library(readr)})
dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))
old <- readRDS(
  ".archive_before_final_fledging_recentring/models/fledging_bb_constant_precision.rds"
)
start <- list(
  beta = unname(fixef(old)$cond), betadisp = unname(fixef(old)$disp),
  theta = unname(getME(old, "theta"))
)
m <- glmmTMB(
  cbind(fledged, failed) ~ rel_LD * age_class + temp_mean_c + year_c +
    (1 | year_f) + (1 | female_id),
  dispformula = ~ 1, ziformula = ~ 0, family = betabinomial(), data = dat,
  start = start,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 4
  )
)
stopifnot(m$fit$convergence == 0L, isTRUE(m$sdr$pdHess))
saveRDS(m, "models/glmmTMB/fledging_bb_constant_precision.rds")
cat("convergence=", m$fit$convergence, " pdHess=", m$sdr$pdHess, "\n", sep = "")
