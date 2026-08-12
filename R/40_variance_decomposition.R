#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(brms); library(dplyr); library(glmmTMB); library(readr); library(tibble)
})
dir.create("models/brms", recursive = TRUE, showWarnings = FALSE)
dir.create("models/glmmTMB", recursive = TRUE, showWarnings = FALSE)

variance_reduction <- function(sd_full, sd_reduced) {
  100 * (1 - sd_full^2 / sd_reduced^2)
}
random_share <- function(year_sd, female_sd) {
  den <- year_sd^2 + female_sd^2
  c(year = 100 * year_sd^2 / den, female = 100 * female_sd^2 / den)
}
glmm_sd <- function(model, group) attr(VarCorr(model)$cond[[group]], "stddev")[[1]]
brms_sd <- function(model, parameter) {
  posterior_summary(model, pars = parameter)[1, "Estimate"]
}

# Clutch size. Each auxiliary model changes only the component whose annual
# random-intercept variance is being assessed. "Temperature alone" compares a
# temperature-only seasonal specification with the no-seasonal-predictor model.
cs_dat <- read_csv("data_derived/clutch_size_model_data.csv", show_col_types = FALSE) |>
  mutate(female_id = factor(female_id), year_f = factor(year_f),
         age_class = factor(age_class, levels = c("YOUNG", "OLD")))
cs_full <- readRDS("models/brms/clutch_size_main_no_tempsd_full.rds")

cs_formula <- function(loc_season, scale_season) bf(
  reformulate(c("rel_LD", loc_season, "age_class",
                "(1 | female_id)", "(1 | year_f)"), response = "clutch"),
  reformulate(c("rel_LD", scale_season, "age_class", "(1 | year_f)"),
              response = "sigma")
)
fit_cs_aux <- function(name, loc_season, scale_season, seed) {
  path <- file.path("models/brms", paste0(name, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  fit <- update(
    cs_full, formula. = cs_formula(loc_season, scale_season), newdata = cs_dat,
    backend = "cmdstanr", chains = 4, cores = 4,
    iter = 4000, warmup = 2000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    seed = seed, init = 0, refresh = 100
  )
  saveRDS(fit, path)
  fit
}
cs_loc_null <- fit_cs_aux("clutch_variance_location_no_seasonal", character(),
                          c("temp_mean_c", "year_c"), 4101)
cs_loc_temp <- fit_cs_aux("clutch_variance_location_temperature_only", "temp_mean_c",
                          c("temp_mean_c", "year_c"), 4102)
cs_scale_null <- fit_cs_aux("clutch_variance_scale_no_seasonal",
                            c("temp_mean_c", "year_c"), character(), 4103)
cs_scale_temp <- fit_cs_aux("clutch_variance_scale_temperature_only",
                            c("temp_mean_c", "year_c"), "temp_mean_c", 4104)

cs_year <- brms_sd(cs_full, "sd_year_f__Intercept")
cs_female <- brms_sd(cs_full, "sd_female_id__Intercept")
cs_scale_year <- brms_sd(cs_full, "sd_year_f__sigma_Intercept")
cs_shares <- random_share(cs_year, cs_female)
cs_out <- tribble(
  ~component, ~quantity, ~estimate, ~unit,
  "Location", "Year random-intercept SD", cs_year, "clutch-size scale",
  "Location", "Female random-intercept SD", cs_female, "clutch-size scale",
  "Location", "Share of modelled random-intercept variance: year", cs_shares["year"], "%",
  "Location", "Share of modelled random-intercept variance: female", cs_shares["female"], "%",
  "Location", "Reduction in year random-intercept variance: temperature and year",
  variance_reduction(cs_year, brms_sd(cs_loc_null, "sd_year_f__Intercept")), "%",
  "Location", "Reduction in year random-intercept variance: temperature alone",
  variance_reduction(brms_sd(cs_loc_temp, "sd_year_f__Intercept"),
                     brms_sd(cs_loc_null, "sd_year_f__Intercept")), "%",
  "Scale", "Year random-intercept SD", cs_scale_year, "log-sigma scale",
  "Scale", "Reduction in year random-intercept variance: temperature and year",
  variance_reduction(cs_scale_year, brms_sd(cs_scale_null, "sd_year_f__sigma_Intercept")), "%",
  "Scale", "Reduction in year random-intercept variance: temperature alone",
  variance_reduction(brms_sd(cs_scale_temp, "sd_year_f__sigma_Intercept"),
                     brms_sd(cs_scale_null, "sd_year_f__sigma_Intercept")), "%"
)
write_csv(cs_out, "tables/clutch_variance_decomposition.csv")

# Fledging success. The precision equation is held fixed in both auxiliary fits.
fs_dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))
fs_full <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")
fs_ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 1)
fs_disp <- ~ rel_LD + temp_mean_c + year_c + age_class
fit_fs_aux <- function(path, seasonal) {
  if (file.exists(path)) return(readRDS(path))
  loc <- reformulate(c("rel_LD * age_class", seasonal,
                       "(1 | year_f)", "(1 | female_id)"),
                     response = "cbind(fledged, failed)")
  fit <- glmmTMB(loc, dispformula = fs_disp, ziformula = ~0,
                 family = betabinomial(), data = fs_dat, control = fs_ctrl)
  stopifnot(fit$fit$convergence == 0L, isTRUE(fit$sdr$pdHess))
  saveRDS(fit, path); fit
}
fs_null <- fit_fs_aux("models/glmmTMB/fledging_variance_no_seasonal.rds", character())
fs_temp <- fit_fs_aux("models/glmmTMB/fledging_variance_temperature_only.rds", "temp_mean_c")
fs_year <- glmm_sd(fs_full, "year_f")
fs_female <- glmm_sd(fs_full, "female_id")
fs_shares <- random_share(fs_year, fs_female)
fs_out <- tribble(
  ~component, ~quantity, ~estimate, ~unit,
  "Location", "Year random-intercept SD", fs_year, "logit scale",
  "Location", "Female random-intercept SD", fs_female, "logit scale",
  "Location", "Share of modelled random-intercept variance: year", fs_shares["year"], "%",
  "Location", "Share of modelled random-intercept variance: female", fs_shares["female"], "%",
  "Location", "Reduction in year random-intercept variance: temperature and year",
  variance_reduction(fs_year, glmm_sd(fs_null, "year_f")), "%",
  "Location", "Reduction in year random-intercept variance: temperature alone",
  variance_reduction(glmm_sd(fs_temp, "year_f"), glmm_sd(fs_null, "year_f")), "%",
  "Precision", "Random effects", NA_real_, "none specified"
)
write_csv(fs_out, "tables/fledging_variance_decomposition.csv")

# Hurdle-Poisson. The zero equation is held fixed; the decomposition concerns
# the positive-count component, which contains a year random intercept. The
# female random intercept was removed after diagnosis confirmed a boundary fit.
hu_full <- readRDS("models/glmmTMB/fledged_count_hurdle_poisson.rds")
hu_zi <- ~ rel_LD + temp_mean_c + year_c + age_class
fit_hu_aux <- function(path, seasonal) {
  if (file.exists(path)) return(readRDS(path))
  cond <- reformulate(c("rel_LD * age_class", seasonal,
                        "(1 | year_f)"), response = "fledged")
  fit <- glmmTMB(cond, ziformula = hu_zi, family = truncated_poisson(),
                 data = fs_dat, control = fs_ctrl)
  stopifnot(fit$fit$convergence == 0L, isTRUE(fit$sdr$pdHess))
  saveRDS(fit, path); fit
}
hu_null <- fit_hu_aux("models/glmmTMB/realised_variance_no_seasonal.rds", character())
hu_temp <- fit_hu_aux("models/glmmTMB/realised_variance_temperature_only.rds", "temp_mean_c")
hu_year <- glmm_sd(hu_full, "year_f")
hu_out <- tribble(
  ~component, ~quantity, ~estimate, ~unit,
  "Positive count", "Year random-intercept SD", hu_year, "log scale",
  "Positive count", "Reduction in year random-intercept variance: temperature and year",
  variance_reduction(hu_year, glmm_sd(hu_null, "year_f")), "%",
  "Positive count", "Reduction in year random-intercept variance: temperature alone",
  variance_reduction(glmm_sd(hu_temp, "year_f"), glmm_sd(hu_null, "year_f")), "%",
  "Complete failure", "Random effects", NA_real_, "none specified"
)
write_csv(hu_out, "tables/realised_fledgling_variance_decomposition.csv")

print(cs_out, n = Inf); print(fs_out, n = Inf); print(hu_out, n = Inf)
