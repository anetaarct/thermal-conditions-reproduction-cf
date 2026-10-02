#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(glmmTMB); library(readr); library(tibble)})
dir.create("models/glmmTMB", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))
stopifnot(nrow(dat) == 9790L, nlevels(dat$year_f) == 39L)

ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 4)
loc_main <- cbind(fledged, failed) ~ rel_LD * age_class + temp_mean_c + year_c +
  (1 | year_f) + (1 | female_id)
disp_main <- ~ rel_LD + temp_mean_c + year_c + age_class
loc_sens <- update(loc_main, . ~ . + temp_sd_c)
disp_sens <- ~ rel_LD + temp_mean_c + temp_sd_c + year_c + age_class

fit <- function(loc, disp, start = NULL) glmmTMB(
  loc, dispformula = disp, ziformula = ~0, family = betabinomial(),
  data = dat, control = ctrl, start = start
)
main_file <- "models/glmmTMB/fledging_bb_main_no_tempsd.rds"
sens_file <- "models/glmmTMB/fledging_bb_sensitivity_tempsd.rds"
m_main <- if (file.exists(main_file)) readRDS(main_file) else fit(loc_main, disp_main)
saveRDS(m_main, main_file)
message("Saved main model; pdHess = ", m_main$sdr$pdHess)
m_sens <- if (file.exists(sens_file)) {
  readRDS(sens_file)
} else {
  old_sens <- readRDS(
    ".archive_before_final_fledging_recentring/models/fledging_bb_sensitivity_tempsd.rds"
  )
  fit(loc_sens, disp_sens, list(
    beta = unname(fixef(old_sens)$cond),
    betadisp = unname(fixef(old_sens)$disp),
    theta = unname(getME(old_sens, "theta"))
  ))
}
saveRDS(m_sens, sens_file)
message("Saved sensitivity model; pdHess = ", m_sens$sdr$pdHess)
stopifnot(m_main$fit$convergence == 0L, m_sens$fit$convergence == 0L,
          isTRUE(m_main$sdr$pdHess), isTRUE(m_sens$sdr$pdHess))

tidy_model <- function(m, label) {
  sm <- summary(m)$coefficients
  bind_rows(
    as.data.frame(sm$cond) |> rownames_to_column("term") |> as_tibble() |>
      transmute(model = label, component = "location_logit", term,
                beta = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`),
    as.data.frame(sm$disp) |> rownames_to_column("term") |> as_tibble() |>
      transmute(model = label, component = "precision_log_phi", term,
                beta = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`)
  )
}
coefs <- bind_rows(tidy_model(m_main, "main_no_temp_sd"),
                   tidy_model(m_sens, "sensitivity_with_temp_sd"))
write_csv(coefs, "tables/fledging_bb_main_sensitivity_coefficients.csv")

comparison <- full_join(
  filter(coefs, model == "main_no_temp_sd") |>
    select(component, term, main_beta = beta, main_SE = SE, main_p = p),
  filter(coefs, model == "sensitivity_with_temp_sd") |>
    select(component, term, sensitivity_beta = beta, sensitivity_SE = SE,
           sensitivity_p = p),
  by = c("component", "term")
)
write_csv(comparison, "tables/fledging_bb_main_vs_tempsd_sensitivity.csv")
write_csv(tibble(
  model = c("main_no_temp_sd", "sensitivity_with_temp_sd"),
  n = c(nobs(m_main), nobs(m_sens)), AIC = c(AIC(m_main), AIC(m_sens)),
  pdHess = c(m_main$sdr$pdHess, m_sens$sdr$pdHess),
  convergence_code = c(m_main$fit$convergence, m_sens$fit$convergence)
), "tables/fledging_bb_main_sensitivity_fit.csv")
message("Wrote FS main and sensitivity comparisons.")
