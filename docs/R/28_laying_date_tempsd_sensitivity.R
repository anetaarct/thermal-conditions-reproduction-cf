#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(glmmTMB); library(readr); library(tibble)})
dir.create("models/glmmTMB", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

raw <- read_csv("data_derived/reproductive_model_data.csv", show_col_types = FALSE) |>
  filter(between(year, 1982, 2025), year != 2020,
         age_class %in% c("YOUNG", "OLD")) |>
  transmute(
    lay_date = as.numeric(laying_date), temp_mean = as.numeric(temp_mean),
    temp_sd = as.numeric(temp_sd), year = as.integer(year),
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = if_else(is.na(female_id) | female_id == "",
                        paste0("UNRINGED_", source_row), female_id),
    year_f = factor(year)
  ) |>
  filter(if_all(everything(), ~ !is.na(.x))) |>
  droplevels()

dat <- raw |>
  mutate(temp_mean_c = temp_mean - mean(temp_mean),
         temp_sd_c = temp_sd - mean(temp_sd), year_c = year - mean(year),
         female_id = factor(female_id))

main <- readRDS("models/glmmTMB/laying_date_age_temperature.rds")
sens <- glmmTMB(
  lay_date ~ temp_mean_c * age_class + temp_sd_c + year_c +
    (1 | year_f) + (1 | female_id),
  dispformula = ~ age_class + temp_mean_c + temp_sd_c,
  family = gaussian(), data = dat,
  control = glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000))
)
stopifnot(sens$fit$convergence == 0L, isTRUE(sens$sdr$pdHess))
saveRDS(sens, "models/glmmTMB/laying_date_tempsd_sensitivity.rds")

tidy_model <- function(m, label) {
  sm <- summary(m)$coefficients
  bind_rows(
    as.data.frame(sm$cond) |> rownames_to_column("term") |> as_tibble() |>
      transmute(model = label, component = "location", term, beta = Estimate,
                SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`),
    as.data.frame(sm$disp) |> rownames_to_column("term") |> as_tibble() |>
      transmute(model = label, component = "scale_log_variance", term, beta = Estimate,
                SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`)
  )
}
coefs <- bind_rows(tidy_model(main, "main_no_temp_sd"),
                   tidy_model(sens, "sensitivity_with_temp_sd"))
comparison <- full_join(
  filter(coefs, model == "main_no_temp_sd") |>
    select(component, term, main_beta = beta, main_SE = SE, main_p = p),
  filter(coefs, model == "sensitivity_with_temp_sd") |>
    select(component, term, sensitivity_beta = beta, sensitivity_SE = SE,
           sensitivity_p = p), by = c("component", "term")
)
write_csv(coefs, "tables/laying_date_main_sensitivity_coefficients.csv")
write_csv(comparison, "tables/laying_date_main_vs_tempsd_sensitivity.csv")
write_csv(tibble(
  model = c("main_no_temp_sd", "sensitivity_with_temp_sd"),
  AIC = c(AIC(main), AIC(sens)), pdHess = c(main$sdr$pdHess, sens$sdr$pdHess),
  convergence_code = c(main$fit$convergence, sens$fit$convergence)
), "tables/laying_date_main_sensitivity_fit.csv")
message("Completed laying-date temp_sd sensitivity model.")
