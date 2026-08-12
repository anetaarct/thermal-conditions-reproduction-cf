#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(glmmTMB)
})

dat <- read_csv("data_derived/reproductive_model_data.csv",
                show_col_types = FALSE) |>
  transmute(
    source_row = as.integer(source_row),
    clutch = as.numeric(clutch), rel_LD = as.numeric(rel_LD),
    temp_mean = as.numeric(temp_mean), temp_sd = as.numeric(temp_sd),
    year_c = as.numeric(year_c), year_f = factor(year_f),
    female_id = if_else(is.na(female_id) | female_id == "",
                        paste0("UNRINGED_", source_row), female_id),
    age_class = factor(age_class, levels = c("YOUNG", "OLD"))
  ) |>
  filter(if_all(everything(), ~ !is.na(.x))) |>
  droplevels() |>
  mutate(
    female_id = factor(female_id),
    temp_mean_c = temp_mean - mean(temp_mean),
    temp_sd_c = temp_sd - mean(temp_sd)
  ) |>
  group_by(year_f) |>
  mutate(n_broods_year = n(), log_n_broods = log(n())) |>
  ungroup() |>
  mutate(log_n_broods_c = log_n_broods -
           mean(log_n_broods[!duplicated(year_f)]))

stopifnot(nrow(dat) > 0, all(between(dat$clutch, 2, 9)),
          nlevels(dat$year_f) == 44)

dispersion <- tibble(
  n = nrow(dat), females = n_distinct(dat$female_id),
  seasons = n_distinct(dat$year_f), mean_clutch = mean(dat$clutch),
  variance_clutch = var(dat$clutch),
  variance_to_mean = variance_clutch / mean_clutch,
  minimum = min(dat$clutch), maximum = max(dat$clutch)
)

clutch_distribution <- dat |> count(clutch, name = "records") |>
  mutate(proportion = records / sum(records))

mean_formula <- clutch ~ rel_LD + temp_mean_c + temp_sd_c + year_c + age_class +
  (1 | female_id) + (1 | year_f)
disp_formula <- ~ rel_LD + temp_mean_c + temp_sd_c + year_c + age_class

model_constant_scale <- glmmTMB(
  mean_formula, dispformula = ~1, family = gaussian(), data = dat,
  REML = FALSE
)
model_location_scale <- glmmTMB(
  mean_formula, dispformula = disp_formula, family = gaussian(), data = dat,
  REML = FALSE
)

model_scale_sample_size <- glmmTMB(
  mean_formula,
  dispformula = ~ rel_LD + temp_mean_c + temp_sd_c + year_c + age_class +
    log_n_broods_c,
  family = gaussian(), data = dat, REML = FALSE
)

comparison_raw <- anova(model_constant_scale, model_location_scale)
comparison <- tibble(
  model = c("Constant scale", "Location-scale"),
  df = comparison_raw$Df,
  AIC = comparison_raw$AIC,
  BIC = comparison_raw$BIC,
  logLik = comparison_raw$logLik,
  deviance = comparison_raw$deviance,
  chi_square = comparison_raw$Chisq,
  chi_df = comparison_raw$`Chi Df`,
  p = comparison_raw$`Pr(>Chisq)`
) |>
  mutate(delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC))

extract_component <- function(model, component) {
  x <- summary(model)$coefficients[[component]]
  as.data.frame(x) |>
    rownames_to_column("term") |>
    as_tibble() |>
    transmute(
      component = component, term,
      estimate = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`,
      CI_low = estimate - 1.96 * SE, CI_high = estimate + 1.96 * SE
    )
}

coefficients <- bind_rows(
  extract_component(model_location_scale, "cond"),
  extract_component(model_location_scale, "disp")
)

female_repeats <- dat |> count(female_id, name = "records")
repeat_summary <- female_repeats |>
  mutate(repeat_class = case_when(records == 1 ~ "1 record",
                                  records == 2 ~ "2 records",
                                  TRUE ~ "3+ records")) |>
  count(repeat_class, name = "females") |>
  mutate(proportion = females / sum(females))

vc <- VarCorr(model_location_scale)$cond
random_effects <- tibble(
  grouping_factor = c("female_id", "year_f"),
  SD = c(attr(vc$female_id, "stddev"), attr(vc$year_f, "stddev"))
)

clutch_two <- read_csv("data_derived/reproductive_model_data.csv",
                       show_col_types = FALSE) |>
  semi_join(select(dat, source_row), by = "source_row") |>
  filter(clutch == 2) |>
  summarise(
    records = n(), seasons = n_distinct(year),
    females = n_distinct(female_id, na.rm = TRUE),
    fledged_observed = sum(!is.na(fledged)),
    zero_fledged = sum(fledged == 0, na.rm = TRUE),
    median_fledged = median(fledged, na.rm = TRUE),
    median_rel_LD = median(rel_LD, na.rm = TRUE),
    min_rel_LD = min(rel_LD, na.rm = TRUE),
    max_rel_LD = max(rel_LD, na.rm = TRUE)
  )

dat_no_two <- dat |> filter(clutch > 2) |> droplevels()
model_constant_no_two <- update(model_constant_scale, data = dat_no_two)
model_scale_no_two <- update(model_location_scale, data = dat_no_two)
sens_raw <- anova(model_constant_no_two, model_scale_no_two)
sensitivity_comparison <- tibble(
  sample = "Clutch > 2", n = nrow(dat_no_two),
  AIC_constant = AIC(model_constant_no_two),
  AIC_location_scale = AIC(model_scale_no_two),
  delta_AIC_support = AIC(model_constant_no_two) - AIC(model_scale_no_two),
  chi_square = sens_raw$Chisq[2], chi_df = sens_raw$`Chi Df`[2],
  p = sens_raw$`Pr(>Chisq)`[2]
)

coef_no_two <- bind_rows(
  extract_component(model_scale_no_two, "cond"),
  extract_component(model_scale_no_two, "disp")
) |>
  select(component, term, estimate_no_two = estimate)

coefficient_sensitivity <- coefficients |>
  select(component, term, estimate_full = estimate) |>
  left_join(coef_no_two, by = c("component", "term")) |>
  mutate(change = estimate_no_two - estimate_full)

annual_clutch <- dat |>
  group_by(year = as.integer(as.character(year_f))) |>
  summarise(n_broods = n(), mean_clutch = mean(clutch),
            sd_clutch = sd(clutch), .groups = "drop")
annual_cor_test <- cor.test(annual_clutch$mean_clutch, annual_clutch$sd_clutch)
annual_mean_sd_correlation <- tibble(
  n_seasons = nrow(annual_clutch), r = unname(annual_cor_test$estimate),
  CI_low = annual_cor_test$conf.int[1], CI_high = annual_cor_test$conf.int[2],
  p = annual_cor_test$p.value
)

sample_size_coefficients <- bind_rows(
  extract_component(model_location_scale, "disp") |>
    transmute(term, estimate_primary = estimate, SE_primary = SE,
              p_primary = p),
  extract_component(model_scale_sample_size, "disp") |>
    transmute(term, estimate_n_adjusted = estimate, SE_n_adjusted = SE,
              p_n_adjusted = p)
) |>
  group_by(term) |>
  summarise(across(everything(), ~ first(.x[!is.na(.x)])), .groups = "drop") |>
  mutate(change = estimate_n_adjusted - estimate_primary)

location_results <- coefficients |>
  filter(component == "cond") |>
  select(term, estimate, SE, CI_low, CI_high, p)

# glmmTMB models Gaussian dispersion on the log-variance scale.
scale_results <- coefficients |>
  filter(component == "disp") |>
  mutate(
    variance_multiplier = exp(estimate),
    variance_change_percent = 100 * (variance_multiplier - 1),
    multiplier_CI_low = exp(CI_low),
    multiplier_CI_high = exp(CI_high)
  ) |>
  select(term, estimate_log_variance = estimate, SE, CI_low, CI_high, p,
         variance_multiplier, variance_change_percent,
         multiplier_CI_low, multiplier_CI_high)

model_list <- list(
  constant_scale = model_constant_scale,
  location_scale = model_location_scale,
  constant_scale_no_two = model_constant_no_two,
  location_scale_no_two = model_scale_no_two,
  location_scale_sample_size = model_scale_sample_size
)
convergence <- bind_rows(lapply(names(model_list), function(nm) {
  m <- model_list[[nm]]
  tibble(model = nm, optimizer_code = m$fit$convergence,
         positive_definite_hessian = isTRUE(m$sdr$pdHess),
         optimizer_message = m$fit$message)
}))

write_csv(dat, "data_derived/clutch_size_model_data.csv")
write_csv(dispersion, "tables/clutch_size_dispersion_diagnostic.csv")
write_csv(clutch_distribution, "tables/clutch_size_distribution.csv")
write_csv(comparison, "tables/clutch_size_scale_model_comparison.csv")
write_csv(coefficients, "tables/clutch_size_location_scale_coefficients.csv")
write_csv(repeat_summary, "tables/clutch_size_female_repeats.csv")
write_csv(random_effects, "tables/clutch_size_random_effect_sd.csv")
write_csv(clutch_two, "tables/clutch_size_two_diagnostic.csv")
write_csv(sensitivity_comparison, "tables/clutch_size_no_two_model_comparison.csv")
write_csv(coefficient_sensitivity, "tables/clutch_size_no_two_coefficient_sensitivity.csv")
write_csv(location_results, "tables/clutch_size_location_results.csv")
write_csv(scale_results, "tables/clutch_size_scale_results.csv")
write_csv(convergence, "tables/clutch_size_convergence.csv")
write_csv(annual_clutch, "data_derived/annual_clutch_summary.csv")
write_csv(annual_mean_sd_correlation,
          "tables/clutch_size_annual_mean_sd_correlation.csv")
write_csv(sample_size_coefficients,
          "tables/clutch_size_sample_size_scale_sensitivity.csv")
saveRDS(
  list(constant_scale = model_constant_scale,
       location_scale = model_location_scale,
       constant_scale_no_two = model_constant_no_two,
       location_scale_no_two = model_scale_no_two,
       location_scale_sample_size = model_scale_sample_size),
  "models/clutch_size_location_scale_glmmTMB.rds"
)

print(as.data.frame(dispersion), row.names = FALSE)
print(as.data.frame(comparison), row.names = FALSE)
print(as.data.frame(repeat_summary), row.names = FALSE)
print(as.data.frame(random_effects), row.names = FALSE)
print(as.data.frame(clutch_two), row.names = FALSE)
print(as.data.frame(sensitivity_comparison), row.names = FALSE)
