#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

dat <- read_csv("data_derived/reproductive_model_data.csv",
                show_col_types = FALSE) |>
  filter(year >= 1981, year <= 2019, !is.na(fledged), fledged > 0) |>
  mutate(missing_age = is.na(age_class), missing_rel_LD = is.na(rel_LD),
         model_complete = !missing_age & !missing_rel_LD)

missing_patterns <- dat |>
  count(missing_age, missing_rel_LD, name = "records") |>
  mutate(proportion = records / nrow(dat))

missing_summary <- tibble(
  quantity = c("All positive fledging records", "Complete model records",
               "Excluded from model", "Missing age class",
               "Missing rel_LD", "Missing both age and rel_LD"),
  value = c(nrow(dat), sum(dat$model_complete), sum(!dat$model_complete),
            sum(dat$missing_age), sum(dat$missing_rel_LD),
            sum(dat$missing_age & dat$missing_rel_LD))
)

timing_by_age_availability <- dat |>
  filter(!is.na(rel_LD)) |>
  group_by(age_status = if_else(missing_age, "Age missing", "Age available")) |>
  summarise(records = n(), mean_rel_LD = mean(rel_LD),
            median_rel_LD = median(rel_LD),
            p05_rel_LD = quantile(rel_LD, 0.05),
            p95_rel_LD = quantile(rel_LD, 0.95), .groups = "drop")

age_missing_model <- glm(missing_age ~ rel_LD, data = filter(dat, !is.na(rel_LD)),
                         family = binomial())
cf <- coef(summary(age_missing_model))["rel_LD", ]
age_missing_timing_test <- tibble(
  term = "rel_LD", log_odds_beta = unname(cf["Estimate"]),
  SE = unname(cf["Std. Error"]), z = unname(cf["z value"]),
  p = unname(cf["Pr(>|z|)"]), odds_ratio = exp(log_odds_beta),
  OR_CI_low = exp(log_odds_beta - 1.96 * SE),
  OR_CI_high = exp(log_odds_beta + 1.96 * SE)
)

missing_rel_by_year <- dat |>
  group_by(year) |>
  summarise(records = n(), missing_rel_LD = sum(missing_rel_LD),
            proportion_missing_rel_LD = mean(missing_rel_LD), .groups = "drop")

write_csv(missing_summary, "tables/fledging_positive_missing_summary.csv")
write_csv(missing_patterns, "tables/fledging_positive_missing_patterns.csv")
write_csv(timing_by_age_availability,
          "tables/fledging_positive_timing_by_age_availability.csv")
write_csv(age_missing_timing_test,
          "tables/fledging_positive_age_missing_timing_test.csv")
write_csv(missing_rel_by_year,
          "tables/fledging_positive_missing_rel_by_year.csv")

print(as.data.frame(missing_summary), row.names = FALSE)
print(as.data.frame(missing_patterns), row.names = FALSE)
print(as.data.frame(timing_by_age_availability), row.names = FALSE)
print(as.data.frame(age_missing_timing_test), row.names = FALSE)
