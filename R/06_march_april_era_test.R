#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

annual <- read_csv("data_derived/annual_raw_median_laying_date.csv",
                   show_col_types = FALSE) |>
  select(year, median_LD)

monthly <- read_csv("data_derived/march_april_temperature_by_year.csv",
                    show_col_types = FALSE) |>
  select(year, month, mean_temperature) |>
  tidyr::pivot_wider(names_from = month, values_from = mean_temperature,
                     names_prefix = "month_") |>
  rename(temp_march = month_3, temp_april = month_4)

annual_model <- annual |>
  inner_join(monthly, by = "year") |>
  mutate(era = factor(if_else(year <= 2002, "1980-2002", "2003-2025"),
                      levels = c("1980-2002", "2003-2025")))

stopifnot(nrow(annual_model) == 45L, !anyNA(annual_model))

model <- lm(median_LD ~ temp_march * era + temp_april * era,
            data = annual_model)
model_no_interactions <- lm(median_LD ~ temp_march + temp_april + era,
                            data = annual_model)

cf <- coef(summary(model))
ci <- confint(model)
coefficients <- tibble(
  term = rownames(cf), estimate = cf[, "Estimate"], SE = cf[, "Std. Error"],
  CI_low = ci[, 1], CI_high = ci[, 2], t = cf[, "t value"],
  p = cf[, "Pr(>|t|)"]
)

b <- coef(model)
V <- vcov(model)
slope_row <- function(variable, era_label, weights) {
  estimate <- sum(weights * b[names(weights)])
  variance <- as.numeric(t(weights) %*% V[names(weights), names(weights)] %*% weights)
  se <- sqrt(variance)
  tibble(variable = variable, era = era_label, estimate = estimate, SE = se,
         CI_low = estimate - 1.96 * se, CI_high = estimate + 1.96 * se)
}

march_early <- c(temp_march = 1)
april_early <- c(temp_april = 1)
march_interaction <- grep("temp_march.*era|era.*temp_march", names(b), value = TRUE)
april_interaction <- grep("temp_april.*era|era.*temp_april", names(b), value = TRUE)
stopifnot(length(march_interaction) == 1L, length(april_interaction) == 1L)
march_late <- setNames(c(1, 1), c("temp_march", march_interaction))
april_late <- setNames(c(1, 1), c("temp_april", april_interaction))

slopes <- bind_rows(
  slope_row("March temperature", "1980-2002", march_early),
  slope_row("March temperature", "2003-2025", march_late),
  slope_row("April temperature", "1980-2002", april_early),
  slope_row("April temperature", "2003-2025", april_late)
)

joint <- anova(model_no_interactions, model)
joint_test <- tibble(
  comparison = "both temperature-by-era interactions",
  numerator_df = joint$Df[2], denominator_df = df.residual(model),
  F = joint$F[2], p = joint$`Pr(>F)`[2],
  R_squared_full = summary(model)$r.squared,
  adjusted_R_squared_full = summary(model)$adj.r.squared
)

write_csv(annual_model, "data_derived/annual_march_april_era_model_data.csv")
write_csv(coefficients, "tables/march_april_era_model_coefficients.csv")
write_csv(slopes, "tables/march_april_era_slopes.csv")
write_csv(joint_test, "tables/march_april_era_joint_interaction_test.csv")
saveRDS(model, "models/climwin/march_april_era_interaction_model.rds")

print(as.data.frame(coefficients), row.names = FALSE)
print(as.data.frame(slopes), row.names = FALSE)
print(as.data.frame(joint_test), row.names = FALSE)
