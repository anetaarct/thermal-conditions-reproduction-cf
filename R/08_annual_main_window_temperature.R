#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

# Environmental series is independent of breeding-field coverage.
expected_years <- 1980:2025

temperature <- read_csv("data_derived/smhi_hoburg_daily_temperature.csv",
                        show_col_types = FALSE) |>
  transmute(date = as.Date(date), temperature = as.numeric(temperature))

annual <- temperature |>
  mutate(year = as.integer(format(date, "%Y")),
         ref = as.Date(sprintf("%d-05-07", year))) |>
  filter(date >= ref - 42, date <= ref) |>
  group_by(year) |>
  summarise(temp_mean = mean(temperature), temp_sd = sd(temperature),
            n_days = n(), window_start = min(date), window_end = max(date),
            .groups = "drop") |>
  filter(year %in% expected_years)

stopifnot(identical(annual$year, expected_years),
          all(annual$n_days == 43L), !anyNA(annual))

mean_trend <- lm(temp_mean ~ year, data = annual)
sd_trend <- lm(temp_sd ~ year, data = annual)

trend_row <- function(response, model) {
  cf <- coef(summary(model))["year", ]
  tibble(
    response = response,
    beta_per_year = unname(cf["Estimate"]),
    SE = unname(cf["Std. Error"]),
    CI_low = beta_per_year - 1.96 * SE,
    CI_high = beta_per_year + 1.96 * SE,
    p = unname(cf["Pr(>|t|)"]),
    R_squared = summary(model)$r.squared
  )
}

trends <- bind_rows(
  trend_row("temp_mean", mean_trend),
  trend_row("temp_sd", sd_trend)
)

correlation_row <- function(variable_x, variable_y, x, y) {
  test <- cor.test(x, y, method = "pearson")
  tibble(
    variable_x = variable_x,
    variable_y = variable_y,
    r = unname(test$estimate),
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2],
    p = test$p.value,
    abs_r_above_0_6 = abs(r) > 0.6
  )
}

correlations <- bind_rows(
  correlation_row("temp_mean", "year", annual$temp_mean, annual$year),
  correlation_row("temp_sd", "year", annual$temp_sd, annual$year),
  correlation_row("temp_mean", "temp_sd", annual$temp_mean, annual$temp_sd)
)

# Check the annual centring used for relative laying date.
relative_laying <- read_csv("data_derived/breeding_reproduction_clean.csv",
                            show_col_types = FALSE) |>
  transmute(year = as.integer(YEAR), LD = as.numeric(LD)) |>
  filter(year %in% expected_years, !is.na(LD)) |>
  group_by(year) |>
  mutate(rel_LD = LD - mean(LD)) |>
  ungroup() |>
  left_join(select(annual, year, temp_mean), by = "year")

correlations <- bind_rows(
  correlations,
  correlation_row("rel_LD", "temp_mean",
                  relative_laying$rel_LD, relative_laying$temp_mean)
)

temperature_range <- tibble(
  coldest_year = annual$year[which.min(annual$temp_mean)],
  coldest_temp_mean = min(annual$temp_mean),
  warmest_year = annual$year[which.max(annual$temp_mean)],
  warmest_temp_mean = max(annual$temp_mean),
  difference_C = warmest_temp_mean - coldest_temp_mean
)

write_csv(annual, "data_derived/annual_main_window_temperature.csv")
write_csv(trends, "tables/main_window_temperature_trends.csv")
write_csv(correlations, "tables/main_window_temperature_correlations.csv")
write_csv(temperature_range, "tables/main_window_temperature_range.csv")
saveRDS(list(mean_trend = mean_trend, sd_trend = sd_trend),
        "models/climwin/main_window_temperature_trends.rds")

print(as.data.frame(trends), row.names = FALSE)
print(as.data.frame(correlations), row.names = FALSE)
