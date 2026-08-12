#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble); library(climwin)})

annual <- read_csv("data_derived/annual_raw_median_laying_date.csv",
                   show_col_types = FALSE) |>
  arrange(year) |>
  mutate(bdate = as.Date(bdate), year_c = year - mean(year))

windows <- read_csv("tables/climwin_windows.csv", show_col_types = FALSE)
main <- windows |> filter(analysis == "Main")
stopifnot(nrow(main) == 1L)

main_object <- readRDS(
  "models/climwin/slidingwin_main_calendar_index_ref07may_range17w.rds"
)
best_dataset_row <- main_object[[1]]$Dataset[1, ]
stopifnot(best_dataset_row$WindowOpen == main$window_open_weeks_before_refday,
          best_dataset_row$WindowClose == main$window_close_weeks_before_refday)
fixed_data <- main_object[[1]]$BestModelData
stopifnot(nrow(fixed_data) == nrow(annual))
stopifnot(isTRUE(all.equal(fixed_data$yvar, annual$median_LD,
                           tolerance = 1e-10)))

model_data <- annual |>
  mutate(temp_main = fixed_data$climate)

model <- lm(median_LD ~ temp_main * year_c, data = model_data)
cf <- coef(summary(model))
ci <- confint(model)
coefficients <- tibble(
  term = rownames(cf), estimate = cf[, "Estimate"], SE = cf[, "Std. Error"],
  CI_low = ci[, 1], CI_high = ci[, 2], t = cf[, "t value"],
  p = cf[, "Pr(>|t|)"]
)

# Conditional temperature slopes at the first, mean, and last observed years.
b <- coef(model); V <- vcov(model)
slope_at_year <- function(target_year) {
  yc <- target_year - mean(model_data$year)
  weights <- c(temp_main = 1, `temp_main:year_c` = yc)
  estimate <- sum(weights * b[names(weights)])
  se <- sqrt(as.numeric(t(weights) %*% V[names(weights), names(weights)] %*% weights))
  tibble(year = target_year, year_c = yc, temperature_slope = estimate, SE = se,
         CI_low = estimate - 1.96 * se, CI_high = estimate + 1.96 * se)
}

slopes <- bind_rows(slope_at_year(min(model_data$year)),
                    slope_at_year(round(mean(model_data$year))),
                    slope_at_year(max(model_data$year)))

summary_table <- tibble(
  n_seasons = nrow(model_data),
  main_window_weeks = paste0(main$window_open_weeks_before_refday, "-",
                             main$window_close_weeks_before_refday),
  R_squared = summary(model)$r.squared,
  adjusted_R_squared = summary(model)$adj.r.squared
)

write_csv(model_data, "data_derived/main_window_time_interaction_data.csv")
write_csv(coefficients, "tables/main_window_time_interaction_coefficients.csv")
write_csv(slopes, "tables/main_window_time_conditional_slopes.csv")
write_csv(summary_table, "tables/main_window_time_interaction_summary.csv")
saveRDS(model, "models/climwin/main_window_time_interaction_model.rds")

print(as.data.frame(coefficients), row.names = FALSE)
print(as.data.frame(slopes), row.names = FALSE)
print(as.data.frame(summary_table), row.names = FALSE)
