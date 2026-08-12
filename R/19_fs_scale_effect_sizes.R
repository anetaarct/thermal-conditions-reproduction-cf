#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

annual <- read_csv("data_derived/annual_main_window_temperature.csv",
                   show_col_types = FALSE) |>
  filter(year >= 1981, year <= 2019)
post <- read_csv("tables/brms_fs_full_posterior_95.csv",
                 show_col_types = FALSE)

temp_range <- max(annual$temp_mean) - min(annual$temp_mean)
year_span <- max(annual$year) - min(annual$year)

transform_effect <- function(parameter, contrast, label) {
  x <- post |> filter(.data$parameter == parameter)
  stopifnot(nrow(x) == 1L)
  tibble(
    contrast = label, contrast_value = contrast,
    log_sigma_beta = x$estimate,
    SD_ratio = exp(x$estimate * contrast),
    SD_change_percent = 100 * (SD_ratio - 1),
    SD_ratio_CrI_low = exp(x$CrI_2.5 * contrast),
    SD_ratio_CrI_high = exp(x$CrI_97.5 * contrast),
    SD_change_CrI_low_percent = 100 * (SD_ratio_CrI_low - 1),
    SD_change_CrI_high_percent = 100 * (SD_ratio_CrI_high - 1)
  )
}

effects <- bind_rows(
  transform_effect("b_sigma_temp_mean_c", temp_range,
                   "Warmest versus coldest observed spring"),
  transform_effect("b_sigma_year_c", year_span,
                   "2019 versus 1981 temporal contrast")
)

range_summary <- tibble(
  coldest_year = annual$year[which.min(annual$temp_mean)],
  coldest_temp_mean = min(annual$temp_mean),
  warmest_year = annual$year[which.max(annual$temp_mean)],
  warmest_temp_mean = max(annual$temp_mean),
  observed_range_C = temp_range,
  temporal_span_years = year_span
)

write_csv(effects, "tables/brms_fs_scale_observed_contrasts.csv")
write_csv(range_summary, "tables/brms_fs_scale_contrast_ranges.csv")
print(as.data.frame(range_summary), row.names = FALSE)
print(as.data.frame(effects), row.names = FALSE)
