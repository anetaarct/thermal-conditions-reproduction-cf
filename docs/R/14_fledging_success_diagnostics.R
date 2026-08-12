#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

dat <- read_csv("data_derived/reproductive_model_data.csv",
                show_col_types = FALSE) |>
  filter(!is.na(fledged))

annual <- dat |>
  group_by(year) |>
  summarise(n_broods = n(), zero_broods = sum(fledged == 0),
            p_zero = mean(fledged == 0), mean_fledged = mean(fledged),
            variance_fledged = var(fledged), .groups = "drop")

zero_trend <- lm(p_zero ~ year, data = annual)
cf <- coef(summary(zero_trend))["year", ]
zero_trend_result <- tibble(
  seasons = nrow(annual), beta_per_year = unname(cf["Estimate"]),
  SE = unname(cf["Std. Error"]),
  CI_low = beta_per_year - 1.96 * SE,
  CI_high = beta_per_year + 1.96 * SE,
  p = unname(cf["Pr(>|t|)"]), R_squared = summary(zero_trend)$r.squared
)

positive <- dat |> filter(fledged > 0)
diagnostics <- tibble(
  quantity = c("Broods with fledged observed", "Broods with fledged = 0",
               "Overall proportion fledged = 0", "Positive fledging records",
               "Mean fledged | fledged > 0", "Variance fledged | fledged > 0",
               "Variance/mean | fledged > 0"),
  value = c(nrow(dat), sum(dat$fledged == 0), mean(dat$fledged == 0),
            nrow(positive), mean(positive$fledged), var(positive$fledged),
            var(positive$fledged) / mean(positive$fledged))
)

write_csv(annual, "data_derived/annual_fledging_zero_summary.csv")
write_csv(diagnostics, "tables/fledging_success_diagnostics.csv")
write_csv(zero_trend_result, "tables/fledging_zero_trend.csv")
saveRDS(zero_trend, "models/fledging_zero_trend_lm.rds")

print(as.data.frame(diagnostics), row.names = FALSE)
print(as.data.frame(zero_trend_result), row.names = FALSE)
