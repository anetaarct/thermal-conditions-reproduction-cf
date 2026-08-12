#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
warm_threshold_c <- 8
threshold_grid_c <- c(4, 5, 6, 8)

temperature <- read_csv("data_derived/smhi_hoburg_daily_temperature.csv",
                        show_col_types = FALSE) |>
  transmute(date = as.Date(date), temperature = as.numeric(temperature)) |>
  filter(!is.na(date), !is.na(temperature)) |>
  mutate(year = as.integer(format(date, "%Y")),
         month = as.integer(format(date, "%m")))

stopifnot(!anyDuplicated(temperature$date))
stopifnot(all(diff(temperature$date) == 1))

monthly <- temperature |>
  filter(month %in% c(3L, 4L)) |>
  group_by(year, month) |>
  summarise(
    mean_temperature = mean(temperature),
    sd_temperature = sd(temperature),
    warm_days_above_8c = sum(temperature > warm_threshold_c),
    maximum_temperature = max(temperature),
    n_days = n(),
    .groups = "drop"
  )

march <- monthly |>
  filter(month == 3L) |>
  transmute(year, mean_mar = mean_temperature, sd_mar = sd_temperature,
            warm_days_above_8c, max_mar = maximum_temperature, n_days_mar = n_days)

april <- monthly |>
  filter(month == 4L) |>
  transmute(year, mean_apr = mean_temperature, sd_apr = sd_temperature,
            max_apr = maximum_temperature, n_days_apr = n_days)

threshold_counts <- temperature |>
  filter(month == 3L) |>
  tidyr::crossing(threshold_c = threshold_grid_c) |>
  group_by(year, threshold_c) |>
  summarise(warm_days = sum(temperature > threshold_c), .groups = "drop")

warm_result <- threshold_counts |>
  group_by(threshold_c) |>
  group_modify(~{
    if (sd(.x$warm_days) == 0) {
      return(tibble(n_years = nrow(.x), total_warm_days = sum(.x$warm_days),
                    beta_days_per_year = NA_real_, SE = NA_real_,
                    CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
                    R_squared = NA_real_))
    }
    fit <- lm(warm_days ~ year, data = .x)
    cf <- coef(summary(fit))["year", ]
    tibble(n_years = nrow(.x), total_warm_days = sum(.x$warm_days),
           beta_days_per_year = unname(cf["Estimate"]),
           SE = unname(cf["Std. Error"]),
           CI_low = beta_days_per_year - 1.96 * SE,
           CI_high = beta_days_per_year + 1.96 * SE,
           p = unname(cf["Pr(>|t|)"]), R_squared = summary(fit)$r.squared)
  }) |>
  ungroup()

joined <- inner_join(march, april, by = "year") |>
  mutate(period = if_else(year <= 2002, "1980-2002", "2003-2025"))

correlation_row <- function(dat, label, matched_laying_years) {
  test <- cor.test(dat$mean_mar, dat$mean_apr, method = "pearson")
  tibble(
    sample = label,
    excludes_2020 = matched_laying_years,
    n_years = nrow(dat),
    r_march_april = unname(test$estimate),
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2],
    p = test$p.value
  )
}

correlations <- bind_rows(
  joined |> group_by(period) |> group_modify(~correlation_row(.x, "all climate years", FALSE)) |>
    ungroup(),
  joined |> filter(year != 2020) |> group_by(period) |>
    group_modify(~correlation_row(.x, "laying-data years", TRUE)) |> ungroup()
)

matched_cor <- correlations |> filter(sample == "laying-data years")
z1 <- atanh(matched_cor$r_march_april[matched_cor$period == "1980-2002"])
z2 <- atanh(matched_cor$r_march_april[matched_cor$period == "2003-2025"])
n1 <- matched_cor$n_years[matched_cor$period == "1980-2002"]
n2 <- matched_cor$n_years[matched_cor$period == "2003-2025"]
z_difference <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
correlation_difference <- tibble(
  comparison = "1980-2002 versus 2003-2025 (laying-data years)",
  fisher_z_difference = z_difference,
  p_two_sided = 2 * pnorm(abs(z_difference), lower.tail = FALSE)
)

write_csv(monthly, "data_derived/march_april_temperature_by_year.csv")
write_csv(march, "tables/march_temperature_diagnostics.csv")
write_csv(warm_result, "tables/march_warm_days_trend.csv")
write_csv(correlations, "tables/march_april_temperature_correlations.csv")
write_csv(correlation_difference, "tables/march_april_correlation_difference.csv")

print(as.data.frame(warm_result), row.names = FALSE)
print(as.data.frame(correlations), row.names = FALSE)
print(as.data.frame(correlation_difference), row.names = FALSE)
