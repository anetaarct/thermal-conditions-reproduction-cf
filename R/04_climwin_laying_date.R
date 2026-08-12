#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(climwin)
})

# Force English month abbreviations/names in all exported calendar labels.
Sys.setlocale("LC_TIME", "C")

dir.create("models/climwin", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# Fixed, preregistered analysis settings.
set.seed(20260808)
search_range <- c(17, 0)

birds_all <- read_csv(
  "data_raw/birds/DB_MISMATCH_current_onset_1980-2025.csv",
  show_col_types = FALSE
) |>
  transmute(
    year = as.integer(YEAR),
    lay_date = as.numeric(LD),
    laying_calendar = as.Date(laying_date, format = "%d/%m/%Y")
  ) |>
  filter(year >= 1980, year <= 2025, !is.na(laying_calendar)) |>
  mutate(lay_date_from_may = as.numeric(
    laying_calendar - as.Date(sprintf("%d-05-01", year))
  ))

stopifnot(all(between(birds_all$lay_date_from_may, -10, 60)))
stopifnot(all(birds_all$lay_date == birds_all$lay_date_from_may,
              na.rm = TRUE))

# The response uses every valid laying date. Age is retained in the raw source
# but is not used to filter or adjust the annual population median.
annual <- birds_all |>
  group_by(year) |>
  summarise(
    n_broods = n(),
    median_LD = median(lay_date_from_may),
    p05_lay_from_may = as.numeric(quantile(lay_date_from_may, 0.05,
                                           names = FALSE, type = 7)),
    .groups = "drop"
  ) |>
  arrange(year) |>
  mutate(
    season_index = row_number(),
    bdate = as.Date(sprintf("%d-05-01", year)) + round(median_LD),
    p05_lay_calendar = as.Date(sprintf("%d-05-01", year)) + p05_lay_from_may
  )

expected_years <- c(1980:2019, 2021:2025)
stopifnot(nrow(annual) == 45L, identical(annual$year, expected_years))
n_first <- 23L

temperature <- read_csv("data_derived/smhi_hoburg_daily_temperature.csv",
                               show_col_types = FALSE) |>
  transmute(date = as.Date(date), temperature = as.numeric(temperature)) |>
  filter(!is.na(date), !is.na(temperature))

stopifnot(!anyDuplicated(temperature$date))
stopifnot(all(diff(temperature$date) == 1))

# Keep the descriptive population median separate from the fixed reference day
# used to anchor all absolute climate windows.
ref_ld <- median(birds_all$lay_date)
refday <- c(7L, 5L)
ref_date_template <- as.Date("2001-05-07")

run_sliding <- function(dat) {
  baseline <- lm(median_LD ~ 1, data = dat)
  slidingwin(
    xvar = list(temperature = temperature$temperature),
    cdate = temperature$date,
    bdate = dat$bdate,
    baseline = baseline,
    type = "absolute",
    refday = refday,
    range = search_range,
    stat = "mean",
    func = "lin",
    cinterval = "week"
  )
}

run_random <- function(dat, repeats = 500L) {
  baseline <- lm(median_LD ~ 1, data = dat)
  randwin(
    repeats = repeats,
    xvar = list(temperature = temperature$temperature),
    cdate = temperature$date,
    bdate = dat$bdate,
    baseline = baseline,
    type = "absolute",
    refday = refday,
    range = search_range,
    stat = "mean",
    func = "lin",
    cinterval = "week"
  )
}

main_win <- run_sliding(annual)
first_win <- run_sliding(filter(annual, season_index <= n_first))
second_win <- run_sliding(filter(annual, season_index > n_first))

saveRDS(main_win, "models/climwin/slidingwin_main_calendar_index_ref07may_range17w.rds")
saveRDS(first_win, "models/climwin/slidingwin_first23_calendar_index_ref07may_range17w.rds")
saveRDS(second_win, "models/climwin/slidingwin_second22_calendar_index_ref07may_range17w.rds")

# This is intentionally expensive and is cached as an RDS object.
rand_path <- "models/climwin/randwin_main_calendar_index_ref07may_range17w_500.rds"
if (file.exists(rand_path)) {
  main_rand <- readRDS(rand_path)
} else {
  main_rand <- run_random(annual, repeats = 500L)
  saveRDS(main_rand, rand_path)
}

stopifnot(nrow(main_rand[[1]]) == 500L)

main_dataset <- main_win[[1]]$Dataset
first_dataset <- first_win[[1]]$Dataset
second_dataset <- second_win[[1]]$Dataset

best_row <- function(x, minimum_length_weeks = 0) {
  x |>
    filter(WindowOpen - WindowClose >= minimum_length_weeks) |>
    slice_min(deltaAICc, n = 1, with_ties = FALSE)
}
main_best <- best_row(main_dataset, minimum_length_weeks = 2)
first_best <- best_row(first_dataset, minimum_length_weeks = 2)
second_best <- best_row(second_dataset, minimum_length_weeks = 2)

# Both false-positive diagnostics required by climwin.
p_delta_aicc <- pvalue(main_dataset, main_rand[[1]], metric = "AIC",
                       sample.size = nrow(annual))
p_c <- pvalue(main_dataset, main_rand[[1]], metric = "C",
              sample.size = nrow(annual))

main_med <- medwin(main_dataset, cw = 0.95)

calendar_label <- function(days_before) {
  format(ref_date_template - as.integer(round(days_before)), "%d %b")
}

window_row <- function(label, best, years) {
  tibble(
    analysis = label,
    years = years,
    window_open_weeks_before_refday = best$WindowOpen,
    window_close_weeks_before_refday = best$WindowClose,
    boundary_difference_weeks = best$WindowOpen - best$WindowClose,
    window_open_days_before_refday = best$WindowOpen * 7,
    window_close_days_before_refday = best$WindowClose * 7,
    window_open_calendar = calendar_label(best$WindowOpen * 7),
    window_close_calendar = calendar_label(best$WindowClose * 7),
    delta_AICc = best$deltaAICc
  )
}

windows <- bind_rows(
  window_row("Main", main_best, "1980-2025 (no 2020)"),
  window_row("First 23 seasons", first_best, "1980-2002"),
  window_row("Remaining 22 seasons", second_best, "2003-2025 (no 2020)")
)

# Out-of-sample evaluation: singlewin applies exactly the same weekly
# aggregation as the training slidingwin, without a new window search.
test_data <- annual |>
  filter(season_index > n_first)

oos_win <- singlewin(
  xvar = list(temperature = temperature$temperature),
  cdate = temperature$date,
  bdate = test_data$bdate,
  baseline = lm(median_LD ~ 1, data = test_data),
  type = "absolute",
  refday = refday,
  range = c(first_best$WindowOpen, first_best$WindowClose),
  stat = "mean",
  func = "lin",
  cinterval = "week"
  )

oos_model <- oos_win$BestModel
oos_summary <- summary(oos_model)
oos_coef <- coef(oos_summary)["climate", ]
oos <- tibble(
  training_years = "1980-2002 (23 seasons)",
  test_years = "2003-2025 (22 seasons; no 2020)",
  selected_open_days_before_refday = first_best$WindowOpen * 7,
  selected_close_days_before_refday = first_best$WindowClose * 7,
  temperature_beta = unname(oos_coef["Estimate"]),
  temperature_SE = unname(oos_coef["Std. Error"]),
  temperature_CI_low = temperature_beta - 1.96 * temperature_SE,
  temperature_CI_high = temperature_beta + 1.96 * temperature_SE,
  test_R_squared = oos_summary$r.squared,
  test_p = coef(oos_summary)["climate", "Pr(>|t|)"]
)

# medwin returns the median opening and closing boundaries of the 95% model set.
med_values <- as.numeric(unlist(main_med)) * 7
med_result <- tibble(
  confidence_set = "95%",
  median_window_open = med_values[1],
  median_window_close = med_values[2],
  median_open_calendar = calendar_label(med_values[1]),
  median_close_calendar = calendar_label(med_values[2])
)

false_positive <- tibble(
  randwin_repeats = 500L,
  P_delta_AICc = as.character(p_delta_aicc),
  P_c = as.numeric(p_c)
)

effect_row_from_model <- function(analysis, model) {
  cf <- coef(summary(model))["climate", ]
  tibble(
    analysis = analysis,
    temperature_beta = unname(cf["Estimate"]),
    temperature_SE = unname(cf["Std. Error"]),
    temperature_CI_low = temperature_beta - 1.96 * temperature_SE,
    temperature_CI_high = temperature_beta + 1.96 * temperature_SE,
    p = unname(cf["Pr(>|t|)"])
  )
}

effect_row_from_dataset <- function(analysis, best) {
  beta <- best$ModelBeta
  se <- best$Std.Error
  t_value <- beta / se
  tibble(
    analysis = analysis,
    temperature_beta = beta,
    temperature_SE = se,
    temperature_CI_low = beta - 1.96 * se,
    temperature_CI_high = beta + 1.96 * se,
    p = 2 * pt(abs(t_value), df = best$sample.size - 2, lower.tail = FALSE)
  )
}

temperature_effects <- bind_rows(
  effect_row_from_dataset("Main: full dataset window", main_best),
  effect_row_from_dataset("First 23 seasons: selected window", first_best),
  effect_row_from_dataset("Remaining 22 seasons: selected window", second_best),
  effect_row_from_model("OOS: first-period window tested in remaining seasons", oos_model)
)

main_close <- as.integer(main_best$WindowClose) * 7
overlap_check <- annual |>
  mutate(
    window_close_calendar = as.Date(sprintf("%d-%02d-%02d", year,
                                             refday[2], refday[1])) - main_close,
    window_reaches_p05_lay = window_close_calendar >= p05_lay_calendar
  ) |>
  select(year, n_broods, p05_lay_from_may, p05_lay_calendar,
         window_close_calendar, window_reaches_p05_lay)

settings <- tibble(
  setting = c("analysis level", "analysis years", "response", "age filtering or adjustment",
              "minimum window length in all window selections",
              "type", "reference laying date (database LD index; LD = 0 is 1 May)",
              "reference calendar day", "range",
              "stat", "func", "cinterval", "baseline", "randwin repeats"),
  value = c("year (season)", "1980-2019 and 2021-2025 (45 seasons)",
            "annual raw median laying date from all valid broods",
            "none", "boundary difference at least 2 weeks", "absolute",
            as.character(ref_ld), format(ref_date_template, "%d %B"),
            "17 to 0 weeks (119 to 0 days)",
            "mean", "lin", "week", "median_LD ~ 1", "500")
)

write_csv(annual, "data_derived/annual_raw_median_laying_date.csv")
write_csv(windows, "tables/climwin_windows.csv")
write_csv(false_positive, "tables/climwin_false_positive_tests.csv")
write_csv(temperature_effects, "tables/climwin_temperature_effects.csv")
write_csv(oos, "tables/climwin_out_of_sample.csv")
write_csv(med_result, "tables/climwin_medwin.csv")
write_csv(overlap_check, "tables/climwin_window_lay_overlap_by_year.csv")
write_csv(settings, "tables/climwin_settings.csv")
saveRDS(oos_model, "models/climwin/out_of_sample_model.rds")
saveRDS(oos_win, "models/climwin/out_of_sample_singlewin.rds")

print(as.data.frame(windows), row.names = FALSE)
print(as.data.frame(false_positive), row.names = FALSE)
print(as.data.frame(oos), row.names = FALSE)
print(as.data.frame(med_result), row.names = FALSE)
