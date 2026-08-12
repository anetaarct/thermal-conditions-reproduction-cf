#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})

# 1980 was a test season and is excluded from reproductive modelling.
expected_years <- c(1981:2019, 2021:2025)

norm_ring <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NAN")] <- NA_character_
  x
}

temperature <- read_csv("data_derived/annual_main_window_temperature.csv",
                        show_col_types = FALSE) |>
  select(year, temp_mean, temp_sd)

raw <- read_csv(
  "data_raw/birds/DB_MISMATCH_current_onset_1980-2025.csv",
  col_types = cols(.default = col_character()), show_col_types = FALSE
) |>
  transmute(
    source_row = row_number(),
    female_id = norm_ring(FRING),
    year = as.integer(YEAR),
    laying_date = as.numeric(LD),
    clutch = as.numeric(CS),
    fledged = as.numeric(FLEDGE),
    age_class_raw = trimws(FAGE_CAT)
  ) |>
  filter(year %in% expected_years) |>
  mutate(
    age_class = case_when(
      tolower(age_class_raw) == "young" ~ "YOUNG",
      tolower(age_class_raw) == "old" ~ "OLD",
      TRUE ~ NA_character_
    ),
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    age_source = if_else(!is.na(age_class), "FAGE_CAT", NA_character_),
    valid_clutch = !is.na(clutch) & between(clutch, 2, 9),
    valid_fledged = is.na(fledged) |
      (!is.na(clutch) & fledged >= 0 & fledged <= clutch),
    zero_fledged = !is.na(fledged) & fledged == 0
  )

excluded <- raw |> filter(!valid_clutch | !valid_fledged)

model_data <- raw |>
  filter(valid_clutch, valid_fledged) |>
  group_by(year) |>
  mutate(
    annual_median_LD = median(laying_date, na.rm = TRUE),
    rel_LD = laying_date - annual_median_LD
  ) |>
  ungroup() |>
  mutate(
    year_f = factor(year, levels = expected_years),
    year_c = year - mean(expected_years),
    age_available = !is.na(age_class)
  ) |>
  left_join(temperature, by = "year") |>
  select(source_row, female_id, year, year_f, year_c, laying_date,
         annual_median_LD, rel_LD, clutch, fledged, zero_fledged,
         age_class, age_source, age_available, temp_mean, temp_sd)

stopifnot(
  setequal(unique(model_data$year), expected_years),
  all(between(model_data$clutch, 2, 9)),
  all(is.na(model_data$fledged) |
        (model_data$fledged >= 0 & model_data$fledged <= model_data$clutch)),
  !anyNA(model_data$temp_mean), !anyNA(model_data$temp_sd),
  all(model_data$year == as.integer(as.character(model_data$year_f))),
  all(is.na(model_data$age_class) | model_data$age_source == "FAGE_CAT"),
  all(is.na(model_data$age_source) == is.na(model_data$age_class))
)

# With median centring, the annual median of rel_LD is exactly zero. The
# individual-level correlation need not be algebraically zero, so both checks
# are retained and reported.
centering_by_year <- model_data |>
  filter(!is.na(rel_LD)) |>
  group_by(year) |>
  summarise(median_rel_LD = median(rel_LD), temp_mean = first(temp_mean),
            .groups = "drop")

centering_checks <- tibble(
  check = c("Individual-level cor(rel_LD, temp_mean)",
            "Maximum absolute annual median(rel_LD)"),
  value = c(
    cor(model_data$rel_LD, model_data$temp_mean, use = "complete.obs"),
    max(abs(centering_by_year$median_rel_LD))
  )
)

audit <- tibble(
  quantity = c(
    "Source records in 44 non-test seasons", "Analysis-ready clutch records",
    "Excluded: clutch outside 2-9 or missing",
    "Excluded: fledged below 0 or above clutch",
    "Records with fledged observed", "Records with zero fledged retained",
    "Records with female age class", "Records without female age class",
    "Unique females among records with a ring"
  ),
  value = c(
    nrow(raw), nrow(model_data), sum(!raw$valid_clutch),
    sum(raw$valid_clutch & !raw$valid_fledged), sum(!is.na(model_data$fledged)),
    sum(model_data$zero_fledged), sum(model_data$age_available),
    sum(!model_data$age_available), n_distinct(model_data$female_id, na.rm = TRUE)
  )
)

age_coverage <- model_data |>
  mutate(age_class = as.character(age_class),
         age_class = if_else(is.na(age_class), "MISSING", age_class)) |>
  count(year, age_class, name = "records") |>
  group_by(year) |>
  mutate(proportion = records / sum(records)) |>
  ungroup()

write_csv(model_data, "data_derived/reproductive_model_data.csv", na = "")
write_csv(excluded, "tables/reproductive_model_exclusions.csv", na = "")
write_csv(audit, "tables/reproductive_model_data_audit.csv")
write_csv(age_coverage, "tables/reproductive_model_age_coverage.csv")
write_csv(centering_checks, "tables/reproductive_model_centering_checks.csv")

print(as.data.frame(audit), row.names = FALSE)
print(as.data.frame(centering_checks), row.names = FALSE)
