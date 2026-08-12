#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

data_file <- "data_derived/fledging_model_data.csv"
dat <- read_csv(data_file, show_col_types = FALSE)

# Re-centre only after every restriction defining the final analytical sample.
# Subtracting the annual median from the existing relative date is equivalent
# to subtracting the annual median from the original laying date.
dat <- dat |>
  group_by(year_f) |>
  mutate(rel_LD = rel_LD - median(rel_LD, na.rm = TRUE)) |>
  ungroup()

annual <- dat |>
  group_by(year_f) |>
  summarise(
    annual_median = median(rel_LD),
    .groups = "drop"
  )

checks <- tibble(
  check = c(
    "Maximum absolute annual median of rel_LD (days)",
    "Correlation of rel_LD with temp_mean_c"
  ),
  value = c(
    max(abs(annual$annual_median)),
    cor(dat$rel_LD, dat$temp_mean_c)
  )
)

stopifnot(max(abs(annual$annual_median)) < 1e-10)
write_csv(dat, data_file)
write_csv(checks, "tables/fledging_final_relLD_centering_checks.csv")
print(checks)
