#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

comparison_table <- function(m0, m1) {
  x <- anova(m0, m1)
  tibble(
    model = c("Constant scale", "Location-scale"),
    df = x$Df, AIC = x$AIC, BIC = x$BIC, logLik = x$logLik,
    deviance = x$deviance, chi_square = x$Chisq,
    chi_df = x$`Chi Df`, p = x$`Pr(>Chisq)`
  ) |>
    mutate(delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC))
}

dat_cs <- read_csv("data_derived/clutch_size_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))

cs0 <- glmmTMB(
  clutch ~ rel_LD + temp_mean_c + year_c + age_class +
    (1 | female_id) + (1 | year_f),
  dispformula = ~1, family = gaussian(), REML = FALSE, data = dat_cs
)
cs1 <- update(
  cs0,
  dispformula = ~ rel_LD + temp_mean_c + year_c + age_class
)
write_csv(comparison_table(cs0, cs1),
          "tables/clutch_size_scale_model_comparison.csv")

message("Clutch-size Gaussian structure test completed without temp_sd_c.")
