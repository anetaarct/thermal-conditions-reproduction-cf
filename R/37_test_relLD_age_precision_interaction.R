#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id),
    year_f = factor(year_f)
  )

m_additive <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")
ctrl <- glmmTMBControl(
  optCtrl = list(iter.max = 20000, eval.max = 20000),
  parallel = 4
)
old_interaction <- readRDS(
  ".archive_before_final_fledging_recentring/models/fledging_bb_precision_relLD_age_interaction.rds"
)
m_interaction <- glmmTMB(
  cbind(fledged, failed) ~ rel_LD * age_class + temp_mean_c + year_c +
    (1 | year_f) + (1 | female_id),
  dispformula = ~ rel_LD * age_class + temp_mean_c + year_c,
  ziformula = ~ 0, family = betabinomial(), data = dat, control = ctrl,
  start = list(
    beta = unname(fixef(old_interaction)$cond),
    betadisp = unname(fixef(old_interaction)$disp),
    theta = unname(getME(old_interaction, "theta"))
  )
)
saveRDS(
  m_interaction,
  "models/glmmTMB/fledging_bb_precision_relLD_age_interaction.rds"
)

cmp <- anova(m_additive, m_interaction)
comparison <- as.data.frame(cmp) |>
  rownames_to_column("model") |>
  as_tibble() |>
  rename(
    df = Df, AIC = AIC, BIC = BIC, logLik = logLik,
    deviance = deviance, chi_square = Chisq,
    chi_df = `Chi Df`, p = `Pr(>Chisq)`
  ) |>
  mutate(
    model = c("Additive precision", "rel_LD x age precision"),
    delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC)
  )
write_csv(comparison, "tables/fledging_precision_relLD_age_anova.csv")

interaction_coef <- summary(m_interaction)$coefficients$disp |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble() |>
  rename(beta = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`) |>
  filter(term == "rel_LD:age_classOLD") |>
  mutate(CI_low = beta - 1.96 * SE, CI_high = beta + 1.96 * SE)
write_csv(interaction_coef, "tables/fledging_precision_relLD_age_interaction.csv")

print(cmp)
print(interaction_coef)
