#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
})

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id), year_f = factor(year_f)
  )
form <- fledged ~ rel_LD * age_class + temp_mean_c + year_c +
  (1 | year_f) + (1 | female_id)
zi_form <- ~ rel_LD + temp_mean_c + year_c + age_class

ctrl_bfgs <- glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS"),
  optCtrl = list(maxit = 20000), parallel = 1
)
ctrl_default <- glmmTMBControl(
  optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 1
)

m_hurdle_nb_bfgs <- glmmTMB(
  form, ziformula = zi_form, family = truncated_nbinom2(),
  data = dat, control = ctrl_bfgs
)
saveRDS(m_hurdle_nb_bfgs,
        "models/glmmTMB/fledged_count_hurdle_nbinom2_bfgs.rds")

m_hurdle_poisson <- glmmTMB(
  form, ziformula = zi_form, family = truncated_poisson(),
  data = dat, control = ctrl_default
)
saveRDS(m_hurdle_poisson,
        "models/glmmTMB/fledged_count_hurdle_poisson.rds")

diagnostics <- data.frame(
  model = c("hurdle_nbinom2_bfgs", "hurdle_poisson"),
  convergence = c(
    m_hurdle_nb_bfgs$fit$convergence,
    m_hurdle_poisson$fit$convergence
  ),
  pdHess = c(m_hurdle_nb_bfgs$sdr$pdHess, m_hurdle_poisson$sdr$pdHess),
  AIC = c(AIC(m_hurdle_nb_bfgs), AIC(m_hurdle_poisson)),
  BIC = c(BIC(m_hurdle_nb_bfgs), BIC(m_hurdle_poisson)),
  logLik = c(logLik(m_hurdle_nb_bfgs), logLik(m_hurdle_poisson))
)
write_csv(diagnostics, "tables/fledged_count_hurdle_optimizer_comparison.csv")
print(diagnostics)
