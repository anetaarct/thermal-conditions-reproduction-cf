#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

set.seed(20260813)
dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id), year_f = factor(year_f)
  )
m_pois <- readRDS("models/glmmTMB/fledged_count_hurdle_poisson.rds")
m_nb <- readRDS("models/glmmTMB/fledged_count_hurdle_nbinom2_bfgs.rds")

# Moments of a zero-truncated Poisson with untruncated rate lambda.
lambda <- exp(as.numeric(predict(m_pois, type = "link")))
p_zero_base <- exp(-lambda)
mean_trunc <- lambda / (1 - p_zero_base)
var_trunc <- mean_trunc * (1 + lambda - mean_trunc)
positive <- dat$fledged > 0
n_positive <- sum(positive)
n_cond_parameters <- nrow(summary(m_pois)$coefficients$cond)
pearson_stat <- sum(
  (dat$fledged[positive] - mean_trunc[positive])^2 / var_trunc[positive]
)
pearson_dispersion <- pearson_stat / (n_positive - n_cond_parameters)

# Simulation reference distribution conditional on the fitted random effects.
# Inverse-CDF sampling from Poisson conditional on Y > 0 keeps the same fitted
# positive broods and isolates the variance assumption of the positive component.
lambda_positive <- lambda[positive]
p0_positive <- exp(-lambda_positive)
sim_dispersion <- replicate(1000, {
  u <- runif(n_positive, min = p0_positive, max = 1)
  y <- qpois(u, lambda_positive)
  sum((y - mean_trunc[positive])^2 / var_trunc[positive]) /
    (n_positive - n_cond_parameters)
})

summary_out <- tibble(
  observed_positive_broods = n_positive,
  observed_Pearson_dispersion = pearson_dispersion,
  simulated_median = median(sim_dispersion),
  simulation_CI_low = unname(quantile(sim_dispersion, .025)),
  simulation_CI_high = unname(quantile(sim_dispersion, .975)),
  proportion_simulated_at_least_observed = mean(sim_dispersion >= pearson_dispersion),
  hurdle_poisson_AIC = AIC(m_pois),
  hurdle_nbinom2_AIC = AIC(m_nb),
  delta_AIC_nbinom2 = AIC(m_nb) - AIC(m_pois),
  hurdle_nbinom2_theta = sigma(m_nb),
  poisson_convergence = m_pois$fit$convergence,
  poisson_pdHess = m_pois$sdr$pdHess,
  nbinom2_convergence = m_nb$fit$convergence,
  nbinom2_pdHess = m_nb$sdr$pdHess,
  nsim = length(sim_dispersion)
)
write_csv(summary_out, "tables/fledged_hurdle_positive_dispersion_summary.csv")
write_csv(
  tibble(simulation = seq_along(sim_dispersion), dispersion = sim_dispersion),
  "tables/fledged_hurdle_positive_dispersion_simulations.csv"
)
print(summary_out, width = Inf)
