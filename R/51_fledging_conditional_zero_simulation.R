#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(glmmTMB); library(readr); library(tibble)})

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))
m <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")

# Unlike simulate.glmmTMB(), this check conditions on the fitted year and
# female random effects. For glmmTMB's beta-binomial parameterisation,
# alpha = mu * phi and beta = (1 - mu) * phi.
mu <- as.numeric(predict(m, type = "conditional"))
phi <- as.numeric(predict(m, type = "disp"))
stopifnot(length(mu) == nrow(dat), length(phi) == nrow(dat),
          all(mu > 0 & mu < 1), all(phi > 0))

set.seed(20260812)
nsim <- 1000L
zero_proportion <- vapply(seq_len(nsim), function(i) {
  latent_p <- rbeta(nrow(dat), mu * phi, (1 - mu) * phi)
  y <- rbinom(nrow(dat), size = dat$clutch, prob = latent_p)
  mean(y == 0)
}, numeric(1))

simulations <- tibble(simulation = seq_len(nsim), zero_proportion)
summary_out <- tibble(
  simulation_type = "Conditional on fitted year and female random effects",
  observed = mean(dat$fledged == 0),
  simulated_mean = mean(zero_proportion),
  simulated_median = median(zero_proportion),
  simulation_CI_low = unname(quantile(zero_proportion, .025)),
  simulation_CI_high = unname(quantile(zero_proportion, .975)),
  nsim = nsim
)
write_csv(simulations, "tables/fledging_zero_frequency_conditional_simulations.csv")
write_csv(summary_out, "tables/fledging_zero_frequency_conditional_summary.csv")
print(summary_out)
