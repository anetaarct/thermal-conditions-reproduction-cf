#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(glmmTMB)
  library(readr); library(tibble); library(tidyr)
})
set.seed(20260815)

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f),
         year = as.numeric(as.character(year_f)))
form <- fledged ~ rel_LD * age_class + temp_mean_c + year_c +
  (1 | year_f)
zi_form <- ~ rel_LD + temp_mean_c + year_c + age_class
ctrl_default <- glmmTMBControl(
  optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 1
)
ctrl_bfgs <- glmmTMBControl(
  optimizer = optim, optArgs = list(method = "BFGS"), parallel = 1
)

m_pois <- glmmTMB(form, ziformula = zi_form, family = truncated_poisson(),
                  data = dat, control = ctrl_default)
m_nb <- glmmTMB(form, ziformula = zi_form, family = truncated_nbinom2(),
                data = dat, control = ctrl_bfgs)
saveRDS(m_pois, "models/glmmTMB/fledged_count_hurdle_poisson.rds")
saveRDS(m_nb, "models/glmmTMB/fledged_count_hurdle_nbinom2_bfgs.rds")

metrics <- tibble(
  model = c("Hurdle NB2", "Hurdle Poisson"),
  convergence = c(m_nb$fit$convergence, m_pois$fit$convergence),
  pdHess = c(m_nb$sdr$pdHess, m_pois$sdr$pdHess),
  AIC = c(AIC(m_nb), AIC(m_pois)), BIC = c(BIC(m_nb), BIC(m_pois)),
  logLik = c(as.numeric(logLik(m_nb)), as.numeric(logLik(m_pois)))
)
write_csv(metrics, "tables/fledged_count_hurdle_optimizer_comparison.csv")

tidy_part <- function(x, component, label) {
  as.data.frame(x) |> rownames_to_column("term") |> as_tibble() |>
    transmute(model = "fledged_count_hurdle_poisson", component = label,
              term, beta = Estimate, SE = `Std. Error`,
              z = `z value`, p = `Pr(>|z|)`)
}
coefs <- bind_rows(
  tidy_part(summary(m_pois)$coefficients$cond, "cond", "positive_count"),
  tidy_part(summary(m_pois)$coefficients$zi, "zi", "zero_probability")
)
write_csv(coefs, "tables/fledged_count_hurdle_coefficients.csv")

# Zero-frequency calibration.
sims <- simulate(m_pois, nsim = 1000, seed = 20260815)
sim_zero <- vapply(sims, function(x) mean(as.numeric(x) == 0), numeric(1))
zero_summary <- tibble(
  observed = mean(dat$fledged == 0), simulated_mean = mean(sim_zero),
  simulated_median = median(sim_zero),
  simulation_CI_low = unname(quantile(sim_zero, .025)),
  simulation_CI_high = unname(quantile(sim_zero, .975)), nsim = length(sim_zero)
)
write_csv(zero_summary, "tables/fledged_count_hurdle_zero_frequency_summary.csv")

# Pearson dispersion of the positive component and its simulation reference.
lambda <- exp(as.numeric(predict(m_pois, type = "link")))
p0 <- exp(-lambda); mean_trunc <- lambda / (1 - p0)
var_trunc <- mean_trunc * (1 + lambda - mean_trunc)
positive <- dat$fledged > 0; n_positive <- sum(positive)
n_cond_parameters <- nrow(summary(m_pois)$coefficients$cond)
pearson <- sum((dat$fledged[positive] - mean_trunc[positive])^2 /
                 var_trunc[positive]) / (n_positive - n_cond_parameters)
lp <- lambda[positive]; p0p <- exp(-lp)
sim_disp <- replicate(1000, {
  u <- runif(n_positive, min = p0p, max = 1); y <- qpois(u, lp)
  sum((y - mean_trunc[positive])^2 / var_trunc[positive]) /
    (n_positive - n_cond_parameters)
})
disp_summary <- tibble(
  observed_positive_broods = n_positive, observed_Pearson_dispersion = pearson,
  simulated_median = median(sim_disp),
  simulation_CI_low = unname(quantile(sim_disp, .025)),
  simulation_CI_high = unname(quantile(sim_disp, .975)),
  proportion_simulated_at_least_observed = mean(sim_disp >= pearson),
  hurdle_poisson_AIC = AIC(m_pois), hurdle_nbinom2_AIC = AIC(m_nb),
  delta_AIC_nbinom2 = AIC(m_nb) - AIC(m_pois),
  hurdle_nbinom2_theta = sigma(m_nb),
  poisson_convergence = m_pois$fit$convergence,
  poisson_pdHess = m_pois$sdr$pdHess,
  nbinom2_convergence = m_nb$fit$convergence,
  nbinom2_pdHess = m_nb$sdr$pdHess, nsim = length(sim_disp)
)
write_csv(disp_summary, "tables/fledged_hurdle_positive_dispersion_summary.csv")
write_csv(tibble(simulation = seq_along(sim_disp), dispersion = sim_disp),
          "tables/fledged_hurdle_positive_dispersion_simulations.csv")

# Update chapter-level numerical diagnostics.
diag_file <- "tables/fledging_chapter_diagnostics.csv"
diagnostics <- read_csv(diag_file, show_col_types = FALSE) |>
  filter(model != "fledged_hurdle_poisson") |>
  bind_rows(tibble(model = "fledged_hurdle_poisson",
                   pdHess = m_pois$sdr$pdHess,
                   convergence_code = m_pois$fit$convergence))
write_csv(diagnostics, diag_file)

pred_df <- function(nd) {
  pr <- predict(m_pois, newdata = nd, type = "response", se.fit = TRUE,
                re.form = NA, allow.new.levels = TRUE)
  bind_cols(nd, tibble(estimate = as.numeric(pr$fit), SE = as.numeric(pr$se.fit))) |>
    mutate(CI_low = pmax(0, estimate - 1.96 * SE), CI_high = estimate + 1.96 * SE)
}
ref <- function(n) tibble(
  rel_LD = rep(0, n), temp_mean_c = rep(0, n), year_c = rep(0, n),
  age_class = factor(rep("YOUNG", n), levels = levels(dat$age_class)),
  female_id = factor(rep(levels(dat$female_id)[1], n), levels = levels(dat$female_id)),
  year_f = factor(rep(levels(dat$year_f)[1], n), levels = levels(dat$year_f))
)
blue <- "#0072B2"; orange <- "#D55E00"
age_cols <- c(YOUNG = blue, OLD = orange)
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))

tq <- quantile(dat$temp_mean_c, c(.05, .95))
temp_center <- read_csv("tables/fledging_bb_centring_constants.csv", show_col_types = FALSE)$temp_mean[1]
nd_temp <- ref(160) |>
  mutate(temp_mean_c = seq(tq[1], tq[2], length.out = n()),
         temp_mean = temp_mean_c + temp_center)
d_temp <- pred_df(nd_temp) |> select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d_temp, "tables/fig5_5_fledged_count_temperature.csv")
p_temp <- ggplot(d_temp, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted number fledged") + theme_pub
ggsave("figures/Fig5_5_fledged_count_temperature.png", p_temp,
       width = 6.4, height = 4.8, dpi = 400)

rq <- quantile(dat$rel_LD, c(.05, .95))
nd_ld <- crossing(rel_LD = seq(rq[1], rq[2], length.out = 160),
                  age_class = factor(levels(dat$age_class), levels = levels(dat$age_class))) |>
  mutate(temp_mean_c = 0, year_c = 0,
         female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
         year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f)))
d_ld <- pred_df(nd_ld) |> select(rel_LD, age_class, estimate, SE, CI_low, CI_high)
write_csv(d_ld, "tables/fig5_6_fledged_count_rellD_age.csv")
p_ld <- ggplot(d_ld, aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) + scale_colour_manual(values = age_cols) +
  scale_fill_manual(values = age_cols) +
  labs(x = "Laying date relative to annual median (days)",
       y = "Predicted number fledged") + theme_pub
ggsave("figures/Fig5_6_fledged_count_rellD_age.png", p_ld,
       width = 6.4, height = 4.8, dpi = 400)

year_center <- mean(dat$year - dat$year_c)
nd_year <- ref(160) |>
  mutate(year = seq(min(dat$year), max(dat$year), length.out = n()),
         year_c = year - year_center)
d_year <- pred_df(nd_year) |> select(year, estimate, SE, CI_low, CI_high)
write_csv(d_year, "tables/fig_realised_fledged_year.csv")
p_year <- ggplot(d_year, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = "Predicted number fledged") + theme_pub
ggsave("figures/Fig_realised_fledged_year.png", p_year,
       width = 6.4, height = 4.8, dpi = 400)

print(metrics); print(zero_summary); print(disp_summary, width = Inf); print(coefs)
