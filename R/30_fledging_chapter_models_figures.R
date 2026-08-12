#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tibble)
  library(tidyr)
})
Sys.setlocale("LC_TIME", "C")
set.seed(20260810)
dir.create("models/glmmTMB", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id), year_f = factor(year_f),
    year = as.integer(as.character(year_f)),
    success = fledged / clutch, zero = as.integer(fledged == 0)
  )
stopifnot(
  nrow(dat) == 9716L, nlevels(dat$year_f) == 38L,
  min(dat$year) == 1982L, max(dat$year) == 2019L,
  sum(dat$zero) == 2311L
)

ctrl <- glmmTMBControl(
  optCtrl = list(iter.max = 20000, eval.max = 20000), parallel = 1
)
location_formula <- cbind(fledged, failed) ~
  rel_LD * age_class + temp_mean_c + year_c +
  (1 | year_f) + (1 | female_id)

# Main beta-binomial structure comparison on the complete dataset.
main_file <- "models/glmmTMB/fledging_bb_main_no_tempsd.rds"
if (file.exists(main_file)) {
  m_bb_scale <- readRDS(main_file)
} else {
  m_bb_scale <- glmmTMB(
    location_formula,
    dispformula = ~ rel_LD + temp_mean_c + year_c + age_class,
    ziformula = ~ 0, family = betabinomial(), data = dat, control = ctrl
  )
  saveRDS(m_bb_scale, main_file)
}

constant_file <- "models/glmmTMB/fledging_bb_constant_precision.rds"
if (file.exists(constant_file)) {
  m_bb_constant <- readRDS(constant_file)
} else {
  m_bb_constant <- glmmTMB(
    location_formula, dispformula = ~ 1, ziformula = ~ 0,
    family = betabinomial(), data = dat, control = ctrl
  )
  saveRDS(m_bb_constant, constant_file)
}

stopifnot(
  isTRUE(m_bb_scale$sdr$pdHess), isTRUE(m_bb_constant$sdr$pdHess),
  m_bb_scale$fit$convergence == 0L, m_bb_constant$fit$convergence == 0L
)

model_metrics <- function(m, label) tibble(
  model = label,
  df = attr(logLik(m), "df"),
  AIC = AIC(m), BIC = BIC(m), logLik = as.numeric(logLik(m))
)
bb_metrics <- bind_rows(
  model_metrics(m_bb_constant, "Constant precision"),
  model_metrics(m_bb_scale, "Predictor-dependent precision")
) |>
  mutate(delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC))
bb_lrt <- anova(m_bb_constant, m_bb_scale)
bb_comparison <- bb_metrics |>
  mutate(
    LR_chisq = c(NA_real_, bb_lrt$Chisq[2]),
    LRT_df = c(NA_real_, bb_lrt$`Chi Df`[2]),
    LRT_p = c(NA_real_, bb_lrt$`Pr(>Chisq)`[2])
  )
write_csv(bb_comparison, "tables/fledging_bb_precision_comparison.csv")

tidy_component <- function(m, model_name) {
  sm <- summary(m)$coefficients
  ans <- as.data.frame(sm$cond) |>
    rownames_to_column("term") |>
    as_tibble() |>
    transmute(
      model = model_name, component = "location", term,
      beta = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`
    )
  if (!is.null(sm$disp)) ans <- bind_rows(
    ans,
    as.data.frame(sm$disp) |>
      rownames_to_column("term") |>
      as_tibble() |>
      transmute(
        model = model_name, component = "precision_log_phi", term,
        beta = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`
      )
  )
  ans
}
bb_coef <- tidy_component(m_bb_scale, "beta_binomial_main")
write_csv(bb_coef, "tables/fledging_bb_main_coefficients.csv")

age_trends <- emtrends(
  m_bb_scale, ~ age_class, var = "rel_LD", component = "cond"
) |>
  confint(level = 0.95) |>
  as.data.frame() |>
  as_tibble()
write_csv(age_trends, "tables/fledging_relLD_age_emtrends.csv")

summary_table <- tibble(
  years = paste0(min(dat$year), "-", max(dat$year)),
  seasons = nlevels(dat$year_f), broods = nrow(dat),
  females = nlevels(dat$female_id), zero_broods = sum(dat$zero),
  zero_percent = 100 * mean(dat$zero)
)
write_csv(summary_table, "tables/fledging_chapter_sample_summary.csv")

# Simulation-based adequacy check for complete brood failures.
bb_sim <- simulate(m_bb_scale, nsim = 1000, seed = 20260810)
sim_zero <- vapply(bb_sim, function(x) mean(as.numeric(x) == 0), numeric(1))
observed_zero <- mean(dat$fledged == 0)
zero_summary <- tibble(
  observed = observed_zero,
  simulated_mean = mean(sim_zero), simulated_median = median(sim_zero),
  simulation_CI_low = unname(quantile(sim_zero, .025)),
  simulation_CI_high = unname(quantile(sim_zero, .975)),
  proportion_as_or_more_extreme = mean(
    abs(sim_zero - mean(sim_zero)) >= abs(observed_zero - mean(sim_zero))
  ),
  nsim = length(sim_zero)
)
write_csv(tibble(simulation = seq_along(sim_zero), zero_proportion = sim_zero),
          "tables/fledging_zero_frequency_simulations.csv")
write_csv(zero_summary, "tables/fledging_zero_frequency_summary.csv")

# Absolute realised fledgling production: Poisson versus negative binomial.
count_formula <- fledged ~ rel_LD * age_class + temp_mean_c + year_c +
  (1 | year_f) + (1 | female_id)
poisson_file <- "models/glmmTMB/fledged_count_poisson.rds"
nb_file <- "models/glmmTMB/fledged_count_nbinom2.rds"
if (file.exists(poisson_file)) {
  m_count_poisson <- readRDS(poisson_file)
} else {
  m_count_poisson <- glmmTMB(
    count_formula, ziformula = ~ 0, family = poisson(), data = dat, control = ctrl
  )
  saveRDS(m_count_poisson, poisson_file)
}
if (file.exists(nb_file)) {
  m_count_nb <- readRDS(nb_file)
} else {
  m_count_nb <- glmmTMB(
    count_formula, ziformula = ~ 0, family = nbinom2(), data = dat, control = ctrl
  )
  saveRDS(m_count_nb, nb_file)
}
stopifnot(
  isTRUE(m_count_poisson$sdr$pdHess), isTRUE(m_count_nb$sdr$pdHess),
  m_count_poisson$fit$convergence == 0L, m_count_nb$fit$convergence == 0L
)

count_metrics <- bind_rows(
  model_metrics(m_count_poisson, "Poisson"),
  model_metrics(m_count_nb, "Negative binomial (NB2)")
) |>
  mutate(
    delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC),
    Pearson_dispersion = c(
      sum(residuals(m_count_poisson, type = "pearson")^2) / df.residual(m_count_poisson),
      sum(residuals(m_count_nb, type = "pearson")^2) / df.residual(m_count_nb)
    )
  )
write_csv(count_metrics, "tables/fledged_count_model_comparison.csv")
count_coef <- tidy_component(m_count_nb, "fledged_count_nbinom2") |>
  filter(component == "location")
write_csv(count_coef, "tables/fledged_count_nbinom2_coefficients.csv")
count_age_trends <- emtrends(m_count_nb, ~ age_class, var = "rel_LD") |>
  confint(level = .95) |>
  as.data.frame() |>
  as_tibble()
write_csv(count_age_trends, "tables/fledged_count_relLD_age_emtrends.csv")

# Equivalent simulation-based count-model diagnostics, including zero frequency.
count_sim <- simulate(m_count_nb, nsim = 1000, seed = 20260811)
count_sim_zero <- vapply(count_sim, function(x) mean(as.numeric(x) == 0), numeric(1))
count_zero_summary <- tibble(
  observed = observed_zero,
  simulated_mean = mean(count_sim_zero),
  simulated_median = median(count_sim_zero),
  simulation_CI_low = unname(quantile(count_sim_zero, .025)),
  simulation_CI_high = unname(quantile(count_sim_zero, .975)),
  nsim = length(count_sim_zero)
)
write_csv(count_zero_summary, "tables/fledged_count_zero_frequency_summary.csv")

# Hurdle model: zero occurrence and positive fledgling production are modelled
# as separate processes. The response-scale prediction integrates both parts.
hurdle_file <- "models/glmmTMB/fledged_count_hurdle_poisson.rds"
stopifnot(file.exists(hurdle_file))
m_count_hurdle <- readRDS(hurdle_file)
stopifnot(
  isTRUE(m_count_hurdle$sdr$pdHess),
  m_count_hurdle$fit$convergence == 0L
)
hurdle_coef <- tidy_component(m_count_hurdle, "fledged_count_hurdle_poisson")
zi_coef <- summary(m_count_hurdle)$coefficients$zi |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble() |>
  transmute(
    model = "fledged_count_hurdle_poisson", component = "zero_probability",
    term, beta = Estimate, SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`
  )
hurdle_coef <- bind_rows(
  hurdle_coef |> filter(component == "location") |>
    mutate(component = "positive_count"),
  zi_coef
)
write_csv(hurdle_coef, "tables/fledged_count_hurdle_coefficients.csv")
write_csv(
  model_metrics(m_count_hurdle, "Hurdle Poisson"),
  "tables/fledged_count_hurdle_metrics.csv"
)
hurdle_sim <- simulate(m_count_hurdle, nsim = 1000, seed = 20260812)
hurdle_sim_zero <- vapply(
  hurdle_sim, function(x) mean(as.numeric(x) == 0), numeric(1)
)
hurdle_zero_summary <- tibble(
  observed = observed_zero,
  simulated_mean = mean(hurdle_sim_zero),
  simulated_median = median(hurdle_sim_zero),
  simulation_CI_low = unname(quantile(hurdle_sim_zero, .025)),
  simulation_CI_high = unname(quantile(hurdle_sim_zero, .975)),
  nsim = length(hurdle_sim_zero)
)
write_csv(
  hurdle_zero_summary,
  "tables/fledged_count_hurdle_zero_frequency_summary.csv"
)

# Targeted simulation check for a mechanical temperature effect on phi.
# The fitted conditional means (including fitted random effects) are held fixed,
# while 500 beta-binomial datasets are generated with constant phi.
p_fixed <- as.numeric(predict(m_bb_constant, type = "conditional"))
phi_constant <- as.numeric(sigma(m_bb_constant))
n_trials <- dat$clutch
X_phi <- model.matrix(~ rel_LD + temp_mean_c + year_c + age_class, data = dat)
bb_phi_nll <- function(beta, y) {
  phi_i <- exp(drop(X_phi %*% beta))
  alpha <- pmax(p_fixed * phi_i, 1e-10)
  beta_shape <- pmax((1 - p_fixed) * phi_i, 1e-10)
  -sum(
    lchoose(n_trials, y) + lbeta(y + alpha, n_trials - y + beta_shape) -
      lbeta(alpha, beta_shape)
  )
}
fit_phi_aux <- function(y) {
  start <- c(log(phi_constant), rep(0, ncol(X_phi) - 1))
  lower <- c(-10, rep(-2, ncol(X_phi) - 1))
  upper <- c(10, rep(2, ncol(X_phi) - 1))
  fit <- optim(
    start, bb_phi_nll, y = y, method = "L-BFGS-B",
    lower = lower, upper = upper,
    control = list(maxit = 500, factr = 1e7)
  )
  at_boundary <- any(abs(fit$par - lower) < 1e-5 | abs(fit$par - upper) < 1e-5)
  c(convergence = fit$convergence, at_boundary = as.integer(at_boundary), fit$par)
}
observed_phi_aux <- fit_phi_aux(dat$fledged)
n_phi_sim <- 200L
phi_sim_results <- matrix(NA_real_, nrow = n_phi_sim, ncol = ncol(X_phi) + 2L)
colnames(phi_sim_results) <- c("convergence", "at_boundary", colnames(X_phi))
for (s in seq_len(n_phi_sim)) {
  latent_p <- rbeta(
    nrow(dat), p_fixed * phi_constant, (1 - p_fixed) * phi_constant
  )
  y_sim <- rbinom(nrow(dat), size = n_trials, prob = latent_p)
  phi_sim_results[s, ] <- fit_phi_aux(y_sim)
}
phi_sim <- as_tibble(phi_sim_results) |>
  mutate(simulation = row_number())
write_csv(phi_sim, "tables/fledging_phi_constant_simulations.csv")
temp_col <- "temp_mean_c"
valid_temp <- phi_sim[[temp_col]][phi_sim$convergence == 0 & phi_sim$at_boundary == 0]
phi_mechanism_summary <- tibble(
  observed_auxiliary_beta = unname(observed_phi_aux[2 + match(temp_col, colnames(X_phi))]),
  simulated_median = median(valid_temp),
  simulation_CI_low = unname(quantile(valid_temp, .025)),
  simulation_CI_high = unname(quantile(valid_temp, .975)),
  proportion_at_least_observed = mean(valid_temp >=
    unname(observed_phi_aux[2 + match(temp_col, colnames(X_phi))])),
  successful_interior_fits = length(valid_temp),
  boundary_fits = sum(phi_sim$at_boundary == 1), simulations = n_phi_sim
)
write_csv(
  phi_mechanism_summary,
  "tables/fledging_phi_temperature_constant_simulation_summary.csv"
)

age_cols <- c(YOUNG = "#0072B2", OLD = "#D55E00")
blue <- "#0072B2"
orange <- "#D55E00"
theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"), legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = .3)
  )
ref_grid <- function(n) tibble(
  rel_LD = rep(0, n), temp_mean_c = rep(0, n), year_c = rep(0, n),
  age_class = factor(rep("YOUNG", n), levels = levels(dat$age_class)),
  female_id = factor(rep(levels(dat$female_id)[1], n), levels = levels(dat$female_id)),
  year_f = factor(rep(levels(dat$year_f)[1], n), levels = levels(dat$year_f))
)
pred_df <- function(m, nd, type = "response") {
  pr <- predict(
    m, newdata = nd, type = type, se.fit = TRUE,
    re.form = NA, allow.new.levels = TRUE
  )
  bind_cols(
    nd,
    tibble(
      estimate = as.numeric(pr$fit), SE = as.numeric(pr$se.fit),
      CI_low = estimate - 1.96 * SE, CI_high = estimate + 1.96 * SE
    )
  )
}
temp_center <- read_csv(
  "tables/fledging_bb_centring_constants.csv", show_col_types = FALSE
)$temp_mean[1]
tq <- quantile(dat$temp_mean_c, c(.05, .95))
rq <- quantile(dat$rel_LD, c(.05, .95))

# Main beta-binomial figures.
nd1 <- ref_grid(160) |>
  mutate(
    temp_mean_c = seq(tq[1], tq[2], length.out = n()),
    temp_mean = temp_mean_c + temp_center
  )
d1 <- pred_df(m_bb_scale, nd1) |>
  select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d1, "tables/fig5_1_success_temperature.csv")
p1 <- ggplot(d1, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted fledging success") +
  theme_pub
ggsave("figures/Fig5_1_success_temperature.png", p1, width = 6.4, height = 4.8, dpi = 400)

# Joint temperature figure: mean success and conditional consistency.
d1_phi <- pred_df(m_bb_scale, nd1, type = "disp") |>
  select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d1_phi, "tables/fig5_1b_phi_temperature.csv")
p1_phi <- ggplot(d1_phi, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(orange, .2)) +
  geom_line(colour = orange, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = expression("Consistency (" * phi * ")")) +
  theme_pub
save_two(
  "figures/Fig5_1_success_phi_temperature.png", p1, p1_phi,
  width = 10, height = 4.6, dpi = 400
)

nd2 <- crossing(
  rel_LD = seq(rq[1], rq[2], length.out = 160),
  age_class = factor(levels(dat$age_class), levels = levels(dat$age_class))
) |>
  mutate(
    temp_mean_c = 0, year_c = 0,
    female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
    year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f))
  )
d2 <- pred_df(m_bb_scale, nd2) |>
  select(rel_LD, age_class, estimate, SE, CI_low, CI_high)
write_csv(d2, "tables/fig5_2_success_rellD_age.csv")
p2 <- ggplot(d2, aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = age_cols) + scale_fill_manual(values = age_cols) +
  labs(
    x = "Laying date relative to annual median (days)",
    y = "Predicted fledging success"
  ) + theme_pub
ggsave("figures/Fig5_2_success_rellD_age.png", p2, width = 6.4, height = 4.8, dpi = 400)

nd3 <- ref_grid(160) |>
  mutate(rel_LD = seq(rq[1], rq[2], length.out = n()))
d3 <- pred_df(m_bb_scale, nd3, type = "disp") |>
  select(rel_LD, estimate, SE, CI_low, CI_high)
write_csv(d3, "tables/fig5_3_phi_rellD.csv")
p3 <- ggplot(d3, aes(rel_LD, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(orange, .2)) +
  geom_line(colour = orange, linewidth = 1) +
  labs(
    x = "Laying date relative to annual median (days)",
    y = "Beta-binomial precision (phi)"
  ) + theme_pub
ggsave("figures/Fig5_3_phi_rellD.png", p3, width = 6.4, height = 4.8, dpi = 400)

# Precision by female age class.
nd_age_phi <- ref_grid(2) |>
  mutate(age_class = factor(levels(dat$age_class), levels = levels(dat$age_class)))
d_age_phi <- pred_df(m_bb_scale, nd_age_phi, type = "disp") |>
  select(age_class, estimate, SE, CI_low, CI_high)
write_csv(d_age_phi, "tables/fig5_3b_phi_age.csv")
p_age_phi <- ggplot(d_age_phi, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .12, linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = age_cols, guide = "none") +
  labs(x = "Female age class", y = expression("Consistency (" * phi * ")")) + theme_pub
ggsave("figures/Fig5_3b_phi_age.png", p_age_phi, width = 5.4, height = 4.6, dpi = 400)

# Constant-phi simulation diagnostic for the Appendix.
phi_hist_data <- phi_sim |>
  filter(convergence == 0, at_boundary == 0) |>
  select(simulation, temp_mean_c)
write_csv(phi_hist_data, "tables/figA_phi_constant_simulation.csv")
p_phi_hist <- ggplot(phi_hist_data, aes(temp_mean_c)) +
  geom_histogram(bins = 30, fill = scales::alpha(blue, .7), colour = "white") +
  geom_vline(
    xintercept = phi_mechanism_summary$observed_auxiliary_beta,
    colour = orange, linewidth = 1.1
  ) +
  labs(
    x = "Temperature coefficient in log-phi under constant phi",
    y = "Simulated datasets"
  ) + theme_pub
ggsave("figures/FigA_phi_constant_simulation.png", p_phi_hist,
       width = 6.4, height = 4.8, dpi = 400)

# Complete-failure diagnostic figure.
zero_plot_data <- tibble(zero_proportion = sim_zero)
p4 <- ggplot(zero_plot_data, aes(zero_proportion)) +
  geom_histogram(bins = 35, fill = scales::alpha(blue, .7), colour = "white") +
  geom_vline(xintercept = observed_zero, colour = orange, linewidth = 1.1) +
  labs(
    x = "Simulated proportion of broods with zero fledglings",
    y = "Simulated datasets"
  ) + theme_pub
ggsave("figures/Fig5_4_zero_frequency_diagnostic.png", p4, width = 6.4, height = 4.8, dpi = 400)

# Absolute fledgling production figures from the selected hurdle-Poisson model.
d5 <- pred_df(m_count_hurdle, nd1) |>
  select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d5, "tables/fig5_5_fledged_count_temperature.csv")
p5 <- ggplot(d5, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted number fledged") +
  theme_pub
ggsave("figures/Fig5_5_fledged_count_temperature.png", p5, width = 6.4, height = 4.8, dpi = 400)

d6 <- pred_df(m_count_hurdle, nd2) |>
  select(rel_LD, age_class, estimate, SE, CI_low, CI_high)
write_csv(d6, "tables/fig5_6_fledged_count_rellD_age.csv")
p6 <- ggplot(d6, aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = age_cols) + scale_fill_manual(values = age_cols) +
  labs(
    x = "Laying date relative to annual median (days)",
    y = "Predicted number fledged"
  ) + theme_pub
ggsave("figures/Fig5_6_fledged_count_rellD_age.png", p6, width = 6.4, height = 4.8, dpi = 400)

diagnostics <- tibble(
  model = c(
    "beta_binomial_constant", "beta_binomial_scale",
    "fledged_poisson", "fledged_nbinom2", "fledged_hurdle_poisson"
  ),
  pdHess = c(
    m_bb_constant$sdr$pdHess, m_bb_scale$sdr$pdHess,
    m_count_poisson$sdr$pdHess, m_count_nb$sdr$pdHess,
    m_count_hurdle$sdr$pdHess
  ),
  convergence_code = c(
    m_bb_constant$fit$convergence, m_bb_scale$fit$convergence,
    m_count_poisson$fit$convergence, m_count_nb$fit$convergence,
    m_count_hurdle$fit$convergence
  )
)
write_csv(diagnostics, "tables/fledging_chapter_diagnostics.csv")
message("Completed revised fledging-success analyses and figures.")
