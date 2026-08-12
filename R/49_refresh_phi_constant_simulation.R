#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(glmmTMB); library(readr); library(tibble)
})

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f))
m_full <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")
m_constant <- readRDS("models/glmmTMB/fledging_bb_constant_precision.rds")
p_fixed <- as.numeric(predict(m_constant, type = "conditional"))
phi_constant <- as.numeric(sigma(m_constant))
n_trials <- dat$clutch
X_phi <- model.matrix(~ rel_LD + temp_mean_c + year_c + age_class, data = dat)

bb_phi_nll <- function(beta, y) {
  phi_i <- exp(drop(X_phi %*% beta))
  alpha <- pmax(p_fixed * phi_i, 1e-10)
  beta_shape <- pmax((1 - p_fixed) * phi_i, 1e-10)
  -sum(lchoose(n_trials, y) +
         lbeta(y + alpha, n_trials - y + beta_shape) -
         lbeta(alpha, beta_shape))
}
fit_phi_aux <- function(y) {
  start <- c(log(phi_constant), rep(0, ncol(X_phi) - 1))
  lower <- c(-10, rep(-2, ncol(X_phi) - 1))
  upper <- c(10, rep(2, ncol(X_phi) - 1))
  fit <- optim(start, bb_phi_nll, y = y, method = "L-BFGS-B",
               lower = lower, upper = upper,
               control = list(maxit = 500, factr = 1e7))
  boundary <- any(abs(fit$par - lower) < 1e-5 | abs(fit$par - upper) < 1e-5)
  c(convergence = fit$convergence, at_boundary = as.integer(boundary), fit$par)
}

observed <- fit_phi_aux(dat$fledged)
set.seed(20260811)
simulated_y <- lapply(seq_len(200L), function(i) {
  latent_p <- rbeta(nrow(dat), p_fixed * phi_constant,
                    (1 - p_fixed) * phi_constant)
  rbinom(nrow(dat), size = n_trials, prob = latent_p)
})
cluster <- parallel::makeCluster(4L)
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterExport(
  cluster,
  c("X_phi", "p_fixed", "phi_constant", "n_trials", "bb_phi_nll", "fit_phi_aux"),
  envir = environment()
)
results <- parallel::parLapply(cluster, simulated_y, fit_phi_aux)
parallel::stopCluster(cluster)
on.exit(NULL, add = FALSE)

mat <- do.call(rbind, results)
colnames(mat) <- c("convergence", "at_boundary", colnames(X_phi))
sim <- as_tibble(mat) |> mutate(simulation = row_number())
write_csv(sim, "tables/fledging_phi_constant_simulations.csv")
valid <- sim$temp_mean_c[sim$convergence == 0 & sim$at_boundary == 0]
observed_beta <- unname(observed[2 + match("temp_mean_c", colnames(X_phi))])
summary_out <- tibble(
  observed_auxiliary_beta = observed_beta,
  simulated_median = median(valid),
  simulation_CI_low = unname(quantile(valid, .025)),
  simulation_CI_high = unname(quantile(valid, .975)),
  proportion_at_least_observed = mean(valid >= observed_beta),
  successful_interior_fits = length(valid),
  boundary_fits = sum(sim$at_boundary == 1),
  simulations = nrow(sim)
)
write_csv(summary_out, "tables/fledging_phi_temperature_constant_simulation_summary.csv")

blue <- "#0072B2"; orange <- "#D55E00"
fig_data <- sim |> filter(convergence == 0, at_boundary == 0) |>
  select(simulation, temp_mean_c)
write_csv(fig_data, "tables/figA_phi_constant_simulation.csv")
p <- ggplot(fig_data, aes(temp_mean_c)) +
  geom_histogram(bins = 30, fill = scales::alpha(blue, .7), colour = "white") +
  geom_vline(xintercept = observed_beta, colour = orange, linewidth = 1.1) +
  labs(x = "Temperature coefficient in log-phi under constant phi",
       y = "Simulated datasets") +
  theme_classic(base_size = 11)
ggsave("figures/FigA_phi_constant_simulation.png", p,
       width = 6.4, height = 4.8, dpi = 400, bg = "white")
print(summary_out)
