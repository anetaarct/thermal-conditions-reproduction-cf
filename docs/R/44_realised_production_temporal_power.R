#!/usr/bin/env Rscript

.libPaths(c(".Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(glmmTMB)
})

set.seed(20260811)
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args)) as.integer(args[1]) else as.integer(Sys.getenv("POWER_NSIM", "500"))
m_count_hurdle <- readRDS("models/glmmTMB/fledged_count_hurdle_poisson.rds")
temperature_models <- readRDS("models/climwin/main_window_temperature_trends.rds")

# Extract every numerical input from the fitted objects.
temperature_trend <- unname(coef(temperature_models$mean_trend)["year"])
dat <- model.frame(m_count_hurdle)
ref <- dat[1, , drop = FALSE]
ref$rel_LD <- 0
ref$age_class <- factor("YOUNG", levels = levels(dat$age_class))
ref$year_c <- 0
ref$year_f <- factor(levels(dat$year_f)[1], levels = levels(dat$year_f))
low <- high <- ref
low$temp_mean_c <- -0.5
high$temp_mean_c <- 0.5
mu_low <- predict(m_count_hurdle, low, type = "response", re.form = NA)
mu_high <- predict(m_count_hurdle, high, type = "response", re.form = NA)
temperature_effect_log <- log(mu_high / mu_low)
temperature_effect_percent <- 100 * (exp(temperature_effect_log) - 1)
beta_expected <- temperature_effect_log * temperature_trend
expected_percent_year <- 100 * (exp(beta_expected) - 1)
sigma_year <- unname(attr(VarCorr(m_count_hurdle)$cond$year_f, "stddev")[1])

power_for_n <- function(n) {
  se <- sigma_year / sqrt(n * (n^2 - 1) / 12)
  z <- abs(beta_expected) / se
  data.frame(se = se, z_score = z, power = pnorm(z - qnorm(0.975)))
}
n_values <- c(39L, 50L, 60L, 80L, 100L, 120L)
analytic <- cbind(data.frame(seasons = n_values),
                  do.call(rbind, lapply(n_values, power_for_n)))
n_80 <- (2:1000)[which(vapply(2:1000, function(n) power_for_n(n)$power >= 0.80,
                              logical(1)))[1]]
inputs <- data.frame(
  quantity = c("Temperature trend", "Marginal temperature effect",
               "Expected temporal trend", "Expected annual percent change",
               "Between-year SD", "Seasons for 80% analytical power"),
  estimate = c(temperature_trend, temperature_effect_log, beta_expected,
               expected_percent_year, sigma_year, n_80),
  unit = c("degrees C/year", "log mean/degree C", "log mean/year", "%/year",
           "log scale", "seasons")
)
write.csv(analytic, "tables/realised_temporal_power_analytic.csv", row.names = FALSE)
write.csv(inputs, "tables/realised_temporal_power_inputs.csv", row.names = FALSE)

# The expected trend is imposed only on the positive-count component.
cond_names <- colnames(model.matrix(m_count_hurdle, component = "cond"))
zi_names <- colnames(model.matrix(m_count_hurdle, component = "zi"))
cond_year <- match("year_c", cond_names)
zi_year <- match("year_c", zi_names)
stopifnot(!is.na(cond_year), !is.na(zi_year))
m_power <- m_count_hurdle
beta_pos <- which(names(m_power$fit$par) == "beta")
beta_zi <- which(names(m_power$fit$par) == "betazi")
m_power$fit$par[beta_pos[cond_year]] <- beta_expected
m_power$fit$par[beta_zi[zi_year]] <- 0
beta_pos_full <- which(names(m_power$fit$parfull) == "beta")
beta_zi_full <- which(names(m_power$fit$parfull) == "betazi")
m_power$fit$parfull[beta_pos_full[cond_year]] <- beta_expected
m_power$fit$parfull[beta_zi_full[zi_year]] <- 0

simulated <- simulate(m_power, nsim = nsim, seed = 20260811)
fit_one <- function(i) {
  if (i %% 10L == 0L) message("Completed ", i, " of ", nsim, " simulations")
  sim_dat <- dat
  sim_dat$fledged <- simulated[[i]]
  fit <- try(glmmTMB(
    fledged ~ rel_LD * age_class + temp_mean_c + year_c + (1 | year_f),
    ziformula = ~ rel_LD + temp_mean_c + year_c + age_class,
    family = truncated_poisson(), data = sim_dat
  ), silent = TRUE)
  if (inherits(fit, "try-error"))
    return(data.frame(simulation = i, estimate = NA_real_, SE = NA_real_, z = NA_real_,
                  p = NA_real_, detected = NA, convergence = NA_integer_, pdHess = NA))
  tab <- summary(fit)$coefficients$cond
  data.frame(simulation = i, estimate = tab["year_c", "Estimate"],
         SE = tab["year_c", "Std. Error"], z = tab["year_c", "z value"],
         p = tab["year_c", "Pr(>|z|)"], detected = tab["year_c", "Pr(>|z|)"] < 0.05,
         convergence = fit$fit$convergence, pdHess = isTRUE(fit$sdr$pdHess))
}
results <- do.call(rbind, lapply(seq_len(nsim), fit_one))
valid <- subset(results, convergence == 0L & pdHess & !is.na(detected))
detected_n <- sum(valid$detected)
binom_ci <- binom.test(detected_n, nrow(valid))$conf.int
simulation_summary <- data.frame(
  simulations_requested = nsim, simulations_valid = nrow(valid), detections = detected_n,
  power = detected_n / nrow(valid), CI_low = binom_ci[1], CI_high = binom_ci[2]
)
write.csv(results, "tables/realised_temporal_power_simulations.csv", row.names = FALSE)
write.csv(simulation_summary, "tables/realised_temporal_power_simulation_summary.csv", row.names = FALSE)
saveRDS(list(inputs = inputs, analytic = analytic, simulations = results,
             simulation_summary = simulation_summary),
        "models/glmmTMB/realised_temporal_power.rds")
