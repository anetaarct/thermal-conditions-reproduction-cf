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
    female_id = factor(female_id), year_f = factor(year_f)
  )
m_full <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")
m_constant <- readRDS("models/glmmTMB/fledging_bb_constant_precision.rds")

# Hold the fitted conditional mean, including fitted year and female effects,
# fixed while estimating the prespecified precision equation.
p_fixed <- as.numeric(predict(m_full, type = "conditional"))
X <- model.matrix(~ rel_LD + temp_mean_c + year_c + age_class, data = dat)
y <- dat$fledged
n <- dat$clutch
positive <- y > 0

log_bb_pmf <- function(y, n, p, phi) {
  alpha <- pmax(p * phi, 1e-10)
  beta_shape <- pmax((1 - p) * phi, 1e-10)
  lchoose(n, y) + lbeta(y + alpha, n - y + beta_shape) -
    lbeta(alpha, beta_shape)
}

nll_full <- function(coef) {
  phi <- exp(drop(X %*% coef))
  -sum(log_bb_pmf(y, n, p_fixed, phi))
}

nll_positive_truncated <- function(coef) {
  phi <- exp(drop(X[positive, , drop = FALSE] %*% coef))
  log_pmf <- log_bb_pmf(
    y[positive], n[positive], p_fixed[positive], phi
  )
  log_p0 <- log_bb_pmf(
    rep(0, sum(positive)), n[positive], p_fixed[positive], phi
  )
  log_survival <- log1p(-pmin(exp(log_p0), 1 - 1e-12))
  -sum(log_pmf - log_survival)
}

fit_precision <- function(objective) {
  start <- c(log(as.numeric(sigma(m_constant))), rep(0, ncol(X) - 1))
  fit <- optim(
    start, objective, method = "L-BFGS-B",
    lower = c(-10, rep(-2, ncol(X) - 1)),
    upper = c(10, rep(2, ncol(X) - 1)),
    hessian = TRUE,
    control = list(maxit = 1000, factr = 1e7)
  )
  covariance <- tryCatch(solve(fit$hessian), error = function(e) {
    matrix(NA_real_, ncol(X), ncol(X))
  })
  se <- sqrt(pmax(diag(covariance), 0))
  tibble(
    term = colnames(X), beta = fit$par, SE = se,
    CI_low = beta - 1.96 * SE, CI_high = beta + 1.96 * SE,
    z = beta / SE, p = 2 * pnorm(abs(z), lower.tail = FALSE),
    convergence = fit$convergence,
    at_boundary = any(abs(fit$par - c(-10, rep(-2, ncol(X) - 1))) < 1e-5 |
      abs(fit$par - c(10, rep(2, ncol(X) - 1))) < 1e-5)
  )
}

result <- bind_rows(
  fit_precision(nll_full) |> mutate(sample = "all_broods", .before = 1),
  fit_precision(nll_positive_truncated) |>
    mutate(sample = "positive_broods_fixed_mean_diagnostic", .before = 1)
)

# Primary conditional-survival check: jointly re-estimate the location and
# precision equations under the correctly zero-truncated beta-binomial
# likelihood, while retaining fitted year/female contributions as offsets.
X_mu <- model.matrix(~ rel_LD * age_class + temp_mean_c + year_c, data = dat)
b_mu_start <- fixef(m_full)$cond[colnames(X_mu)]
eta_full <- as.numeric(predict(m_full, type = "link"))
random_offset <- eta_full - drop(X_mu %*% b_mu_start)
b_phi_start <- fixef(m_full)$disp[colnames(X)]

nll_positive_joint <- function(par) {
  b_mu <- par[seq_len(ncol(X_mu))]
  b_phi <- par[ncol(X_mu) + seq_len(ncol(X))]
  eta <- random_offset[positive] +
    drop(X_mu[positive, , drop = FALSE] %*% b_mu)
  p <- plogis(eta)
  phi <- exp(drop(X[positive, , drop = FALSE] %*% b_phi))
  log_pmf <- log_bb_pmf(y[positive], n[positive], p, phi)
  log_p0 <- log_bb_pmf(rep(0, sum(positive)), n[positive], p, phi)
  -sum(log_pmf - log1p(-pmin(exp(log_p0), 1 - 1e-12)))
}

joint_start <- c(b_mu_start, b_phi_start)
joint_lower <- c(rep(-10, ncol(X_mu)), -10, rep(-3, ncol(X) - 1))
joint_upper <- c(rep(10, ncol(X_mu)), 10, rep(3, ncol(X) - 1))
joint_fit <- optim(
  joint_start, nll_positive_joint, method = "L-BFGS-B",
  lower = joint_lower, upper = joint_upper, hessian = TRUE,
  control = list(maxit = 3000, factr = 1e7)
)
joint_cov <- tryCatch(solve(joint_fit$hessian), error = function(e) {
  matrix(NA_real_, length(joint_start), length(joint_start))
})
joint_se <- sqrt(pmax(diag(joint_cov), 0))
joint_names <- c(paste0("location:", colnames(X_mu)),
                 paste0("precision:", colnames(X)))
joint_result <- tibble(
  sample = "positive_broods_zero_truncated_joint",
  term = joint_names, beta = joint_fit$par, SE = joint_se,
  CI_low = beta - 1.96 * SE, CI_high = beta + 1.96 * SE,
  z = beta / SE, p = 2 * pnorm(abs(z), lower.tail = FALSE),
  convergence = joint_fit$convergence,
  at_boundary = abs(beta - joint_lower) < 1e-5 |
    abs(beta - joint_upper) < 1e-5
)
result <- bind_rows(result, joint_result)
write_csv(result, "tables/fledging_phi_temperature_conditional_survival.csv")
write_csv(
  tibble(
    all_broods = nrow(dat), zero_broods = sum(!positive),
    positive_broods = sum(positive)
  ),
  "tables/fledging_phi_temperature_conditional_survival_sample.csv"
)
print(filter(result, grepl("temp_mean_c$", term)), width = Inf)
