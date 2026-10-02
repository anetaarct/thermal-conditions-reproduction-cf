#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

dir.create("models/glmmTMB", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

raw <- read_csv("data_derived/reproductive_model_data.csv", show_col_types = FALSE) |>
  filter(between(year, 1981, 2025), year != 2020,
         age_class %in% c("YOUNG", "OLD")) |>
  transmute(
    lay_date = as.numeric(laying_date),
    temp_mean = as.numeric(temp_mean),
    year = as.integer(year),
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = if_else(is.na(female_id) | female_id == "",
                        paste0("UNRINGED_", source_row), female_id),
    year_f = factor(year)
  ) |>
  filter(if_all(everything(), ~ !is.na(.x))) |>
  droplevels()

centres <- raw |>
  summarise(temp_mean = mean(temp_mean), year = mean(year))

dat <- raw |>
  mutate(
    temp_mean_c = temp_mean - centres$temp_mean,
    year_c = year - centres$year,
    female_id = factor(female_id)
  ) |>
  select(lay_date, temp_mean_c, year_c, age_class, female_id, year_f)

observed_years <- sort(unique(as.integer(as.character(dat$year_f))))
stopifnot(
  length(observed_years) == 44L,
  identical(setdiff(1981:2025, observed_years), 2020L),
  abs(mean(dat$temp_mean_c)) < 1e-10,
  abs(mean(dat$year_c)) < 1e-10
)

write_csv(dat, "data_derived/laying_date_model_data.csv")
write_csv(centres, "tables/laying_date_age_temperature_centres.csv")

model <- glmmTMB(
  lay_date ~ temp_mean_c * age_class + year_c +
    (1 | year_f) + (1 | female_id),
  dispformula = ~ age_class + temp_mean_c,
  family = gaussian(),
  data = dat,
  control = glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000))
)

saveRDS(model, "models/glmmTMB/laying_date_age_temperature.rds")

stopifnot(model$fit$convergence == 0L, isTRUE(model$sdr$pdHess))

coefficient_table <- function(model) {
  sm <- summary(model)$coefficients
  bind_rows(
    as.data.frame(sm$cond) |>
      rownames_to_column("term") |>
      as_tibble() |>
      transmute(component = "location", term, beta = Estimate,
                SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`),
    as.data.frame(sm$disp) |>
      rownames_to_column("term") |>
      as_tibble() |>
      transmute(component = "scale_log_variance", term, beta = Estimate,
                SE = `Std. Error`, z = `z value`, p = `Pr(>|z|)`)
  )
}

coefficients <- coefficient_table(model)
write_csv(coefficients, "tables/laying_date_age_temperature_coefficients.csv")

# Conditional temperature slopes from the joint fixed-effect covariance matrix.
b <- fixef(model)$cond
V <- vcov(model)$cond
main_name <- "temp_mean_c"
int_name <- grep("temp_mean_c:age_classOLD|age_classOLD:temp_mean_c",
                 names(b), value = TRUE)
stopifnot(length(int_name) == 1L)

linear_combo <- function(weights, label) {
  est <- sum(weights * b[names(weights)])
  se <- sqrt(drop(t(weights) %*% V[names(weights), names(weights), drop = FALSE] %*% weights))
  tibble(contrast = label, estimate = est, SE = se,
         CI_low = est - 1.96 * se, CI_high = est + 1.96 * se,
         z = est / se, p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
}

slopes <- bind_rows(
  linear_combo(setNames(1, main_name), "YOUNG slope"),
  linear_combo(setNames(c(1, 1), c(main_name, int_name)), "OLD slope"),
  linear_combo(setNames(1, int_name), "OLD - YOUNG slope difference")
)
write_csv(slopes, "tables/laying_date_temperature_slopes_by_age.csv")

# Plot A: temperature slopes within the observed range, fixed effects only.
temp_range_c <- range(dat$temp_mean_c)
grid_temp <- expand.grid(
  temp_mean_c = seq(temp_range_c[1], temp_range_c[2], length.out = 160),
  age_class = levels(dat$age_class),
  year_c = 0,
  female_id = levels(dat$female_id)[1],
  year_f = levels(dat$year_f)[1]
) |>
  mutate(age_class = factor(age_class, levels = levels(dat$age_class)))
pred_temp <- predict(model, newdata = grid_temp, type = "response", se.fit = TRUE,
                     re.form = NA, allow.new.levels = TRUE)
plot_a_data <- grid_temp |>
  mutate(
    temp_mean = temp_mean_c + centres$temp_mean,
    estimate = as.numeric(pred_temp$fit),
    SE = as.numeric(pred_temp$se.fit),
    CI_low = estimate - 1.96 * SE,
    CI_high = estimate + 1.96 * SE
  ) |>
  select(temp_mean, age_class, estimate, SE, CI_low, CI_high)
write_csv(plot_a_data, "tables/laying_date_figure_A_temperature.csv")

# Plot B: predicted location at centred temperature and year.
grid_age <- tibble(
  temp_mean_c = 0, year_c = 0,
  age_class = factor(levels(dat$age_class), levels = levels(dat$age_class)),
  female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
  year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f))
)
pred_age <- predict(model, newdata = grid_age, type = "response", se.fit = TRUE,
                    re.form = NA, allow.new.levels = TRUE)
plot_b_data <- grid_age |>
  mutate(estimate = as.numeric(pred_age$fit), SE = as.numeric(pred_age$se.fit),
         CI_low = estimate - 1.96 * SE, CI_high = estimate + 1.96 * SE) |>
  select(age_class, estimate, SE, CI_low, CI_high)
write_csv(plot_b_data, "tables/laying_date_figure_B_age_location.csv")

# Plot C: Gaussian dispersion is variance on exp(eta); residual SD = exp(eta/2).
bd <- fixef(model)$disp
Vd <- vcov(model)$disp
scale_for_age <- function(age) {
  w <- setNames(rep(0, length(bd)), names(bd)); w["(Intercept)"] <- 1
  age_term <- grep("age_classOLD", names(bd), value = TRUE)
  if (age == "OLD" && length(age_term) == 1L) w[age_term] <- 1
  eta <- sum(w * bd)
  se_eta <- sqrt(drop(t(w) %*% Vd %*% w))
  tibble(
    age_class = age,
    log_variance = eta,
    SE_log_variance = se_eta,
    sigma = exp(eta / 2),
    CI_low = exp((eta - 1.96 * se_eta) / 2),
    CI_high = exp((eta + 1.96 * se_eta) / 2)
  )
}
plot_c_data <- bind_rows(scale_for_age("YOUNG"), scale_for_age("OLD")) |>
  mutate(age_class = factor(age_class, levels = levels(dat$age_class)))
write_csv(plot_c_data, "tables/laying_date_figure_C_age_scale.csv")

age_cols <- c(YOUNG = "#0072B2", OLD = "#D55E00")
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        legend.position = "top", legend.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

p_a <- ggplot(plot_a_data, aes(temp_mean, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.17, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = age_cols, labels = c("Young (1 year)", "Old (>=2 years)")) +
  scale_fill_manual(values = age_cols, labels = c("Young (1 year)", "Old (>=2 years)")) +
  labs(x = "Mean spring temperature (°C)", y = "Predicted laying date\n(days after 30 April)") +
  theme_pub

p_b <- ggplot(plot_b_data, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.12, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = age_cols, guide = "none") +
  scale_x_discrete(labels = c("Young\n(1 year)", "Old\n(>=2 years)")) +
  labs(x = "Female age class", y = "Predicted laying date\n(days after 30 April)") +
  theme_pub

p_c <- ggplot(plot_c_data, aes(age_class, sigma, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.12, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = age_cols, guide = "none") +
  scale_x_discrete(labels = c("Young\n(1 year)", "Old\n(>=2 years)")) +
  labs(x = "Female age class", y = "Inconsistency in laying date (σ)") +
  theme_pub

ggsave("figures/laying_date_A_temperature_by_age.png", p_a, width = 6.2, height = 4.8, dpi = 400)
ggsave("figures/laying_date_B_location_by_age.png", p_b, width = 4.5, height = 4.8, dpi = 400)
ggsave("figures/laying_date_C_scale_by_age.png", p_c, width = 4.5, height = 4.8, dpi = 400)
ggsave("figures/laying_date_A_temperature_by_age.pdf", p_a, width = 6.2, height = 4.8, device = cairo_pdf)
ggsave("figures/laying_date_B_location_by_age.pdf", p_b, width = 4.5, height = 4.8, device = cairo_pdf)
ggsave("figures/laying_date_C_scale_by_age.pdf", p_c, width = 4.5, height = 4.8, device = cairo_pdf)

write_csv(tibble(
  n = nrow(dat), seasons = nlevels(dat$year_f), first_year = min(observed_years),
  last_year = max(observed_years), missing_years = paste(setdiff(1981:2025, observed_years), collapse = ","),
  pdHess = model$sdr$pdHess, convergence_code = model$fit$convergence
), "tables/laying_date_age_temperature_diagnostics.csv")

message("Completed laying-date age-temperature location-scale analysis.")
