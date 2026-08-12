#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

blue <- "#0072B2"; orange <- "#D55E00"
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))

# Temperature: location and scale panels, including the unsupported scale effect.
d_temp_mean <- read_csv("tables/fig4_1_clutch_temperature.csv", show_col_types = FALSE)
coef_all <- read_csv("tables/brms_cs_main_no_tempsd_coefficients.csv", show_col_types = FALSE)
b0_temp_sigma <- coef_all |>
  filter(component == "scale_log_sigma", term == "Intercept") |> slice(1)
b_temp_sigma <- coef_all |>
  filter(component == "scale_log_sigma", term == "temp_mean_c") |> slice(1)
temp_center <- mean(d_temp_mean$temp_mean)
d_temp_sigma <- d_temp_mean |>
  transmute(temp_mean, x = temp_mean - temp_center,
            eta = b0_temp_sigma$median + b_temp_sigma$median * x,
            SE = sqrt(b0_temp_sigma$SE^2 + x^2 * b_temp_sigma$SE^2),
            estimate = exp(eta), CI_low = exp(eta - 1.96 * SE),
            CI_high = exp(eta + 1.96 * SE)) |>
  select(temp_mean, estimate, CI_low, CI_high)
write_csv(d_temp_sigma, "tables/fig4_1b_sigma_temperature.csv")
p_temp_mean <- ggplot(d_temp_mean, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted clutch size") + theme_pub
p_temp_sigma <- ggplot(d_temp_sigma, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(orange, .2)) +
  geom_line(colour = orange, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)",
       y = expression("Inconsistency (" * sigma * ")")) + theme_pub
png("figures/Fig4_1_clutch_temperature.png", width = 10, height = 4.6, units = "in", res = 400)
grid::grid.newpage(); lay <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = lay))
print(p_temp_mean, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_temp_sigma, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

d_mean <- read_csv("tables/fig4_cs_mean_timing_predictions.csv", show_col_types = FALSE) |>
  transmute(rel_LD, estimate, CI_low = lower, CI_high = upper)
d_sigma <- read_csv("tables/fig4_cs_sd_timing_predictions.csv", show_col_types = FALSE) |>
  transmute(rel_LD, estimate, CI_low = lower, CI_high = upper)
write_csv(d_mean, "tables/fig4_2a_clutch_rellD.csv")

p_mean <- ggplot(d_mean, aes(rel_LD, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Laying date relative to annual median (days)",
       y = "Predicted clutch size") + theme_pub
p_sigma <- ggplot(d_sigma, aes(rel_LD, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(orange, .2)) +
  geom_line(colour = orange, linewidth = 1) +
  labs(x = "Laying date relative to annual median (days)",
       y = expression("Inconsistency (" * sigma * ")")) + theme_pub
write_csv(d_sigma, "tables/fig4_2b_sigma_rellD.csv")
png("figures/Fig4_2_clutch_rellD.png", width = 10, height = 4.6, units = "in", res = 400)
grid::grid.newpage(); lay <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = lay))
print(p_mean, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sigma, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

coef <- coef_all |>
  filter(component == "location")
b0 <- coef |> filter(term == "Intercept") |> slice(1)
b_age <- coef |> filter(term == "age_classOLD") |> slice(1)
age <- tibble(
  age_class = factor(c("YOUNG", "OLD"), levels = c("YOUNG", "OLD")),
  estimate = c(b0$median, b0$median + b_age$median),
  SE = c(b0$SE, sqrt(b0$SE^2 + b_age$SE^2))
) |>
  mutate(CI_low = estimate - 1.96 * SE, CI_high = estimate + 1.96 * SE)
write_csv(age, "tables/fig4_3a_clutch_age.csv")
p_age_mean <- ggplot(age, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .10,
                linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(YOUNG = blue, OLD = orange), guide = "none") +
  labs(x = "Female age class", y = "Predicted clutch size") +
  theme_pub

b0_sigma <- read_csv("tables/brms_cs_main_no_tempsd_coefficients.csv", show_col_types = FALSE) |>
  filter(component == "scale_log_sigma", term == "Intercept") |> slice(1)
b_age_sigma <- read_csv("tables/brms_cs_main_no_tempsd_coefficients.csv", show_col_types = FALSE) |>
  filter(component == "scale_log_sigma", term == "age_classOLD") |> slice(1)
age_sigma <- tibble(
  age_class = factor(c("YOUNG", "OLD"), levels = c("YOUNG", "OLD")),
  eta = c(b0_sigma$median, b0_sigma$median + b_age_sigma$median),
  SE = c(b0_sigma$SE, sqrt(b0_sigma$SE^2 + b_age_sigma$SE^2))
) |>
  mutate(estimate = exp(eta), CI_low = exp(eta - 1.96 * SE),
         CI_high = exp(eta + 1.96 * SE))
write_csv(age_sigma, "tables/fig4_3b_sigma_age.csv")
p_age_sigma <- ggplot(age_sigma, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .10, linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(YOUNG = blue, OLD = orange), guide = "none") +
  labs(x = "Female age class", y = expression("Inconsistency (" * sigma * ")")) +
  theme_pub

png("figures/Fig4_3_clutch_age.png", width = 10, height = 4.6, units = "in", res = 400)
grid::grid.newpage(); lay <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = lay))
print(p_age_mean, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_age_sigma, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
