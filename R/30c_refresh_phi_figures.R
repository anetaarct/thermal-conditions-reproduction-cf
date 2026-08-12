#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

m <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")
dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id), year_f = factor(year_f),
    year = as.numeric(as.character(year_f))
  )
blue <- "#0072B2"
orange <- "#D55E00"
age_cols <- c(YOUNG = blue, OLD = orange)
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
pred_disp <- function(nd) {
  pr <- predict(
    m, newdata = nd, type = "disp", se.fit = TRUE,
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
pred_response <- function(nd) {
  pr <- predict(
    m, newdata = nd, type = "response", se.fit = TRUE,
    re.form = NA, allow.new.levels = TRUE
  )
  bind_cols(
    nd,
    tibble(
      estimate = as.numeric(pr$fit), SE = as.numeric(pr$se.fit),
      CI_low = pmax(0, estimate - 1.96 * SE),
      CI_high = pmin(1, estimate + 1.96 * SE)
    )
  )
}
temp_center <- read_csv(
  "tables/fledging_bb_centring_constants.csv", show_col_types = FALSE
)$temp_mean[1]
tq <- quantile(dat$temp_mean_c, c(.05, .95))
nd_temp <- ref_grid(160) |>
  mutate(
    temp_mean_c = seq(tq[1], tq[2], length.out = n()),
    temp_mean = temp_mean_c + temp_center
  )
d_phi_temp <- pred_disp(nd_temp) |>
  select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d_phi_temp, "tables/fig5_1b_phi_temperature.csv")

d_success <- pred_response(nd_temp) |>
  select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d_success, "tables/fig5_1_success_temperature.csv")
p_success <- ggplot(d_success, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted fledging success") +
  theme_pub
p_phi_temp <- ggplot(d_phi_temp, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(orange, .2)) +
  geom_line(colour = orange, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = expression("Consistency (" * phi * ")")) +
  theme_pub
png("figures/Fig5_1_success_phi_temperature.png", width = 10, height = 4.6,
    units = "in", res = 400)
grid::grid.newpage()
layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p_success, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_phi_temp, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

nd_age <- ref_grid(2) |>
  mutate(age_class = factor(levels(dat$age_class), levels = levels(dat$age_class)))
d_age_success <- pred_response(nd_age) |>
  select(age_class, estimate, SE, CI_low, CI_high)
write_csv(d_age_success, "tables/fig5_3a_success_age.csv")
d_age <- pred_disp(nd_age) |>
  select(age_class, estimate, SE, CI_low, CI_high)
write_csv(d_age, "tables/fig5_3b_phi_age.csv")
p_age_success <- ggplot(d_age_success, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .12, linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = age_cols, guide = "none") +
  labs(x = "Female age class", y = "Predicted fledging success") + theme_pub
p_age <- ggplot(d_age, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .12, linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = age_cols, guide = "none") +
  labs(x = "Female age class", y = expression("Consistency (" * phi * ")")) + theme_pub
png("figures/Fig5_3b_phi_age.png", width = 10, height = 4.6,
    units = "in", res = 400)
grid::grid.newpage()
layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p_age_success, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_age, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Relative laying date: location (by age) and precision in one figure.
ld_limits <- quantile(dat$rel_LD, c(.05, .95))
ld_seq <- seq(ld_limits[1], ld_limits[2], length.out = 160)
nd_ld_phi <- tidyr::crossing(
  rel_LD = ld_seq,
  age_class = factor(levels(dat$age_class), levels = levels(dat$age_class))
) |>
  mutate(
    temp_mean_c = 0, year_c = 0,
    female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
    year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f))
  )
d_ld_success <- pred_response(nd_ld_phi) |>
  select(rel_LD, age_class, estimate, SE, CI_low, CI_high)
write_csv(d_ld_success, "tables/fig5_2_success_rellD_age.csv")
d_ld_phi <- pred_disp(nd_ld_phi) |>
  select(rel_LD, age_class, estimate, SE, CI_low, CI_high)
write_csv(d_ld_phi, "tables/fig5_3_phi_rellD.csv")
p_ld_success <- ggplot(d_ld_success,
                       aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = age_cols) + scale_fill_manual(values = age_cols) +
  labs(x = "Laying date relative to annual median (days)",
       y = "Predicted fledging success") + theme_pub
p_ld_phi <- ggplot(d_ld_phi,
                   aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = age_cols) + scale_fill_manual(values = age_cols) +
  labs(x = "Laying date relative to annual median (days)",
       y = expression("Consistency (" * phi * ")")) + theme_pub
png("figures/Fig5_2_success_rellD_age.png", width = 10, height = 4.6,
    units = "in", res = 400)
grid::grid.newpage(); layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p_ld_success, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_ld_phi, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Time: conditional predictions with temperature held at its study-wide mean.
annual_temperature <- dat |>
  distinct(year, year_c, temp_mean_c) |>
  arrange(year)
nd_year <- ref_grid(nrow(annual_temperature)) |>
  mutate(
    year = annual_temperature$year,
    year_c = annual_temperature$year_c,
    temp_mean_c = 0,
    prediction = "Conditional: mean temperature"
  )
d_year_success <- pred_response(nd_year) |>
  select(year, temp_mean_c, prediction, estimate, SE, CI_low, CI_high)
d_year_phi <- pred_disp(nd_year) |>
  select(year, temp_mean_c, prediction, estimate, SE, CI_low, CI_high)
write_csv(d_year_success, "tables/fig5_4a_success_year.csv")
write_csv(d_year_phi, "tables/fig5_4b_phi_year.csv")
p_year_success <- ggplot(d_year_success, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = blue, alpha = .16, colour = NA) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = "Predicted fledging success") +
  theme_pub + theme(legend.position = "none")
p_year_phi <- ggplot(d_year_phi, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = blue, alpha = .16, colour = NA) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = expression("Consistency (" * phi * ")")) +
  theme_pub + theme(legend.position = "none")
png("figures/Fig5_4_time_location_scale.png", width = 10, height = 4.6,
    units = "in", res = 400, bg = "white")
grid::grid.newpage(); layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p_year_success, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_year_phi, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

phi_sim <- read_csv("tables/fledging_phi_constant_simulations.csv", show_col_types = FALSE) |>
  filter(convergence == 0, at_boundary == 0) |>
  select(simulation, temp_mean_c)
phi_summary <- read_csv(
  "tables/fledging_phi_temperature_constant_simulation_summary.csv",
  show_col_types = FALSE
)
write_csv(phi_sim, "tables/figA_phi_constant_simulation.csv")
p_hist <- ggplot(phi_sim, aes(temp_mean_c)) +
  geom_histogram(bins = 30, fill = scales::alpha(blue, .7), colour = "white") +
  geom_vline(
    xintercept = phi_summary$observed_auxiliary_beta,
    colour = orange, linewidth = 1.1
  ) +
  labs(
    x = "Temperature coefficient in log-phi under constant phi",
    y = "Simulated datasets"
  ) + theme_pub
ggsave("figures/FigA_phi_constant_simulation.png", p_hist,
       width = 6.4, height = 4.8, dpi = 400)
