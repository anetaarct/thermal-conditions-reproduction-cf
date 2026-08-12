#!/usr/bin/env Rscript
.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(brms); library(dplyr); library(ggplot2); library(posterior); library(readr); library(tibble)
})

m <- readRDS("models/brms/clutch_size_main_no_tempsd_full.rds")
dat <- read_csv("data_derived/clutch_size_model_data.csv", show_col_types = FALSE) |>
  mutate(year = as.integer(as.character(year_f)))
annual <- dat |>
  distinct(year, year_c, temp_mean_c) |>
  left_join(dat |> group_by(year) |>
              summarise(mean_clutch = mean(clutch), sd_clutch = sd(clutch), .groups = "drop"),
            by = "year") |>
  arrange(year)

dr <- as_draws_df(m, variable = c(
  "b_Intercept", "b_temp_mean_c", "b_year_c",
  "b_sigma_Intercept", "b_sigma_temp_mean_c", "b_sigma_year_c"))

curve_draws <- function(intercept, beta_temp, beta_year, temperatures) {
  outer(intercept, rep(1, nrow(annual))) +
    outer(beta_temp, temperatures) + outer(beta_year, annual$year_c)
}
mean_cond <- curve_draws(dr$b_Intercept, dr$b_temp_mean_c, dr$b_year_c,
                         rep(0, nrow(annual)))
mean_obs <- curve_draws(dr$b_Intercept, dr$b_temp_mean_c, dr$b_year_c,
                        annual$temp_mean_c)
sigma_cond <- exp(curve_draws(dr$b_sigma_Intercept, dr$b_sigma_temp_mean_c,
                              dr$b_sigma_year_c, rep(0, nrow(annual))))
sigma_obs <- exp(curve_draws(dr$b_sigma_Intercept, dr$b_sigma_temp_mean_c,
                             dr$b_sigma_year_c, annual$temp_mean_c))

summarise_curve <- function(x, label) tibble(
  year = annual$year, prediction = label,
  estimate = apply(x, 2, median),
  CI_low = apply(x, 2, quantile, .025),
  CI_high = apply(x, 2, quantile, .975))
mean_data <- bind_rows(
  summarise_curve(mean_cond, "Conditional: mean temperature"),
  summarise_curve(mean_obs, "Marginal: observed annual temperature")) |>
  left_join(annual, by = "year")
sigma_data <- bind_rows(
  summarise_curve(sigma_cond, "Conditional: mean temperature"),
  summarise_curve(sigma_obs, "Marginal: observed annual temperature")) |>
  left_join(annual, by = "year")
write_csv(mean_data, "tables/fig4_4a_annual_clutch.csv")
write_csv(sigma_data, "tables/fig4_4b_annual_sigma.csv")

summarise_change <- function(x, component, prediction, unit, relative = FALSE) {
  z <- if (relative) 100 * (x[, ncol(x)] / x[, 1] - 1) else x[, ncol(x)] - x[, 1]
  tibble(component, prediction, estimate = median(z),
         CI_low = quantile(z, .025), CI_high = quantile(z, .975), unit)
}
write_csv(bind_rows(
  summarise_change(mean_cond, "Mean clutch size", "Conditional: mean temperature", "eggs"),
  summarise_change(mean_obs, "Mean clutch size", "Marginal: observed annual temperature", "eggs"),
  summarise_change(sigma_cond, "Residual scale", "Conditional: mean temperature", "%", TRUE),
  summarise_change(sigma_obs, "Residual scale", "Marginal: observed annual temperature", "%", TRUE)),
  "tables/clutch_temporal_change_summary.csv")

blue <- "#0072B2"; orange <- "#D55E00"
cols <- c("Conditional: mean temperature" = blue,
          "Marginal: observed annual temperature" = orange)
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), legend.position = "top",
        legend.title = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))
p_mean <- ggplot(mean_data, aes(year, estimate, colour = prediction, fill = prediction)) +
  geom_point(data = annual, aes(year, mean_clutch), inherit.aes = FALSE,
             colour = "grey45", alpha = .65) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .13, colour = NA) +
  geom_line(linewidth = 1) + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) + labs(x = "Year", y = "Mean clutch size") + theme_pub
p_sigma <- ggplot(sigma_data, aes(year, estimate, colour = prediction, fill = prediction)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .13, colour = NA) +
  geom_line(linewidth = 1) + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) + labs(x = "Year", y = "Residual inconsistency (sigma)") + theme_pub
png("figures/Fig4_4_annual_clutch_mean_sd.png", width = 10, height = 4.8, units = "in", res = 400)
grid::grid.newpage(); layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p_mean, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sigma, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
