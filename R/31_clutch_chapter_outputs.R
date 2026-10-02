#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggplot2)
  library(posterior)
  library(readr)
  library(tidyr)
})
Sys.setlocale("LC_TIME", "C")
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

model_file <- "models/brms/clutch_size_main_no_tempsd_full.rds"
stopifnot(file.exists(model_file))
m <- readRDS(model_file)
dat <- read_csv("data_derived/clutch_size_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f),
         year = as.integer(as.character(year_f)))

stopifnot(nrow(dat) > 12000L, n_distinct(dat$year_f) == 44L)

ps <- as.data.frame(posterior_summary(m, probs = c(.025, .975))) |>
  tibble::rownames_to_column("parameter") |>
  as_tibble() |>
  rename(mean = Estimate, SE = Est.Error, CI_low = Q2.5, CI_high = Q97.5)

keep <- ps |>
  filter(grepl("^b_", parameter) | parameter == "sd_year_f__sigma_Intercept") |>
  pull(parameter)
draw_diagnostics <- as_draws(m, variable = keep) |>
  summarise_draws("median", "rhat", "ess_bulk", "ess_tail") |>
  as_tibble() |>
  rename(parameter = variable, Rhat = rhat, Bulk_ESS = ess_bulk,
         Tail_ESS = ess_tail)

coef_table <- ps |>
  filter(parameter %in% keep) |>
  left_join(draw_diagnostics, by = "parameter") |>
  mutate(component = case_when(
    grepl("^b_sigma_", parameter) ~ "scale_log_sigma",
    parameter == "sd_year_f__sigma_Intercept" ~ "scale_random_SD",
    TRUE ~ "location"
  ), term = case_when(
    component == "scale_log_sigma" ~ sub("^b_sigma_", "", parameter),
    component == "location" ~ sub("^b_", "", parameter),
    TRUE ~ parameter
  )) |>
  select(component, term, parameter, median, mean, SE, CI_low, CI_high,
         Rhat, Bulk_ESS, Tail_ESS)
write_csv(coef_table, "tables/brms_cs_main_no_tempsd_coefficients.csv")

sample_table <- tibble(
  first_year = min(dat$year), last_year = max(dat$year),
  seasons = n_distinct(dat$year), broods = nrow(dat),
  females = n_distinct(dat$female_id), clutch_min = min(dat$clutch),
  clutch_max = max(dat$clutch), temp_mean_center = mean(dat$temp_mean),
  year_center = mean(dat$year)
)
write_csv(sample_table, "tables/clutch_chapter_sample_summary.csv")

summarise_draws_matrix <- function(x) {
  tibble(
    estimate = apply(x, 2, median),
    CI_low = apply(x, 2, quantile, .025),
    CI_high = apply(x, 2, quantile, .975)
  )
}

reference_grid <- function(n) tibble(
  rel_LD = rep(0, n), temp_mean_c = rep(0, n), year_c = rep(0, n),
  age_class = factor(rep("YOUNG", n), levels = levels(dat$age_class)),
  female_id = factor(rep(levels(dat$female_id)[1], n), levels = levels(dat$female_id)),
  year_f = factor(rep(levels(dat$year_f)[1], n), levels = levels(dat$year_f))
)

blue <- "#0072B2"; orange <- "#D55E00"
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))

# Figure 4.1: location response to mean temperature.
tq <- quantile(dat$temp_mean_c, c(.05, .95))
nd1 <- reference_grid(160) |>
  mutate(temp_mean_c = seq(tq[[1]], tq[[2]], length.out = n()),
         temp_mean = temp_mean_c + mean(dat$temp_mean))
d1 <- bind_cols(nd1, summarise_draws_matrix(
  posterior_epred(m, newdata = nd1, re_formula = NA)
)) |> select(temp_mean, estimate, CI_low, CI_high)
write_csv(d1, "tables/fig4_1_clutch_temperature.csv")
p1 <- ggplot(d1, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (°C)", y = "Predicted clutch size") + theme_pub
ggsave("figures/Fig4_1_clutch_temperature.png", p1, width = 6.4, height = 4.8, dpi = 400)

# Figure 4.2: location and scale responses to relative laying date.
rq <- quantile(dat$rel_LD, c(.05, .95))
nd2 <- tidyr::crossing(
  rel_LD = seq(rq[[1]], rq[[2]], length.out = 160),
  age_class = factor(c("YOUNG", "OLD"), levels = levels(dat$age_class))
) |>
  mutate(
    temp_mean_c = 0,
    year_c = 0,
    female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
    year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f))
  )
d2_mean <- bind_cols(nd2, summarise_draws_matrix(
  posterior_epred(m, newdata = nd2, re_formula = NA)
)) |> select(rel_LD, age_class, estimate, CI_low, CI_high)
write_csv(d2_mean, "tables/fig4_2a_clutch_rellD.csv")
d2 <- bind_cols(nd2, summarise_draws_matrix(
  posterior_linpred(m, newdata = nd2, dpar = "sigma", transform = TRUE,
                    re_formula = NA)
)) |> select(rel_LD, age_class, estimate, CI_low, CI_high)
write_csv(d2, "tables/fig4_2b_sigma_rellD.csv")
p2_mean <- ggplot(d2_mean, aes(rel_LD, estimate, colour = age_class,
                               fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16,
              colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(YOUNG = blue, OLD = orange), name = NULL) +
  scale_fill_manual(values = c(YOUNG = blue, OLD = orange), name = NULL) +
  labs(x = "Laying date relative to annual median (days)",
       y = "Predicted clutch size") +
  theme_pub + theme(legend.position = "top")
p2 <- ggplot(d2, aes(rel_LD, estimate, colour = age_class,
                     fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16,
              colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(YOUNG = blue, OLD = orange), name = NULL) +
  scale_fill_manual(values = c(YOUNG = blue, OLD = orange), name = NULL) +
  labs(x = "Laying date relative to annual median (days)",
       y = "Inconsistency (σ)") +
  theme_pub + theme(legend.position = "top")
png("figures/Fig4_2_clutch_rellD.png", width = 10, height = 4.6,
    units = "in", res = 400)
grid::grid.newpage(); layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p2_mean, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Figure 4.4: conditional predictions over time at study-wide mean temperature.
annual <- dat |> group_by(year, year_c, temp_mean_c) |>
  summarise(mean_clutch = mean(clutch), sd_clutch = sd(clutch), .groups = "drop")
nd_cond <- reference_grid(nrow(annual)) |>
  mutate(year_c = annual$year_c, year = annual$year, temp_mean_c = 0)
nd_obs <- nd_cond |> mutate(temp_mean_c = annual$temp_mean_c)

mean_cond_draws <- posterior_epred(m, newdata = nd_cond, re_formula = NA)
mean_obs_draws <- posterior_epred(m, newdata = nd_obs, re_formula = NA)
sigma_cond_draws <- posterior_linpred(
  m, newdata = nd_cond, dpar = "sigma", transform = TRUE, re_formula = NA)
sigma_obs_draws <- posterior_linpred(
  m, newdata = nd_obs, dpar = "sigma", transform = TRUE, re_formula = NA)

summarise_curve <- function(draws, label) {
  bind_cols(tibble(year = annual$year), summarise_draws_matrix(draws)) |>
    mutate(prediction = label)
}
d3_mean <- summarise_curve(mean_cond_draws, "Conditional: mean temperature")
d3_sd <- summarise_curve(sigma_cond_draws, "Conditional: mean temperature")
write_csv(left_join(d3_mean, annual, by = "year"), "tables/fig4_4a_annual_clutch.csv")
write_csv(left_join(d3_sd, annual, by = "year"), "tables/fig4_4b_annual_sigma.csv")

summarise_change <- function(x, component, prediction, unit, relative = FALSE) {
  change <- if (relative) 100 * (x[, ncol(x)] / x[, 1] - 1) else
    x[, ncol(x)] - x[, 1]
  tibble(component, prediction, estimate = median(change),
         CI_low = quantile(change, .025), CI_high = quantile(change, .975), unit)
}
time_change <- bind_rows(
  summarise_change(mean_cond_draws, "Mean clutch size", "Conditional: mean temperature", "eggs"),
  summarise_change(mean_obs_draws, "Mean clutch size", "Marginal: observed annual temperature", "eggs"),
  summarise_change(sigma_cond_draws, "Residual scale", "Conditional: mean temperature", "%", TRUE),
  summarise_change(sigma_obs_draws, "Residual scale", "Marginal: observed annual temperature", "%", TRUE))
write_csv(time_change, "tables/clutch_temporal_change_summary.csv")

p3a <- ggplot(d3_mean, aes(year, estimate)) +
  geom_point(data = annual, aes(year, mean_clutch), inherit.aes = FALSE,
             colour = "grey45", alpha = .65) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = blue, alpha = .13, colour = NA) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = "Mean clutch size") +
  theme_pub + theme(legend.position = "none")
p3b <- ggplot(d3_sd, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = blue, alpha = .13, colour = NA) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = expression("Conditional residual SD (" * sigma * ")")) +
  theme_pub + theme(legend.position = "none")

png("figures/Fig4_4_annual_clutch_mean_sd.png", width = 10, height = 4.6,
    units = "in", res = 400)
grid::grid.newpage(); layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p3a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p3b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Figure 4.4: significant age-class difference in mean clutch size.
nd4 <- reference_grid(2) |>
  mutate(age_class = factor(levels(dat$age_class), levels = levels(dat$age_class)))
d4 <- bind_cols(nd4, summarise_draws_matrix(
  posterior_epred(m, newdata = nd4, re_formula = NA)
)) |> select(age_class, estimate, CI_low, CI_high)
write_csv(d4, "tables/fig4_3a_clutch_age.csv")
age_cols <- c(YOUNG = blue, OLD = orange)
p4 <- ggplot(d4, aes(age_class, estimate, colour = age_class)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .12, linewidth = .7) +
  geom_point(size = 3) +
  scale_colour_manual(values = age_cols, guide = "none") +
  labs(x = "Female age class", y = "Predicted clutch size") + theme_pub
ggsave("figures/Fig4_3_clutch_age.png", p4, width = 5.4, height = 4.6, dpi = 400)

message("Clutch-size chapter outputs created.")
