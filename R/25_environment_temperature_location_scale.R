#!/usr/bin/env Rscript
Sys.setlocale('LC_CTYPE', 'English_United States.utf8')

.libPaths(c(".Rlib", "../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tibble)
  library(patchwork)
})

dir.create("models/glmmTMB", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

annual <- read_csv("data_derived/annual_main_window_temperature.csv", show_col_types = FALSE) |>
  filter(between(year, 1980, 2025)) |>
  transmute(year = as.integer(year), temp_mean = as.numeric(temp_mean)) |>
  filter(if_all(everything(), ~ !is.na(.x))) |>
  arrange(year) |>
  mutate(year_c = year - mean(year))

stopifnot(nrow(annual) == 46L, identical(annual$year, 1980:2025),
          abs(mean(annual$year_c)) < 1e-10)

model <- glmmTMB(
  temp_mean ~ year_c,
  dispformula = ~ year_c,
  family = gaussian(),
  data = annual,
  control = glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000))
)

stopifnot(model$fit$convergence == 0L, isTRUE(model$sdr$pdHess))
saveRDS(model, "models/glmmTMB/environment_temperature_location_scale.rds")

sm <- summary(model)$coefficients
coefficients <- bind_rows(
  as.data.frame(sm$cond) |>
    rownames_to_column("term") |>
    as_tibble() |>
    transmute(component = "location", term, beta = Estimate, SE = `Std. Error`,
              z = `z value`, p = `Pr(>|z|)`, interpretation = case_when(
                term == "year_c" & beta > 0 ~ "mean spring temperature increases through time",
                term == "year_c" & beta < 0 ~ "mean spring temperature decreases through time",
                TRUE ~ "intercept"
              )),
  as.data.frame(sm$disp) |>
    rownames_to_column("term") |>
    as_tibble() |>
    transmute(component = "scale_log_sigma", term, beta = Estimate, SE = `Std. Error`,
              z = `z value`, p = `Pr(>|z|)`, interpretation = case_when(
                term == "year_c" & `Pr(>|z|)` >= 0.05 ~ "no clear temporal trend in residual interannual variability",
                term == "year_c" & beta > 0 ~ "residual interannual variability around the mean trend increases",
                term == "year_c" & beta < 0 ~ "residual interannual variability around the mean trend decreases",
                TRUE ~ "baseline log residual SD"
              ))
)
coefficients <- coefficients |> mutate(CI_low = beta - 1.96 * SE,
                                        CI_high = beta + 1.96 * SE)
write_csv(coefficients, "tables/environment_temperature_location_scale_coefficients.csv")

grid <- tibble(
  year = seq(min(annual$year), max(annual$year), length.out = 240),
  year_c = year - mean(annual$year)
)
pred <- predict(model, newdata = grid, type = "response", se.fit = TRUE,
                re.form = NA)
plot_data <- grid |>
  mutate(estimate = as.numeric(pred$fit), SE = as.numeric(pred$se.fit),
         CI_low = estimate - 1.96 * SE, CI_high = estimate + 1.96 * SE)
write_csv(plot_data, "tables/environment_temperature_time_predictions.csv")

# glmmTMB 1.1.14 models Gaussian dispersion as log residual SD directly.
bd <- fixef(model)$disp
Vd <- vcov(model)$disp
scale_values <- lapply(grid$year_c, function(yc) {
  w <- c("(Intercept)" = 1, "year_c" = yc)
  eta <- sum(w[names(bd)] * bd)
  se_eta <- sqrt(drop(t(w[names(bd)]) %*% Vd %*% w[names(bd)]))
  tibble(log_sigma = eta, SE_log_sigma = se_eta,
         sigma = exp(eta),
         CI_low = exp(eta - 1.96 * se_eta),
         CI_high = exp(eta + 1.96 * se_eta))
}) |>
  bind_rows()
scale_data <- bind_cols(select(grid, year, year_c), scale_values)
write_csv(scale_data, "tables/environment_temperature_scale_predictions.csv")

blue <- "#0072B2"
p <- ggplot(annual, aes(year, temp_mean)) +
  geom_ribbon(data = plot_data, aes(x = year, ymin = CI_low, ymax = CI_high),
              inherit.aes = FALSE, fill = scales::alpha(blue, 0.20)) +
  geom_line(data = plot_data, aes(x = year, y = estimate), inherit.aes = FALSE,
            colour = blue, linewidth = 1) +
  geom_point(shape = 21, fill = "white", colour = "grey25", size = 2.2,
             stroke = 0.55) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  labs(x = "Year", y = "Mean spring temperature (°C)",
       subtitle = "Thermal window: 26 March-7 May") +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

ggsave("figures/environment_temperature_time.png", p, width = 7.2, height = 4.8, dpi = 400)
ggsave("figures/environment_temperature_time.pdf", p, width = 7.2, height = 4.8,
       device = cairo_pdf)

p_scale <- ggplot(scale_data, aes(year, sigma)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = orange <- "#D55E00", alpha = 0.2) +
  geom_line(colour = orange, linewidth = 1) +
  labs(x = "Year", y = "Interannual residual SD (σ)") + theme_classic(base_size = 11)
combined <- (p + labs(subtitle = NULL)) + p_scale +
  plot_annotation(tag_levels = "A")
ggsave("figures/Fig2_temperature_location_scale_time.png", combined,
       width = 8.4, height = 4.2, dpi = 400)
ggsave("figures/environment_temperature_scale.png", p_scale,
       width = 6.4, height = 4.8, dpi = 400)

write_csv(tibble(
  n = nrow(annual), first_year = min(annual$year), last_year = max(annual$year),
  missing_years = paste(setdiff(1980:2025, annual$year), collapse = ","),
  year_centre = mean(annual$year), pdHess = model$sdr$pdHess,
  convergence_code = model$fit$convergence
), "tables/environment_temperature_location_scale_diagnostics.csv")

message("Completed annual temperature location-scale model.")
