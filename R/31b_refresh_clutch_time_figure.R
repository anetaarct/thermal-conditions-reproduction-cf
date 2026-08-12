#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

mean_data <- read_csv("tables/fig4_4a_annual_clutch.csv", show_col_types = FALSE)
sigma_data <- read_csv("tables/fig4_4b_annual_sigma.csv", show_col_types = FALSE)

blue <- "#0072B2"
orange <- "#D55E00"
theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = .3)
  )

p_mean <- ggplot(mean_data, aes(year, mean_clutch)) +
  geom_point(colour = "grey45", alpha = .65) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high, fill = prediction),
              alpha = .13, colour = NA) +
  geom_line(aes(y = estimate, colour = prediction), linewidth = 1) +
  scale_colour_manual(values = c("Conditional: mean temperature" = blue,
                                 "Marginal: observed annual temperature" = orange)) +
  scale_fill_manual(values = c("Conditional: mean temperature" = blue,
                               "Marginal: observed annual temperature" = orange)) +
  labs(x = "Year", y = "Annual mean clutch size", colour = NULL, fill = NULL) +
  theme_pub

p_sigma <- ggplot(sigma_data, aes(year, estimate, colour = prediction, fill = prediction)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .13, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Conditional: mean temperature" = blue,
                                 "Marginal: observed annual temperature" = orange)) +
  scale_fill_manual(values = c("Conditional: mean temperature" = blue,
                               "Marginal: observed annual temperature" = orange)) +
  labs(x = "Year", y = "Conditional residual SD", colour = NULL, fill = NULL) +
  theme_pub

png(
  "figures/Fig4_4_annual_clutch_mean_sd.png",
  width = 10, height = 4.6, units = "in", res = 400
)
grid::grid.newpage()
layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p_mean, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_sigma, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
