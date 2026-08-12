#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

blue <- "#0072B2"
orange <- "#D55E00"
age_cols <- c(YOUNG = blue, OLD = orange)
theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "top", legend.title = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = .3)
  )

d1 <- read_csv("tables/fig5_1_success_temperature.csv", show_col_types = FALSE)
p1 <- ggplot(d1, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted fledging success") +
  theme_pub
ggsave("figures/Fig5_1_success_temperature.png", p1, width = 6.4, height = 4.8, dpi = 400)

d2 <- read_csv("tables/fig5_2_success_rellD_age.csv", show_col_types = FALSE)
p2 <- ggplot(d2, aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = age_cols) + scale_fill_manual(values = age_cols) +
  labs(
    x = "Laying date relative to annual median (days)",
    y = "Predicted fledging success"
  ) + theme_pub
ggsave("figures/Fig5_2_success_rellD_age.png", p2, width = 6.4, height = 4.8, dpi = 400)

