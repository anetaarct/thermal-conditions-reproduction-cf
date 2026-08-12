#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(readr)})

blue <- "#0072B2"
orange <- "#D55E00"
theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"), legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = .3)
  )
trajectory_cols <- c(
  "Conditional (mean temperature)" = blue,
  "Observed annual temperature" = orange
)
a <- read_csv("tables/fig5_4a_success_year.csv", show_col_types = FALSE)
b <- read_csv("tables/fig5_4b_phi_year.csv", show_col_types = FALSE)
make_plot <- function(d, y_label) {
  ggplot(d, aes(year, estimate, colour = trajectory, fill = trajectory,
                group = trajectory)) +
    geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = trajectory_cols) +
    scale_fill_manual(values = trajectory_cols) +
    labs(x = "Year", y = y_label) + theme_pub
}
p1 <- make_plot(a, "Predicted fledging success")
p2 <- make_plot(b, expression("Consistency (" * phi * ")"))
legend <- cowplot::get_legend(p1 + theme(legend.position = "top"))
panels <- gridExtra::arrangeGrob(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), nrow = 1
)
combined <- gridExtra::arrangeGrob(
  legend, panels, ncol = 1, heights = c(.12, 1)
)
ggsave(
  "figures/Fig5_4_time_location_scale.png", combined,
  width = 10, height = 4.6, dpi = 400, bg = "white"
)
