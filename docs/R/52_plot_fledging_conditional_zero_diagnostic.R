#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(readr)})
d <- read_csv("tables/fledging_zero_frequency_conditional_simulations.csv",
              show_col_types = FALSE)
s <- read_csv("tables/fledging_zero_frequency_conditional_summary.csv",
              show_col_types = FALSE)
p <- ggplot(d, aes(zero_proportion)) +
  geom_histogram(bins = 35, fill = "#0072B2", alpha = .7, colour = "white") +
  geom_vline(xintercept = s$observed, colour = "#D55E00", linewidth = 1.1) +
  labs(x = "Simulated proportion of broods with zero fledglings",
       y = "Simulated datasets") +
  theme_classic(base_size = 11)
ggsave("figures/Fig5_4_zero_frequency_diagnostic.png", p,
       width = 6.4, height = 4.8, dpi = 400, bg = "white")
