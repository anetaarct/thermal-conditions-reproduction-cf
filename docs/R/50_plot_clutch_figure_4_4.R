#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(readr)})

blue <- "#0072B2"
theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = .3),
    legend.position = "none"
  )
a <- read_csv("tables/fig4_4a_annual_clutch.csv", show_col_types = FALSE)
b <- read_csv("tables/fig4_4b_annual_sigma.csv", show_col_types = FALSE)
a <- a[grepl("^Conditional", a$prediction), ]
b <- b[grepl("^Conditional", b$prediction), ]
write_csv(a, "tables/fig4_4a_annual_clutch.csv")
write_csv(b, "tables/fig4_4b_annual_sigma.csv")
stopifnot(
  length(unique(a$prediction)) == 1L,
  length(unique(b$prediction)) == 1L,
  !any(grepl("Marginal|observed annual", a$prediction, ignore.case = TRUE)),
  !any(grepl("Marginal|observed annual", b$prediction, ignore.case = TRUE))
)
p1 <- ggplot(a, aes(year, estimate)) +
  geom_point(aes(y = mean_clutch), colour = "grey45", alpha = .65) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = blue, alpha = .13, colour = NA) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = "Mean clutch size") + theme_pub
p2 <- ggplot(b, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = blue, alpha = .13, colour = NA) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = expression("Conditional residual SD (" * sigma * ")")) +
  theme_pub

png("figures/Fig4_4_annual_clutch_mean_sd.png", width = 10, height = 4.6,
    units = "in", res = 400, bg = "white")
grid::grid.newpage()
layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
