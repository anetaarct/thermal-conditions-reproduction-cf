#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

m <- readRDS("models/glmmTMB/fledged_count_hurdle_poisson.rds")
dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id), year_f = factor(year_f),
    year = as.numeric(as.character(year_f))
  )

year_center <- mean(dat$year - dat$year_c)
nd <- tibble(
  year = seq(min(dat$year), max(dat$year), length.out = 160),
  rel_LD = 0, temp_mean_c = 0,
  age_class = factor("YOUNG", levels = levels(dat$age_class)),
  female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
  year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f))
) |>
  mutate(year_c = year - year_center)

pr <- predict(
  m, newdata = nd, type = "response", se.fit = TRUE,
  re.form = NA, allow.new.levels = TRUE
)
d <- nd |>
  transmute(
    year, estimate = as.numeric(pr$fit), SE = as.numeric(pr$se.fit),
    CI_low = pmax(0, estimate - 1.96 * SE),
    CI_high = estimate + 1.96 * SE
  )
write_csv(d, "tables/fig_realised_fledged_year.csv")

blue <- "#0072B2"
p <- ggplot(d, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
              fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = "Predicted number fledged") +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))
ggsave("figures/Fig_realised_fledged_year.png", p,
       width = 6.4, height = 4.8, dpi = 400)
