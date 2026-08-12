#!/usr/bin/env Rscript

.libPaths(c(".Rlib", "../CF phenology/gotland-flycatcher-phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(glmmTMB)
  library(readr); library(tibble); library(tidyr)
})
m <- readRDS("models/glmmTMB/fledged_count_hurdle_poisson.rds")
dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(age_class = factor(age_class, levels = c("YOUNG", "OLD")),
         female_id = factor(female_id), year_f = factor(year_f),
         year = as.numeric(as.character(year_f)))
pred_df <- function(nd) {
  pr <- predict(m, newdata = nd, type = "response", se.fit = TRUE,
                re.form = NA, allow.new.levels = TRUE)
  bind_cols(nd, tibble(estimate = as.numeric(pr$fit), SE = as.numeric(pr$se.fit))) |>
    mutate(CI_low = pmax(0, estimate - 1.96 * SE), CI_high = estimate + 1.96 * SE)
}
ref <- function(n) tibble(
  rel_LD = rep(0, n), temp_mean_c = rep(0, n), year_c = rep(0, n),
  age_class = factor(rep("YOUNG", n), levels = levels(dat$age_class)),
  female_id = factor(rep(levels(dat$female_id)[1], n), levels = levels(dat$female_id)),
  year_f = factor(rep(levels(dat$year_f)[1], n), levels = levels(dat$year_f)))
blue <- "#0072B2"; orange <- "#D55E00"; age_cols <- c(YOUNG = blue, OLD = orange)
theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title = element_text(face = "bold"), legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))

tq <- quantile(dat$temp_mean_c, c(.05, .95))
tc <- read_csv("tables/fledging_bb_centring_constants.csv", show_col_types = FALSE)$temp_mean[1]
nd <- ref(160) |> mutate(temp_mean_c = seq(tq[1], tq[2], length.out = n()),
                        temp_mean = temp_mean_c + tc)
d <- pred_df(nd) |> select(temp_mean, estimate, SE, CI_low, CI_high)
write_csv(d, "tables/fig5_5_fledged_count_temperature.csv")
p <- ggplot(d, aes(temp_mean, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Mean spring temperature (deg C)", y = "Predicted number fledged") + theme_pub
ggsave("figures/Fig5_5_fledged_count_temperature.png", p, width = 6.4, height = 4.8, dpi = 400)

rq <- quantile(dat$rel_LD, c(.05, .95))
nd <- crossing(rel_LD = seq(rq[1], rq[2], length.out = 160),
               age_class = factor(levels(dat$age_class), levels = levels(dat$age_class))) |>
  mutate(temp_mean_c = 0, year_c = 0,
         female_id = factor(levels(dat$female_id)[1], levels = levels(dat$female_id)),
         year_f = factor(levels(dat$year_f)[1], levels = levels(dat$year_f)))
d <- pred_df(nd) |> select(rel_LD, age_class, estimate, SE, CI_low, CI_high)
write_csv(d, "tables/fig5_6_fledged_count_rellD_age.csv")
p <- ggplot(d, aes(rel_LD, estimate, colour = age_class, fill = age_class)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .16, colour = NA) +
  geom_line(linewidth = 1) + scale_colour_manual(values = age_cols) +
  scale_fill_manual(values = age_cols) +
  labs(x = "Laying date relative to annual median (days)", y = "Predicted number fledged") + theme_pub
ggsave("figures/Fig5_6_fledged_count_rellD_age.png", p, width = 6.4, height = 4.8, dpi = 400)

yc <- mean(dat$year - dat$year_c)
nd <- ref(160) |> mutate(year = seq(min(dat$year), max(dat$year), length.out = n()),
                         year_c = year - yc)
d <- pred_df(nd) |> select(year, estimate, SE, CI_low, CI_high)
write_csv(d, "tables/fig_realised_fledged_year.csv")
p <- ggplot(d, aes(year, estimate)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = scales::alpha(blue, .2)) +
  geom_line(colour = blue, linewidth = 1) +
  labs(x = "Year", y = "Predicted number fledged") + theme_pub
ggsave("figures/Fig_realised_fledged_year.png", p, width = 6.4, height = 4.8, dpi = 400)
