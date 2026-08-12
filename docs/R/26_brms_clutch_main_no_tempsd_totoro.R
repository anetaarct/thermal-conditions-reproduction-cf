#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms); library(readr); library(dplyr); library(ggplot2); library(posterior); library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args)) tolower(args[1]) else "pilot"
stopifnot(run_mode %in% c("pilot", "full"))
dir.create("models/brms", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

dat <- read_csv("data_derived/clutch_size_model_data.csv", show_col_types = FALSE) |>
  mutate(female_id = factor(female_id), year_f = factor(year_f),
         age_class = factor(age_class, levels = c("YOUNG", "OLD")))
stopifnot(nrow(dat) == 12627L, nlevels(dat$year_f) == 44L,
          !anyNA(dat[c("clutch", "rel_LD", "temp_mean_c", "year_c",
                       "age_class", "female_id", "year_f")]))

formula_cs <- bf(
  clutch ~ rel_LD + temp_mean_c + year_c + age_class +
    (1 | female_id) + (1 | year_f),
  sigma ~ rel_LD + temp_mean_c + year_c + age_class + (1 | year_f)
)

priors_cs <- c(
  prior(normal(6, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(normal(-0.2, 0.5), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.3), class = "b", dpar = "sigma"),
  prior(exponential(5), class = "sd", dpar = "sigma")
)

settings <- if (run_mode == "pilot") {
  list(iter = 500, warmup = 250, chains = 2, cores = 2,
       adapt_delta = 0.95, max_treedepth = 12)
} else {
  list(iter = 4000, warmup = 2000, chains = 4, cores = 4,
       adapt_delta = 0.99, max_treedepth = 15)
}

fit_file <- file.path("models/brms", paste0("clutch_size_main_no_tempsd_", run_mode))
start_time <- Sys.time()
m_cs_main <- brm(
  formula_cs, data = dat, family = gaussian(), prior = priors_cs,
  backend = "cmdstanr", chains = settings$chains, cores = settings$cores,
  iter = settings$iter, warmup = settings$warmup,
  control = list(adapt_delta = settings$adapt_delta,
                 max_treedepth = settings$max_treedepth),
  seed = 1981, init = 0, save_pars = save_pars(all = TRUE),
  file = fit_file, file_refit = "on_change", refresh = 50
)

elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
np <- nuts_params(m_cs_main)
all_summary <- as_draws(m_cs_main) |> summarise_draws()
diagnostics <- tibble(
  run_mode, elapsed_minutes,
  divergences = sum(np$Parameter == "divergent__" & np$Value == 1),
  treedepth_hits = sum(np$Parameter == "treedepth__" & np$Value >= settings$max_treedepth),
  max_rhat = max(all_summary$rhat, na.rm = TRUE),
  min_bulk_ess = min(all_summary$ess_bulk, na.rm = TRUE)
)
write_csv(diagnostics, paste0("tables/brms_cs_main_no_tempsd_", run_mode, "_diagnostics.csv"))
write_csv(all_summary, paste0("tables/brms_cs_main_no_tempsd_", run_mode, "_draw_summary.csv"))

post95 <- posterior_summary(m_cs_main, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rownames_to_column("variable") |>
  as_tibble()
write_csv(post95, paste0("tables/brms_cs_main_no_tempsd_", run_mode, "_posterior_95.csv"))

p <- pp_check(m_cs_main, type = "stat_grouped", stat = "sd", group = "year_f")
ggsave(paste0("figures/brms_cs_main_no_tempsd_", run_mode, "_pp_sd_by_year.png"),
       p, width = 10, height = 7, dpi = 300)
print(diagnostics)
