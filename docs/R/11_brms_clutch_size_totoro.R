#!/usr/bin/env Rscript

# Usage on TOTORO:
#   Rscript R/11_brms_clutch_size_totoro.R pilot
#   Rscript R/11_brms_clutch_size_totoro.R full

suppressPackageStartupMessages({
  library(brms)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args)) tolower(args[1]) else "pilot"
stopifnot(run_mode %in% c("pilot", "full"))

dir.create("models/brms", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

dat <- read_csv("data_derived/clutch_size_model_data.csv",
                show_col_types = FALSE) |>
  mutate(
    female_id = factor(female_id),
    year_f = factor(year_f),
    age_class = factor(age_class, levels = c("YOUNG", "OLD"))
  )

expected_n <- 12627L
if (nrow(dat) != expected_n) {
  stop("Model-data mismatch: expected ", expected_n,
       " rows (the exact glmmTMB sample), but found ", nrow(dat),
       ". Re-run and synchronise R/09, R/10, and the derived CSV before TOTORO.")
}
stopifnot(nlevels(dat$year_f) == 44L,
          !anyNA(dat[c("clutch", "rel_LD", "temp_mean_c", "temp_sd_c",
                       "year_c", "age_class", "female_id", "year_f")]))

formula_cs <- bf(
  clutch ~ rel_LD + temp_mean_c + temp_sd_c + year_c + age_class +
    (1 | female_id) + (1 | year_f),
  sigma ~ rel_LD + temp_mean_c + temp_sd_c + year_c + age_class +
    (1 | year_f)
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

fit_file <- file.path("models/brms", paste0("clutch_size_", run_mode))
start_time <- Sys.time()

m_cs <- brm(
  formula = formula_cs,
  data = dat,
  family = gaussian(),
  prior = priors_cs,
  backend = "cmdstanr",
  chains = settings$chains,
  cores = settings$cores,
  iter = settings$iter,
  warmup = settings$warmup,
  control = list(adapt_delta = settings$adapt_delta,
                 max_treedepth = settings$max_treedepth),
  seed = 1981,
  init = 0,
  save_pars = save_pars(all = TRUE),
  file = fit_file,
  file_refit = "on_change",
  refresh = 50
)

elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time,
                                       units = "mins"))
np <- nuts_params(m_cs)
divergences <- sum(np$Parameter == "divergent__" & np$Value == 1)
treedepth_hits <- sum(np$Parameter == "treedepth__" &
                        np$Value >= settings$max_treedepth)

all_draw_summary <- as_draws(m_cs) |> summarise_draws()
draw_summary <- all_draw_summary |>
  filter(grepl("sigma|sd_", variable))

run_diagnostics <- tibble(
  run_mode, chains = settings$chains, iter = settings$iter,
  warmup = settings$warmup, adapt_delta = settings$adapt_delta,
  max_treedepth = settings$max_treedepth,
  elapsed_minutes, divergences, treedepth_hits,
  max_rhat = max(all_draw_summary$rhat, na.rm = TRUE),
  min_bulk_ess = min(all_draw_summary$ess_bulk, na.rm = TRUE)
)

write_csv(run_diagnostics,
          file.path("tables", paste0("brms_cs_", run_mode,
                                      "_diagnostics.csv")))
write_csv(draw_summary,
          file.path("tables", paste0("brms_cs_", run_mode,
                                      "_sigma_draw_summary.csv")))
write_csv(all_draw_summary,
          file.path("tables", paste0("brms_cs_", run_mode,
                                      "_all_draw_summary.csv")))

# Primary posterior predictive check for the annual variance structure.
p_sd_year <- pp_check(m_cs, type = "stat_grouped", stat = "sd",
                      group = "year_f")
ggsave(file.path("figures", paste0("brms_cs_", run_mode,
                                    "_pp_sd_by_year.png")),
       p_sd_year, width = 10, height = 7, dpi = 300)

print(run_diagnostics)
print(draw_summary)

if (run_mode == "pilot" && divergences > 10) {
  warning("Pilot has >10 divergences. Keep adapt_delta >= 0.99 for the full run and inspect geometry before interpretation.")
}
if (treedepth_hits > 0) {
  warning("Transitions reached max_treedepth: ", treedepth_hits,
          ". Inspect geometry and runtime before the full run.")
}
