#!/usr/bin/env Rscript

.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(readr)
  library(tibble)
})

dir.create("models/glmmTMB", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

dat <- read_csv("data_derived/fledging_model_data.csv", show_col_types = FALSE) |>
  mutate(
    age_class = factor(age_class, levels = c("YOUNG", "OLD")),
    female_id = factor(female_id),
    year_f = factor(year_f)
  )
stopifnot(nrow(dat) == 9716L, nlevels(dat$year_f) == 38L)

ctrl <- glmmTMBControl(
  optCtrl = list(iter.max = 30000, eval.max = 30000), parallel = 1
)

# Current thermal model: temperature is present in location and precision.
m_thermal <- readRDS("models/glmmTMB/fledging_bb_main_no_tempsd.rds")

# Attribution test 1: remove temperature only from log(phi). This directly
# tests whether temperature absorbs the temporal signal in the precision part.
m_no_temp_phi <- update(
  m_thermal,
  dispformula = ~ rel_LD + year_c + age_class,
  data = dat,
  control = ctrl
)

# Attribution test 2: remove temperature from both model components, while
# retaining exactly the same sample, response, random effects, and age terms.
m_no_temp_anywhere <- update(
  m_thermal,
  formula = cbind(fledged, clutch - fledged) ~
    rel_LD * age_class + year_c + (1 | year_f) + (1 | female_id),
  dispformula = ~ rel_LD + year_c + age_class,
  data = dat,
  control = ctrl
)

stopifnot(
  isTRUE(m_thermal$sdr$pdHess),
  isTRUE(m_no_temp_phi$sdr$pdHess),
  isTRUE(m_no_temp_anywhere$sdr$pdHess),
  m_no_temp_phi$fit$convergence == 0L,
  m_no_temp_anywhere$fit$convergence == 0L
)

saveRDS(m_no_temp_phi,
        "models/glmmTMB/fledging_bb_no_temperature_in_phi.rds")
saveRDS(m_no_temp_anywhere,
        "models/glmmTMB/fledging_bb_no_temperature_anywhere.rds")

extract_phi_year <- function(model, specification) {
  tab <- as.data.frame(summary(model)$coefficients$disp) |>
    rownames_to_column("term") |>
    as_tibble() |>
    filter(term == "year_c")
  tibble(
    specification,
    beta_year_log_phi = tab$Estimate,
    SE = tab$`Std. Error`,
    z = tab$`z value`,
    p = tab$`Pr(>|z|)`,
    direction = if_else(tab$Estimate > 0,
                        "higher phi: greater consistency",
                        "lower phi: lower consistency"),
    AIC = AIC(model),
    n = nobs(model),
    pdHess = model$sdr$pdHess,
    convergence_code = model$fit$convergence
  )
}

result <- bind_rows(
  extract_phi_year(m_thermal,
                   "temperature in location and phi (current main model)"),
  extract_phi_year(m_no_temp_phi,
                   "temperature in location only"),
  extract_phi_year(m_no_temp_anywhere,
                   "temperature absent from location and phi")
)

write_csv(result, "tables/fledging_phi_time_temperature_attribution.csv")
print(result, n = Inf)
