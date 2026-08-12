#!/usr/bin/env Rscript
.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(lme4)})

d <- read_csv("data_derived/breeding_reproduction_clean.csv",
              col_types=cols(.default=col_character()),show_col_types=FALSE) |>
  transmute(
    female_id=factor(toupper(trimws(FRING))),
    year=factor(as.integer(YEAR)),
    age=as.numeric(F_EXACT_AGE_YEARS),
    LD=as.numeric(LD),
    CS=as.numeric(CS)
  ) |>
  filter(!is.na(female_id),!is.na(year),!is.na(age),age>=1) |>
  mutate(age_binary=factor(if_else(age==1,"YOUNG","OLD"),levels=c("YOUNG","OLD")),
         age_capped=factor(if_else(age>=5,"5+",as.character(as.integer(age))),
                           levels=c("1","2","3","4","5+")))

fit_response <- function(response) {
  z <- d |> filter(!is.na(.data[[response]])) |> droplevels()
  f_bin <- as.formula(paste(response,"~ age_binary + (1|female_id) + (1|year)"))
  f_cat <- as.formula(paste(response,"~ age_capped + (1|female_id) + (1|year)"))
  f_quad <- as.formula(paste(response,"~ poly(age,2) + (1|female_id) + (1|year)"))
  mods <- list(binary=lmer(f_bin,data=z,REML=FALSE),
               categorical_capped=lmer(f_cat,data=z,REML=FALSE),
               quadratic=lmer(f_quad,data=z,REML=FALSE))
  a <- do.call(AIC,unname(mods)); a$model<-names(mods);rownames(a)<-NULL
  a$delta_AIC<-a$AIC-min(a$AIC)
  a <- a |> transmute(response=response,model,df,AIC,delta_AIC,
                     records=nrow(z),females=n_distinct(z$female_id))
  b <- coef(summary(mods$binary))
  effect <- tibble(response=response,contrast="OLD minus YOUNG",
                   estimate=b["age_binaryOLD","Estimate"],SE=b["age_binaryOLD","Std. Error"],
                   t=b["age_binaryOLD","t value"],records=nrow(z),females=n_distinct(z$female_id))
  list(aic=a,effect=effect,models=mods)
}

ld <- fit_response("LD")
cs <- fit_response("CS")
aic <- bind_rows(ld$aic,cs$aic)
effects <- bind_rows(ld$effect,cs$effect)
saveRDS(list(LD=ld$models,CS=cs$models),"data_derived/age_classification_validation_models.rds")
write_csv(aic,"tables/age_classification_model_comparison.csv")
write_csv(effects,"tables/age_binary_effects.csv")

coverage <- d |> count(age,name="records") |> arrange(age)
write_csv(coverage,"tables/exact_age_validation_coverage.csv")

summary <- aic |> group_by(response) |> arrange(delta_AIC,.by_group=TRUE) |>
  summarise(best_model=first(model),binary_delta_AIC=delta_AIC[model=="binary"],
            categorical_delta_AIC=delta_AIC[model=="categorical_capped"],
            quadratic_delta_AIC=delta_AIC[model=="quadratic"],.groups="drop")
write_csv(summary,"tables/age_classification_validation_summary.csv")
print(as.data.frame(aic),row.names=FALSE)
print(as.data.frame(effects),row.names=FALSE)
