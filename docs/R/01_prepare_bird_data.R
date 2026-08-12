#!/usr/bin/env Rscript
.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble)})

norm_ring <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NAN")] <- NA_character_
  x
}

adult <- read_csv("data_raw/birds/CF_adult_birds_1980-2019.csv", col_types=cols(.default=col_character()), show_col_types=FALSE) |>
  mutate(source_row = row_number(), FRING = norm_ring(FRING), MRING = norm_ring(MRING),
         across(any_of(c("YEAR","LD","CS","HS","FLEDGE","RECRUIT","NEWRECRUIT")), as.numeric))
nest <- read_csv("data_raw/birds/CF_nestlings_1981-2019.csv", col_types=cols(.default=col_character()), show_col_types=FALSE) |>
  mutate(Ring = norm_ring(Ring), Year = as.integer(Year))

ring_year <- nest |>
  filter(!is.na(Ring), !is.na(Year)) |>
  distinct(Ring, Year) |>
  count(Ring, name="n_hatch_years") |>
  left_join(nest |> filter(!is.na(Ring), !is.na(Year)) |> group_by(Ring) |>
              summarise(F_EXACT_HATCH_YEAR=min(Year), .groups="drop"), by="Ring") |>
  mutate(F_EXACT_HATCH_YEAR=if_else(n_hatch_years==1L,F_EXACT_HATCH_YEAR,NA_integer_))

prepared <- adult |>
  left_join(ring_year |> select(FRING=Ring,F_EXACT_HATCH_YEAR,n_hatch_years), by="FRING") |>
  mutate(
    F_EXACT_AGE_YEARS = if_else(!is.na(F_EXACT_HATCH_YEAR), YEAR-F_EXACT_HATCH_YEAR, NA_real_),
    F_AGE_EXACT_KNOWN = !is.na(F_EXACT_AGE_YEARS),
    F_AGE_ZERO_FLAG = F_EXACT_AGE_YEARS==0,
    invalid_HS_above_CS = !is.na(HS)&!is.na(CS)&HS>CS,
    invalid_FLEDGE_above_HS = !is.na(FLEDGE)&!is.na(HS)&FLEDGE>HS,
    logical_record = !invalid_HS_above_CS & !invalid_FLEDGE_above_HS
  )

flags <- prepared |> filter(!logical_record | F_AGE_ZERO_FLAG | n_hatch_years>1)
clean <- prepared |> filter(logical_record)

summary <- tibble(
  quantity=c("source breeding rows","logically consistent rows","HS above CS",
             "FLEDGE above HS","rows with confirmed female age","unique females with confirmed age",
             "breeding records with exact age zero","rings linked to multiple nestling years"),
  value=c(nrow(prepared),nrow(clean),sum(prepared$invalid_HS_above_CS,na.rm=TRUE),
          sum(prepared$invalid_FLEDGE_above_HS,na.rm=TRUE),sum(prepared$F_AGE_EXACT_KNOWN,na.rm=TRUE),
          n_distinct(prepared$FRING[prepared$F_AGE_EXACT_KNOWN]),sum(prepared$F_AGE_ZERO_FLAG,na.rm=TRUE),
          n_distinct(prepared$FRING[prepared$n_hatch_years>1],na.rm=TRUE))
)

write_csv(prepared,"data_derived/breeding_reproduction_flagged.csv",na="")
write_csv(clean,"data_derived/breeding_reproduction_clean.csv",na="")
write_csv(flags,"tables/bird_data_quality_flags.csv",na="")
write_csv(summary,"tables/bird_data_preparation_summary.csv")
age_summary <- tibble(
  quantity=c("breeding records with confirmed exact age","unique females with confirmed exact age",
             "minimum confirmed age","maximum confirmed age","records with age zero",
             "records with age at least one"),
  value=c(sum(prepared$F_AGE_EXACT_KNOWN,na.rm=TRUE),
          n_distinct(prepared$FRING[prepared$F_AGE_EXACT_KNOWN]),
          min(prepared$F_EXACT_AGE_YEARS,na.rm=TRUE),max(prepared$F_EXACT_AGE_YEARS,na.rm=TRUE),
          sum(prepared$F_AGE_ZERO_FLAG,na.rm=TRUE),sum(prepared$F_EXACT_AGE_YEARS>=1,na.rm=TRUE))
)
write_csv(age_summary,"tables/exact_female_age_summary.csv")
print(as.data.frame(summary), row.names=FALSE)
