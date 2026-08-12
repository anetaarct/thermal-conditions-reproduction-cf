#!/usr/bin/env Rscript
.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib",
            "../../repo_gotland_flycatcher_phenology/Rlib", .libPaths()))
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble)})

read_smhi <- function(path, station) {
  x <- read.delim(path, sep=";", skip=10, fileEncoding="UTF-8-BOM",
                  check.names=FALSE, stringsAsFactors=FALSE)
  tibble(date=as.Date(x[[3]]), temperature=as.numeric(gsub(",",".",x[[4]])),
         quality=x[[5]], station=station) |>
    filter(!is.na(date), !is.na(temperature))
}

old <- read_smhi("data_raw/smhi/hoburg_68550_daily_mean_temperature.csv","Hoburg 68550")
new <- read_smhi("data_raw/smhi/hoburg_a_68560_daily_mean_temperature.csv","Hoburg A 68560")

daily <- full_join(old |> select(date,temp_old=temperature,q_old=quality),
                   new |> select(date,temp_new=temperature,q_new=quality),by="date") |>
  mutate(temperature=rowMeans(cbind(temp_old,temp_new),na.rm=TRUE),
         source=case_when(!is.na(temp_old)&!is.na(temp_new)~"mean of both stations",
                          !is.na(temp_old)~"Hoburg 68550",TRUE~"Hoburg A 68560"),
         quality=case_when(q_old=="G"|q_new=="G"~"G",TRUE~coalesce(q_new,q_old))) |>
  arrange(date)

overlap <- daily |>
  filter(!is.na(temp_old), !is.na(temp_new)) |>
  mutate(difference_old_minus_new = temp_old - temp_new)

overlap_summary <- overlap |>
  summarise(
    overlap_days = n(),
    first_overlap = as.character(min(date)),
    last_overlap = as.character(max(date)),
    mean_difference_old_minus_new = mean(difference_old_minus_new),
    sd_difference = sd(difference_old_minus_new),
    rmse = sqrt(mean(difference_old_minus_new^2)),
    correlation = cor(temp_old, temp_new)
  )

calendar <- tibble(date=seq(as.Date("1980-01-01"),as.Date("2025-12-31"),by="day")) |>
  left_join(daily,by="date") |>
  mutate(year=as.integer(format(date,"%Y")),doy=as.integer(format(date,"%j")))

coverage <- calendar |> group_by(year) |>
  summarise(days=n(),observed=sum(!is.na(temperature)),coverage=observed/days,.groups="drop")
summary <- tibble(quantity=c("first date","last date","daily records","observed daily temperatures",
                             "minimum annual coverage","median annual coverage"),
                  value=c(as.character(min(calendar$date)),as.character(max(calendar$date)),nrow(calendar),
                          sum(!is.na(calendar$temperature)),round(min(coverage$coverage),4),round(median(coverage$coverage),4)))

write_csv(calendar,"data_derived/smhi_hoburg_daily_temperature.csv",na="")
write_csv(coverage,"tables/temperature_coverage_by_year.csv")
write_csv(summary,"tables/temperature_data_summary.csv")
write_csv(overlap_summary,"tables/smhi_hoburg_station_overlap_summary.csv")
print(as.data.frame(summary), row.names=FALSE)
print(as.data.frame(overlap_summary), row.names=FALSE)
