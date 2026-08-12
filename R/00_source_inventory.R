#!/usr/bin/env Rscript
.libPaths(c("../../repo_gotland_flycatcher_phenology/.r-lib", .libPaths()))
suppressPackageStartupMessages({library(readr);library(tibble)})
files <- list.files("data_raw",recursive=TRUE,full.names=TRUE)
info <- file.info(files)
out <- tibble(file=gsub("\\\\","/",files),bytes=info$size,
              modified=format(info$mtime,"%Y-%m-%d %H:%M:%S"))
write_csv(out,"tables/source_file_inventory.csv")
print(as.data.frame(out), row.names=FALSE)
