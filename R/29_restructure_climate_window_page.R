#!/usr/bin/env Rscript

path <- "climate_window.qmd"
x <- readLines(path, encoding = "UTF-8", warn = FALSE)

if (any(x == "## Appendix: diagnostics and rejected hypotheses")) {
  stop("The climate-window page has already been restructured.")
}

at <- function(heading) {
  z <- which(x == heading)
  stopifnot(length(z) == 1L)
  z
}

i_packages <- at("## Packages and reproducible files")
i_assumptions <- at("## Prespecified assumptions")
i_march <- at("## March temperature structure")
i_annual <- at("## Annual mean and variability in the main window")
i_overlap <- at("## Does the window extend into laying?")
i_cached <- at("## Cached model objects")

packages_body <- x[(i_packages + 1L):(i_assumptions - 1L)]
march_body <- x[(i_march + 1L):(i_annual - 1L)]
overlap_body <- x[(i_overlap + 1L):(i_cached - 1L)]
cached_body <- x[(i_cached + 1L):length(x)]

before_packages <- x[1L:(i_packages - 1L)]
assumptions_to_march <- x[i_assumptions:(i_march - 1L)]
annual_to_overlap <- x[i_annual:(i_overlap - 1L)]

march_pointer <- c(
  "Three candidate explanations for the failed out-of-sample test were examined — an increase in warm March episodes, a decline in the predictive value of March relative to April, and a linear change in thermal sensitivity through time. None was supported (see Appendix). The window reported here is unchanged by these diagnostics.",
  ""
)

overlap_pointer <- c(
  "In all 45 seasons the main thermal window closed before the annual fifth percentile of laying dates, so the window never overlaps the laying period (per-season table in Appendix).",
  ""
)

callout <- function(title, body) c(
  "::: {.callout-note collapse=\"true\"}",
  paste0("## ", title),
  body,
  ":::",
  ""
)

appendix <- c(
  "## Appendix: diagnostics and rejected hypotheses",
  "",
  callout("Packages and reproducible files", packages_body),
  callout("March temperature structure: three hypotheses tested and rejected", march_body),
  callout("Annual check that the window does not extend into laying", overlap_body),
  callout("Cached model objects", cached_body)
)

out <- c(
  before_packages,
  assumptions_to_march,
  march_pointer,
  annual_to_overlap,
  overlap_pointer,
  appendix
)

# The modelling data use annual-median centring; correct the inconsistent prose.
out <- sub(
  "`rel_LD = LD - mean\\(LD\\)` within each year",
  "`rel_LD = LD - median(LD)` within each year",
  out
)

writeLines(out, path, useBytes = TRUE)
message("Restructured climate_window.qmd and harmonised the rel_LD definition.")
