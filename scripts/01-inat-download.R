### ============================================================
### iNaturalist Data Download
### ============================================================
### Downloads research-grade observations from iNaturalist
### for specified taxa, saves results as CSV/RDS.
###
### Usage:
###   source("scripts/01_inat_download.R")
###   df <- inat_download("Grona", maxresults = 10000)
###   df <- inat_download("Medicago lupulina", query = "flower")
### ============================================================

library(rinat)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

#' Download iNaturalist observations for a taxon
#'
#' @param taxon_name  Character. Taxon to query (e.g. "Grona", "Medicago lupulina").
#' @param maxresults  Integer. Maximum observations to retrieve (default 10000).
#' @param quality     Character. Quality grade filter (default "research").
#' @param query       Character or NULL. Optional free-text query (e.g. "flower").
#' @param save        Logical. If TRUE, saves CSV + RDS to data/.
#' @return A tibble of iNaturalist observations.
inat_download <- function(taxon_name,
                          maxresults = 10000,
                          quality    = "research",
                          query      = NULL,
                          save       = TRUE) {

    message("Querying iNaturalist for: ", taxon_name,
            if (!is.null(query)) paste0(" (query: ", query, ")") else "")

    args <- list(
        taxon_name = taxon_name,
        maxresults = maxresults,
        quality    = quality
    )
    if (!is.null(query)) args$query <- query

    df <- do.call(get_inat_obs, args) |> as_tibble()

    message("  Retrieved ", nrow(df), " observations")

    if (save && nrow(df) > 0) {
        tag <- gsub(" ", "_", tolower(taxon_name))
        if (!is.null(query)) tag <- paste0(tag, "_", tolower(query))
        dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)

        csv_path <- file.path("data/raw", paste0(tag, ".csv"))
        rds_path <- file.path("data/raw", paste0(tag, ".rds"))
        write_csv(df, csv_path)
        saveRDS(df, rds_path)
        message("  Saved to ", csv_path, " and ", rds_path)
    }

    df
}

#' Subset observations to a geographic bounding box
#'
#' @param df   Tibble with latitude/longitude columns.
#' @param lat  Numeric vector of length 2: c(min_lat, max_lat).
#' @param lon  Numeric vector of length 2: c(min_lon, max_lon).
#' @return Filtered tibble.
subset_area <- function(df, lat, lon) {
    df |>
        filter(between(latitude, lat[1], lat[2]),
               between(longitude, lon[1], lon[2]))
}

### ============================================================
### Quick example (uncomment to run)
### ============================================================
# grona   <- inat_download("Grona", maxresults = 10000)
# grona_tw <- subset_area(grona, c(21, 26), c(120, 123))
#
# lupulina <- inat_download("Medicago lupulina", maxresults = 5000)
# lupulina_flower <- inat_download("Medicago lupulina", query = "flower", maxresults = 1000)
