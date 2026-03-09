### ============================================================
### TAI Herbarium Specimen Search
### ============================================================
### Searches the Plants of Taiwan database (tai2.ntu.edu.tw)
### for digitised herbarium specimens and returns a tidy table.
###
### The TAI system aggregates ~187k specimens from TAI, NTUF,
### KYO, TI, PH, and Kanagawa herbaria.
###
### Usage:
###   source("scripts/02-tai-specimen.R")
###   df <- tai_search("Grona triflora")
###   df <- tai_search("Desmodium triflorum")
###   df <- tai_search("Grona triflora", herbarium = "TAI")
### ============================================================

library(httr2)
library(jsonlite)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

TAI_BASE <- "https://tai2.ntu.edu.tw"

#' Search TAI herbarium specimens
#'
#' @param species_name  Character. Scientific name to search (e.g. "Grona triflora").
#' @param collector     Character. Collector name filter (default "").
#' @param location      Character. Collection locality filter (default "").
#' @param year_from     Character. Start year (default "").
#' @param year_to       Character. End year (default "").
#' @param coll_no       Character. Collection number (default "").
#' @param type_only     Logical. If TRUE, restrict to type specimens (default FALSE).
#' @param herbarium     Character. Herbarium code: "", "TAI", "NTUF", "Kanagawa",
#'                      "KYO", "TI", "Upenn" (default "" = all).
#' @param tai_id        Character. TAI accession number (default "").
#' @param save          Logical. If TRUE, saves CSV + RDS to data/.
#' @return A tibble of specimen records.
tai_search <- function(species_name,
                       collector  = "",
                       location   = "",
                       year_from  = "",
                       year_to    = "",
                       coll_no    = "",
                       type_only  = FALSE,
                       herbarium  = "",
                       tai_id     = "",
                       save       = TRUE) {

    message("Searching TAI for: ", species_name)

    # Build the query string the same way the site's JS does:
    #   /search_specimen/<form serialized params with & encoded as %26 in path>
    type2_val <- if (type_only) "%type%" else ""

    params <- paste0(
        "species=",   utils::URLencode(species_name, reserved = TRUE),
        "&mb=",       utils::URLencode(collector, reserved = TRUE),
        "&loc=",      utils::URLencode(location, reserved = TRUE),
        "&year1=",    year_from,
        "&year2=",    year_to,
        "&collno=",   utils::URLencode(coll_no, reserved = TRUE),
        "&type2=",    utils::URLencode(type2_val, reserved = TRUE),
        "&herbarium=", utils::URLencode(herbarium, reserved = TRUE),
        "&taiid=",    utils::URLencode(tai_id, reserved = TRUE)
    )

    # The TAI server expects the query serialized into the URL path with
    # ampersands encoded as %26 (the JS $.ajax sends it that way).
    path_params <- gsub("&", "%26", params)

    # Step 1: get a session page to obtain a valid CSRF token + cookies
    encoded_name <- utils::URLencode(species_name, reserved = TRUE)
    session_resp <- request(paste0(TAI_BASE, "/search/2/", encoded_name)) |>
        req_perform()

    cookies <- resp_headers(session_resp)
    cookie_values <- cookies[names(cookies) == "set-cookie"]
    page_body <- resp_body_string(session_resp)
    csrf <- str_match(page_body, 'csrf-token"\\s+content="([^"]+)"')[, 2]

    if (is.na(csrf)) {
        stop("Could not obtain CSRF token from TAI. The site may be down.")
    }

    # Extract the Laravel session cookie
    cookie_str <- paste(
        str_extract(unlist(cookie_values), "^[^;]+"),
        collapse = "; "
    )

    # Step 2: call the specimen search API
    api_url <- paste0(TAI_BASE, "/search_specimen/", path_params)

    api_resp <- request(api_url) |>
        req_headers(
            `X-CSRF-TOKEN`     = csrf,
            `X-Requested-With` = "XMLHttpRequest",
            `Accept`           = "application/json",
            `Cookie`           = cookie_str
        ) |>
        req_perform()

    body <- resp_body_json(api_resp, simplifyVector = FALSE)

    records <- body$result
    if (length(records) == 0) {
        message("  No specimens found.")
        return(tibble())
    }

    message("  Found ", length(records), " specimens")

    # Flatten into a tidy tibble
    df <- map_dfr(records, flatten_specimen)

    if (save && nrow(df) > 0) {
        tag <- gsub(" ", "_", tolower(species_name))
        dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
        csv_path <- file.path("data/raw", paste0("tai_", tag, ".csv"))
        rds_path <- file.path("data/raw", paste0("tai_", tag, ".rds"))
        write_csv(df, csv_path)
        saveRDS(df, rds_path)
        message("  Saved to ", csv_path, " and ", rds_path)
    }

    df
}


### ============================================================
### Internal helper
### ============================================================

#' Flatten one specimen JSON record into a single-row tibble
flatten_specimen <- function(rec) {
    loc <- rec$locinfo %||% list()
    tibble(
        tai_id        = rec$taiid        %||% NA_character_,
        species       = rec$species      %||% NA_character_,
        spcm_species  = rec$spcmspecies  %||% NA_character_,
        year          = rec$year         %||% NA_character_,
        collector     = rec$collinfo     %||% NA_character_,
        locality      = loc$loc          %||% NA_character_,
        locality_en   = loc$locE         %||% NA_character_,
        district      = loc$district     %||% NA_character_,
        country       = loc$country      %||% NA_character_,
        longitude     = loc$X            %||% NA_real_,
        latitude      = loc$Y            %||% NA_real_,
        name          = rec$name         %||% NA_character_,
        type          = rec$type         %||% NA_character_,
        family        = rec$apgfamily    %||% NA_character_,
        family_zh     = rec$chapgfamily  %||% NA_character_,
        genus         = rec$genus        %||% NA_character_,
        genus_zh      = rec$chgenus      %||% NA_character_,
        name_zh       = rec$chname       %||% NA_character_,
        code          = rec$code         %||% NA_character_,
        result_type   = rec$resultType   %||% NA_character_
    )
}


### ============================================================
### Quick example (uncomment to run)
### ============================================================
# grona_triflora  <- tai_search("Grona triflora")
# grona_tw        <- grona_triflora |> filter(country == "Taiwan")
# desmodium_trif  <- tai_search("Desmodium triflorum")
