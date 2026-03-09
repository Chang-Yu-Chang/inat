### ============================================================
### Clean Smithsonian Specimen Data
### ============================================================
### Reads the raw Smithsonian CSV for Grona triflora, selects
### and cleans relevant columns, and writes a tidy table to
### data/processed/.
###
### Usage:
###   source("scripts/03-clean-smithsonian.R")
### ============================================================

library(dplyr)
library(readr)
library(stringr)
library(janitor)
library(lubridate)

### ============================================================
### 1. Read raw data
### ============================================================
raw <- read_csv("data/raw/smithsonian.csv", show_col_types = FALSE) |>
    clean_names()

### ============================================================
### 2. Select & rename useful columns
### ============================================================
df <- raw |>
    transmute(
        barcode,
        catalog_number,
        family,
        subfamily,
        tribe,
        # Parse the compound "Taxonomic Name" field:
        #   "Name : Identifier : Institution : Date"
        taxonomic_name_raw = taxonomic_name_filed_as_identified_by_identification_date,
        scientific_name    = str_extract(
            taxonomic_name_filed_as_identified_by_identification_date,
            "^[^:]+") |> str_trim(),
        identified_by      = str_match(
            taxonomic_name_filed_as_identified_by_identification_date,
            ":\\s*([^:]+?)\\s*,")[, 2],
        synonyms = other_taxonomic_names_identification_identified_by_date_identified,
        collector   = collector_s,
        coll_number = collection_number,
        date_collected,
        date_parsed = parse_date_time(date_collected,
                                      orders = c("d b Y", "b Y", "Y"),
                                      quiet = TRUE),
        year  = year(date_parsed),
        month = month(date_parsed),
        biogeographic_region,
        country,
        state   = province_state,
        county  = district_county,
        locality = precise_locality,
        island   = island_name,
        latitude  = centroid_latitude,
        longitude = centroid_longitude,
        elevation_m = elevation_m,
        habit,
        microhabitat = microhabitat_description,
        type_citations = type_citations_scientific_name_type_status_verification_degree_citation,
        notes,
        cultivated,
        ezid = ezid
    )

### ============================================================
### 3. Derive extra fields
### ============================================================
df <- df |>
    mutate(
        # Flag whether coordinates are available
        has_coords = !is.na(latitude) & !is.na(longitude),
        # Clean country names
        country = str_trim(country),
        # Extract type status from type_citations if present
        is_type = !is.na(type_citations) & type_citations != ""
    )

### ============================================================
### 4. Summary
### ============================================================
message("Cleaned Smithsonian data: ", nrow(df), " specimens")
message("  With coordinates: ", sum(df$has_coords))
message("  Countries: ", n_distinct(df$country, na.rm = TRUE))
message("  Year range: ",
        min(df$year, na.rm = TRUE), " - ",
        max(df$year, na.rm = TRUE))

### ============================================================
### 5. Save processed data
### ============================================================
dir.create("data/processed", showWarnings = FALSE)

write_csv(df, "data/processed/smithsonian_grona_triflora.csv")
saveRDS(df, "data/processed/smithsonian_grona_triflora.rds")

message("Saved to data/processed/smithsonian_grona_triflora.csv")
