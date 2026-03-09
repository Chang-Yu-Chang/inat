### ============================================================
### Clean TAI Specimen Data
### ============================================================
### Reads a raw TAI RDS file, adds a standardised `county`
### column using the current Taiwan administrative divisions,
### and writes to data/processed/.
###
### Usage:
###   source("scripts/04-clean-tai.R")
###   clean_tai("grona_triflora")
###   clean_tai("grona_heterophylla")
###   clean_tai("grona_heterocarpa")
### ============================================================

library(dplyr)
library(readr)
library(stringr)

### ============================================================
### 1. County cross-reference table
### ============================================================
### Maps the TAI `district` field (e.g. "台北 Taipei") to the
### current administrative division names (county/city).

county_lookup <- tribble(
    ~district_pattern,   ~county_zh,    ~county_en,
    "台北",              "臺北市",      "Taipei",
    "新北",              "新北市",      "New Taipei",
    "桃園",              "桃園市",      "Taoyuan",
    "新竹",              "新竹縣",      "Hsinchu",
    "苗栗",              "苗栗縣",      "Miaoli",
    "台中",              "臺中市",      "Taichung",
    "彰化",              "彰化縣",      "Changhua",
    "南投",              "南投縣",      "Nantou",
    "雲林",              "雲林縣",      "Yunlin",
    "嘉義",              "嘉義縣",      "Chiayi",
    "台南",              "臺南市",      "Tainan",
    "高雄",              "高雄市",      "Kaohsiung",
    "屏東",              "屏東縣",      "Pingtung",
    "台東",              "臺東縣",      "Taitung",
    "花蓮",              "花蓮縣",      "Hualien",
    "宜蘭",              "宜蘭縣",      "Yilan",
    "澎湖",              "澎湖縣",      "Penghu",
    "金門",              "金門縣",      "Kinmen",
    "連江",              "連江縣",      "Lienchiang",
    "基隆",              "基隆市",      "Keelung"
)

### ============================================================
### 2. Helper: map district → county
### ============================================================
match_county <- function(district, lookup = county_lookup) {
    matched_zh <- NA_character_
    matched_en <- NA_character_
    for (i in seq_len(nrow(lookup))) {
        if (!is.na(district) &&
            str_detect(district, fixed(lookup$district_pattern[i]))) {
            matched_zh <- lookup$county_zh[i]
            matched_en <- lookup$county_en[i]
            break
        }
    }
    list(county_zh = matched_zh, county_en = matched_en)
}

### ============================================================
### 3. Main cleaning function
### ============================================================
#' Clean a raw TAI RDS and add county columns
#'
#' @param tag Character. File tag, e.g. "grona_triflora".
#'            Expects data/raw/tai_<tag>.rds to exist.
#' @return Invisibly, the Taiwan-only tibble with county columns.
clean_tai <- function(tag) {
    rds_in <- file.path("data/raw", paste0("tai_", tag, ".rds"))
    if (!file.exists(rds_in)) stop("File not found: ", rds_in)

    raw <- readRDS(rds_in)
    message("\n=== Cleaning: ", tag, " ===")
    message("Raw records: ", nrow(raw))

    # Taiwan-only with county
    df <- raw |>
        filter(country == "Taiwan") |>
        mutate(
            .match    = purrr::map(district, match_county),
            county_zh = purrr::map_chr(.match, "county_zh"),
            county_en = purrr::map_chr(.match, "county_en"),
            .match    = NULL,
            is_taipei   = county_en == "Taipei",
            is_pingtung = county_en == "Pingtung"
        ) |>
        arrange(county_en, as.integer(year))

    message("Taiwan specimens: ", nrow(df))
    message("Counties represented: ", n_distinct(df$county_zh, na.rm = TRUE))

    county_counts <- df |> count(county_zh, county_en, sort = TRUE)
    message("\nSpecimens per county:")
    for (i in seq_len(nrow(county_counts))) {
        message("  ", county_counts$county_zh[i], " (",
                county_counts$county_en[i], "): ", county_counts$n[i])
    }

    taipei_n   <- sum(df$is_taipei,   na.rm = TRUE)
    pingtung_n <- sum(df$is_pingtung, na.rm = TRUE)
    message("\nFocus areas:")
    message("  Taipei (臺北): ", taipei_n, " specimens")
    message("  Pingtung (屏東): ", pingtung_n, " specimens")

    # Full dataset (all countries) with county for Taiwan rows
    df_all <- raw |>
        mutate(
            .match    = purrr::map(district, match_county),
            county_zh = purrr::map_chr(.match, "county_zh"),
            county_en = purrr::map_chr(.match, "county_en"),
            .match    = NULL
        )

    # Save
    dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
    write_csv(df,     file.path("data/processed", paste0("tai_", tag, "_taiwan.csv")))
    saveRDS(df,       file.path("data/processed", paste0("tai_", tag, "_taiwan.rds")))
    write_csv(df_all, file.path("data/processed", paste0("tai_", tag, "_all.csv")))
    saveRDS(df_all,   file.path("data/processed", paste0("tai_", tag, "_all.rds")))
    message("Saved to data/processed/tai_", tag, "_taiwan.csv  (Taiwan only)")
    message("Saved to data/processed/tai_", tag, "_all.csv  (all countries)")

    invisible(df)
}

### ============================================================
### 4. Run for all species
### ============================================================
clean_tai("grona_triflora")
clean_tai("grona_heterophylla")
clean_tai("grona_heterocarpa")
