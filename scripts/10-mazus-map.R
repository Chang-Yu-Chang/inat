### ============================================================
### Geographic Map of Mazus Specimens in Taiwan
### ============================================================
### Produces a map of Taiwan with specimen points for two species:
###   M. pumilus and M. fauriei
### Points colored by collection period, shaped by species.
###
### Usage:
###   source("scripts/10-mazus-map.R")
### ============================================================

library(tidyverse)
library(sf)
library(cowplot)

### ============================================================
### 1. Load and combine Mazus specimen data
### ============================================================
pumilus <- read_csv("data/processed/tai_mazus_pumilus_taiwan.csv", show_col_types = FALSE)
fauriei <- read_csv("data/processed/tai_mazus_fauriei_taiwan.csv", show_col_types = FALSE)

all_specimens <- bind_rows(
    pumilus |> mutate(taxon = "M. pumilus"),
    fauriei |> mutate(taxon = "M. fauriei")
) |>
    filter(!is.na(latitude), !is.na(longitude), !is.na(year), year > 0) |>
    mutate(
        year_period = case_when(
            year < 1945 ~ "Pre-1945",
            year < 1965 ~ "1945\u20131965",
            year < 1985 ~ "1965\u20131985",
            year < 2005 ~ "1985\u20132005",
            TRUE        ~ "2005\u2013Present"
        ),
        year_period = factor(year_period, levels = c(
            "Pre-1945", "1945\u20131965", "1965\u20131985",
            "1985\u20132005", "2005\u2013Present"
        ))
    )

message("Total Mazus specimens with coords: ", nrow(all_specimens))
message("Specimens per species:")
print(table(all_specimens$taxon))
message("Specimens per period:")
print(table(all_specimens$year_period))

### ============================================================
### 2. Taiwan county map (GADM)
### ============================================================
gadm_cache <- "data/gadm41_TWN_2.json"
if (!file.exists(gadm_cache)) {
    message("Downloading Taiwan county boundaries from GADM...")
    dir.create(dirname(gadm_cache), showWarnings = FALSE, recursive = TRUE)
    download.file(
        "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_TWN_2.json",
        gadm_cache, mode = "wb", quiet = TRUE
    )
}
map_tw <- read_sf(gadm_cache) |>
    st_crop(xmin = 119.8, xmax = 122.2, ymin = 21.8, ymax = 25.4)

### ============================================================
### 3. Color palette and shapes
### ============================================================
year_colors <- c(
    "Pre-1945"          = "#2ca02c",
    "1945\u20131965"    = "#1f77b4",
    "1965\u20131985"    = "#9467bd",
    "1985\u20132005"    = "#d62728",
    "2005\u2013Present" = "#8b0000"
)

### ============================================================
### 4. Map: Whole Taiwan
### ============================================================
p_taiwan <- ggplot() +
    geom_sf(data = map_tw, fill = "gray97", color = "gray40", linewidth = 0.3) +
    geom_point(
        data = all_specimens,
        aes(x = longitude, y = latitude, color = year_period, shape = taxon),
        alpha = 0.7, size = 2, stroke = 0.7, fill = NA
    ) +
    scale_color_manual(values = year_colors, name = "Collection period") +
    scale_shape_manual(
        values = c("M. pumilus" = 1, "M. fauriei" = 5),
        name = "Species"
    ) +
    coord_sf(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 12)
    ) +
    labs(
        title = "Mazus specimen collections in Taiwan",
        subtitle = "TAI herbarium records colored by collection period",
        x = "Longitude", y = "Latitude"
    )

### ============================================================
### 5. Save
### ============================================================
dir.create("plots", showWarnings = FALSE)

ggsave("plots/10-mazus-map-taiwan.png", p_taiwan,
       width = 7, height = 9, dpi = 300)

message("\nSaved: plots/10-mazus-map-taiwan.png")
