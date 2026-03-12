### ============================================================
### Geographic Map of Grona Specimens by Collection Year
### ============================================================
### Produces two maps with district contours and sample points
### colored by year period:
###   1) Whole Taiwan
###   2) Greater Taipei area (Taipei Basin)
###
### Year groups:
###   Pre-1945, 1945–1965, 1965–1985, 1985–2005, 2005–Present
###
### Usage:
###   source("scripts/08-map-by-year.R")
### ============================================================

library(tidyverse)
library(sf)
library(terra)
library(elevatr)
library(cowplot)

### ============================================================
### 1. Load and combine all TAI specimen data
### ============================================================
triflora     <- read_csv("data/processed/tai_grona_triflora_taiwan.csv",     show_col_types = FALSE)
heterocarpa  <- read_csv("data/processed/tai_grona_heterocarpa_taiwan.csv",  show_col_types = FALSE)
heterophylla <- read_csv("data/processed/tai_grona_heterophylla_taiwan.csv", show_col_types = FALSE)

all_specimens <- bind_rows(
    triflora     |> mutate(taxon = "G. triflora"),
    heterocarpa  |> mutate(taxon = "G. heterocarpa"),
    heterophylla |> mutate(taxon = "G. heterophylla")
) |>
    filter(!is.na(latitude), !is.na(longitude), year > 0) |>
    mutate(
        year_period = case_when(
            year < 1945 ~ "Pre-1945",
            year <= 1965 ~ "1945\u20131965",
            year <= 1985 ~ "1965\u20131985",
            year <= 2005 ~ "1985\u20132005",
            TRUE         ~ "2005\u2013Present"
        ),
        year_period = factor(year_period, levels = c(
            "Pre-1945", "1945\u20131965", "1965\u20131985",
            "1985\u20132005", "2005\u2013Present"
        ))
    )

message("Total specimens with coords and valid year: ", nrow(all_specimens))
message("Specimens per period:")
print(table(all_specimens$year_period))

### ============================================================
### 2. Taiwan county map (district contours from GADM)
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
### 2b. Geographic features: rivers, mountains, basin contour
### ============================================================

# --- Major rivers (simplified paths) ---
river_list <- list(
    "Tamsui R."   = matrix(c(121.432,25.179, 121.462,25.115, 121.505,25.069, 121.508,25.040), ncol = 2, byrow = TRUE),
    "Dahan R."    = matrix(c(121.505,25.069, 121.480,25.012, 121.420,24.953, 121.356,24.900), ncol = 2, byrow = TRUE),
    "Xindian R."  = matrix(c(121.505,25.069, 121.520,25.023, 121.537,24.985, 121.555,24.950), ncol = 2, byrow = TRUE),
    "Keelung R."  = matrix(c(121.505,25.069, 121.560,25.055, 121.620,25.065, 121.690,25.065, 121.740,25.060), ncol = 2, byrow = TRUE),
    "Zhuoshui R." = matrix(c(120.215,23.835, 120.420,23.815, 120.640,23.773, 120.850,23.810, 120.990,23.850), ncol = 2, byrow = TRUE),
    "Kaoping R."  = matrix(c(120.430,22.480, 120.510,22.575, 120.600,22.700, 120.685,22.835, 120.745,22.975), ncol = 2, byrow = TRUE),
    "Lanyang R."  = matrix(c(121.840,24.685, 121.720,24.668, 121.590,24.660, 121.480,24.665), ncol = 2, byrow = TRUE)
)

rivers_sf <- st_sf(
    name = names(river_list),
    geometry = st_sfc(lapply(river_list, st_linestring), crs = 4326)
)

# River name labels (whole Taiwan)
river_labels_tw <- tibble(
    name = c("Tamsui R.", "Zhuoshui R.", "Kaoping R.", "Lanyang R."),
    lon  = c(121.32, 120.30, 120.44, 121.72),
    lat  = c(25.19, 23.85, 22.54, 24.72)
)

# --- Mountain peaks ---
peaks_taiwan <- tibble(
    name = c("Yushan\n(3952 m)", "Xueshan\n(3886 m)"),
    lon  = c(120.957, 121.231),
    lat  = c(23.470, 24.383)
)

peaks_taipei <- tibble(
    name = c("Qixing Mt.\n(1120 m)", "Datun Mt.\n(1092 m)"),
    lon  = c(121.560, 121.522),
    lat  = c(25.171, 25.177)
)

# --- Taipei Basin approximate contour ---
taipei_basin <- st_sf(
    name = "Taipei Basin",
    geometry = st_sfc(st_polygon(list(matrix(c(
        121.440, 25.015,
        121.445, 25.070,
        121.460, 25.110,
        121.490, 25.140,
        121.530, 25.140,
        121.560, 25.115,
        121.580, 25.075,
        121.582, 25.040,
        121.575, 25.010,
        121.555, 24.990,
        121.520, 24.978,
        121.485, 24.980,
        121.455, 24.995,
        121.440, 25.015
    ), ncol = 2, byrow = TRUE))), crs = 4326)
)

# Crop rivers to Taipei area
taipei_bbox <- c(xmin = 121.35, xmax = 121.85, ymin = 24.85, ymax = 25.25)
rivers_taipei <- st_crop(rivers_sf, taipei_bbox)

# River name labels (Taipei area)
river_labels_tp <- tibble(
    name = c("Tamsui R.", "Dahan R.", "Xindian R.", "Keelung R."),
    lon  = c(121.42, 121.38, 121.56, 121.67),
    lat  = c(25.19, 24.92, 24.94, 25.08)
)

### ============================================================
### 3. Color palette for year periods (blue = early, red = present)
### ============================================================
year_colors <- c(
    "Pre-1945"                   = "#2ca02c",
    "1945\u20131965"             = "#1f77b4",
    "1965\u20131985"             = "#9467bd",
    "1985\u20132005"             = "#d62728",
    "2005\u2013Present"          = "#8b0000"
)

### ============================================================
### 4. Map 1: Whole Taiwan
### ============================================================
p_taiwan <- ggplot() +
    geom_sf(data = map_tw, fill = "gray97", color = "gray40", linewidth = 0.3) +
    geom_sf(data = rivers_sf, color = "steelblue", linewidth = 0.45, alpha = 0.6) +
    geom_text(
        data = river_labels_tw, aes(x = lon, y = lat, label = name),
        size = 1.8, color = "steelblue4", fontface = "italic",
        inherit.aes = FALSE
    ) +
    geom_point(
        data = peaks_taiwan, aes(x = lon, y = lat),
        shape = 17, size = 2, color = "gray30", inherit.aes = FALSE
    ) +
    geom_text(
        data = peaks_taiwan, aes(x = lon, y = lat, label = name),
        size = 2, color = "gray30", vjust = -0.7, fontface = "italic",
        lineheight = 0.8, inherit.aes = FALSE
    ) +
    geom_point(
        data = all_specimens,
        aes(x = longitude, y = latitude, color = year_period, shape = taxon),
        alpha = 0.75, size = 1.8, stroke = 0.7, fill = NA
    ) +
    # Inset rectangle for Greater Taipei
    geom_rect(
        aes(xmin = 121.35, xmax = 121.85, ymin = 24.85, ymax = 25.25),
        color = "black", fill = NA, linewidth = 0.5, linetype = "dashed"
    ) +
    scale_color_manual(values = year_colors, name = "Collection period") +
    scale_shape_manual(
        values = c("G. triflora" = 1, "G. heterocarpa" = 5, "G. heterophylla" = 0),
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
        title = "Grona specimen collections in Taiwan",
        subtitle = "TAI herbarium records colored by collection period",
        x = "Longitude", y = "Latitude"
    )

### ============================================================
### 5. Map 2: Greater Taipei area (Taipei Basin)
### ============================================================
# Crop map to Greater Taipei
taipei_bbox <- c(xmin = 121.35, xmax = 121.85, ymin = 24.85, ymax = 25.25)
map_taipei <- st_crop(map_tw, taipei_bbox)

taipei_specimens <- all_specimens |>
    filter(
        longitude >= taipei_bbox["xmin"], longitude <= taipei_bbox["xmax"],
        latitude  >= taipei_bbox["ymin"], latitude  <= taipei_bbox["ymax"]
    )

message("Specimens in Greater Taipei area: ", nrow(taipei_specimens))

# --- Elevation raster for Taipei area ---
message("Downloading elevation data for Taipei area...")
taipei_elev_bbox <- st_as_sf(
    st_as_sfc(st_bbox(c(xmin = 121.35, xmax = 121.85, ymin = 24.85, ymax = 25.25))),
    crs = 4326
)
dem_raw <- get_elev_raster(taipei_elev_bbox, z = 10, clip = "bbox")
dem_terra <- rast(dem_raw)

# Compute hillshade for terrain relief
slope  <- terrain(dem_terra, v = "slope",  unit = "radians")
aspect <- terrain(dem_terra, v = "aspect", unit = "radians")
hillshade <- shade(slope, aspect, angle = 40, direction = 315)

# Convert to data frame for ggplot
dem_df <- as.data.frame(dem_terra, xy = TRUE) |>
    rename(elevation = 3) |>
    filter(!is.na(elevation))
hs_df <- as.data.frame(hillshade, xy = TRUE) |>
    rename(hillshade = 3) |>
    filter(!is.na(hillshade))

p_taipei <- ggplot() +
    # Hillshade terrain relief
    geom_raster(data = hs_df, aes(x = x, y = y, fill = hillshade),
                show.legend = FALSE) +
    scale_fill_gradient(low = "gray55", high = "white") +
    # Elevation contour lines
    geom_contour(data = dem_df, aes(x = x, y = y, z = elevation),
                 color = "gray50", linewidth = 0.15, alpha = 0.5,
                 breaks = seq(0, 1200, by = 100)) +
    # District borders
    geom_sf(data = map_taipei, fill = NA, color = "gray40", linewidth = 0.4) +
    # Rivers
    geom_sf(data = rivers_taipei, color = "steelblue", linewidth = 0.55, alpha = 0.7) +
    geom_text(
        data = river_labels_tp, aes(x = lon, y = lat, label = name),
        size = 2, color = "steelblue4", fontface = "italic",
        inherit.aes = FALSE
    ) +
    # Mountain peaks
    geom_point(
        data = peaks_taipei, aes(x = lon, y = lat),
        shape = 17, size = 2.5, color = "gray20", inherit.aes = FALSE
    ) +
    geom_text(
        data = peaks_taipei, aes(x = lon, y = lat, label = name),
        size = 2.2, color = "gray20", vjust = -0.7, fontface = "bold.italic",
        lineheight = 0.8, inherit.aes = FALSE
    ) +
    # Specimen points (open shapes)
    geom_point(
        data = taipei_specimens,
        aes(x = longitude, y = latitude, color = year_period, shape = taxon),
        alpha = 0.3, size = 1.3, stroke = 0.5, fill = NA
    ) +
    scale_color_manual(values = year_colors, name = "Collection period") +
    scale_shape_manual(
        values = c("G. triflora" = 1, "G. heterocarpa" = 5, "G. heterophylla" = 0),
        name = "Species"
    ) +
    coord_sf(xlim = c(121.35, 121.85), ylim = c(24.85, 25.25), clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 12)
    ) +
    labs(
        title = "Greater Taipei area (Taipei Basin)",
        subtitle = "TAI herbarium records with elevation relief",
        x = "Longitude", y = "Latitude"
    )

### ============================================================
### 6. Combined panel figure
### ============================================================
p_combined <- plot_grid(
    p_taiwan  + theme(legend.position = "none"),
    p_taipei  + theme(legend.position = "none"),
    ncol = 2, labels = c("A", "B"), rel_widths = c(1, 1.1)
)

# Extract shared legend
legend <- get_legend(
    p_taiwan + theme(legend.position = "bottom", legend.box = "horizontal")
)

p_final <- plot_grid(
    p_combined, legend,
    ncol = 1, rel_heights = c(1, 0.15)
) +
    theme(plot.background = element_rect(fill = "white", color = NA))

### ============================================================
### 7. Save
### ============================================================
dir.create("plots", showWarnings = FALSE)

ggsave("plots/08-grona-map-by-year-taiwan.png",   p_taiwan,  width = 7, height = 9, dpi = 300)
ggsave("plots/08-grona-map-by-year-taipei.png",   p_taipei,  width = 7, height = 6, dpi = 300)
ggsave("plots/08-grona-map-by-year-combined.png",  p_final,  width = 12, height = 8, dpi = 300)

message("\nSaved plots to plots/08-grona-map-by-year-*.png")
