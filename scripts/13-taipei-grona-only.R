### ============================================================
### Taipei Map + Chronological Figure - Grona Species Only
### ============================================================
### Three Grona species: G. triflora, G. heterocarpa, G. heterophylla
###
### Panel A: Taiwan-wide overview map
### Panel B: Greater Taipei zoom map with elevation hillshade
### Panel C: Temperature + population dual-axis
### Panel D: Specimen histogram over time
###
### Usage:
###   source("scripts/13-taipei-grona-only.R")
### ============================================================

library(tidyverse)
library(sf)
library(terra)
library(elevatr)
library(cowplot)

### ============================================================
### 1. Shared definitions
### ============================================================
x_lim <- c(1900, 2026)

period_bands <- tibble(
    label = c("(I) Pre-1945", "(II) 1945–1965", "(III) 1965–1985",
              "(IV) 1985–2005", "(V) 2005–Present"),
    xmin  = c(1900, 1945, 1965, 1985, 2005),
    xmax  = c(1945, 1965, 1985, 2005, 2026)
)

year_colors <- c(
    "(I) Pre-1945"          = "#2ca02c",
    "(II) 1945–1965"        = "#1f77b4",
    "(III) 1965–1985"       = "#9467bd",
    "(IV) 1985–2005"        = "#d62728",
    "(V) 2005–Present"      = "#8b0000"
)

period_bands$fill <- year_colors[period_bands$label]
period_bands$full_label <- paste0(
    period_bands$label, "\n",
    c("Japanese Colonial", "Post-War Transition",
      "Industrial Expansion", "Metropolitan",
      "Contemporary")
)
period_boundaries <- c(1945, 1965, 1985, 2005)

# 3 species: shapes
species_shapes <- c(
    "G. triflora"      = 1,   # circle
    "G. heterocarpa"   = 5,   # diamond
    "G. heterophylla"  = 0    # square
)

### ============================================================
### 2. Load specimen data (3 Grona species)
### ============================================================
grona_triflora     <- read_csv("data/processed/tai_grona_triflora_taiwan.csv",     show_col_types = FALSE)
grona_heterocarpa  <- read_csv("data/processed/tai_grona_heterocarpa_taiwan.csv",  show_col_types = FALSE)
grona_heterophylla <- read_csv("data/processed/tai_grona_heterophylla_taiwan.csv", show_col_types = FALSE)

all_specimens <- bind_rows(
    grona_triflora     |> mutate(taxon = "G. triflora"),
    grona_heterocarpa  |> mutate(taxon = "G. heterocarpa"),
    grona_heterophylla |> mutate(taxon = "G. heterophylla")
) |>
    filter(!is.na(latitude), !is.na(longitude), !is.na(year),
           year > 1899, year < 2100) |>
    mutate(
        year_period = case_when(
            year < 1945 ~ "(I) Pre-1945",
            year < 1965 ~ "(II) 1945–1965",
            year < 1985 ~ "(III) 1965–1985",
            year < 2005 ~ "(IV) 1985–2005",
            TRUE        ~ "(V) 2005–Present"
        ),
        year_period = factor(year_period, levels = c(
            "(I) Pre-1945", "(II) 1945–1965", "(III) 1965–1985",
            "(IV) 1985–2005", "(V) 2005–Present"
        )),
        taxon = factor(taxon, levels = names(species_shapes))
    )

# Greater Taipei: Taipei + New Taipei
taipei_specimens <- all_specimens |>
    filter(county_en %in% c("Taipei", "New Taipei"))

message("All 3 Grona species — Greater Taipei: ", nrow(taipei_specimens))
message("By species:")
print(table(taipei_specimens$taxon))

### ============================================================
### 3. Taiwan overview map (Panel A)
### ============================================================
gadm_cache <- "data/gadm41_TWN_2.json"
map_tw <- read_sf(gadm_cache)
map_tw_crop <- st_crop(map_tw, xmin = 119.8, xmax = 122.2, ymin = 21.8, ymax = 25.4)

taipei_bbox <- c(xmin = 121.35, xmax = 121.85, ymin = 24.85, ymax = 25.25)

p_taiwan <- ggplot() +
    geom_sf(data = map_tw_crop, fill = "gray97", color = "gray40", linewidth = 0.3) +
    geom_point(
        data = filter(all_specimens, longitude >= 120, latitude < 25.5),
        aes(x = longitude, y = latitude, color = year_period, shape = taxon),
        alpha = 0.7, size = 1.3, stroke = 0.5, fill = NA
    ) +
    geom_rect(
        aes(xmin = taipei_bbox["xmin"], xmax = taipei_bbox["xmax"],
            ymin = taipei_bbox["ymin"], ymax = taipei_bbox["ymax"]),
        color = "black", fill = NA, linewidth = 0.5, linetype = 1
    ) +
    scale_color_manual(values = year_colors, name = "Collection period") +
    scale_shape_manual(values = species_shapes, name = "Species") +
    scale_x_continuous(limits = c(118, 122.5), breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(21.8, 26.5), breaks = scales::pretty_breaks(n = 5)) +
    coord_sf(clip = "off") +
    theme_minimal(base_size = 10) +
    theme(
        legend.position  = "none",
        plot.margin      = margin(5, 5, 5, 5),
        axis.title       = element_blank(),
        axis.ticks       = element_blank()
    )

### ============================================================
### 3b. Greater Taipei zoom map with elevation (Panel B)
### ============================================================

# Elevation raster
message("Downloading elevation data for Taipei area...")
taipei_elev_bbox <- st_as_sf(
    st_as_sfc(st_bbox(c(xmin = 121.35, xmax = 121.85, ymin = 24.85, ymax = 25.25))),
    crs = 4326
)
dem_raw <- get_elev_raster(taipei_elev_bbox, z = 10, clip = "bbox")
dem_terra <- rast(dem_raw)

slope  <- terrain(dem_terra, v = "slope",  unit = "radians")
aspect <- terrain(dem_terra, v = "aspect", unit = "radians")
hillshade <- shade(slope, aspect, angle = 40, direction = 315)

dem_df <- as.data.frame(dem_terra, xy = TRUE) |> rename(elevation = 3) |> filter(!is.na(elevation))
hs_df  <- as.data.frame(hillshade, xy = TRUE) |> rename(hillshade = 3) |> filter(!is.na(hillshade))

# Filter specimens within map bbox
map_specimens <- taipei_specimens |>
    filter(
        longitude >= taipei_bbox["xmin"], longitude <= taipei_bbox["xmax"],
        latitude  >= taipei_bbox["ymin"], latitude  <= taipei_bbox["ymax"]
    )

p_map <- ggplot() +
    geom_raster(data = hs_df, aes(x = x, y = y, fill = hillshade),
                show.legend = FALSE) +
    scale_fill_gradient(low = "gray55", high = "white") +
    geom_sf(data = map_tw, fill = NA, color = "gray40", linewidth = 0.4) +
    geom_point(
        data = map_specimens,
        aes(x = longitude, y = latitude, color = year_period, shape = taxon),
        alpha = 0.7, size = .8, stroke = 0.5, fill = NA
    ) +
    scale_color_manual(values = year_colors, name = "Collection period") +
    scale_shape_manual(values = species_shapes, name = "Species") +
    coord_sf(xlim = c(121.35, 121.85), ylim = c(24.85, 25.25)) +
    scale_x_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n=2)) +
    scale_y_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n=2)) +
    theme_minimal(base_size = 10) +
    theme(
        legend.position  = "right",
        legend.key.size  = unit(0.35, "cm"),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7, face = "italic"),
        plot.margin      = margin(13, 6, 5, 5), 
        plot.background = element_rect(color = "black", fill = "snow", linewidth = .8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .6),
        axis.title = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 5)
    )

### ============================================================
### 4. Load temperature + population
### ============================================================
parse_cwb <- function(path) {
    raw <- read_csv(path, col_names = FALSE, skip = 1, show_col_types = FALSE)
    raw <- raw |> filter(X1 != "極端值")
    years <- suppressWarnings(as.integer(raw[[1]]))
    month_mat <- raw[, 2:13] |>
        mutate(across(everything(), ~ {
            v <- str_extract(.x, "^[^/]+") |> str_trim()
            suppressWarnings(as.numeric(v))
        }))
    n_valid <- rowSums(!is.na(month_mat))
    annual_mean <- ifelse(n_valid >= 6, rowMeans(month_mat, na.rm = TRUE), NA_real_)
    tibble(year = years, annual_mean = annual_mean) |> filter(!is.na(year))
}

tmax_annual <- parse_cwb("data/raw/466920-2025-MaxAirTemperature-month.csv") |>
    rename(mean_monthly_max = annual_mean)
tmin_annual <- parse_cwb("data/raw/466920-2025-MinAirTemperature-month.csv") |>
    rename(mean_monthly_min = annual_mean)
temp_annual <- inner_join(tmax_annual, tmin_annual, by = "year")

pop_data <- read_csv("data/raw/taipei_city_population.csv", show_col_types = FALSE)

T_MIN <- 0; T_RANGE <- 40; P_MAX <- 3.2
pop_to_temp <- function(p) p * (T_RANGE / P_MAX) + T_MIN
temp_to_pop <- function(t) (t - T_MIN) * (P_MAX / T_RANGE)

temp_long <- temp_annual |>
    filter(year >= 1900) |>
    pivot_longer(cols = c(mean_monthly_max, mean_monthly_min),
                 names_to = "metric", values_to = "temp") |>
    mutate(metric = recode(metric,
        mean_monthly_max = "Mean monthly max",
        mean_monthly_min = "Mean monthly min"
    ))

### ============================================================
### 5. Panel C: Temperature + Population
### ============================================================
p_dual <- ggplot() +
    geom_rect(
        data = period_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = period_bands$fill, alpha = 0.08, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = period_boundaries,
               linetype = "dotted", color = "gray50", linewidth = 0.3) +
    geom_line(
        data = pop_data,
        aes(x = year, y = pop_to_temp(population / 1e6)),
        color = "black", linewidth = 0.7, alpha = 0.7, linetype = 1
    ) +
    geom_point(
        data = pop_data,
        aes(x = year, y = pop_to_temp(population / 1e6)),
        color = "black", size = 0.8, alpha = 0.7
    ) +
    geom_line(
        data = temp_long,
        aes(x = year, y = temp, color = metric),
        linewidth = 0.6
    ) +
    scale_y_continuous(
        name = expression("Temperature (" * degree * "C)"),
        sec.axis = sec_axis(~ temp_to_pop(.), name = "Population (millions)")
    ) +
    scale_color_manual(
        values = c("Mean monthly max" = "grey30", "Mean monthly min" = "grey60"),
        name = NULL
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    annotate("text", x = 1912, y = 34,
             label = "Mean monthly max", color = "grey30",
             size = 2.5, fontface = "italic", hjust = 0) +
    annotate("text", x = 1912, y = 16,
             label = "Mean monthly min", color = "grey60",
             size = 2.5, fontface = "italic", hjust = 0) +
    annotate("text", x = 1960, y = pop_to_temp(2),
             label = "Population", color = "black",
             size = 2.5, fontface = "italic") +
    theme_minimal(base_size = 10) +
    theme(
        axis.title.x       = element_blank(),
        axis.text.x        = element_blank(),
        axis.text          = element_text(size = 6),
        axis.ticks.x       = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = .5),
        axis.title.y = element_text(size = 8),
        legend.position    = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(0, 10, 5, 10)
    )

### ============================================================
### 6. Panel D: Specimen histogram (3 Grona species)
### ============================================================
bin_width <- 5
taipei_specimens <- taipei_specimens |>
    mutate(
        bin_start  = floor((year - 1900) / bin_width) * bin_width + 1900,
        bin_center = bin_start + bin_width / 2
    )

bar_data <- taipei_specimens |>
    count(bin_start, bin_center, year_period, name = "count")

species_dots <- taipei_specimens |>
    count(bin_center, taxon, name = "n") |>
    mutate(x_offset = case_when(
        taxon == "G. triflora"      ~ -1.2,
        taxon == "G. heterocarpa"   ~  0.0,
        taxon == "G. heterophylla"  ~  1.2
    ))

p_hist <- ggplot() +
    geom_rect(
        data = period_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = period_bands$fill, alpha = 0.08, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = period_boundaries,
               linetype = "dotted", color = "gray50", linewidth = 0.3) +
    geom_col(
        data = bar_data,
        aes(x = bin_center, y = count, fill = year_period),
        width = bin_width - 0.5, alpha = 0.5, color = "white", linewidth = 0.2
    ) +
    geom_point(
        data = species_dots,
        aes(x = bin_center + x_offset, y = n, shape = taxon),
        size = 1, stroke = 0.6, color = "gray20", fill = NA, alpha = 0.8
    ) +
    scale_fill_manual(values = year_colors, guide = "none") +
    scale_shape_manual(values = species_shapes, name = "Species") +
    scale_x_continuous(breaks = seq(1905, 2025, 20), expand = c(0, 0), name = "Year") +
    coord_cartesian(xlim = x_lim) +
    scale_y_continuous(name = "No. of specimens") +
    theme_minimal(base_size = 10) +
    theme(
        legend.position    = c(0.85, 0.68),
        legend.background  = element_rect(fill = "snow", color = "black", linewidth = 0.3),
        legend.key.size    = unit(0.35, "cm"),
        legend.title       = element_blank(),
        legend.text        = element_text(face = "italic", size = 5),
        legend.margin = margin(1,1,1,1),
        legend.spacing.y = unit(2, "mm"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(0, 10, 5, 10),
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1.2),
        axis.title = element_text(size = 8),
        axis.ticks.x = element_line(color = "black", linewidth = .5)
    )

### ============================================================
### 7. Strip panel with period labels
### ============================================================
strip_data <- period_bands |> mutate(x_mid = (xmin + xmax) / 2)

p_strip <- ggplot() +
    geom_text(
        data = strip_data,
        aes(x = xmin, y = 0, label = full_label, color = full_label),
        angle = 30, hjust = 0, vjust = 2, size = 2.7, lineheight = 0.7
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = setNames(period_bands$fill, period_bands$full_label)) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    guides(color = "none")

### ============================================================
### 8. Combine: maps on left, chronological panels on right
### ============================================================
# Left: Taiwan overview (A) with Taipei zoom inset (B)
# Inset placement (ggdraw 0-1 coords)

# Trapezoid shade: connect rect corners to inset corners
zoom_shade <- data.frame(
    x = c(0.58, 0.74, 0.74, 0.58),
    y = c(0.94, 0.7, 0.63, 0.64)
)

p_left <- ggdraw() +
    draw_plot(p_taiwan) +
    geom_polygon(
        data = zoom_shade, aes(x = x, y = y),
        fill = "grey70", alpha = 0.25, color = NA
    ) +
    draw_plot(
        p_map + theme(legend.position = "none"),
        x = 0.18, y = 0.59, width = 0.4, height = 0.4
    ) +
    draw_plot_label(label = c("A", ""), x = c(0, 0.22), y = c(1, 0.94), size = 12)

# Right column: strip + temp/pop (C) + histogram (D)
p_right <- plot_grid(
    p_strip,
    p_dual,
    p_hist,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(0.4, 1, 1)
) +
    draw_plot_label(label = c("B", "C"), x = c(0, 0), y = c(1, 0.47), size = 12)

# Full figure: maps (A+B) | chronological (C+D)
p_final <- plot_grid(
    p_left,
    p_right,
    ncol = 2,
    rel_widths = c(1, 1)
)

# Extract legends
legend_period <- get_legend(
    p_map +
        guides(shape = "none", color = guide_legend(override.aes = list(size = 5))) +
        theme(legend.position = "bottom",
              legend.title = element_text(size = 8),
              legend.text  = element_text(size = 8))
)

legend_row <- legend_period

p_out <- plot_grid(
    p_final, legend_row,
    ncol = 1, rel_heights = c(1, 0.08)
) +
    theme(plot.background = element_rect(fill = "white", color = NA))

### ============================================================
### 9. Save
### ============================================================
dir.create("plots", showWarnings = FALSE)

ggsave("plots/13-taipei-grona-only.png", p_out,
       width = 8, height = 5, dpi = 800)

message("Saved: plots/13-taipei-grona-only.png")
