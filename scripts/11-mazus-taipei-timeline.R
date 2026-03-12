### ============================================================
### Mazus Specimen Counts Over Time — Greater Taipei
### ============================================================
### Histogram of Mazus pumilus and M. fauriei specimen counts
### in greater Taipei (Taipei + New Taipei) over 5-year bins.
###
### Usage:
###   source("scripts/11-mazus-taipei-timeline.R")
### ============================================================

library(tidyverse)
library(cowplot)

### ============================================================
### 1. Load and filter Mazus specimen data
### ============================================================
pumilus <- read_csv("data/processed/tai_mazus_pumilus_taiwan.csv", show_col_types = FALSE)
fauriei <- read_csv("data/processed/tai_mazus_fauriei_taiwan.csv", show_col_types = FALSE)

all_specimens <- bind_rows(
    pumilus |> mutate(taxon = "M. pumilus"),
    fauriei |> mutate(taxon = "M. fauriei")
) |>
    filter(!is.na(year), year > 1899, year < 2100,
           county_en %in% c("Taipei", "New Taipei")) |>
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

message("Greater Taipei Mazus specimens: ", nrow(all_specimens))
message("By species: ")
print(table(all_specimens$taxon))

### ============================================================
### 2. Shared definitions
### ============================================================
x_lim <- c(1900, 2026)

period_bands <- tibble(
    label = c("Pre-1945", "1945\u20131965", "1965\u20131985",
              "1985\u20132005", "2005\u2013Present"),
    xmin  = c(1900, 1945, 1965, 1985, 2005),
    xmax  = c(1945, 1965, 1985, 2005, 2026)
)

year_colors <- c(
    "Pre-1945"          = "#2ca02c",
    "1945\u20131965"    = "#1f77b4",
    "1965\u20131985"    = "#9467bd",
    "1985\u20132005"    = "#d62728",
    "2005\u2013Present" = "#8b0000"
)

period_bands$fill <- year_colors[period_bands$label]
period_boundaries <- c(1945, 1965, 1985, 2005)

### ============================================================
### 3. Bin specimens into 5-year intervals
### ============================================================
bin_width <- 5
all_specimens <- all_specimens |>
    mutate(
        bin_start  = floor((year - 1900) / bin_width) * bin_width + 1900,
        bin_center = bin_start + bin_width / 2
    )

# Bar heights (total per bin)
bar_data <- all_specimens |>
    count(bin_start, bin_center, year_period, name = "count")

# Species dots per bin
species_dots <- all_specimens |>
    count(bin_center, taxon, name = "n") |>
    mutate(x_offset = case_when(
        taxon == "M. pumilus" ~ -1.0,
        taxon == "M. fauriei" ~  1.0
    ))

### ============================================================
### 4. Plot
### ============================================================
p_hist <- ggplot() +
    # Period background bands
    geom_rect(
        data = period_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = period_bands$fill, alpha = 0.08,
        inherit.aes = FALSE
    ) +
    # Period boundary lines
    geom_vline(xintercept = period_boundaries,
               linetype = "dotted", color = "gray50", linewidth = 0.3) +
    # Histogram bars
    geom_col(
        data = bar_data,
        aes(x = bin_center, y = count, fill = year_period),
        width = bin_width - 0.5, alpha = 0.5, color = "white", linewidth = 0.2
    ) +
    # Species dots
    geom_point(
        data = species_dots,
        aes(x = bin_center + x_offset, y = n, shape = taxon),
        size = 2, stroke = 0.6, color = "gray20", fill = NA, alpha = 0.8
    ) +
    scale_fill_manual(values = year_colors, guide = "none") +
    scale_shape_manual(
        values = c("M. pumilus" = 1, "M. fauriei" = 5),
        name = "Species"
    ) +
    scale_x_continuous(breaks = seq(1905, 2025, 10), expand = c(0, 0), name = "Year") +
    coord_cartesian(xlim = x_lim) +
    scale_y_continuous(name = "No. of specimens") +
    theme_minimal(base_size = 10) +
    theme(
        legend.position    = c(0.90, 0.85),
        legend.background  = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key.size    = unit(0.4, "cm"),
        legend.title       = element_text(size = 9),
        legend.text        = element_text(face = "italic", size = 8),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", size = 12),
        plot.margin        = margin(10, 10, 5, 10)
    ) +
    labs(title = "Mazus specimen collections in Greater Taipei",
         subtitle = "TAI herbarium \u2014 Taipei + New Taipei")

### ============================================================
### 5. Save
### ============================================================
dir.create("plots", showWarnings = FALSE)

p_final <- p_hist +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("plots/11-mazus-taipei-timeline.png", p_final,
       width = 8, height = 4, dpi = 300)

message("Saved: plots/11-mazus-taipei-timeline.png")
