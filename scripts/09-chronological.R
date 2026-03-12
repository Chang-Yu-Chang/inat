### ============================================================
### Chronological Changes in the Taipei Basin (1900–Present)
### ============================================================
### Composite figure with shared x-axis (year) and two panels:
###   A) Dual y-axis: mean monthly max/min temperature
###      (CWB Station 466920, Daan) + Taipei City
###      population (census/registration data)
###   B) Histogram of Grona specimen collections with individual
###      species dots overlaid
###
### Data sources:
###   Temperature — CWB CODIS Station 466920 (Daan, Taipei),
###     monthly max/min air temperature, 1897–2025.
###   Population — Japanese colonial censuses, ROC household
###     registration, DGBAS census data (data/raw/).
###
### Usage:
###   source("scripts/09-chronological.R")
### ============================================================

library(tidyverse)
library(cowplot)

### ============================================================
### 1. Shared definitions: periods, colors, x-axis limits
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
period_bands$full_label <- paste0(period_bands$label, "\n", c("Japanese Colonial", "Post-War Transition", "Industrial Expansion", "Metropolitan Consolidate", "Anthropocene Resilience"))
period_boundaries <- c(1945, 1965, 1985, 2005)

### ============================================================
### 2. Load CWB station temperature data (Station 466920, Daan)
### ============================================================
# CWB monthly extreme CSVs: each cell is "value / datetime" or "-- / --"
parse_cwb <- function(path) {
    raw <- read_csv(path, col_names = FALSE, skip = 1, show_col_types = FALSE)
    # Drop last row if it's "極端值"
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
message("CWB temperature: ", nrow(temp_annual), " years (",
        min(temp_annual$year), "\u2013", max(temp_annual$year), ")")

### ============================================================
### 3. Load population data (census / registration)
### ============================================================
pop_data <- read_csv("data/raw/taipei_city_population.csv", show_col_types = FALSE)
message("Population data: ", nrow(pop_data), " records (",
        min(pop_data$year), "\u2013", max(pop_data$year), ")")

### ============================================================
### 4. Load specimen data
### ============================================================
triflora     <- read_csv("data/processed/tai_grona_triflora_taiwan.csv",     show_col_types = FALSE)
heterocarpa  <- read_csv("data/processed/tai_grona_heterocarpa_taiwan.csv",  show_col_types = FALSE)
heterophylla <- read_csv("data/processed/tai_grona_heterophylla_taiwan.csv", show_col_types = FALSE)

all_specimens <- bind_rows(
    triflora     |> mutate(taxon = "G. triflora"),
    heterocarpa  |> mutate(taxon = "G. heterocarpa"),
    heterophylla |> mutate(taxon = "G. heterophylla")
) |>
    filter(!is.na(latitude), !is.na(longitude), year > 1899, year < 2100,
           county_en == "Taipei") |>
    mutate(
        year_period = case_when(
            year < 1945 ~ "Pre-1945",
            year < 1965 ~ "1945\u20131965",
            year < 1985 ~ "1965\u20131985",
            year < 2005 ~ "1985\u20132005",
            TRUE         ~ "2005\u2013Present"
        ),
        year_period = factor(year_period, levels = c(
            "Pre-1945", "1945\u20131965", "1965\u20131985",
            "1985\u20132005", "2005\u2013Present"
        ))
    )

### ============================================================
### 5. Panel A: dual y-axis — temperature + population
### ============================================================
# Axis transformation constants
# Left y: temperature (°C),  range 0–40 (monthly extremes)
# Right y: population (millions), range 0–3.2
T_MIN <- 0; T_RANGE <- 40     # temperature axis: 0–40
P_MAX <- 3.2                   # population axis: 0–3.2
pop_to_temp <- function(p) p * (T_RANGE / P_MAX) + T_MIN
temp_to_pop <- function(t) (t - T_MIN) * (P_MAX / T_RANGE)

# Reshape temperature for two-line legend
temp_long <- temp_annual |>
    filter(year >= 1900) |>
    pivot_longer(cols = c(mean_monthly_max, mean_monthly_min),
                 names_to = "metric", values_to = "temp") |>
    mutate(metric = recode(metric,
        mean_monthly_max = "Mean monthly max",
        mean_monthly_min = "Mean monthly min"
    ))

p_dual <- ggplot() +
    # Period background bands (reduced opacity)
    geom_rect(
        data = period_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = period_bands$fill, alpha = 0.08,
        inherit.aes = FALSE
    ) +
    # Period boundary lines
    geom_vline(xintercept = period_boundaries,
               linetype = "dotted", color = "gray50", linewidth = 0.3) +
    # Population line (mapped to temperature scale)
    geom_line(
        data = pop_data,
        aes(x = year, y = pop_to_temp(population / 1e6)),
        color = "black", linewidth = 0.7, alpha = 0.7, linetype = 2
    ) +
    geom_point(
        data = pop_data,
        aes(x = year, y = pop_to_temp(population / 1e6)),
        color = "black", size = 0.8, alpha = 0.7
    ) +
    # Temperature lines
    geom_line(
        data = temp_long,
        aes(x = year, y = temp, color = metric),
        linewidth = 0.6
    ) +
    # Dual y-axes
    scale_y_continuous(
        name = expression("Temperature (" * degree * "C)"),
        sec.axis = sec_axis(~ temp_to_pop(.), name = "Population (millions)")
    ) +
    scale_color_manual(
        values = c("Mean monthly max" = "grey10", "Mean monthly min" = "grey70"),
        name = NULL
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    # Annotate lines directly (no legend)
    annotate("text", x = 1912, y = 34,
             label = "Mean monthly max", color = "grey10",
             size = 3, fontface = "italic", hjust = 0) +
    annotate("text", x = 1912, y = 16,
             label = "Mean monthly min", color = "grey70",
             size = 3, fontface = "italic", hjust = 0) +
    annotate("text", x = 1960, y = pop_to_temp(2),
             label = "Population", color = "black",
             size = 3, fontface = "italic") +
    theme_minimal(base_size = 10) +
    theme(
        axis.title.x       = element_blank(),
        axis.text.x        = element_blank(),
        axis.text = element_text(size = 6),
        axis.ticks.x       = element_blank(),
        legend.position    = "none", 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(0, 10, 5, 10)
    )

### ============================================================
### 6. Panel B: specimen histogram + species dots
### ============================================================
# Bin specimens into 5-year intervals
bin_width <- 5
all_specimens <- all_specimens |>
    mutate(
        bin_start  = floor((year - 1900) / bin_width) * bin_width + 1900,
        bin_center = bin_start + bin_width / 2
    )

# Compute bar heights (total per bin)
bar_data <- all_specimens |>
    count(bin_start, bin_center, year_period, name = "count")

# Compute species counts per bin for dot overlay
# Position dots at actual species count with horizontal offset per species
species_dots <- all_specimens |>
    count(bin_center, taxon, name = "n") |>
    mutate(x_offset = case_when(
        taxon == "G. heterocarpa"  ~ -1.4,
        taxon == "G. heterophylla" ~  0.0,
        taxon == "G. triflora"     ~  1.4
    ))

p_hist <- ggplot() +
    # Period background bands (reduced opacity)
    geom_rect(
        data = period_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = period_bands$fill, alpha = 0.08,
        inherit.aes = FALSE
    ) +
    # Period boundary lines
    geom_vline(xintercept = period_boundaries,
               linetype = "dotted", color = "gray50", linewidth = 0.3) +
    # Histogram bars (total per bin)
    geom_col(
        data = bar_data,
        aes(x = bin_center, y = count, fill = year_period),
        width = bin_width - 0.5, alpha = 0.5, color = "white", linewidth = 0.2
    ) +
    # Species dots: y = actual count per species, x = offset within bin
    geom_point(
        data = species_dots,
        aes(x = bin_center + x_offset, y = n, shape = taxon),
        size = 1, stroke = 0.6, color = "gray20", fill = NA, alpha = 0.8
    ) +
    scale_fill_manual(values = year_colors, guide = "none") +
    scale_shape_manual(
        values = c("G. triflora" = 1, "G. heterocarpa" = 5, "G. heterophylla" = 0),
        name = "Species"
    ) +
    scale_x_continuous(breaks = seq(1905, 2025, 10), expand = c(0, 0), name = "Year") +
    coord_cartesian(xlim = x_lim) +
    scale_y_continuous(name = "No. of specimens") +
    theme_minimal(base_size = 10) +
    theme(
        legend.position    = c(0.85, .75),
        legend.background  = element_rect(fill = "grey99", color = "black", linewidth = .3),
        legend.key.size    = unit(0.35, "cm"),
        legend.title       = element_blank(),
        legend.text        = element_text(face = "italic", size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(0, 10, 5, 10)
    )

### ============================================================
### 7. Strip panel with period labels
### ============================================================
strip_data <- period_bands |>
    mutate(x_mid = (xmin + xmax) / 2)

p_strip <- ggplot() +
    geom_text(
        data = strip_data,
        aes(x = xmin, y = 0, label = full_label, color = full_label),
        angle = 30, hjust = 0, vjust = 2,
        size = 3.2
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = setNames(period_bands$fill, period_bands$full_label)) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    theme_void() +
    theme(
        plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(color = "none")

### ============================================================
### 8. Combine panels
### ============================================================
p_panels <- plot_grid(
    p_strip,
    p_dual,
    p_hist,
    ncol = 1, align = "v", axis = "lr",
    labels = c("", "A", "B"), scale = 0.95,
    rel_heights = c(0.5, 1, 1)
) + 
    theme(plot.background = element_rect(fill = "white", color = NA))

### ============================================================
### 8. Save
### ============================================================
dir.create("plots", showWarnings = FALSE)

ggsave("plots/09-chronological-taipei.png", p_panels,
       width = 6, height = 5, dpi = 800)

message("Saved: plots/09-chronological-taipei.png")


