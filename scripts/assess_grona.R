### ============================================================
### 0. Libraries
### ============================================================
library(tidyverse)
library(janitor)
library(rinat)
library(cowplot)
library(sf)
library(MASS)   # kde2d
library(FNN)    # nearest neighbor
library(ggthemes)


### ============================================================
### 1. Download iNaturalist data (research grade)
### ============================================================
gro <- get_inat_obs(
    taxon_name = "Grona",
    maxresults = 10000,
    quality = "research"
) %>% as_tibble()

# Clean duplicated epithet if needed
gro <- gro %>%
    mutate(scientific_name = str_replace(scientific_name, "heterocarpos heterocarpos","heterocarpos")) %>%
    mutate(scientific_name = str_replace(scientific_name, "heterocarpos strigosa","heterocarpos"))


### ============================================================
### 2. Subset to Taiwan only
### ============================================================
subset_area <- function(x, lat, lon){
    x %>%
        filter(between(latitude, lat[1], lat[2]),
               between(longitude, lon[1], lon[2]))
}

gro_tw <- subset_area(gro, c(21, 26), c(120, 123)) %>%
    mutate(scientific_name = factor(scientific_name, c("Grona triflora", "Grona heterophylla", "Grona heterocarpos")))

### ============================================================
### 3. Taiwan map
### ============================================================
map_tw <- map_data("world") %>%
    as_tibble() %>%
    filter(region == "Taiwan") %>%
    filter(lat > 21, long > 120)

#tab_use <- read_csv("~/Downloads/98136_95-104年國土利用調查成果鄉鎮市區統計資料（95年版土地使用分類系統表，3級分類）.csv")



### ============================================================
### PANEL A — Taiwan wide-map (hexbin density + raw points)
### ============================================================
pA <- ggplot() +
    geom_polygon(data = map_tw,
                 aes(x = long, y = lat, group = group),
                 fill = "gray98", color = "gray50", linewidth = 0.3) +

    # density layer
    stat_bin_hex(data = gro_tw,
                 aes(x = longitude, y = latitude),
                 bins = 40, alpha = 0.7) +
    scale_fill_gradient(low = "grey99", high = "grey40", name = "Density") +
    # raw points
    geom_point(data = gro_tw, aes(longitude, latitude, color = scientific_name), alpha = 0.6, size = 0.3) +
    scale_color_brewer(palette = "Set2", name = "Species") +
    coord_fixed(clip = "off") +
    theme_minimal() +
    theme(
        legend.position = "right",
        legend.background = element_blank()
    ) +
    guides(
        #fill = guide_legend(direction = "horizontal")
        color = guide_legend(override.aes = list(size = 2))
    ) +
    labs(x = "Longitude", y = "Latitude")


### ============================================================
### PANEL B — Species counts barplot
### ============================================================
species_counts <- gro_tw %>%
    count(scientific_name)

pB <- ggplot(species_counts) +
    geom_col(aes(x = reorder(scientific_name, -n), y = n, fill = scientific_name), color = "black") +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    labs(y = "Count") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none"
    )


### ============================================================
### PANEL C — Zoom-in panel (automatic detection of densest region)
### ============================================================
#
# # KDE for density
# dens <- kde2d(gro_tw$longitude, gro_tw$latitude, n = 300)
# peak_idx <- which(dens$z == max(dens$z), arr.ind = TRUE)
# peak_lon <- dens$x[peak_idx[1]]
# peak_lat <- dens$y[peak_idx[2]]
#
# # Define a ~2 km radius bounding box
# zoom_radius <- 0.02  # in degrees (~2 km)
# zoom_box <- gro_tw %>%
#     filter(between(longitude, peak_lon - zoom_radius, peak_lon + zoom_radius),
#            between(latitude,  peak_lat - zoom_radius, peak_lat + zoom_radius))
#
# pC <- ggplot(zoom_box) +
#     geom_polygon(data = map_tw,
#                  aes(long, lat, group = group),
#                  fill = NA, color = "gray60", linewidth = 0.2) +
#     geom_point(aes(longitude, latitude, color = scientific_name),
#                size = 2, alpha = 0.85) +
#     scale_color_brewer(palette = "Set2", name = NULL) +
#     coord_fixed(clip = "off") +
#     theme_bw() +
#     labs(
# #        subtitle = sprintf("Zoomed-in region near (%.3f, %.3f)", peak_lat, peak_lon),
#         x = "Longitude", y = "Latitude"
#     ) +
#     theme(legend.position = "none")
#

### ============================================================
### PANEL D — Interspecific nearest neighbor distances
### ============================================================

# helper: degree to meters
deg_to_m <- function(d) d * 111000

dist_list <- list()
species_vec <- unique(gro_tw$scientific_name)

for (sp in species_vec) {
    focal  <- gro_tw %>% filter(scientific_name == sp)
    others <- gro_tw %>% filter(scientific_name != sp)

    if (nrow(focal) == 0 | nrow(others) == 0) next

    focal_mat  <- as.matrix(focal[, c("longitude", "latitude")])
    others_mat <- as.matrix(others[, c("longitude", "latitude")])

    nn <- get.knnx(others_mat, focal_mat, k = 1)

    focal$nn_dist_m <- deg_to_m(nn$nn.dist)
    dist_list[[sp]] <- focal
}

distances_all <- bind_rows(dist_list)

pD <- distances_all %>%
    ggplot(aes(nn_dist_m)) +
    geom_histogram(bins = 40,
                   fill = "#9ecae1", color = "black", alpha = 0.85) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    labs(
        #title = "D. Interspecific Nearest-Neighbor Distance",
        #subtitle = "Distances < 50–100 m indicate frequent fine-scale co-occurrence",
        x = "Distance to nearest other Grona species (meters)",
        y = "Number of observations"
    )


### ============================================================
### Combine all panels
### ============================================================
p_full <- plot_grid(
    pA,
    plot_grid(pB, pD, ncol = 1, labels = c("B", "C"), align = "v", axis = "lr"),
    ncol = 2, labels = c("A", ""),
    rel_widths = c(1.2, 1), scale = c(0.95, 1),
    align = "hv", axis = "tl"
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave("plots/grona-four-panel.png", p_full, width = 10, height = 6, dpi = 300)

