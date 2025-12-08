### ============================================================
### 0. Libraries
### ============================================================

library(tidyverse)
library(twmap) # for taiwan map
library(janitor)
library(rinat)
library(cowplot)
library(sf)
library(MASS)   # kde2d
library(FNN)    # nearest neighbor
library(ggthemes)
library(geosphere) # for computing distance


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
map_tw <- st_as_sf(twmap::tw_county) %>%
    st_crop(
        xmin = 119.8, xmax = 123,
        ymin = 21,  ymax = 25.4
    )

# map_tw <- map_data("world") %>%
#     as_tibble() %>%
#     filter(region == "Taiwan") %>%
#     filter(lat > 21, long > 120)


ggplot(map_tw) +
    geom_sf()
### ============================================================
### PANEL A — Taiwan wide-map (hexbin density + raw points)
### ============================================================
pA <- ggplot() +
    geom_sf(data = map_tw, fill = "gray98", color = "gray50", linewidth = 0.3) +
    # density layer
    stat_bin_2d(data = gro_tw, aes(x = longitude, y = latitude), bins = 40, alpha = 0.7) +
    scale_fill_gradient(low = "grey99", high = "grey40", name = "Density") +
    # raw points
    geom_point(data = gro_tw, aes(longitude, latitude, color = scientific_name), alpha = 0.6, size = 0.3) +
    # zoom in
    geom_rect(aes(xmin = 121.4, xmax = 121.6, ymin = 24.9, ymax = 25.1), color = "black", fill = "NA", linewidth = .5) +
    scale_color_brewer(palette = "Set2", name = "Species") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    coord_sf(clip = "off") +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        legend.background = element_blank()
    ) +
    guides(
        color = "none"
    ) +
    labs(x = "Longitude", y = "Latitude")

### ============================================================
### PANEL A inset — Species counts barplot
### ============================================================
species_counts <- gro_tw %>%
    count(scientific_name) %>%
    mutate(sp_name = str_replace(scientific_name, "Grona", "G."))

pA2 <- ggplot(species_counts) +
    geom_col(aes(x = reorder(sp_name, -n), y = n, fill = scientific_name), color = "black") +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(position = "right") +
    theme_bw() +
    labs(x = "", y = "Count") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5, face = "bold.italic"),
        axis.text.y = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_rect(color = "black", fill = "white", linewidth = .5)
    )


### ============================================================
### Panel B - zoom in map
### ============================================================
map_tw_sub <- st_crop(
    map_tw,
    xmin = 121.4, xmax = 121.6,
    ymin = 24.9,  ymax = 25.1
)

gro_tw_sub <- filter(
    gro_tw,
    longitude > 121.4, longitude < 121.6,
    latitude > 24.9, latitude < 25.1
)

pB <- ggplot() +
    geom_sf(data = map_tw_sub, fill = "gray98", color = "gray50", linewidth = 0.3) +
    # raw points
    geom_point(data = gro_tw_sub, aes(longitude, latitude, color = scientific_name), alpha = 0.6, size = 1) +
    scale_color_brewer(palette = "Set2", name = "Species") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0,0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0,0)) +
    coord_sf(clip = "off") +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(angle = 45, size = 5, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
    ) +
    guides(
        color = "none"
    ) +
    labs(x = "Longitude", y = "Latitude")

### ============================================================
### Panel C - between species distance
### ============================================================
library(tidyverse)
library(FNN)

# Helper to convert degrees to meters
deg_to_m <- function(d) d * 111000

# Function: nearest distance from species sp1 to sp2
nearest_between_species <- function(df, sp1, sp2){
    focal  <- df %>% filter(scientific_name == sp1)
    others <- df %>% filter(scientific_name == sp2)

    if(nrow(focal) == 0 | nrow(others) == 0){
        return(NULL)
    }

    nn <- get.knnx(
        data = as.matrix(others[, c("longitude", "latitude")]),
        query = as.matrix(focal[, c("longitude", "latitude")]),
        k = 1
    )

    tibble(
        sp_pair = paste(sp1, "→", sp2),
        sp1 = sp1,
        sp2 = sp2,
        nn_dist_m = deg_to_m(nn$nn.dist[,1])
    )
}

# Species present
species <- gro_tw$scientific_name %>% unique() %>% sort()

# All heterospecific pairs
pairs <- expand_grid(sp1 = species, sp2 = species) %>%
    filter(sp1 != sp2)

# Compute distances for all 6 directional pairs
dist_all <- map2_dfr(pairs$sp1, pairs$sp2, ~nearest_between_species(gro_tw, .x, .y))

pC <- dist_all %>%
    filter(nn_dist_m < 200) %>%
    #filter(sp1 %in% c("Grona triflora") | sp1 %in% c("Grona heterophylla")) %>%
    ggplot(aes(nn_dist_m, fill = sp_pair)) +
    geom_histogram(alpha = 0.3, bins = 40, position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
    facet_grid(sp1~.) +
    coord_cartesian(clip = "off") +
    labs(
        x = "Nearest heterospecific distance (meters)",
        y = "Count",
        fill = "Species Pair"
    ) +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_blank()
    )

### ============================================================
### Combine all panels
### ============================================================
pB2 <- ggdraw() + draw_image(here::here("plots/grona.jpg"), scale = 1)

p <- plot_grid(
    pA,
    plot_grid(
        plot_grid(pB, pB2, nrow = 1, rel_widths = c(1, 1.2), labels = c("B", "C"), label_x = c(-0.1, 0)),
        pC,
        rel_heights = c(1,2),
        ncol = 1, scale = .95, labels = c("", "D")),
    ncol = 2, labels = c("A", ""),
    rel_widths = c(1, 1.2), scale = c(1, 1)
) +
    theme(plot.background = element_rect(color = NA, fill = "white")) +
    draw_plot(pA2, .03, .75, .15, .23)

ggsave("plots/grona.png", p, width = 8, height = 6, dpi = 300)

