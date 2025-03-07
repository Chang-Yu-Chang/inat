#' Assess iNaturalist data

library(tidyverse)
library(janitor)
library(rinat) # for assessing iNaturalist data

# Query iNaturalist
## Lup
lup <- get_inat_obs(
    taxon_name  = "Medicago lupulina",
    maxresults = 5000,
    quality = "research"
) %>% as_tibble

## Lup flower
lup_f <- get_inat_obs(
    taxon_name  = "Medicago lupulina",
    query = "flower",
    maxresults = 1000,
    quality = "research"
) %>% as_tibble


# Filter for North america based observations
subset_na <- function (inat_data) {
    inat_data %>%
        filter(latitude > 14 & latitude < 71) %>%
        filter(longitude > -125 & longitude < -65) %>%
        # Extract the season
        mutate(month = month(datetime) %>% ordered(1:12))
}

lup_na <- subset_na(lup)
lup_f_na <- subset_na(lup_f)


# Map for north america
map_na <- map_data("world") %>%
    as_tibble() %>%
    filter(region %in% c("USA", "Canada", "Mexico")) %>%
    filter(!subregion %in% c("Alaska", "Hawaii"))

# Color palette
seasonal_palette <- c(
    "#0000FF",  # January (cold - blue)
    "#0099CC",  # February (cold - same as January)
    "#3399FF",  # March
    "#99CCFF",  # April
    "#FFDD99",  # May (neutral transition)
    "#FFAA33",  # June
    "#FF3333",  # July
    "#FFAA33",  # August (warm - red)
    "#FFDD99",  # September (still warm)
    "#99CCFF",  # October (cooling down)
    "#3399FF",  # November (cool tone)
    "#0099CC"   # December (cold - same as January and February)
)

names(seasonal_palette) <- 1:12

# Lupulina ob
table(lup_na$month)
p <- lup_na %>%
    ggplot() +
    geom_polygon(data = map_na, aes(x = long, y = lat, group = group), fill = "gray99", color = "gray40", linewidth = 0.1) +
    geom_point(aes(x = longitude, y = latitude, color = month), size = 1, stroke = .5, shape = 21, alpha = 0.8) +
    coord_fixed(clip = "off") +
    scale_color_manual(values = seasonal_palette) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs()
ggsave(here::here("plot/01-lupulina_na.png"), p, width = 8, height = 5)

# Lupulina flower
table(lup_f_na$month)
p <- lup_f_na %>%
    ggplot() +
    geom_polygon(data = map_na, aes(x = long, y = lat, group = group), fill = "gray99", color = "gray40", linewidth = 0.1) +
    geom_point(aes(x = longitude, y = latitude, color = month), size = 1, stroke = .5, shape = 21, alpha = 0.8) +
    coord_fixed(clip = "off") +
    scale_color_manual(values = seasonal_palette) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs()
ggsave(here::here("plot/02-lupulina_flower_na.png"), p, width = 8, height = 5)


#


if (F) {
    p <- lup_f %>%
        mutate(month = month(datetime) %>% ordered(1:12)) %>%
        ggplot() +
        geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), fill = "gray99", color = "gray40", linewidth = 0.1) +
        geom_point(aes(x = longitude, y = latitude, color = month), size = 1, stroke = .5, shape = 21, alpha = 0.8) +
        coord_fixed(clip = "off") +
        scale_color_manual(values = seasonal_palette) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs(title = "M. lupulina flower")
    ggsave(here::here("plot/02-lupulina_flower_world.png"), p, width = 10, height = 8)



}
