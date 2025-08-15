#'

library(tidyverse)
library(janitor)
library(rinat) # for assessing iNaturalist data

list_sp <- c("littoralis", "lupulina", "orbicularis", "sativa", "laciniata", "truncatula", "polymorpha", "italica", "rigiduloides")

# Query iNaturalist
med <- get_inat_obs(
    taxon_name  = "Medicago",
    maxresults = 5000,
    quality = "research"
) %>% as_tibble

meds <- med %>%
    filter(str_detect(scientific_name, paste(list_sp, collapse = "|")))



table(meds$scientific_name)

meds %>%
    ggplot(aes(x = scientific_name, y = abs(latitude))) +
    #geom_violin(draw_quantiles = c(.05, .5, .95)) +
    geom_jitter(size = 1, shape = 21, width = .2, alpha = .4) +
    #geom_text(data = mwc, aes(x = species_name, label = n), y = -Inf, vjust = -1) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides(color = 'none') +
    labs()
ggsave(paste0(folder_data, "scraping/medicago/05a-latitude.png"), plot = p, width = 6, height = 5)



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
