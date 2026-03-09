# inat

Tools for downloading biodiversity observation and specimen data.

## Setup

```r
# renv manages all package versions — just restore:
renv::restore()
```

## Data layout

```
data/
  raw/          # untouched downloads (gitignored)
  processed/    # cleaned tables (gitignored)
```

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/01-inat-download.R` | Download iNaturalist research-grade observations for any taxon |
| `scripts/02-tai-specimen.R` | Search herbarium specimens from TAI (tai2.ntu.edu.tw) |
| `scripts/03-clean-smithsonian.R` | Clean raw Smithsonian *Grona triflora* CSV → tidy table |
| `scripts/04-clean-tai.R` | Clean TAI specimens & add county column for *G. triflora*, *G. heterophylla*, *G. heterocarpa* |
| `scripts/05-assess-grona.R` | *Grona* species analysis & mapping in Taiwan |
| `scripts/06-assess-lupulina.R` | *Medicago lupulina* seasonal mapping (North America) |
| `scripts/07-medicago.R` | Multi-species *Medicago* latitude analysis |

## Quick start

```r
# 1. iNaturalist data
source("scripts/01-inat-download.R")
grona <- inat_download("Grona", maxresults = 10000)
grona_tw <- subset_area(grona, lat = c(21, 26), lon = c(120, 123))

# 2. TAI herbarium specimens (searches both Grona and Desmodium names)
source("scripts/02-tai-specimen.R")
specimens <- tai_search("Grona triflora")
specimens_tw <- specimens |> dplyr::filter(country == "Taiwan")

# 3. Clean Smithsonian data (place raw CSV in data/raw/ first)
source("scripts/03-clean-smithsonian.R")

# 4. Clean TAI data & add county columns (runs all three species)
source("scripts/04-clean-tai.R")
# Or clean one species at a time:
# clean_tai("grona_heterophylla")
```

Downloaded data is saved to `data/` as both CSV and RDS.
