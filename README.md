# California Ocean & Coastal Monitoring Inventory

An interactive Leaflet map of California ocean and coastal monitoring programs. Displays monitoring coverage as hex grid cells, survey transects, and point-based discharger/WWTP monitoring stations.

---

## Repository Structure

```
ca-ocean-monitoring-map/
├── R/
│   ├── build_program_layer.R      # Process one monitoring program → hex GeoJSON
│   ├── build_discharger_layer.R   # Process discharger CSVs → point GeoJSON
│   └── build_combine_map.R        # Combine all layers → Master_Inventory.geojson
├── web/
│   └── index.html                 # Interactive map (copy to Monitoring_Outputs/)
├── README.md
└── .gitignore
```

---

## Output Folder Structure

Running the R scripts populates a `Monitoring_Outputs/` folder (not committed to GitHub):

```
Monitoring_Outputs/
├── index.html                    ← Copy from web/
├── Master_Inventory.geojson      ← Combined monitoring hex grid
├── transects.csv                 ← Combined survey transects
├── gebco_2025_*.tif              ← GEBCO bathymetry raster (download separately)
├── Dischargers/
│   └── Dischargers.geojson
├── CalCOFI/
│   └── CalCOFI.geojson
├── NOAA THL/
│   └── NOAA THL.geojson
└── SCCWRP/
    └── SCCWRP.geojson
```

---

## Prerequisites

**R packages:**
```r
install.packages(c("readr", "dplyr", "tidyr", "stringr", "purrr",
                   "sf", "terra", "janitor", "tidyverse"))
```

**Input data not included in this repo:**
- Monitoring program CSVs in per-program folders
- Discharger monitoring CSVs
- CA boundary shapefile (`CA_State.shp`)
- Attribute lookup table (`Attribute_Table.csv`)
- GEBCO 2025 bathymetry GeoTIFF (download from [GEBCO](https://www.gebco.net/data_and_products/gridded_bathymetry_data/))

---

## How to Run

### Step 1 — Build each monitoring program layer

Edit USER SETTINGS at the top of `R/build_program_layer.R` and run once per program folder (or chunk). Outputs a GeoJSON and `transects.csv` per program.

### Step 2 — Build discharger layer

Edit USER SETTINGS in `R/build_discharger_layer.R` and run. Outputs `Dischargers/Dischargers.geojson`.

### Step 3 — Combine everything

Edit paths at the top of `R/build_combine_map.R` and run. Outputs `Master_Inventory.geojson` and `transects.csv`.

### Step 4 — Serve the map

Copy `web/index.html` into `Monitoring_Outputs/` then serve with a local web server:

**VS Code:** Right-click `index.html` → *Open with Live Server*

**Python:**
```bash
cd path/to/Monitoring_Outputs
python -m http.server 8000
# Open http://localhost:8000
```

---

## Map Features

- **Hex grid** — monitoring coverage by program, colored and patterned per program
- **Filters** — filter by program, parameter, or EOV group
- **Transects** — survey cruise tracks loaded from `transects.csv`
- **Dischargers** — WWTP/ocean discharger stations, colored per facility
- **Bathymetry** — optional GEBCO 2025 seafloor depth layer
- **Popups** — parameters, depth range, GEBCO seafloor depth, coordinates

---

## Adding a New Monitoring Program

1. Place CSVs in a subfolder of `Monitoring_Outputs/`
2. Run `build_program_layer.R` pointing to that folder
3. Run `build_combine_map.R` to regenerate `Master_Inventory.geojson`

## Adding a New Discharger Layer

1. Run `build_discharger_layer.R` with updated folder/name settings
2. Add folder name to `discharger_folder_names` in `build_combine_map.R`
3. Add entry to `DISCHARGER_SOURCES` in `web/index.html`:
```javascript
const DISCHARGER_SOURCES = [
  { path: 'Dischargers/Dischargers.geojson', label: 'Dischargers' },
  { path: 'NewLayer/NewLayer.geojson',        label: 'New Layer'  }
];
```
