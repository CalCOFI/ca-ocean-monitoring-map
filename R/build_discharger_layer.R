############################################
#V.b

# build_discharger_layer.R
#
# Reads one or more discharger/WWTP monitoring
# CSVs and outputs a point-marker GeoJSON.
# Map display and styling is handled by index.html.
#
# Run this BEFORE build_combined_map.R
############################################

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(sf)
library(janitor)

############################################
# USER SETTINGS — edit these
############################################

discharger_folder     <- "C:/Users/bhuan/Downloads/Monitoring Data/Dischargers"
discharger_layer_name <- "Dischargers"

output_root   <- "C:/Users/bhuan/Downloads/Monitoring_Outputs"
output_folder <- file.path(output_root, discharger_layer_name)
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

############################################
# HELPER FUNCTIONS
############################################

clean_cols <- function(df) {
  df %>% mutate(across(where(is.character), str_squish))
}

collapse_unique <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- x[!is.na(x) & str_trim(x) != ""]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

is_placeholder <- function(x) {
  str_trim(str_to_lower(as.character(x))) %in%
    c("", "na", "n/a", "nan", "null", "none", "nd", "n.d.",
      "-999", "-9999", "999", "9999", "see figure", "tbd", "unknown")
}

# Detect latitude column from common naming patterns
detect_lat_col <- function(df) {
  patterns <- c("lat.*mid|mid.*lat", "^latitude$|^lat$|lat_",
                "latitude.*n|lat.*n", "^lat")
  for (pat in patterns) {
    hits <- names(df)[str_detect(names(df), regex(pat, ignore_case = TRUE))]
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

# Detect longitude column from common naming patterns
detect_lon_col <- function(df) {
  patterns <- c("lon.*mid|long.*mid|mid.*lon", "^longitude$|^lon$|^long$|lon_",
                "longitude.*w|lon.*w", "^lon|^long")
  for (pat in patterns) {
    hits <- names(df)[str_detect(names(df), regex(pat, ignore_case = TRUE))]
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

# CA longitudes should be negative (~-117 to -124); flip if stored as positive
fix_longitude_sign <- function(vals, file_label = "") {
  numeric_vals <- suppressWarnings(as.numeric(vals))
  med <- median(numeric_vals, na.rm = TRUE)
  if (!is.na(med) && med > 100 && med < 130) {
    message("  [", file_label, "] Flipping longitude sign to negative.")
    return(-abs(numeric_vals))
  }
  numeric_vals
}

############################################
# READ + CLEAN ONE CSV
############################################

read_discharger_file <- function(file_path) {
  message("Reading: ", basename(file_path))

  df_raw <- read_csv(file_path, show_col_types = FALSE,
                     col_types = cols(.default = col_character()))

  # Sanitize non-ASCII characters in column names
  raw_names <- names(df_raw) %>%
    str_replace_all("\u00b5", "u") %>%
    str_replace_all("[^[:ascii:]]", "u") %>%
    str_replace_all("/", "_per_")
  names(df_raw) <- raw_names

  df <- df_raw %>%
    clean_names() %>%
    clean_cols() %>%
    mutate(source_file = basename(file_path),
           source_path = file_path)

  lat_col <- detect_lat_col(df)
  lon_col <- detect_lon_col(df)

  if (is.na(lat_col) || is.na(lon_col)) {
    warning("Could not detect lat/lon in ", basename(file_path), " — skipping.")
    return(NULL)
  }

  cat("  lat:", lat_col, "| lon:", lon_col, "\n")

  lat_vals <- suppressWarnings(as.numeric(df[[lat_col]]))
  lon_vals <- fix_longitude_sign(df[[lon_col]], basename(file_path))

  bad <- is_placeholder(df[[lat_col]]) | is_placeholder(df[[lon_col]]) |
         is.na(lat_vals) | is.na(lon_vals)

  if (any(bad)) message("  Dropping ", sum(bad), " row(s) with bad coordinates.")

  df %>%
    mutate(latitude_std  = lat_vals,
           longitude_std = lon_vals,
           coord_ok      = !bad) %>%
    filter(coord_ok) %>%
    select(-coord_ok)
}

############################################
# READ ALL CSVs IN FOLDER
############################################

csv_files <- list.files(discharger_folder, pattern = "\\.csv$",
                        full.names = TRUE, recursive = TRUE) %>%
  # Exclude any output subfolders
  .[!str_detect(., regex("(/|\\\\)(output[^/\\\\]*)(/|\\\\)", ignore_case = TRUE))]

if (length(csv_files) == 0) stop("No CSV files found in discharger_folder.")

cat("\nFound", length(csv_files), "CSV file(s):\n")
print(basename(csv_files))

all_df <- map(csv_files, read_discharger_file) %>%
  compact() %>%
  bind_rows()

if (nrow(all_df) == 0) stop("No plottable rows after reading all files.")
cat("\nTotal plottable rows:", nrow(all_df), "\n")

############################################
# COLLAPSE TO ONE ROW PER UNIQUE LOCATION
############################################

COORD_COLS <- c("latitude_std", "longitude_std", "source_file", "source_path")

# Re-scan files to identify original lat/lon column names to exclude
lat_col_names <- map_chr(csv_files, ~ {
  tryCatch(read_csv(.x, n_max = 1, show_col_types = FALSE,
                    col_types = cols(.default = col_character())) %>%
             clean_names() %>% detect_lat_col(),
           error = function(e) NA_character_)
}) %>% na.omit() %>% unique()

lon_col_names <- map_chr(csv_files, ~ {
  tryCatch(read_csv(.x, n_max = 1, show_col_types = FALSE,
                    col_types = cols(.default = col_character())) %>%
             clean_names() %>% detect_lon_col(),
           error = function(e) NA_character_)
}) %>% na.omit() %>% unique()

exclude_cols  <- unique(c(lat_col_names, lon_col_names, COORD_COLS))
attr_cols     <- setdiff(names(all_df), exclude_cols)

# Determine which attribute columns have multiple values per location
loc_df <- all_df %>%
  mutate(.loc_key = paste(round(latitude_std, 5), round(longitude_std, 5), sep = "_"))

multi_val_cols  <- attr_cols[map_lgl(attr_cols, ~ {
  if (!.x %in% names(loc_df)) return(FALSE)
  loc_df %>% group_by(.loc_key) %>%
    summarise(n = n_distinct(.data[[.x]], na.rm = TRUE), .groups = "drop") %>%
    pull(n) %>% max(na.rm = TRUE) > 1
})]
single_val_cols <- setdiff(attr_cols, multi_val_cols)

grouped <- loc_df %>%
  group_by(.loc_key, latitude_std, longitude_std) %>%
  summarise(
    across(all_of(multi_val_cols),  ~ collapse_unique(.x)),
    across(all_of(single_val_cols), ~ { v <- unique(na.omit(as.character(.x))); if (length(v) == 0) NA_character_ else v[1] }),
    source_files = collapse_unique(source_file),
    n_rows       = n(),
    .groups      = "drop"
  ) %>%
  select(-.loc_key) %>%
  mutate(point_id = paste0("DISCH_", row_number()))

cat("\nUnique station points:", nrow(grouped), "\n")

############################################
# WRITE GEOJSON
############################################

points_sf <- grouped %>%
  st_as_sf(coords = c("longitude_std", "latitude_std"), crs = 4326, remove = FALSE)

geojson_path <- file.path(output_folder, paste0(discharger_layer_name, ".geojson"))
suppressWarnings(st_write(points_sf, geojson_path, delete_dsn = TRUE, quiet = TRUE))

cat("\n===== BUILD COMPLETE =====\n")
cat("CSV files read:        ", length(csv_files), "\n")
cat("Total rows parsed:     ", nrow(all_df), "\n")
cat("Unique station points: ", nrow(grouped), "\n")
cat("GeoJSON:               ", geojson_path, "\n")
cat("==========================\n")
