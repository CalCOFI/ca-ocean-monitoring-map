############################################
# build_master_map.R
#
# Step 1 — Chunk combine: merges split program
#           GeoJSONs (e.g. CalCOFI_1, CalCOFI_2)
#           into one GeoJSON per program.
#
# Step 2 — Master combine: merges all program
#           GeoJSONs into Master_Inventory.geojson
#           (excludes discharger point layers).
#
# Outputs: Master_Inventory.geojson, transects.csv
#
# Run build_discharger_layer.R first.
############################################

library(tidyverse)
library(sf)

output_root   <- "C:/Users/bhuan/Downloads/Monitoring_Outputs"
combined_name <- "Master_Inventory"

# Names of discharger output folders — excluded from hex combine
discharger_folder_names <- c(
  "Dischargers"
)

# ── Helper ─────────────────────────────────
collapse_unique <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- x[!is.na(x) & str_trim(x) != "" & x != "NA"]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

############################################
# STEP 1 — CHUNK COMBINE
# Merges split program folders (e.g. CalCOFI_1,
# CalCOFI_2) into one GeoJSON per program.
############################################

cat("=== STEP 1: Chunk combine ===\n")

chunk_folders <- list.dirs(output_root, full.names = TRUE, recursive = TRUE) %>%
  { gsub("\\\\", "/", .) } %>%
  .[str_detect(basename(.), "[_ ]?\\d+$")] %>%
  .[!str_detect(., regex("Master_", ignore_case = TRUE))]

if (length(chunk_folders) == 0) {
  cat("No chunk folders found — skipping.\n")
} else {
  
  folder_map <- tibble(
    full_path  = chunk_folders,
    chunk_name = basename(chunk_folders),
    base_name  = str_remove(chunk_name, "[_ ]?\\d+$")
  )
  
  for (prog in unique(folder_map$base_name)) {
    prog_chunks <- folder_map %>% filter(base_name == prog) %>% arrange(chunk_name)
    chunk_names <- prog_chunks$chunk_name
    cat("\nCombining", prog, ":", paste(chunk_names, collapse = ", "), "\n")
    
    chunk_files <- file.path(prog_chunks$full_path, paste0(chunk_names, ".geojson"))
    missing     <- chunk_files[!file.exists(chunk_files)]
    if (length(missing) > 0) { warning("Missing: ", paste(missing, collapse = ", ")); next }
    
    # Single chunk — just copy
    if (length(chunk_names) == 1) {
      out_folder <- file.path(output_root, prog)
      dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
      file.copy(chunk_files, file.path(out_folder, paste0(prog, ".geojson")), overwrite = TRUE)
      cat("  Copied directly.\n")
      transect_src <- file.path(prog_chunks$full_path[1], "transects.csv")
      if (file.exists(transect_src))
        file.copy(transect_src, file.path(out_folder, "transects.csv"), overwrite = TRUE)
      next
    }
    
    chunk_sf_list <- map(chunk_files, st_read, quiet = TRUE)
    
    # Drop empty chunks
    empty <- map_lgl(chunk_sf_list, ~ nrow(.x) == 0)
    if (any(empty)) {
      cat("  Skipping empty chunks:", paste(chunk_names[empty], collapse = ", "), "\n")
      chunk_sf_list <- chunk_sf_list[!empty]
      prog_chunks   <- prog_chunks[!empty, ]
      chunk_names   <- chunk_names[!empty]
    }
    if (length(chunk_sf_list) == 0) { warning("All chunks empty for ", prog); next }
    
    # Align param columns across chunks
    all_param_cols <- chunk_sf_list %>%
      map(~ names(.x)[str_detect(names(.x), "^param_")]) %>%
      reduce(union)
    
    chunk_sf_list <- map(chunk_sf_list, function(sf_obj) {
      gebco_vals <- sf_obj[["Gebco.Mean.Depth"]]  # save before coercion
      for (col in setdiff(all_param_cols, names(sf_obj))) sf_obj[[col]] <- 0L
      out <- sf_obj %>% mutate(across(-geometry, as.character))
      out[["Gebco.Mean.Depth"]] <- gebco_vals      # restore numeric
      out
    })
    
    combined_sf <- bind_rows(chunk_sf_list) %>%
      rename_with(~ str_replace_all(.x, "\\.", " "), everything())
    
    combined_sf <- combined_sf %>%
      mutate(across(where(is.character), ~ na_if(str_trim(.x), "NA"))) %>%

      mutate(`Gebco Mean Depth` = suppressWarnings(as.numeric(`Gebco Mean Depth`)))
    
    names(combined_sf) <- names(combined_sf) %>%
      str_replace("Depth Range  m ", "Depth Range (m)")
    
    lat_col <- names(combined_sf)[str_detect(names(combined_sf), regex("centroid.*lat", ignore_case = TRUE))][1]
    lon_col <- names(combined_sf)[str_detect(names(combined_sf), regex("centroid.*lon", ignore_case = TRUE))][1]
    
    if (is.na(lat_col) | is.na(lon_col)) { warning("No centroid columns for ", prog); next }
    
    merged_df <- combined_sf %>%
      st_drop_geometry() %>%
      group_by(.data[[lat_col]], .data[[lon_col]]) %>%
      summarise(
        `Program Name`      = first(na.omit(`Program Name`)),
        `Full Program Name` = first(na.omit(`Full Program Name`)),
        `Monitoring Program`= first(na.omit(`Monitoring Program`)),
        `Frequency`         = collapse_unique(Frequency),
        `Platform`          = collapse_unique(Platform),
        `First Year`        = suppressWarnings(min(as.numeric(`First Year`),    na.rm = TRUE)),
        `Last Year`         = suppressWarnings(max(as.numeric(`Last Year`),     na.rm = TRUE)),
        `Gebco Mean Depth` = suppressWarnings(mean(as.numeric(`Gebco Mean Depth`), na.rm = TRUE)),
        `Years Sampled`     = suppressWarnings(max(as.numeric(`Years Sampled`), na.rm = TRUE)),
        `Depth Range (m)`   = collapse_unique(`Depth Range (m)`),
        `Sample Locations`  = suppressWarnings(sum(as.numeric(`Sample Locations`), na.rm = TRUE)),
        `Parameters`        = collapse_unique(Parameters),
        `EOV Groups`        = collapse_unique(`EOV Groups`),
        `Parameter Count`   = suppressWarnings(max(as.numeric(`Parameter Count`), na.rm = TRUE)),
        across(all_of(all_param_cols), ~ max(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        `First Year`    = na_if(as.numeric(`First Year`),     Inf),
        `Gebco Mean Depth` = na_if(`Gebco Mean Depth`,          NaN),
        `Last Year`     = na_if(as.numeric(`Last Year`),     -Inf),
        `Years Sampled` = na_if(as.numeric(`Years Sampled`), -Inf)
      )
    
    merged_sf <- combined_sf %>%
      select(all_of(c(lat_col, lon_col)), geometry) %>%
      distinct(.data[[lat_col]], .data[[lon_col]], .keep_all = TRUE) %>%
      left_join(merged_df, by = c(lat_col, lon_col)) %>%
      st_as_sf()
    
    out_folder <- file.path(output_root, prog)
    dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
    out_path <- file.path(out_folder, paste0(prog, ".geojson"))
    suppressWarnings(st_write(merged_sf, out_path, delete_dsn = TRUE, quiet = TRUE))
    cat("  Written:", out_path, "(", nrow(merged_sf), "rows)\n")
    
    # Combine transects for this program
    transect_files <- file.path(prog_chunks$full_path, "transects.csv") %>% .[file.exists(.)]
    if (length(transect_files) > 0) {
      map_dfr(transect_files, read_csv, show_col_types = FALSE,
              col_types = cols(.default = col_character())) %>%
        distinct() %>%
        write_csv(file.path(out_folder, "transects.csv"))
      cat("  Transects written.\n")
    }
  }
}

############################################
# STEP 2 — MASTER COMBINE
# Merges all monitoring hex GeoJSONs.
# Excludes discharger point layers.
############################################

cat("\n=== STEP 2: Master combine ===\n")

geojson_files <- list.files(output_root, pattern = "\\.geojson$",
                            full.names = TRUE, recursive = TRUE) %>%
  { gsub("\\\\", "/", .) } %>%
  .[!str_detect(., regex("Master_Inventory|Master_Parameter_Long", ignore_case = TRUE))] %>%
  .[!str_detect(., regex("/chunks/", ignore_case = TRUE))] %>%
  .[!str_detect(basename(dirname(.)), "^.*[_ ]?\\d+$")] %>%
  # Exclude discharger folders by folder name
  .[!basename(dirname(.)) %in% discharger_folder_names] %>%
  # Exclude root-level discharger GeoJSONs
  .[!tools::file_path_sans_ext(basename(.)) %in% discharger_folder_names]

cat("Monitoring GeoJSON files found:", length(geojson_files), "\n")
print(basename(geojson_files))

combined_sf <- map(geojson_files, st_read, quiet = TRUE) %>%
  map(~ {
    # Column has dots when read from GeoJSON ("Gebco.Mean.Depth")
    gebco_col  <- names(.x)[grep("gebco", names(.x), ignore.case = TRUE)][1]
    gebco_vals <- if (!is.na(gebco_col)) .x[[gebco_col]] else NULL
    out <- .x %>% mutate(across(-geometry, as.character))
    if (!is.na(gebco_col) && !is.null(gebco_vals)) out[[gebco_col]] <- gebco_vals
    out
  }) %>%
  bind_rows() %>%
  st_make_valid()

# Restore readable column names (GeoJSON driver converts spaces to dots)
names(combined_sf) <- names(combined_sf) %>%
  str_replace_all("\\.", " ") %>%
  str_replace("Depth Range  m ", "Depth Range (m)") %>%
  str_squish()

# Drop internal GeoJSON metadata column
combined_sf <- combined_sf %>%
  select(-any_of(c("Geometry Types", "Geometry.Types")))

# Replace Inf/NaN from numeric aggregations
combined_sf <- combined_sf %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x) | is.nan(.x), NA, .x)))

cat("Total rows:", nrow(combined_sf), "\n")
cat("Programs:", paste(unique(combined_sf$`Program Name`), collapse = ", "), "\n")

# Write directly to final path (avoids Windows file-lock issues with temp rename)
final_path <- file.path(output_root, paste0(combined_name, ".geojson"))
if (file.exists(final_path)) { Sys.sleep(0.5); file.remove(final_path) }

suppressWarnings(
  st_write(combined_sf, final_path, delete_dsn = TRUE, quiet = TRUE,
           layer_options = "RFC7946=NO")
)
cat("Written:", final_path, "\n")

############################################
# STEP 2b — MASTER TRANSECTS
############################################

cat("\nCombining transects...\n")

transect_files <- list.files(output_root, pattern = "^transects\\.csv$",
                             full.names = TRUE, recursive = TRUE) %>%
  { gsub("\\\\", "/", .) } %>%
  .[!str_detect(., regex("/chunks/", ignore_case = TRUE))] %>%
  .[!str_detect(dirname(.), "[_ ]?\\d+$")] %>%
  .[!str_detect(., regex("Master_Inventory", ignore_case = TRUE))]

if (length(transect_files) > 0) {
  map_dfr(transect_files, read_csv, show_col_types = FALSE,
          col_types = cols(.default = col_character())) %>%
    distinct() %>%
    write_csv(file.path(output_root, "transects.csv"))
  cat("Master transects written.\n")
} else {
  cat("No transect files found.\n")
}

cat("\n=== PIPELINE COMPLETE ===\n")
cat("Outputs: Master_Inventory.geojson, transects.csv\n")