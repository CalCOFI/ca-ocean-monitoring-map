###################################################################################
# build_program_layer.R
#
# Processes one monitoring program's CSV files into a hex grid GeoJSON.
#
# Run once per program
# After all programs are processed, run build_combined_map.R to combine everything.
#
# Outputs (written to output_folder):
#   - <program>.geojson       hex grid map layer
#   - <program>_map.html      standalone Leaflet map
#   - transects.csv           survey transect lines
#   - diagnostic CSVs         parameter + QA tables
###################################################################################

###########
# LIBRARIES
###########

library(tidyverse)
library(janitor)
library(sf)
library(terra)
library(here)
library(glue)

############################
# USER SETTINGS — edit these
############################

# Temp directory (R uses this for large file operations)
dir_tmp  <- here("tmp")
dir_data <- here("data")
dir_out  <- here("outputs")
dir.create(dir_tmp, showWarnings = F)
Sys.setenv(VROOM_TEMP_PATH = dir_tmp)
Sys.setenv(TMPDIR          = dir_tmp)
Sys.setenv(TMP             = dir_tmp)
Sys.setenv(TEMP            = dir_tmp)
# dir.create("C:/Users/bhuan/Documents/R_temp", showWarnings = FALSE, recursive = TRUE)

# Folder containing this program's CSV files.
# Chunk suffix (e.g. _1, _2) is stripped automatically for output naming.
# program_folder <- "C:/Users/bhuan/Downloads/Monitoring Data/CalCOFI_5"
program_folder <- glue("{dir_data}/CalCOFI_5")

# output_root   <- "C:/Users/bhuan/Downloads/Monitoring_Outputs"
# ca_boundary_path    <- "C:/Users/bhuan/Downloads/Monitoring Data/ca_state/CA_State.shp"
# attribute_table_path <- "C:/Users/bhuan/Downloads/Monitoring Data/Attribute_Table.csv"
# gebco_raster_path <- "C:/Users/bhuan/Downloads/Monitoring_Outputs/gebco_2025_n48.0_s30.0_w-130.0_e-110.0_geotiff.tif"
output_root          <- dir_out
ca_boundary_path     <- glue("{dir_out}/ca_state/CA_State.shp")
attribute_table_path <- glue("{dir_out}/Attribute_Table.csv")
gebco_raster_path    <- glue("{dir_out}/gebco_2025_n48.0_s30.0_w-130.0_e-110.0_geotiff.tif")

start_year         <- 2000    # Exclude records before this year
active_cutoff_year <- 2015    # Programs with no data after this year are flagged inactive

buffer_miles  <- 13.4         # Coastal buffer for spatial clip
buffer_meters <- buffer_miles * 1609.34

apply_coastal_clip       <- TRUE
hex_cellsize_m           <- 3500
min_real_values_required <- 5

##############################################
# DERIVED PATHS (auto-set from program_folder)
##############################################

program_name <- basename(program_folder) %>% str_remove("[_ ]?\\d+$")
chunk_name   <- basename(program_folder)
output_folder <- file.path(output_root, program_name, chunk_name)
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

############################################
# LOAD ATTRIBUTE TABLE
# Provides program-level metadata (frequency,
# platform) and per-parameter cross-references.
############################################

attr_table_raw <- read_csv(attribute_table_path, show_col_types = FALSE) %>%
  clean_names() %>%
  filter(!is.na(acronym), str_trim(acronym) != "",
         !is.na(standard_parameter), str_trim(standard_parameter) != "") %>%
  mutate(across(everything(), str_trim))

attr_param_lookup <- attr_table_raw %>%
  select(acronym, standard_parameter,
         attr_frequency = frequency,
         attr_platform  = sampling_platform_ship_buoy_etc) %>%
  filter(!is.na(attr_frequency) | !is.na(attr_platform)) %>%
  distinct()

# Program-level metadata (full name, frequency, platform)
attr_table_programs <- read_csv(attribute_table_path, show_col_types = FALSE) %>%
  clean_names() %>%
  filter(!is.na(acronym), str_trim(acronym) != "") %>%
  mutate(across(everything(), str_trim))

program_metadata <- attr_table_programs %>%
  group_by(acronym) %>%
  summarise(
    full_name = first(na.omit(program)),
    frequency = first(na.omit(frequency)),
    platform  = first(na.omit(sampling_platform_ship_buoy_etc)),
    .groups = "drop"
  ) %>%
  rename(program_name = acronym) %>%
  mutate(across(everything(), ~ if_else(is.na(.x) | str_trim(.x) == "", "Unknown", str_trim(.x))))

program_meta <- program_metadata %>% filter(program_name == !!program_name)

if (nrow(program_meta) == 0) {
  warning("'", program_name, "' not found in Attribute_Table.csv — ",
          "add a row to set full_name, frequency, and platform.")
  display_name       <- chunk_name
  program_full_name  <- "Unknown"
  sampling_frequency <- "Unknown"
  program_platform   <- "Unknown"
} else {
  display_name       <- chunk_name
  program_full_name  <- program_meta$full_name
  sampling_frequency <- program_meta$frequency
  program_platform   <- program_meta$platform
}

cat("Program:", program_name, "| Frequency:", sampling_frequency, "| Platform:", program_platform, "\n")

platform_lookup <- tribble(
  ~pattern,                          ~platform,
  "acoustic",                        "Vessel - Towed Acoustic Array",
  "continuous.*fish.*egg|underway",  "Vessel - Continuous Underway Sampler",
  "glider",                          "Glider",
  "buoy|mooring",                    "Mooring / Buoy",
  "aerial",                          "Aerial",
  "shore|land|beach",                "Shore-based",
  "satellite|model",                 "Satellite / Model",
  ".*",                              program_platform
)

######################################################################
# WORKFLOW-INJECTED COLUMN NAMES
# These are added by the script and excluded from parameter detection.
######################################################################

WORKFLOW_INJECTED_COLS <- c(
  "source_row_id", "source_file", "source_path", "file_stub", "program",
  "detected_lat_col", "detected_lon_col", "detected_year_col",
  "detected_date_col", "detected_program_col", "detected_depth_col",
  "detected_coord_role",
  "year_detected", "latitude_std", "longitude_std", "depth_std",
  "sample_point_key", "hex_id", "station_key",
  "geometry_type", "classification_reason",
  ".row_id_temp", ".join_key", "last_year", "activity_status"
)

##################
# HELPER FUNCTIONS
##################

get_file_stub <- function(path) {
  tools::file_path_sans_ext(basename(path)) %>%
    str_replace_all("[^A-Za-z0-9]+", "_")
}

clean_character_cols <- function(df) {
  df %>% mutate(across(where(is.character), str_squish))
}

collapse_unique_text <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- x[!is.na(x) & str_trim(x) != ""]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

detect_col <- function(df, patterns) {
  nm   <- names(df)
  hits <- nm[map_lgl(nm, ~ any(str_detect(.x, regex(patterns, ignore_case = TRUE))))]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

detect_col_priority <- function(df, pattern_vector, fallback_type = NULL) {
  nm <- names(df)
  for (pat in pattern_vector) {
    hits <- nm[str_detect(nm, regex(pat, ignore_case = TRUE))]
    if (length(hits) > 0) return(hits[1])
  }
  if (!is.null(fallback_type) && fallback_type == "lat") {
    hits <- nm[str_detect(nm, regex("latitude|lat|y", ignore_case = TRUE))]
    if (length(hits) > 0) return(hits[1])
  }
  if (!is.null(fallback_type) && fallback_type == "lon") {
    hits <- nm[str_detect(nm, regex("longitude|lon|long|x", ignore_case = TRUE))]
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

# Parses seasonal strings ("Spring 2018", "Q2 2019") into approximate dates
parse_season_to_approx_date <- function(x) {
  x_clean    <- str_squish(str_to_lower(as.character(x)))
  year_match <- str_extract(x_clean, "\\d{4}")
  year_val   <- suppressWarnings(as.integer(year_match))
  season_md  <- case_when(
    str_detect(x_clean, "spring")        ~ "04-15",
    str_detect(x_clean, "summer")        ~ "07-15",
    str_detect(x_clean, "fall|autumn")   ~ "10-15",
    str_detect(x_clean, "winter")        ~ "01-15",
    str_detect(x_clean, "q1|quarter.?1") ~ "02-15",
    str_detect(x_clean, "q2|quarter.?2") ~ "05-15",
    str_detect(x_clean, "q3|quarter.?3") ~ "08-15",
    str_detect(x_clean, "q4|quarter.?4") ~ "11-15",
    TRUE                                  ~ NA_character_
  )
  suppressWarnings(as.Date(
    if_else(!is.na(season_md) & !is.na(year_val),
            paste(year_val, season_md, sep = "-"),
            NA_character_)
  ))
}

contains_seasonal_language <- function(x) {
  x_clean <- str_squish(str_to_lower(as.character(x)))
  any(str_detect(x_clean, "spring|summer|fall|autumn|winter|\\bq[1-4]\\b|quarter"), na.rm = TRUE)
}

parse_date_time_safe <- function(x) {
  x <- as.character(x)
  formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d",
               "%Y-%m-%d %H:%M:%S", "%m/%d/%Y %H:%M:%S", "%Y%m%d")
  result <- rep(as.Date(NA), length(x))
  for (fmt in formats) {
    idx    <- which(is.na(result))
    if (length(idx) == 0) break
    result[idx] <- suppressWarnings(as.Date(x[idx], format = fmt))
  }
  still_na <- which(is.na(result))
  if (length(still_na) > 0)
    result[still_na] <- parse_season_to_approx_date(x[still_na])
  result
}

extract_year <- function(df, year_col = NA, date_col = NA) {
  out <- rep(NA_integer_, nrow(df))
  if (!is.na(year_col) && year_col %in% names(df))
    out <- suppressWarnings(as.integer(df[[year_col]]))
  if (!is.na(date_col) && date_col %in% names(df)) {
    yr <- suppressWarnings(as.integer(format(parse_date_time_safe(df[[date_col]]), "%Y")))
    out[is.na(out)] <- yr[is.na(out)]
  }
  out
}

is_placeholder_value <- function(x) {
  str_trim(str_to_lower(as.character(x))) %in%
    c("", "na", "n/a", "nan", "null", "none", "nd", "n.d.",
      "-999", "-9999", "-99999", "999", "9999", "99999")
}

count_real_values <- function(x) {
  sum(!is.na(as.character(x)) & !is_placeholder_value(x))
}

normalize_col_name <- function(x) {
  x %>% str_to_lower() %>% str_replace_all("[^a-z0-9]+", "")
}

clean_text_value <- function(x) {
  x %>% tolower() %>% str_replace_all("[^a-z0-9]+", " ") %>% str_squish()
}

# Parses coordinate strings including DMS, hemisphere letters, degree symbols
parse_coord_vector <- function(vals, is_lon = FALSE, file_label = "", col_label = "") {
  x_chr   <- str_trim(as.character(vals))
  is_ph   <- is_placeholder_value(x_chr) | x_chr == ""
  hemi    <- str_extract(x_chr, regex("[NSEWnsew]", ignore_case = TRUE))
  x_clean <- str_trim(str_remove_all(
    x_chr,
    paste0("[NSEWnsew\u00b0\u00ba\u2032\u2033\u2019\u201d'\"", "]")
  ))
  result  <- suppressWarnings(as.numeric(x_clean))
  sw_mask <- !is.na(hemi) & str_detect(hemi, regex("[SWsw]"))
  result  <- if_else(!is.na(result) & sw_mask, -abs(result), result)

  need_dms <- which(is.na(result) & !is_ph)
  if (length(need_dms) > 0) {
    dms <- vapply(need_dms, function(i) {
      parts <- suppressWarnings(
        as.numeric(str_split(x_clean[i], "[:\\s]+|(?<=[0-9])-(?=[0-9])")[[1]])
      )
      parts <- parts[!is.na(parts)]
      if      (length(parts) == 2) parts[1] + parts[2] / 60
      else if (length(parts) == 3) parts[1] + parts[2] / 60 + parts[3] / 3600
      else                         NA_real_
    }, numeric(1))
    neg  <- is.na(hemi[need_dms]) & str_starts(x_chr[need_dms], "-")
    sw2  <- !is.na(hemi[need_dms]) & str_detect(hemi[need_dms], regex("[SWsw]"))
    dms  <- if_else(neg | sw2, -abs(dms), dms)
    result[need_dms] <- dms
  }

  if (is_lon) result <- if_else(!is.na(result) & (result < -180 | result > 180), NA_real_, result)
  else        result <- if_else(!is.na(result) & (result < -90  | result > 90),  NA_real_, result)

  n_failed <- sum(is.na(result) & !is_ph & !is.na(x_chr))
  if (n_failed > 0)
    warning(n_failed, " coordinate(s) could not be parsed in ", file_label,
            " [", col_label, "] — those rows will be dropped.")
  result
}

validate_coord_col <- function(df, detected_col, all_patterns) {
  if (!is.na(detected_col) && detected_col %in% names(df) &&
      count_real_values(suppressWarnings(as.numeric(df[[detected_col]]))) > 0)
    return(detected_col)
  for (pat in all_patterns) {
    for (candidate in names(df)[str_detect(names(df), regex(pat, ignore_case = TRUE))]) {
      if (!identical(candidate, detected_col) &&
          count_real_values(suppressWarnings(as.numeric(df[[candidate]]))) > 0)
        return(candidate)
    }
  }
  detected_col
}

interval_to_frequency_label <- function(median_days) {
  case_when(
    is.na(median_days)    ~ "Unknown",
    median_days <= 1      ~ "Daily or sub-daily",
    median_days <= 10     ~ "Weekly",
    median_days <= 45     ~ "Monthly",
    median_days <= 100    ~ "Quarterly",
    median_days <= 200    ~ "Semi-annual",
    median_days <= 400    ~ "Annual",
    TRUE                  ~ "Multi-year / Episodic"
  )
}

# Estimates per-parameter sampling frequency from dates or year counts
compute_parameter_frequency <- function(param_hits_df, raw_df, fallback_frequency) {

  file_date_map <- raw_df %>%
    distinct(source_file, detected_date_col) %>%
    filter(!is.na(detected_date_col))

  parsed_date_lookup <- pmap_dfr(file_date_map, function(source_file, detected_date_col) {
    raw_slice <- raw_df %>%
      filter(source_file == !!source_file) %>%
      select(source_file, source_row_id, any_of(detected_date_col))
    if (!detected_date_col %in% names(raw_slice)) return(tibble())
    raw_vals <- raw_slice[[detected_date_col]]
    raw_slice %>%
      mutate(
        parsed_date           = parse_date_time_safe(.data[[detected_date_col]]),
        frequency_source_hint = if_else(
          contains_seasonal_language(raw_vals),
          "estimated_from_seasonal_labels", "computed_from_dates"
        )
      ) %>%
      select(source_file, source_row_id, parsed_date, frequency_source_hint)
  })

  if (!all(c("source_file","source_row_id","parsed_date","frequency_source_hint") %in% names(parsed_date_lookup)))
    parsed_date_lookup <- tibble(source_file = character(), source_row_id = integer(),
                                 parsed_date = as.Date(character()), frequency_source_hint = character())

  hits_with_dates <- param_hits_df %>%
    left_join(parsed_date_lookup, by = c("source_file", "source_row_id"))

  freq_from_dates <- hits_with_dates %>%
    filter(!is.na(parsed_date)) %>%
    group_by(source_file, standard_parameter) %>%
    summarise(
      n_obs                = n_distinct(parsed_date),
      frequency_source     = first(frequency_source_hint),
      median_interval_days = {
        d <- sort(unique(parsed_date))
        if (length(d) < 2) NA_real_ else {
          ivs <- as.numeric(diff(d))
          lg  <- ivs[ivs >= 14]
          if (length(lg) >= 2) median(lg) else median(ivs)
        }
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(median_interval_days)) %>%
    mutate(median_interval_days = if_else(
      frequency_source == "estimated_from_seasonal_labels" & between(median_interval_days, 75, 110),
      91.25, median_interval_days
    )) %>%
    group_by(standard_parameter) %>%
    slice_min(order_by = median_interval_days, n = 1, with_ties = FALSE) %>%
    ungroup()

  needs_fallback <- param_hits_df %>% distinct(standard_parameter) %>%
    anti_join(freq_from_dates, by = "standard_parameter")

  freq_from_years <- hits_with_dates %>%
    filter(standard_parameter %in% needs_fallback$standard_parameter,
           is.na(parsed_date), !is.na(year_detected)) %>%
    group_by(standard_parameter) %>%
    summarise(
      n_obs        = n(),
      obs_per_year = n() / max(n_distinct(year_detected), 1),
      .groups = "drop"
    ) %>%
    filter(obs_per_year <= 5) %>%
    mutate(
      median_interval_days = case_when(
        obs_per_year >= 3   ~ 91.25,
        obs_per_year >= 1.5 ~ 182,
        TRUE                ~ 365
      ),
      frequency_source = "estimated_from_year_counts"
    )

  bind_rows(
    freq_from_dates %>% select(standard_parameter, n_obs, median_interval_days, frequency_source),
    freq_from_years %>% select(standard_parameter, n_obs, median_interval_days, frequency_source)
  ) %>%
    mutate(frequency_label = interval_to_frequency_label(median_interval_days)) %>%
    right_join(param_hits_df %>% distinct(standard_parameter), by = "standard_parameter") %>%
    mutate(
      frequency_label  = coalesce(frequency_label, fallback_frequency, "Unknown"),
      frequency_source = coalesce(frequency_source, "program_metadata_fallback")
    ) %>%
    select(standard_parameter, frequency_label, frequency_source, median_interval_days, n_obs) %>%
    arrange(standard_parameter)
}

############################################
# COLUMN DETECTION PATTERNS
############################################

lat_priority_patterns  <- c("lat.*mid|mid.*lat", "^latitude$|^lat$|^lat_|_lat$|_lat_|latitude",
                             "^cc_lat", "lat.*start|start.*lat", "lat.*stop|stop.*lat")
lon_priority_patterns  <- c("lon.*mid|long.*mid|mid.*lon|mid.*long",
                             "^longitude$|^lon$|^long$|^lon_|_lon$|_lon_|^long_|_long$|longitude",
                             "^cc_lon", "lon.*start|long.*start|start.*lon|start.*long",
                             "lon.*stop|long.*stop|stop.*lon|stop.*long")
year_patterns  <- "^year$|^year_|_year$|sample_year|sampling_year|obs_year|collection_year|survey_year"
date_patterns  <- "^date$|^date_|_date$|sample_date|cast_date|collection_date|date_time|datetime|julian|^time$|timestamp"
depth_patterns <- "^depth$|^depth_|_depth$|depth_m|depthm|bottom_depth|sample_depth|station_depth|mean_depth|max_depth|pressure"
program_patterns <- "^program$|^program_|_program$|monitoring_program|program_name|project|agency|dataset"

############################################
# METADATA COLUMN EXCLUSION PATTERNS
# Columns matching these are not treated as
# measurement parameters during detection.
############################################

metadata_patterns <- c(
  "^source_file$","^source_path$","^file_stub$","^program$","^source_row_id$","^detected_",
  "^year_detected$","^latitude_std$","^longitude_std$","^depth_std$",
  "^station_key$","^sample_point_key$","^detected_coord_role$",
  "^year$","^month$","^day$","^date$","^time$","^datetime$",
  "^lat$","^latitude$","^lon$","^long$","^longitude$",
  "^depth$","^depth_m$","^depthm$","^pressure$",
  "^sample_id$","^cast_id$","^event_id$","^record_id$","^id$","^index$",
  "^flag","^qc","^quality","^comment$","^comments$","^notes$","^note$",
  "^cruise","^agency$","^platform$","^vessel$","^ship$",
  "^hex_id$","^geometry_type$","^classification_reason$"
)

exclude_patterns_broad <- c(
  "_qual$","_qc$","_flag$","_prec$","_cnt$",
  "^btl_cnt$","^btl_num$","^btlnum$",
  "^station_lat$","^station_lon$","^lat$","^lon$","^long$","^latitude$","^longitude$",
  "^depth_id$","^station_id$","^sample_id$","^cast_id$","^event_id$","^record_id$","^taxon_id$",
  "^collection$","^female$","^male$","^sex$","^mesh_size$","^mesh$",
  "^scale_factor$","^conversion_factor$","^proportion$","^percent_sampled$",
  "^area$","^area_m2$","^bin_number$","^julian_date$","^julian_day$",
  "^latitude_start$","^latitude_stop$","^latitude_mid$",
  "^longitude_start$","^longitude_stop$","^longitude_mid$",
  "^time_sec$","^transect_number$","^length_m$","^width_m$",
  "^effort$","^observer$","^on_effort$","^off_effort$",
  "^segment$","^course$","^speed$","^bearing$","^distance$","^strip_width$",
  "^start_time$","^stop_time$","^end_time$",
  "^line$","^line_number$","^station$","^station_number$",
  "^cst_cnt$","^sta_id$","^depthm$","^recind$",
  "^t_prec$","^t_qual$","^s_prec$","^s_qual$","^p_qual$","^o_qual$",
  "^sthtaq$","^o2satq$","^chlqua$","^phaqua$","^po4q$","^sio3qu$",
  "^no2q$","^no3q$","^nh3q$","^dic_quality_comment$",
  "^lat_deg$","^lon_deg$","^lon_min$","^bottom_d$","^secchi$",
  "^ship_code$","^data_type$","^order_occ$","^cruz_leg$","^inc_str$","^inc_end$",
  "^pst_lan$","^civil_t$","^timezone$","^time_zone$",
  "^wave_dir$","^wave_ht$","^wave_prd$","^wind_dir$","^wind_spd$","^barometer$",
  "^dry_t$","^wet_t$","^wea$","^cloud_typ$","^cloud_amt$","^visibility$",
  "^int_chl$","^depth_strata$","^depth_stratum$","^strata$","^stratum$",
  "^c14as","^c14a","^darkas$","^darkap$","^darkaq$",
  "^meanas$","^meanap$","^meanaq$","^inctim$","^r_depth$","^r_temp$",
  "^r_sal$","^r_dynht$","^r_nuts$","^r_oxy","^lightp$"
)

is_probably_metadata_column <- function(col_name) {
  str_detect(col_name, regex(paste(c(
    "^x$","^y$","^latitude$","^longitude$","^lat$","^lon$","^long$",
    "^date$","^time$","^datetime$","^year$","^month$","^day$","^timestamp$",
    "^station$","^station id$","^station number$","^cast$","^cast id$",
    "^bottle$","^btl$","^line$","^cruise$","^survey id$",
    "^sample id$","^record id$","^event id$","^id$","^index$","^depth id$",
    "^source file$","^source path$","^source row$","^file stub$","^program$",
    "^latitude std$","^longitude std$","^year detected$","^hex id$",
    "^detected lat$","^detected lon$","^sample point key$","^station key$",
    "^geometry type$","^classification reason$","^activity status$","^last year$",
    "^qual$","qual$","^prec$","^flag$"," flag$","^qc$"," qc$",
    "^r depth$","^r dynht$","^r nuts$","^r temp$","^r sal$","^r oxy$",
    "^darkas$","^darkap$","^darkaq$","^meanas$","^meanap$","^meanaq$",
    "^inctim$","^lightp$"
  ), collapse = "|"), ignore_case = TRUE))
}

############################################
# CONTEXT + PARAMETER DICTIONARIES
# These map file/column names → standard
# parameter names and EOV groups.
############################################

file_context_dictionary <- tribble(
  ~pattern, ~context,
  "\\btemp\\b|temperature",                                   "temperature",
  "\\bsalin|salinity",                                        "salinity",
  "\\boxygen\\b|\\bo2\\b|dissolved.*oxygen",                 "oxygen",
  "\\bdensity\\b|\\bsigma\\b|\\btheta\\b",                  "density",
  "\\bnitrate\\b|\\bno3\\b",                                 "nitrate",
  "\\bnitrite\\b|\\bno2\\b",                                 "nitrite",
  "\\bphosphate\\b|\\bpo4\\b",                               "phosphate",
  "\\bsilicate\\b|\\bsio3\\b",                               "silicate",
  "\\bammonium\\b|\\bnh3\\b|\\bnh4\\b",                     "ammonium",
  "\\bchlorophyll\\b|\\bchl\\b|chlora",                     "chlorophyll_a",
  "\\bphaeo\\b|phaeop",                                      "phaeopigments",
  "\\bproductivity\\b|\\b14c\\b|primary.*prod",             "primary_production",
  "\\bdic\\b|dissolved.*inorganic.*carbon|\\btco2\\b",      "dic",
  "\\balkalinity\\b|\\btalk\\b|total.*alk",                 "alkalinity",
  "\\bpco2\\b|\\bfco2\\b",                                  "pco2",
  "(^|[_\\-\\s])ph([_\\-\\s]|$)|\\bseawater_ph\\b",        "ph",
  "transmissiv",                                             "transmissivity",
  "radiative.*flux|\\birradiance\\b|\\bpar\\b",             "radiative_flux",
  "marine.*mammal|\\bcetacean\\b|\\bwhale\\b|\\bdolphin\\b","marine_mammal_visual",
  "acoustic|bioacoustic",                                    "whale_acoustic",
  "\\bseabird\\b|sea.*bird",                                 "seabird",
  "fish.*egg|egg.*count|egg.*stage|continuous.*fish.*egg",  "fish_egg",
  "fish.*larvae|fish.*larval|\\blarval\\b|\\bichthyoplankton\\b|rockfish.*recruit", "fish_larvae",
  "\\binvertebrate\\b|\\binvertabrate\\b",                  "invertebrate",
  "\\bzooplankton\\b|\\bholoplankton\\b",                   "zooplankton",
  "\\bkrill\\b|\\beuphausii",                               "krill",
  "\\btrawl\\b|\\bhaul\\b|length.*frequency|\\bspecimen\\b","trawl_biota",
  "\\bpigment\\b|\\bhplc\\b|\\bfucoxanthin\\b",            "phytoplankton_pigments",
  "\\bphytoplankton\\b",                                    "phytoplankton_abundance",
  "\\bcdom\\b",                                             "cdom",
  "\\bphycoerythrin\\b|(^|[_\\-])pe([_\\-]|$)",            "pe_fluorescence",
  "\\bfluorescence\\b|\\bfluorometer\\b",                   "chlorophyll_fluorescence",
  "\\bmicrob|\\bbacteria\\b|\\b16s\\b|\\b18s\\b|\\botu\\b|\\basv\\b|\\bmrna\\b", "microbial",
  "picoplankton.*abundance|\\bprochlorococcus\\b|\\bsynechococcus\\b", "picoplankton_abundance",
  "bacterial.*abundance|bacteria.*abundance",               "bacterial_abundance"
)

column_pattern_dictionary <- tribble(
  ~pattern, ~measurement_type,
  "\\begg_count\\b|\\beggs_10m2\\b|\\beggs_100m3\\b",               "egg_count",
  "\\begg_stage\\b|\\begg_stages\\b",                                "egg_stage",
  "\\blarval_length\\b|\\blarval_size\\b",                           "larvae_size",
  "\\blarval_count\\b|\\blarval_abundance\\b",                       "larvae_count",
  "\\bdisplacement_volume\\b|\\bzooplankton_volume\\b",              "volume",
  "\\bstandard_length\\b|\\bfork_length\\b|\\btotal_length\\b",     "size",
  "\\babundance\\b|\\bcount\\b|\\bnumber\\b",                        "abundance",
  "\\bscientific_name\\b|\\btaxon\\b|\\bspecies\\b",                "observation_presence",
  "^sex$|^gender$",                                                  "sex_structure",
  "^weight$|^wet_weight$|^dry_weight$",                              "carbon_biomass",
  "\\bpigment\\b|\\bfucoxanthin\\b|\\bperidinin\\b",                "pigment",
  "\\bcdom\\b",                                                      "cdom",
  "\\bphycoerythrin\\b|pe_fluorescence",                            "pe_fluorescence",
  "\\bfluorescence\\b|\\bchl_fluorescence\\b",                      "chlorophyll_fluorescence",
  "\\bmrna\\b|\\bgene_expression\\b|\\btranscript\\b",              "mrna",
  "\\b16s\\b|\\b18s\\b|\\botu\\b|\\basv\\b|\\btaxonomy\\b",        "community_composition",
  "\\bcarbon\\b|\\bbiomass\\b",                                     "carbon_biomass",
  "\\btemperature\\b|\\bt_deg_c\\b|\\bctd_temp\\b",                "direct_value",
  "\\bsalinity\\b|\\bsalnty\\b|\\bctd_sal\\b",                     "direct_value",
  "\\boxygen\\b|\\bo2ml_l\\b|\\bdo_mgl\\b|oxy_mol_kg|dissolved_oxygen", "direct_value",
  "\\bsigma_theta\\b|\\bstheta\\b|\\bdensity\\b",                  "direct_value",
  "\\bnitrate\\b|\\bno3\\b|\\bno3um\\b",                           "direct_value",
  "\\bnitrite\\b|\\bno2\\b|\\bno2um\\b",                           "direct_value",
  "\\bphosphate\\b|\\bpo4\\b|\\bpo4um\\b",                         "direct_value",
  "\\bsilicate\\b|\\bsio3\\b|\\bsio3um\\b",                        "direct_value",
  "\\bammonium\\b|\\bnh3\\b|\\bnh4\\b|\\bnh3um\\b",               "direct_value",
  "\\bchlorophyll\\b|\\bchlora\\b|\\bchl_a\\b|\\bchlorophylla\\b", "direct_value",
  "\\bphaeo\\b|\\bphaeop\\b",                                       "direct_value",
  "\\bproductivity\\b|\\b14c\\b|primary.*prod|\\bintc14\\b",       "direct_value",
  "\\bdic\\b|dissolved.*inorganic.*carbon|\\bdic1\\b|\\bdic2\\b",  "direct_value",
  "\\btalk\\b|\\balkalinity\\b|\\bta1\\b|\\bta2\\b",              "direct_value",
  "^ph$|^p_h$|^ph_|\\bph1\\b|\\bph2\\b|sea.*ph|ph.*total",        "direct_value",
  "\\bpco2\\b|\\bfco2\\b",                                         "direct_value",
  "\\btransmissiv\\b|beam.*attenuation",                           "direct_value",
  "radiative.*flux|\\birradiance\\b|\\bpar\\b",                    "direct_value",
  "\\bo2sat\\b|oxygen.*sat|\\bdo_sat\\b|do.*percent",             "direct_value"
)

context_parameter_map <- tribble(
  ~context,               ~measurement_type,       ~standard_parameter,                                         ~eov_group,
  "temperature",          "direct_value",          "Temperature",                                               "Physical",
  "salinity",             "direct_value",          "Salinity",                                                  "Physical",
  "oxygen",               "direct_value",          "Dissolved Oxygen",                                          "Biogeochemical",
  "density",              "direct_value",          "Density (Sigma Theta)",                                     "Physical",
  "nitrate",              "direct_value",          "Nitrate",                                                   "Biogeochemical",
  "nitrite",              "direct_value",          "Nitrite",                                                   "Biogeochemical",
  "phosphate",            "direct_value",          "Phosphate",                                                 "Biogeochemical",
  "silicate",             "direct_value",          "Silicate",                                                  "Biogeochemical",
  "ammonium",             "direct_value",          "Ammonium",                                                  "Biogeochemical",
  "chlorophyll_a",        "direct_value",          "Chlorophyll-a",                                             "Biogeochemical",
  "phaeopigments",        "direct_value",          "Phaeopigments",                                             "Biogeochemical",
  "primary_production",   "direct_value",          "Primary Production (14C uptake)",                          "Biogeochemical",
  "dic",                  "direct_value",          "Dissolved Inorganic Carbon (DIC)",                         "Biogeochemical",
  "alkalinity",           "direct_value",          "Total Alkalinity",                                          "Biogeochemical",
  "ph",                   "direct_value",          "pH",                                                        "Biogeochemical",
  "pco2",                 "direct_value",          "pCO2",                                                      "Biogeochemical",
  "transmissivity",       "direct_value",          "Transmissivity",                                            "Physical",
  "radiative_flux",       "direct_value",          "Radiative Flux",                                            "Physical",
  "marine_mammal_visual", "observation_presence",  "Marine Mammal Abundance",                                   "Biological",
  "whale_acoustic",       "observation_presence",  "Whale Acoustic",                                            "Biological",
  "seabird",              "observation_presence",  "Seabird Count",                                             "Biological",
  "fish_egg",             "egg_count",             "Fish Egg Counts",                                           "Biological",
  "fish_egg",             "egg_stage",             "Fish Egg Stages",                                           "Biological",
  "fish_egg",             "observation_presence",  "Fish Egg Presence",                                         "Biological",
  "fish_larvae",          "larvae_count",          "Fish Larvae Counts",                                        "Biological",
  "fish_larvae",          "larvae_size",           "Fish Larvae Sizes",                                         "Biological",
  "fish_larvae",          "observation_presence",  "Fish Larvae Presence",                                      "Biological",
  "trawl_biota",          "abundance",             "Fish Larvae Counts",                                        "Biological",
  "trawl_biota",          "size",                  "Fish Larvae Sizes",                                         "Biological",
  "trawl_biota",          "observation_presence",  "Fish Larvae Presence",                                      "Biological",
  "invertebrate",         "size",                  "Invertebrate Size",                                         "Biological",
  "invertebrate",         "abundance",             "Invertebrate Abundance",                                    "Biological",
  "invertebrate",         "observation_presence",  "Invertebrate Abundance",                                    "Biological",
  "zooplankton",          "abundance",             "Zooplankton Abundance",                                     "Biological",
  "zooplankton",          "observation_presence",  "Zooplankton Abundance",                                     "Biological",
  "krill",                "abundance",             "Krill (Euphausiid) Abundance",                             "Biological",
  "krill",                "observation_presence",  "Krill (Euphausiid) Abundance",                             "Biological",
  "phytoplankton_pigments","pigment",              "Phytoplankton Taxon-Specific Pigments",                    "Biological",
  "phytoplankton_abundance","abundance",           "Phytoplankton Abundance",                                   "Biological",
  "cdom",                 "cdom",                  "CDOM Fluorescence",                                         "Biogeochemical",
  "pe_fluorescence",      "pe_fluorescence",       "Phycoerythrin (PE) Fluorescence",                          "Biogeochemical",
  "chlorophyll_fluorescence","chlorophyll_fluorescence","Chlorophyll fluorescence",                            "Biogeochemical",
  "microbial",            "community_composition", "Microbial community composition",                           "Biological",
  "microbial",            "mrna",                  "Microbial Genomics (mRNA)",                                 "Biological",
  "picoplankton_abundance","abundance",            "Picoplankton abundance",                                    "Biological",
  "bacterial_abundance",  "abundance",             "Bacterial abundance",                                       "Biological"
)

parameter_dictionary <- tribble(
  ~standard_parameter,                     ~eov_group,       ~raw_pattern,                                                                                                     ~detection_type, ~assignment_scope,
  "Temperature",                           "Physical",       "\\btemp\\b|\\btemperature\\b|t degc|ctd temp|ctdtemp|cc temp|water temp|\\bsst\\b",                            "measurement",   "column",
  "Salinity",                              "Physical",       "\\bsalinity\\b|\\bsalnty\\b|cc sal|ctd sal|ctdsal|practical salinity",                                          "measurement",   "column",
  "Density (Sigma Theta)",                 "Physical",       "sigma theta|\\bsigma\\b|\\bstheta\\b|s theta|\\bdynht\\b|\\bdensity\\b",                                       "measurement",   "column",
  "Transmissivity",                        "Physical",       "transmissiv|beam.*attenuation|optical.*attenuation",                                                             "measurement",   "column",
  "Radiative Flux",                        "Physical",       "radiative.*flux|solar.*flux|\\birradiance\\b|\\bpar\\b",                                                        "measurement",   "column",
  "Dissolved Oxygen",                      "Biogeochemical", "\\boxygen\\b|\\bo2ml\\b|o2ml l|oxy umol|dissolved.*oxygen|do_mg|do_mgl",                                       "measurement",   "column",
  "Oxygen Saturation",                     "Biogeochemical", "o2sat|oxygen.*sat|sat.*oxygen|\\bdo_sat\\b|do.*percent",                                                        "measurement",   "column",
  "Nitrate",                               "Biogeochemical", "\\bno3\\b|\\bnitrate\\b|no3_um|no3um",                                                                          "measurement",   "column",
  "Nitrite",                               "Biogeochemical", "\\bno2\\b|\\bnitrite\\b|no2_um|no2um",                                                                          "measurement",   "column",
  "Phosphate",                             "Biogeochemical", "\\bpo4\\b|\\bphosphate\\b|po4_um|po4um",                                                                        "measurement",   "column",
  "Silicate",                              "Biogeochemical", "\\bsio3\\b|\\bsilicate\\b|sio3_um|sio3um",                                                                      "measurement",   "column",
  "Ammonium",                              "Biogeochemical", "\\bnh3\\b|\\bnh4\\b|\\bammonium\\b|\\bammonia\\b|nh3_um|nh4_um|nh3um|nh4um",                                   "measurement",   "column",
  "Chlorophyll-a",                         "Biogeochemical", "\\bchlorophyll\\b|cc chl|\\bchlora\\b|\\bchl_a\\b|\\bt_chla\\b|\\bchla\\b|\\blogchl\\b",                      "measurement",   "column",
  "Phaeopigments",                         "Biogeochemical", "\\bphaeop\\b|\\bphaeo\\b|\\bpheopigment\\b",                                                                   "measurement",   "column",
  "Primary Production (14C uptake)",       "Biogeochemical", "\\b14c\\b|\\bc14\\b|14c uptake|primary.*prod|\\bproductivity\\b|\\blightp\\b|\\bmeanap\\b",                    "measurement",   "column",
  "Dissolved Inorganic Carbon (DIC)",      "Biogeochemical", "\\bdic\\b|dissolved inorganic carbon|\\bdic1\\b|\\bdic2\\b",                                                   "measurement",   "column",
  "Total Alkalinity",                      "Biogeochemical", "\\btalk\\b|\\balkalinity\\b|total alkalinity|\\bta1\\b|\\bta2\\b",                                             "measurement",   "column",
  "pH",                                    "Biogeochemical", "^ph$|^p_h$|^ph_|\\bph1\\b|\\bph2\\b|sea.*ph|ph.*total|seawater.*ph",                                          "measurement",   "column",
  "pCO2",                                  "Biogeochemical", "\\bpco2\\b|\\bfco2\\b|co2 partial pressure",                                                                   "measurement",   "column",
  "CDOM Fluorescence",                     "Biogeochemical", "\\bcdom\\b|chromophoric dissolved organic matter",                                                              "measurement",   "column",
  "Chlorophyll fluorescence",              "Biogeochemical", "^cf_rb|^cf_rg|^fv_fm|\\bfluorescence\\b|\\bfluorometer\\b",                                                   "measurement",   "column",
  "Phycoerythrin (PE) Fluorescence",       "Biogeochemical", "\\bphycoerythrin\\b|pe fluorescence|^pe1_|^pe2_|^pe3_",                                                       "measurement",   "column",
  "Phytoplankton Taxon-Specific Pigments", "Biological",     "\\bpigment\\b|\\bfucoxanthin\\b|\\bperidinin\\b|\\bprasinoxanthin\\b|\\bviolaxanthin\\b",                      "measurement",   "column",
  "Picoplankton abundance",                "Biological",     "\\bprochlorococcus\\b|\\bsynechococcus\\b|\\bpicoeukaryote\\b",                                                 "measurement",   "column",
  "Bacterial abundance",                   "Biological",     "heterotrophic.*bacteria|bacterial.*abundance|bacteria.*ug",                                                     "measurement",   "column",
  "Phytoplankton Abundance",               "Biological",     "\\bdiatom\\b|\\bdinoflag\\b|auto.*euk|total.*phyto|phytoplankton.*abundance",                                  "measurement",   "column",
  "Seabird Species",                       "Biological",     "seabird.*species|sea.*bird.*species",                                                                           "observation",   "file_or_column",
  "Seabird Behavior",                      "Biological",     "seabird.*behavior",                                                                                              "observation",   "file_or_column",
  "Seabird Count",                         "Biological",     "seabird.*count|seabird.*abundance",                                                                             "observation",   "file_or_column",
  "Marine Mammal Species",                 "Biological",     "marine mammal.*species|mammal.*species",                                                                        "observation",   "file_or_column",
  "Marine Mammal Behavior",                "Biological",     "marine mammal.*behavior|mammal.*behavior",                                                                      "observation",   "file_or_column",
  "Marine Mammal Abundance",               "Biological",     "marine mammal.*count|marine mammal.*abundance|mammal.*count",                                                   "observation",   "file_or_column",
  "Fish Larvae Counts",                    "Biological",     "fish larvae.*count|larvae.*count|\\bichthyoplankton\\b|rockfish.*recruit",                                      "observation",   "file_or_column",
  "Fish Larvae Sizes",                     "Biological",     "fish larvae.*size|larvae.*size|larval.*length",                                                                  "observation",   "file_or_column",
  "Fish Egg Counts",                       "Biological",     "fish egg.*count|\\begg.*count\\b|continuous.*fish.*egg",                                                        "observation",   "file_or_column",
  "Fish Egg Stages",                       "Biological",     "fish egg.*stage|\\begg.*stage\\b",                                                                              "observation",   "file_or_column",
  "Fish Abundance and Distribution",       "Biological",     "scientific name|\\btaxon\\b|\\bspecies\\b|\\bspecimen\\b",                                                      "observation",   "file_or_column",
  "Fish Size Structure",                   "Biological",     "standard length|total length|fork length",                                                                       "observation",   "file_or_column",
  "Fish Biomass",                          "Biological",     "^weight$|^wet_weight$|^dry_weight$|\\bbiomass\\b",                                                              "observation",   "file_or_column",
  "Fish Population Structure",             "Biological",     "^sex$|^gender$|^male$|^female$",                                                                                "observation",   "file_or_column",
  "Zooplankton Volume",                    "Biological",     "zooplankton.*volume|\\bdisplacement.*volume\\b",                                                                "observation",   "file_or_column",
  "Zooplankton Abundance",                 "Biological",     "zooplankton.*abundance|zooplankton.*count",                                                                      "observation",   "file_or_column",
  "Krill (Euphausiid) Abundance",          "Biological",     "\\bkrill\\b|\\beuphausii",                                                                                     "observation",   "file_or_column",
  "Invertebrate Abundance",                "Biological",     "invertebrate.*abundance|invertebrate.*count|\\bbycatch\\b",                                                     "observation",   "file_or_column"
)

observation_context_rules <- tribble(
  ~context_name,               ~file_pattern,                                                          ~column_pattern,                                          ~standard_parameter,                ~eov_group,
  "seabird_species",           "\\bseabird\\b|sea.*bird",                                              "species|common_name|scientific_name|taxon",             "Seabird Species",                  "Biological",
  "seabird_behavior",          "\\bseabird\\b|sea.*bird",                                              "behavior|behaviour",                                    "Seabird Behavior",                 "Biological",
  "seabird_count",             "\\bseabird\\b|sea.*bird",                                              "\\bcount\\b|\\babundance\\b|\\bnumber\\b",              "Seabird Count",                    "Biological",
  "seabird_presence_generic",  "\\bseabird\\b|sea.*bird",                                              "\\btransect\\b|\\bcruise\\b|\\bsurvey\\b",              "Seabird Count",                    "Biological",
  "seabird_species_generic",   "\\bseabird\\b|sea.*bird",                                              "\\btransect\\b|\\bcruise\\b|\\bsurvey\\b",              "Seabird Species",                  "Biological",
  "seabird_behavior_generic",  "\\bseabird\\b|sea.*bird",                                              "\\btransect\\b|\\bcruise\\b|\\bsurvey\\b",              "Seabird Behavior",                 "Biological",
  "edna_genomics",             "\\bedna\\b|\\bncog\\b|\\b16s\\b|\\b18s\\b|\\bamplicon\\b",            "\\bdate\\b|\\blat\\b|\\blon\\b|\\bsample_id\\b",        "Microbial Genomics (mRNA)",        "Biological",
  "mm_species",                "marine.?mammal|\\bcetacean\\b|\\bwhale\\b|\\bdolphin\\b",              "species|common_name|scientific_name|\\btaxon\\b",       "Marine Mammal Species",            "Biological",
  "mm_behavior",               "marine.?mammal|\\bcetacean\\b|\\bwhale\\b|\\bdolphin\\b",              "behavior|behaviour",                                    "Marine Mammal Behavior",           "Biological",
  "mm_abundance",              "marine.?mammal|\\bcetacean\\b|\\bwhale\\b|\\bdolphin\\b",              "\\bcount\\b|\\babundance\\b|\\bnumber\\b|\\bpod\\b",    "Marine Mammal Abundance",          "Biological",
  "mm_presence_generic",       "marine.?mammal|\\bcetacean\\b|\\bwhale\\b|\\bdolphin\\b",              "\\btransect\\b|\\bcruise\\b|\\bsurvey\\b",              "Marine Mammal Abundance",          "Biological",
  "mm_species_generic",        "marine.?mammal|\\bcetacean\\b|\\bwhale\\b|\\bdolphin\\b",              "\\btransect\\b|\\bcruise\\b|\\bsurvey\\b",              "Marine Mammal Species",            "Biological",
  "larvae_count",              "fish.*larvae|\\blarval\\b|\\bichthyoplankton\\b|rockfish.*recruit",    "\\bcount\\b|\\babundance\\b|\\bnumber\\b",              "Fish Larvae Counts",               "Biological",
  "larvae_size",               "fish.*larvae|\\blarval\\b|\\bichthyoplankton\\b",                      "\\bsize\\b|\\blength\\b|\\bfrequency\\b",               "Fish Larvae Sizes",                "Biological",
  "egg_count",                 "fish.*egg|continuous.*fish.*egg",                                       "\\bcount\\b|\\babundance\\b|\\bnumber\\b",              "Fish Egg Counts",                  "Biological",
  "egg_stage",                 "fish.*egg|continuous.*fish.*egg",                                       "\\bstage\\b",                                           "Fish Egg Stages",                  "Biological",
  "zoo_volume",                "\\bzooplankton\\b",                                                     "\\bvolume\\b|displacement",                             "Zooplankton Volume",               "Biological",
  "zoo_abundance",             "\\bzooplankton\\b",                                                     "\\bcount\\b|\\babundance\\b|\\bnumber\\b",              "Zooplankton Abundance",            "Biological",
  "invert_size",               "\\binvertebrate\\b|\\binvertabrate\\b",                                 "\\bsize\\b|\\blength\\b",                               "Invertebrate Size",                "Biological",
  "invert_abund",              "\\binvertebrate\\b|\\binvertabrate\\b|\\bbycatch\\b|\\btrawl\\b",       "\\bcount\\b|\\babundance\\b|\\bnumber\\b",              "Invertebrate Abundance",           "Biological",
  "krill_abund",               "\\bkrill\\b|\\beuphausii",                                             "\\bcount\\b|\\babundance\\b|\\bnumber\\b",              "Krill (Euphausiid) Abundance",     "Biological",
  "fish_specimen_presence",    "\\bcps\\b|\\bfish\\b|\\btrawl\\b|\\bhaul\\b|\\bspecimen\\b",           "scientific_name|\\btaxon\\b|\\bspecies\\b",             "Fish Abundance and Distribution",  "Biological",
  "fish_specimen_size",        "\\bcps\\b|\\bfish\\b|\\btrawl\\b|\\bhaul\\b|\\bspecimen\\b",           "standard_length|fork_length|total_length",              "Fish Size Structure",              "Biological",
  "fish_specimen_biomass",     "\\bcps\\b|\\bfish\\b|\\btrawl\\b|\\bhaul\\b|\\bspecimen\\b",           "^weight$|^wet_weight$|^dry_weight$",                    "Fish Biomass",                     "Biological",
  "fish_specimen_sex",         "\\bcps\\b|\\bfish\\b|\\btrawl\\b|\\bhaul\\b|\\bspecimen\\b",           "^sex$|^gender$",                                        "Fish Population Structure (Sex)",  "Biological"
)

filename_parameter_dictionary <- tibble(
  pattern = c(
    "\\btemp\\b|\\btemperature\\b", "\\bsalinity\\b|\\bsalin",
    "\\boxygen\\b|\\bo2\\b|dissolved.*oxygen",
    "\\bchl\\b|\\bchlorophyll\\b|\\bphaeo\\b",
    "\\bnitrate\\b|\\bnitrite\\b|\\bammonium\\b|\\bphosphate\\b|\\bsilicate\\b|\\bnutrients\\b",
    "\\bdic\\b|dissolved_inorganic_carbon", "\\balkalinity\\b|\\btalk\\b",
    "(^|[_\\-])ph([_\\-]|$)|seawater_ph", "transmissiv",
    "radiative.*flux|\\birradiance\\b",
    "\\bproductivity\\b|\\b14c\\b|primary.*prod",
    "\\bfish.*egg\\b|\\begg.*count\\b", "\\blarva|\\blarval\\b|\\bichthyoplankton\\b",
    "\\bmammal\\b|\\bcetacean\\b|\\bwhale\\b", "\\bseabird\\b|sea.*bird",
    "\\bacoustic\\b|bioacoustic", "\\bcdom\\b",
    "\\bphycoerythrin\\b|(^|[_\\-])pe([_\\-]|$)",
    "\\bmicrob|\\bbacteria\\b|\\b16s\\b|\\b18s\\b|\\bmrna\\b",
    "\\binvertebrate\\b|\\bspecimen\\b|\\bbycatch\\b|\\btrawl\\b"
  ),
  filename_parameter = c(
    "Temperature", "Salinity", "Dissolved Oxygen", "Chlorophyll-a",
    "Nutrients", "Dissolved Inorganic Carbon (DIC)", "Total Alkalinity",
    "pH", "Transmissivity", "Radiative Flux",
    "Primary Production (14C uptake)", "Fish Egg Counts", "Fish Larvae Counts",
    "Marine Mammal Abundance", "Seabird Count", "Whale Acoustic",
    "CDOM Fluorescence", "Phycoerythrin (PE) Fluorescence",
    "Microbial Genomics (mRNA)", "Invertebrate Abundance"
  )
)

detect_parameter_from_filename <- function(file_name, dictionary_tbl) {
  matches <- dictionary_tbl %>%
    filter(str_detect(file_name, regex(pattern, ignore_case = TRUE)))
  if (nrow(matches) == 0) return(tibble(filename_parameter = NA_character_))
  matches %>% distinct(filename_parameter)
}

################
# FIND CSV FILES
################

csv_files <- list.files(program_folder, pattern = "\\.csv$",
                         full.names = TRUE, recursive = TRUE) %>%
  .[!str_detect(., regex("(/|\\\\)(output[^/\\\\]*)(/|\\\\)", ignore_case = TRUE))]

if (length(csv_files) == 0)        stop("No CSV files found in program_folder.")
if (!file.exists(ca_boundary_path)) stop("Boundary shapefile not found. Check ca_boundary_path.")

cat("CSV files found:", length(csv_files), "\n")
print(basename(csv_files))

##############################
# STEP 1-2: READ + CLEAN FILES
##############################

cat("\n=== STEP 1-2: Read and clean ===\n")

read_and_clean_file <- function(file, default_program) {
  message("Reading: ", basename(file))

  df_raw <- read_csv(file, show_col_types = FALSE,
                     col_types = cols(.default = col_character()))

  raw_names <- names(df_raw) %>%
    str_replace_all("\u00b5", "u") %>%
    str_replace_all("[^[:ascii:]]", "u") %>%
    str_replace_all("/", "_per_")
  names(df_raw) <- raw_names

  df <- df_raw %>%
    clean_names() %>% clean_character_cols() %>%
    mutate(source_row_id = row_number())

  lat_col     <- validate_coord_col(df, detect_col_priority(df, lat_priority_patterns),  lat_priority_patterns)
  lon_col     <- validate_coord_col(df, detect_col_priority(df, lon_priority_patterns),  lon_priority_patterns)
  year_col    <- detect_col(df, year_patterns)
  date_col    <- detect_col(df, date_patterns)
  program_col <- detect_col(df, program_patterns)
  depth_col   <- detect_col(df, depth_patterns)

  # Reject depth columns with implausible values
  if (!is.na(depth_col) && depth_col %in% names(df)) {
    raw_check  <- suppressWarnings(as.numeric(df[[depth_col]]))
    valid_vals <- raw_check[!is.na(raw_check)]
    if (length(valid_vals) == 0 || median(valid_vals) < 0 || median(valid_vals) > 12000)
      depth_col <- NA_character_
  }

  coord_role <- case_when(
    !is.na(lat_col) & str_detect(lat_col, regex("mid",   ignore_case = TRUE)) ~ "mid",
    !is.na(lat_col) & str_detect(lat_col, regex("start", ignore_case = TRUE)) ~ "start",
    !is.na(lat_col) & str_detect(lat_col, regex("stop",  ignore_case = TRUE)) ~ "stop",
    TRUE ~ "standard"
  )

  df %>%
    mutate(
      source_file          = basename(file),
      source_path          = file,
      file_stub            = get_file_stub(file),
      program              = if (!is.na(program_col) && program_col %in% names(df))
                               as.character(.data[[program_col]])
                             else rep(as.character(default_program), nrow(df)),
      program              = if_else(is.na(program) | str_trim(program) == "",
                                     as.character(default_program), as.character(program)),
      detected_lat_col     = lat_col,
      detected_lon_col     = lon_col,
      detected_year_col    = year_col,
      detected_date_col    = date_col,
      detected_program_col = program_col,
      detected_depth_col   = depth_col,
      detected_coord_role  = coord_role,
      year_detected        = extract_year(df, year_col, date_col),
      latitude_std         = if (!is.na(lat_col) && lat_col %in% names(df))
                               parse_coord_vector(df[[lat_col]], FALSE, basename(file), lat_col)
                             else rep(NA_real_, nrow(df)),
      longitude_std        = if (!is.na(lon_col) && lon_col %in% names(df))
                               parse_coord_vector(df[[lon_col]], TRUE, basename(file), lon_col)
                             else rep(NA_real_, nrow(df)),
      depth_std            = if (!is.na(depth_col) && depth_col %in% names(df)) {
                               raw <- suppressWarnings(as.numeric(df[[depth_col]]))
                               if_else(raw >= 0 & raw <= 12000, raw, NA_real_)
                             } else rep(NA_real_, nrow(df)),
      sample_point_key     = paste0(get_file_stub(file), "_ROW_", source_row_id)
    )
}

all_raw_df <- map_dfr(csv_files, read_and_clean_file, default_program = program_name)
cat("Total rows read:", nrow(all_raw_df), "\n")

detection_summary <- all_raw_df %>%
  distinct(source_file, source_path, detected_lat_col, detected_lon_col,
           detected_year_col, detected_date_col, detected_program_col,
           detected_depth_col, detected_coord_role)

############################################
# STEP 3: YEAR FILTER
############################################

cat("\n=== STEP 3: Year filter (>= ", start_year, ") ===\n", sep = "")

all_filtered_df <- all_raw_df %>%
  filter(is.na(year_detected) | year_detected >= start_year)

all_filtered_coords_df_raw <- all_filtered_df %>%
  filter(!is.na(latitude_std), !is.na(longitude_std))

cat("Rows after year filter:", nrow(all_filtered_df),
    "| rows with coords:", nrow(all_filtered_coords_df_raw), "\n")

############################################
# STEP 4: OPTIONAL COASTAL CLIP
############################################

cat("\n=== STEP 4: Coastal clip ===\n")

ca_boundary_proj  <- st_read(ca_boundary_path, quiet = TRUE) %>%
  st_transform(3310) %>% st_union()
ca_coastal_buffer <- st_buffer(ca_boundary_proj, dist = buffer_meters)

if (nrow(all_filtered_coords_df_raw) == 0) {
  all_filtered_coords_df <- all_filtered_coords_df_raw
  clip_summary <- tibble(clip_applied = apply_coastal_clip,
                         rows_before_clip = 0L, rows_after_clip = 0L,
                         rows_removed_by_clip = 0L)
} else {
  coords_sf <- all_filtered_coords_df_raw %>%
    st_as_sf(coords = c("longitude_std", "latitude_std"), crs = 4326, remove = FALSE) %>%
    st_transform(3310)

  if (apply_coastal_clip) {
    keep              <- st_intersects(coords_sf, ca_coastal_buffer, sparse = FALSE)[, 1]
    all_filtered_coords_df <- st_drop_geometry(coords_sf[keep, ])
  } else {
    all_filtered_coords_df <- st_drop_geometry(coords_sf)
  }

  clip_summary <- tibble(
    clip_applied         = apply_coastal_clip,
    rows_before_clip     = nrow(all_filtered_coords_df_raw),
    rows_after_clip      = nrow(all_filtered_coords_df),
    rows_removed_by_clip = nrow(all_filtered_coords_df_raw) - nrow(all_filtered_coords_df)
  )
  cat("Rows after clip:", nrow(all_filtered_coords_df), "\n")
}

##############################################
# GEBCO BATHYMETRY EXTRACTION
# Extracts seafloor depth at each sample point
##############################################

cat("\n=== Extracting GEBCO bathymetry ===\n")

if (file.exists(gebco_raster_path) && nrow(all_filtered_coords_df) > 0) {
  gebco_raster <- terra::rast(gebco_raster_path)
  bathy_pts    <- all_filtered_coords_df %>%
    filter(!is.na(latitude_std), !is.na(longitude_std)) %>%
    select(sample_point_key, longitude_std, latitude_std) %>% distinct()
  pts_vect   <- terra::vect(bathy_pts, geom = c("longitude_std","latitude_std"), crs = "EPSG:4326")
  bathy_vals <- terra::extract(gebco_raster, pts_vect)
  bathy_pts$gebco_depth_m <- bathy_vals[, 2]
  all_filtered_coords_df  <- all_filtered_coords_df %>%
    left_join(bathy_pts %>% select(sample_point_key, gebco_depth_m), by = "sample_point_key")
  cat("Bathymetry extracted for", sum(!is.na(all_filtered_coords_df$gebco_depth_m)), "points\n")
} else {
  if (!file.exists(gebco_raster_path)) cat("WARNING: GEBCO raster not found — depth will be NA\n")
  all_filtered_coords_df$gebco_depth_m <- NA_real_
}

#####################################
# STEP 5: HEX GRID + ASSIGN HEX CELLS
#####################################

cat("\n=== STEP 5: Hex grid ===\n")

if (nrow(all_filtered_coords_df) == 0) {
  hex_grid_sf        <- st_sf(hex_id = character(), geometry = st_sfc(crs = 3310))
  points_with_hex_df <- all_filtered_coords_df %>%
    mutate(hex_id = NA_character_, station_key = NA_character_)
} else {
  points_sf <- all_filtered_coords_df %>%
    st_as_sf(coords = c("longitude_std","latitude_std"), crs = 4326, remove = FALSE) %>%
    st_transform(3310)

  grid_extent <- if (apply_coastal_clip) ca_coastal_buffer else
    st_buffer(st_as_sfc(st_bbox(points_sf)), dist = hex_cellsize_m * 2)

  hex_grid_sf <- st_make_grid(grid_extent, cellsize = hex_cellsize_m, square = FALSE) %>%
    st_sf() %>% mutate(hex_id = paste0("HEX_", row_number())) %>%
    st_set_crs(3310)

  idx       <- st_intersects(points_sf, hex_grid_sf, sparse = TRUE)
  hex_ids   <- vapply(idx, function(i) if (length(i) > 0) hex_grid_sf$hex_id[i[1]] else NA_character_, character(1))
  na_mask   <- is.na(hex_ids)
  if (any(na_mask))
    hex_ids[na_mask] <- hex_grid_sf$hex_id[st_nearest_feature(points_sf[na_mask,], hex_grid_sf)]

  points_with_hex_df <- st_drop_geometry(points_sf) %>%
    mutate(hex_id      = hex_ids,
           station_key = if_else(!is.na(hex_ids), paste(program, hex_ids, sep=" | "), NA_character_))

  hex_grid_sf <- hex_grid_sf %>% filter(hex_id %in% unique(na.omit(hex_ids)))
  cat("Hex cells used:", nrow(hex_grid_sf), "\n")
}

if (nrow(points_with_hex_df) == 0) stop("No rows remain after Step 5.")

#################################################################
# STEP 6: CLASSIFY GEOMETRY TYPE
# Determines whether each file represents
# fixed stations, transect tows, or underway continuous sampling.
#################################################################

cat("\n=== STEP 6: Classify geometry type ===\n")

dataset_classification <- points_with_hex_df %>%
  distinct(source_file, program) %>%
  mutate(
    src_lower     = str_to_lower(source_file),
    geometry_type = case_when(
      str_detect(src_lower, "underway|continuous|flow.?through|flowthrough") ~ "underway",
      str_detect(src_lower, "tow|trawl|transect|haul|catch|length.?frequency|recruitment|survey") ~ "line",
      TRUE ~ "point"
    ),
    classification_reason = case_when(
      geometry_type == "underway" ~ "Continuous underway sampling",
      geometry_type == "line"     ~ "Mobile transect/tow sampling",
      TRUE                        ~ "Fixed station"
    )
  ) %>%
  select(source_file, program, geometry_type, classification_reason)

points_with_hex_df <- points_with_hex_df %>%
  left_join(dataset_classification, by = c("source_file","program"))

############################################
# STEP 7: PARAMETER DETECTION
# Detects which oceanographic parameters are
# present in each file via column names,
# file names, and contextual rules.
############################################

cat("\n=== STEP 7: Parameter detection ===\n")

empty_param_schema <- tibble(
  source_file = character(), raw_parameter_name = character(),
  detected_context = character(), detected_measurement_type = character(),
  standard_parameter = character(), eov_group = character(), mapping_status = character()
)

WORKFLOW_COL_REGEX <- regex(
  paste(c("^detected_","^source_","^file_stub$","^sample_point_key$","^station_key$",
          "^hex_id$","^\\.row_id_temp$","^\\.join_key$","^geometry_type$",
          "^classification_reason$","^program$","^year_detected$",
          "^latitude_std$","^longitude_std$","^depth_std$","^activity_status$","^last_year$"),
        collapse = "|"), ignore_case = TRUE)

# Detects parameters via column name matching against parameter_dictionary
detect_measurement_parameter_columns <- function(df, source_file, dictionary) {
  if (nrow(df) == 0) return(empty_param_schema)
  df_names <- names(df)[!str_detect(names(df), WORKFLOW_COL_REGEX)]
  cleaned  <- clean_text_value(df_names)
  map2_dfr(df_names, cleaned, function(orig, clean) {
    if (is_probably_metadata_column(clean)) return(tibble())
    if (orig %in% names(df) && count_real_values(df[[orig]]) == 0) return(tibble())
    hits <- dictionary %>%
      filter(detection_type == "measurement",
             str_detect(clean, regex(raw_pattern, ignore_case = TRUE)))
    if (nrow(hits) == 0) return(tibble())
    hits %>% transmute(source_file = source_file, raw_parameter_name = orig,
                       detected_context = "column_name", detected_measurement_type = "measurement_column",
                       standard_parameter, eov_group, mapping_status = "matched_from_column")
  }) %>% distinct()
}

# Detects parameters when both file name AND column names match a rule
detect_observation_parameters_from_file_and_columns <- function(df, source_file, rules_table) {
  if (nrow(df) == 0) return(empty_param_schema)
  file_clean <- clean_text_value(source_file)
  col_clean  <- clean_text_value(names(df)[!str_detect(names(df), WORKFLOW_COL_REGEX)])
  pmap_dfr(rules_table, function(context_name, file_pattern, column_pattern, standard_parameter, eov_group) {
    if (str_detect(file_clean, regex(file_pattern, ignore_case = TRUE)) &&
        any(str_detect(col_clean, regex(column_pattern, ignore_case = TRUE)))) {
      tibble(source_file = source_file, raw_parameter_name = paste0("[context] ", context_name),
             detected_context = "file_and_column_context", detected_measurement_type = "observation_context",
             standard_parameter = standard_parameter, eov_group = eov_group,
             mapping_status = "matched_from_file_and_columns")
    } else tibble()
  }) %>% distinct()
}

# Detects specimen-level biological parameters from files with taxon + supporting columns
detect_fish_specimen_observation_parameters <- function(df, source_file) {
  if (nrow(df) == 0) return(empty_param_schema)
  col_clean <- clean_text_value(names(df)[!str_detect(names(df), WORKFLOW_COL_REGEX)])
  has_taxon <- any(str_detect(col_clean, regex("scientific name|taxon|species", ignore_case = TRUE)))
  has_size  <- any(str_detect(col_clean, regex("length|size", ignore_case = TRUE)))
  has_mass  <- any(str_detect(col_clean, regex("weight|biomass", ignore_case = TRUE)))
  has_sex   <- any(str_detect(col_clean, regex("^sex$|gender", ignore_case = TRUE)))
  if (!has_taxon || !any(c(has_size, has_mass, has_sex))) return(empty_param_schema)
  out <- tibble(
    source_file = source_file,
    raw_parameter_name = c("[specimen] taxon","[specimen] size","[specimen] biomass","[specimen] sex"),
    detected_context = "specimen_context", detected_measurement_type = "observation_context",
    standard_parameter = c("Organism Abundance and Distribution","Organism Size Structure",
                           "Organism Biomass","Organism Population Structure"),
    eov_group = "Biological", mapping_status = "matched_from_specimen_context"
  )
  out[c(TRUE, has_size, has_mass, has_sex), , drop = FALSE] %>% distinct()
}

expand_file_level_parameters_to_rows <- function(df, file_level_hits) {
  if (nrow(df) == 0 || nrow(file_level_hits) == 0) return(tibble())
  src <- unique(df$source_file)[1]
  df %>% mutate(.row_id_temp = row_number(), .join_key = 1L) %>%
    left_join(file_level_hits %>%
                select(raw_parameter_name, detected_context, detected_measurement_type,
                       standard_parameter, eov_group, mapping_status) %>%
                distinct() %>% mutate(.join_key = 1L),
              by = ".join_key", relationship = "many-to-many") %>%
    select(-.join_key) %>%
    mutate(source_file = src)
}

split_file_list              <- split(points_with_hex_df, points_with_hex_df$source_file)
row_level_parameter_hits_list <- list()
filename_parameter_hits_list  <- list()
ambiguous_filename_cases_list <- list()
likely_parameter_columns      <- list()

for (this_file in names(split_file_list)) {
  df_file <- split_file_list[[this_file]]
  cat("  Processing:", basename(this_file), "- rows:", nrow(df_file), "\n")

  meas_hits     <- detect_measurement_parameter_columns(df_file, this_file, parameter_dictionary)
  obs_hits      <- detect_observation_parameters_from_file_and_columns(df_file, this_file, observation_context_rules)
  specimen_hits <- detect_fish_specimen_observation_parameters(df_file, this_file)

  fn_raw  <- detect_parameter_from_filename(basename(this_file), filename_parameter_dictionary)
  fn_hits <- if (nrow(fn_raw) > 0 && !all(is.na(fn_raw$filename_parameter)))
    fn_raw %>% filter(!is.na(filename_parameter)) %>%
      transmute(source_file = this_file, raw_parameter_name = paste0("[filename] ", filename_parameter),
                detected_context = "filename_only", detected_measurement_type = "filename_inference",
                standard_parameter = filename_parameter, eov_group = NA_character_,
                mapping_status = "filename_only_candidate")
  else empty_param_schema

  if (nrow(meas_hits) > 0) likely_parameter_columns[[length(likely_parameter_columns)+1]] <- meas_hits
  if (nrow(fn_hits) > 0)  filename_parameter_hits_list[[length(filename_parameter_hits_list)+1]] <- fn_hits

  trusted <- bind_rows(meas_hits, obs_hits, specimen_hits) %>% distinct()

  if (nrow(trusted) == 0 && nrow(fn_hits) > 1)
    ambiguous_filename_cases_list[[length(ambiguous_filename_cases_list)+1]] <-
      fn_hits %>% mutate(mapping_status = "ambiguous_filename_only_review")

  if (nrow(trusted) > 0) {
    expanded <- expand_file_level_parameters_to_rows(df_file, trusted)
    if (nrow(expanded) > 0) row_level_parameter_hits_list[[length(row_level_parameter_hits_list)+1]] <- expanded
  }
}

likely_parameter_columns <- bind_rows(c(list(empty_param_schema), likely_parameter_columns)) %>% distinct()
row_level_parameter_hits <- bind_rows(c(list(empty_param_schema), row_level_parameter_hits_list)) %>% distinct()
filename_parameter_hits  <- bind_rows(c(list(empty_param_schema), filename_parameter_hits_list)) %>% distinct()
ambiguous_filename_cases <- bind_rows(c(list(empty_param_schema), ambiguous_filename_cases_list)) %>% distinct()

cat("Parameter rows detected:", nrow(row_level_parameter_hits), "\n")

###################################
# STEP 8: PARAMETER STANDARDIZATION
###################################

cat("\n=== STEP 8: Standardize parameters ===\n")

enforce_param_schema <- function(df) {
  for (col in c("source_file","raw_parameter_name","detected_context",
                "detected_measurement_type","standard_parameter","eov_group","mapping_status"))
    if (!col %in% names(df)) df[[col]] <- NA_character_
  df
}

standardized_hits <- row_level_parameter_hits %>%
  mutate(across(c(standard_parameter, eov_group, mapping_status), as.character))

standardized_hits_mapped <- standardized_hits %>%
  filter(!is.na(standard_parameter), !is.na(eov_group))

parameter_review_table <- bind_rows(
  enforce_param_schema(standardized_hits) %>%
    select(source_file, raw_parameter_name, detected_context,
           detected_measurement_type, standard_parameter, eov_group, mapping_status) %>%
    distinct() %>% mutate(review_type = "assigned"),
  enforce_param_schema(ambiguous_filename_cases) %>%
    select(source_file, raw_parameter_name, detected_context,
           detected_measurement_type, standard_parameter, eov_group, mapping_status) %>%
    distinct() %>% mutate(review_type = "needs_review"),
  enforce_param_schema(filename_parameter_hits) %>%
    filter(!source_file %in% standardized_hits$source_file) %>%
    select(source_file, raw_parameter_name, detected_context,
           detected_measurement_type, standard_parameter, eov_group, mapping_status) %>%
    distinct() %>% mutate(review_type = "filename_only_candidate")
) %>% distinct()

if (nrow(standardized_hits_mapped) == 0)
  stop("No parameters detected for: ", program_name,
       "\nCheck: coords, year filter, column names vs parameter_dictionary.")

######################################
# STEP 8B: PRELIMINARY ACTIVITY STATUS
######################################

latest_year_by_file <- points_with_hex_df %>%
  filter(!is.na(year_detected)) %>%
  group_by(source_file) %>%
  summarise(last_year = max(year_detected, na.rm = TRUE), .groups = "drop")

standardized_hits_prelim_active <- standardized_hits_mapped %>%
  left_join(latest_year_by_file, by = "source_file") %>%
  mutate(activity_status = case_when(
    !is.na(last_year) & last_year >= active_cutoff_year ~ "Active",
    !is.na(last_year) & last_year < active_cutoff_year  ~ "Inactive",
    TRUE ~ "Unknown"
  ))

##################################
# STEP 8C: PER-PARAMETER FREQUENCY
##################################

cat("\n=== STEP 8C: Per-parameter frequency ===\n")

parameter_frequency_table <- compute_parameter_frequency(
  param_hits_df      = standardized_hits_prelim_active,
  raw_df             = all_raw_df,
  fallback_frequency = sampling_frequency
) %>%
  left_join(
    attr_param_lookup %>% filter(acronym == program_name) %>%
      select(standard_parameter, attr_frequency, attr_platform),
    by = "standard_parameter"
  ) %>%
  mutate(
    frequency_label = case_when(
      frequency_source %in% c("computed_from_dates","estimated_from_seasonal_labels","estimated_from_year_counts") ~ frequency_label,
      !is.na(attr_frequency) & attr_frequency != "" ~ attr_frequency,
      TRUE ~ frequency_label
    ),
    frequency_source = case_when(
      frequency_source %in% c("computed_from_dates","estimated_from_seasonal_labels","estimated_from_year_counts") ~ frequency_source,
      !is.na(attr_frequency) & attr_frequency != "" ~ "attribute_table",
      TRUE ~ frequency_source
    )
  )

print(parameter_frequency_table %>% select(standard_parameter, frequency_label))

################################
# STEP 9: HEX PARAMETER PRESENCE
################################

cat("\n=== STEP 9: Hex parameter presence ===\n")

program_hex_parameter_table <- standardized_hits_mapped %>%
  group_by(program, hex_id, station_key, standard_parameter, eov_group) %>%
  summarise(
    parameter_present    = 1L,
    first_year           = if (all(is.na(year_detected))) NA_integer_ else min(year_detected, na.rm = TRUE),
    last_year            = if (all(is.na(year_detected))) NA_integer_ else max(year_detected, na.rm = TRUE),
    min_depth            = if (all(is.na(depth_std))) NA_real_ else min(depth_std, na.rm = TRUE),
    max_depth            = if (all(is.na(depth_std))) NA_real_ else max(depth_std, na.rm = TRUE),
    geometry_types_found = collapse_unique_text(geometry_type),
    source_files         = collapse_unique_text(source_file),
    .groups = "drop"
  ) %>%
  mutate(program_hex_key = paste(program, hex_id, sep = " | "))

############################################
# STEP 9B: ACTIVITY STATUS FILTER
# Removes hex cells where the most recent
# data predates active_cutoff_year.
############################################

active_program_hex <- standardized_hits_mapped %>%
  group_by(program, hex_id) %>%
  summarise(
    last_year = if (all(is.na(year_detected))) NA_integer_ else max(year_detected, na.rm = TRUE),
    activity_status = case_when(
      !is.na(last_year) & last_year >= active_cutoff_year ~ "active_since_2015",
      !is.na(last_year) & last_year < active_cutoff_year  ~ "inactive_pre_2015",
      TRUE ~ "unknown_year"
    ),
    .groups = "drop"
  )

program_hex_parameter_table_active <- program_hex_parameter_table %>%
  left_join(active_program_hex, by = c("program","hex_id")) %>%
  filter(activity_status != "inactive_pre_2015") %>%
  left_join(parameter_frequency_table %>% select(standard_parameter, frequency_label, frequency_source),
            by = "standard_parameter")

standardized_hits_active <- standardized_hits_mapped %>%
  left_join(active_program_hex, by = c("program","hex_id")) %>%
  filter(activity_status != "inactive_pre_2015")

cat("Active hex rows:", nrow(program_hex_parameter_table_active), "\n")

############################################
# STEP 9C: PER-HEX DEPTH RANGES BY PARAMETER
############################################

hex_parameter_depths <- program_hex_parameter_table_active %>%
  group_by(hex_id, standard_parameter) %>%
  summarise(
    param_min = if (all(is.na(min_depth))) NA_real_ else min(min_depth, na.rm = TRUE),
    param_max = if (all(is.na(max_depth))) NA_real_ else max(max_depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(param_min) | !is.na(param_max)) %>%
  mutate(depth_label = case_when(
    !is.na(param_min) & !is.na(param_max) & round(param_min) != round(param_max) ~
      paste0(standard_parameter, ": ", round(param_min), "\u2013", round(param_max), " m"),
    !is.na(param_min) ~ paste0(standard_parameter, ": ", round(param_min), " m"),
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(depth_label)) %>%
  group_by(hex_id) %>%
  summarise(parameter_depths = collapse_unique_text(depth_label), .groups = "drop")

############################################
# STEP 10: SUMMARIES + WIDE PARAMETER FLAGS
############################################

cat("\n=== STEP 10: Summaries ===\n")

program_hex_parameter_wide <- program_hex_parameter_table_active %>%
  select(program, hex_id, standard_parameter, parameter_present) %>%
  distinct() %>%
  mutate(yes_no_field = paste0("param_", str_replace_all(str_to_lower(standard_parameter), "[^a-z0-9]+", "_"))) %>%
  select(-standard_parameter) %>%
  pivot_wider(names_from = yes_no_field, values_from = parameter_present, values_fill = 0)

program_hex_summary <- standardized_hits_active %>%
  group_by(program, hex_id, activity_status) %>%
  summarise(
    first_year           = if (all(is.na(year_detected))) NA_integer_ else min(year_detected, na.rm = TRUE),
    last_year            = if (all(is.na(year_detected))) NA_integer_ else max(year_detected, na.rm = TRUE),
    year_count           = if (all(is.na(year_detected))) NA_integer_ else
      as.integer(max(year_detected, na.rm = TRUE) - min(year_detected, na.rm = TRUE) + 1L),
    years_with_data      = n_distinct(year_detected, na.rm = TRUE),
    min_depth            = if (all(is.na(depth_std))) NA_real_ else min(depth_std, na.rm = TRUE),
    max_depth            = if (all(is.na(depth_std))) NA_real_ else max(depth_std, na.rm = TRUE),
    parameter_count      = n_distinct(standard_parameter),
    eov_group_count      = n_distinct(eov_group),
    parameters           = collapse_unique_text(standard_parameter),
    eov_groups           = collapse_unique_text(eov_group),
    geometry_types_found = collapse_unique_text(geometry_type),
    source_files         = collapse_unique_text(source_file),
    gebco_mean_depth     = if (all(is.na(gebco_depth_m))) NA_real_ else mean(gebco_depth_m, na.rm = TRUE),
    gebco_min_depth      = if (all(is.na(gebco_depth_m))) NA_real_ else min(gebco_depth_m,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(station_key     = paste(program, hex_id, sep = " | "),
         program_hex_key = paste(program, hex_id, sep = " | "))

# Hex centroids in WGS84 for popup coordinates
hex_centroids_ll <- st_centroid(st_geometry(hex_grid_sf)) %>%
  st_as_sf() %>% st_set_crs(st_crs(hex_grid_sf)) %>%
  mutate(hex_id = hex_grid_sf$hex_id) %>%
  st_transform(4326) %>%
  mutate(centroid_longitude = st_coordinates(.)[,1],
         centroid_latitude  = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  select(hex_id, centroid_longitude, centroid_latitude)

program_hex_inventory <- program_hex_summary %>%
  left_join(program_hex_parameter_wide, by = c("program","hex_id")) %>%
  left_join(hex_centroids_ll, by = "hex_id") %>%
  arrange(program, hex_id)

hex_point_count <- standardized_hits_active %>%
  filter(!is.na(latitude_std), !is.na(longitude_std)) %>%
  group_by(hex_id) %>%
  summarise(sample_location_count = n_distinct(paste(round(latitude_std,5), round(longitude_std,5))),
            .groups = "drop")

hex_inventory_for_map <- program_hex_inventory %>%
  group_by(hex_id) %>%
  summarise(
    programs             = collapse_unique_text(program),
    activity_status      = collapse_unique_text(activity_status),
    geometry_types_found = collapse_unique_text(geometry_types_found),
    first_year           = suppressWarnings(ifelse(all(is.na(first_year)), NA_integer_, min(first_year, na.rm = TRUE))),
    last_year            = suppressWarnings(ifelse(all(is.na(last_year)),  NA_integer_, max(last_year,  na.rm = TRUE))),
    year_count           = suppressWarnings(ifelse(all(is.na(year_count)), NA_integer_, max(year_count, na.rm = TRUE))),
    parameter_count      = max(parameter_count, na.rm = TRUE),
    eov_group_count      = max(eov_group_count, na.rm = TRUE),
    gebco_mean_depth     = suppressWarnings(ifelse(all(is.na(gebco_mean_depth)), NA_real_, mean(gebco_mean_depth, na.rm = TRUE))),
    parameters           = collapse_unique_text(parameters),
    eov_groups           = collapse_unique_text(eov_groups),
    source_files         = collapse_unique_text(source_files),
    centroid_longitude   = first(na.omit(centroid_longitude)),
    centroid_latitude    = first(na.omit(centroid_latitude)),
    min_depth            = suppressWarnings(ifelse(all(is.na(min_depth)), NA_real_, min(min_depth, na.rm = TRUE))),
    max_depth            = suppressWarnings(ifelse(all(is.na(max_depth)), NA_real_, max(max_depth, na.rm = TRUE))),
    across(starts_with("param_"), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(hex_point_count,      by = "hex_id") %>%
  left_join(hex_parameter_depths, by = "hex_id")

############################################
# STEP 11: OVERLAP TABLE + SAMPLED POINTS
############################################

hex_overlap_table <- standardized_hits_active %>%
  group_by(hex_id) %>%
  summarise(
    program_count   = n_distinct(program),
    programs        = collapse_unique_text(program),
    parameter_count = n_distinct(standard_parameter),
    parameters      = collapse_unique_text(standard_parameter),
    eov_groups      = collapse_unique_text(eov_group),
    first_year      = if (all(is.na(year_detected))) NA_integer_ else min(year_detected, na.rm = TRUE),
    last_year       = if (all(is.na(year_detected))) NA_integer_ else max(year_detected, na.rm = TRUE),
    overlap_flag    = ifelse(n_distinct(program) > 1, 1L, 0L),
    .groups = "drop"
  ) %>%
  left_join(hex_centroids_ll, by = "hex_id")

all_active_sampled_points <- standardized_hits_active %>%
  transmute(source_file, source_row_id, sample_point_key, station_key, program,
            year = year_detected, latitude = latitude_std, longitude = longitude_std,
            depth = depth_std, hex_id, geometry_type,
            raw_parameter_name, standard_parameter, eov_group,
            detected_context, detected_measurement_type) %>%
  arrange(program, hex_id, year, source_file)

############################################
# TRANSECT EXPORT
# Produces transects.csv for the HTML map.
# Filters to corridors visited >= 5 years
# and deduplicates to one line per corridor.
############################################

cat("\n=== Transect export ===\n")

mobile_sampling_points <- all_active_sampled_points %>% filter(geometry_type == "line")

# Also catch point-classified files that have start/stop columns
files_with_start_stop <- all_raw_df %>%
  distinct(source_file) %>%
  filter(purrr::map_lgl(source_file, function(sf) {
    nm <- names(all_raw_df %>% filter(source_file == sf) %>%
                  select(where(~ any(!is.na(.) & . != ""))))
    any(str_detect(nm, regex("lat.*(start|begin)|start.*lat", ignore_case = TRUE))) &&
      any(str_detect(nm, regex("lat.*(stop|end)|stop.*lat",   ignore_case = TRUE)))
  })) %>%
  pull(source_file)

transect_source_rows <- bind_rows(
  mobile_sampling_points,
  standardized_hits_active %>%
    filter(source_file %in% files_with_start_stop,
           !source_file %in% mobile_sampling_points$source_file) %>%
    transmute(source_file, source_row_id, sample_point_key, station_key, program,
              year = year_detected, latitude = latitude_std, longitude = longitude_std,
              depth = depth_std, hex_id, geometry_type = "line",
              raw_parameter_name, standard_parameter, eov_group,
              detected_context, detected_measurement_type)
) %>% distinct(source_file, source_row_id, .keep_all = TRUE)

if (nrow(transect_source_rows) > 0) {

  row_param_lookup <- standardized_hits_active %>%
    group_by(source_file, source_row_id) %>%
    summarise(Parameters = collapse_unique_text(standard_parameter),
              `EOV Groups` = collapse_unique_text(eov_group), .groups = "drop")

  file_param_lookup <- standardized_hits_active %>%
    group_by(program) %>%
    summarise(file_parameters = collapse_unique_text(standard_parameter),
              file_eov_groups = collapse_unique_text(eov_group), .groups = "drop")

  # Build transects from files that have explicit start/stop coordinate columns
  expanded <- purrr::map_dfr(files_with_start_stop, function(sf) {
    raw_slice <- all_raw_df %>% filter(source_file == sf)
    nm        <- names(raw_slice)
    lat_s <- nm[str_detect(nm, regex("lat.*(start|begin)|start.*lat", ignore_case = TRUE))][1]
    lon_s <- nm[str_detect(nm, regex("lon.*(start|begin)|start.*lon", ignore_case = TRUE))][1]
    lat_m <- nm[str_detect(nm, regex("lat.*mid|mid.*lat",              ignore_case = TRUE))][1]
    lon_m <- nm[str_detect(nm, regex("lon.*mid|mid.*lon",              ignore_case = TRUE))][1]
    lat_e <- nm[str_detect(nm, regex("lat.*(stop|end)|stop.*lat",      ignore_case = TRUE))][1]
    lon_e <- nm[str_detect(nm, regex("lon.*(stop|end)|stop.*lon",      ignore_case = TRUE))][1]
    if (is.na(lat_s) || is.na(lat_e)) return(tibble())
    d_col   <- unique(na.omit(raw_slice$detected_date_col))[1]
    dep_col <- unique(na.omit(raw_slice$detected_depth_col))[1]
    svy_col <- nm[str_detect(nm, regex("^svy$|^survey$|^cruise$", ignore_case = TRUE))][1]
    raw_slice %>%
      left_join(row_param_lookup, by = c("source_file","source_row_id")) %>%
      transmute(
        Date              = if (!is.na(d_col)   && d_col   %in% nm) .data[[d_col]]   else NA_character_,
        `Latitude Start`  = suppressWarnings(as.numeric(.data[[lat_s]])),
        `Longitude Start` = suppressWarnings(as.numeric(.data[[lon_s]])),
        `Latitude Mid`    = if (!is.na(lat_m) && lat_m %in% nm) suppressWarnings(as.numeric(.data[[lat_m]])) else NA_real_,
        `Longitude Mid`   = if (!is.na(lon_m) && lon_m %in% nm) suppressWarnings(as.numeric(.data[[lon_m]])) else NA_real_,
        `Latitude Stop`   = suppressWarnings(as.numeric(.data[[lat_e]])),
        `Longitude Stop`  = suppressWarnings(as.numeric(.data[[lon_e]])),
        `Depth (m)`       = if (!is.na(dep_col) && dep_col %in% nm) suppressWarnings(as.numeric(.data[[dep_col]])) else NA_real_,
        SVY               = if (!is.na(svy_col) && svy_col %in% nm) .data[[svy_col]] else program_name,
        Parameters        = coalesce(Parameters, NA_character_),
        `EOV Groups`      = coalesce(`EOV Groups`, NA_character_),
        source_file, program = program_name
      ) %>%
      filter(!is.na(`Latitude Start`), !is.na(`Latitude Stop`))
  })

  # Single-point fallback for mobile files without explicit start/stop columns
  single_pt_transects <- transect_source_rows %>%
    filter(!source_file %in% files_with_start_stop, !is.na(latitude), !is.na(longitude)) %>%
    left_join(row_param_lookup, by = c("source_file","source_row_id")) %>%
    transmute(
      Date = as.character(year), `Latitude Start` = latitude, `Longitude Start` = longitude,
      `Latitude Mid` = NA_real_, `Longitude Mid` = NA_real_,
      `Latitude Stop` = latitude, `Longitude Stop` = longitude,
      `Depth (m)` = depth, SVY = program,
      Parameters  = coalesce(Parameters, NA_character_),
      `EOV Groups`= coalesce(`EOV Groups`, NA_character_),
      source_file, program
    )

  all_transects <- bind_rows(expanded, single_pt_transects) %>%
    # Year filter
    filter(is.na(Date) | suppressWarnings(as.integer(str_extract(Date, "\\d{4}"))) >= start_year) %>%
    mutate(
      lat_bin_s = round(suppressWarnings(as.numeric(`Latitude Start`)),  1),
      lon_bin_s = round(suppressWarnings(as.numeric(`Longitude Start`)), 1),
      lat_bin_e = round(suppressWarnings(as.numeric(`Latitude Stop`)),   1),
      lon_bin_e = round(suppressWarnings(as.numeric(`Longitude Stop`)),  1),
      year_val  = suppressWarnings(as.integer(str_extract(Date, "\\d{4}")))
    ) %>%
    # Only keep repeatedly sailed corridors (>= 5 years)
    group_by(lat_bin_s, lon_bin_s, lat_bin_e, lon_bin_e) %>%
    mutate(years_visited = n_distinct(year_val, na.rm = TRUE)) %>%
    filter(is.na(lat_bin_s) | years_visited >= 5) %>%
    # One representative line per corridor (most recent)
    slice_max(order_by = year_val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-lat_bin_s, -lon_bin_s, -lat_bin_e, -lon_bin_e, -years_visited, -year_val)

  if (nrow(all_transects) > 0) {
    write_csv(all_transects, file.path(output_folder, "transects.csv"))
    cat("transects.csv written:", nrow(all_transects), "rows\n")
  } else {
    cat("No transect rows passed filters — transects.csv not written.\n")
  }
} else {
  cat("No transect source rows found — transects.csv not written.\n")
}

############################################
# WRITE DIAGNOSTIC CSV OUTPUTS
############################################

cat("\n=== Writing diagnostic CSVs ===\n")

write_csv(detection_summary,                   file.path(output_folder, paste0(program_name, "_detection_summary.csv")))
write_csv(clip_summary,                        file.path(output_folder, paste0(program_name, "_clip_summary.csv")))
write_csv(dataset_classification,              file.path(output_folder, paste0(program_name, "_geometry_classification.csv")))
write_csv(likely_parameter_columns,            file.path(output_folder, paste0(program_name, "_parameter_columns.csv")))
write_csv(filename_parameter_hits,             file.path(output_folder, paste0(program_name, "_filename_parameter_hits.csv")))
write_csv(ambiguous_filename_cases,            file.path(output_folder, paste0(program_name, "_ambiguous_filename_cases.csv")))
write_csv(parameter_review_table,              file.path(output_folder, paste0(program_name, "_parameter_review.csv")))
write_csv(parameter_frequency_table,           file.path(output_folder, paste0(program_name, "_parameter_frequency.csv")))
write_csv(program_hex_parameter_table_active,  file.path(output_folder, paste0(program_name, "_hex_parameter_table.csv")))
write_csv(program_hex_inventory,               file.path(output_folder, paste0(program_name, "_hex_inventory.csv")))
write_csv(hex_overlap_table,                   file.path(output_folder, paste0(program_name, "_hex_overlap.csv")))
write_csv(all_active_sampled_points,           file.path(output_folder, paste0(program_name, "_active_sampled_points.csv")))

############################################
# GEOJSON MAP EXPORT
############################################

cat("\n=== Writing GeoJSON ===\n")

hex_map_sf <- hex_grid_sf %>%
  left_join(hex_inventory_for_map, by = "hex_id") %>%
  filter(!is.na(programs)) %>%
  st_transform(4326) %>% st_make_valid()

hex_map_sf_export <- hex_map_sf %>%
  mutate(
    `Station Key`       = paste(programs, hex_id, sep = " | "),
    `Program Name`      = program_name,
    `Full Program Name` = program_full_name,
    `Geometry Types`    = case_when(
      str_detect(geometry_types_found, "underway") & str_detect(geometry_types_found, "point") ~ "Underway + Station",
      str_detect(geometry_types_found, "underway") & str_detect(geometry_types_found, "line")  ~ "Underway + Transect",
      str_detect(geometry_types_found, "underway") ~ "Continuous Underway",
      str_detect(geometry_types_found, "line")     ~ "Transect / Tow",
      TRUE                                          ~ "Fixed Station"
    ),
    `Frequency` = map_chr(parameters, function(param_str) {
      hex_params  <- str_split(param_str, "; ")[[1]]
      freq_labels <- parameter_frequency_table %>%
        filter(standard_parameter %in% hex_params) %>%
        pull(frequency_label) %>% unique() %>% sort()
      if (length(freq_labels) == 0) sampling_frequency else paste(freq_labels, collapse = "; ")
    }),
    `Platform` = map2_chr(parameters, source_files, function(param_str, src_str) {
      hex_params <- str_split(param_str, "; ")[[1]]
      attr_p     <- parameter_frequency_table %>%
        filter(standard_parameter %in% hex_params, !is.na(attr_platform), attr_platform != "") %>%
        pull(attr_platform) %>% unique() %>% sort()
      if (length(attr_p) > 0) paste(attr_p, collapse = "; ") else {
        hit <- platform_lookup %>%
          filter(str_detect(str_to_lower(src_str), regex(pattern, ignore_case = TRUE))) %>% slice(1)
        if (nrow(hit) == 0) program_platform else hit$platform
      }
    }),
    `Depth Range (m)` = case_when(
      !is.na(min_depth) & !is.na(max_depth) & min_depth != max_depth ~
        paste0(round(min_depth), "\u2013", round(max_depth), " m"),
      !is.na(min_depth) ~ paste0(round(min_depth), " m"),
      TRUE ~ NA_character_
    )
  ) %>%
  rename(
    `Monitoring Program` = programs,
    `Parameters`         = parameters,
    `EOV Groups`         = eov_groups,
    `Parameter Count`    = parameter_count,
    `First Year`         = first_year,
    `Last Year`          = last_year,
    `Years Sampled`      = year_count,
    `Source Files`       = source_files,
    `Centroid Latitude`  = centroid_latitude,
    `Centroid Longitude` = centroid_longitude,
    `Sample Locations`   = sample_location_count,
    `Gebco Mean Depth`   = gebco_mean_depth
  ) %>%
  select(
    `Program Name`, `Full Program Name`, `Monitoring Program`,
    `Parameters`, `EOV Groups`, `Parameter Count`,
    `First Year`, `Last Year`, `Years Sampled`,
    `Frequency`, `Platform`, `Geometry Types`,
    `Depth Range (m)`, `Sample Locations`,
    `Centroid Latitude`, `Centroid Longitude`,
    `Gebco Mean Depth`= gebco_mean_depth,
    starts_with("param_"),
    geometry
  )

cat("GEBCO values sample:\n")
print(head(hex_map_sf_export$`Gebco Mean Depth`, 10))

suppressWarnings(
  st_write(hex_map_sf_export,
           file.path(output_folder, paste0(display_name, ".geojson")),
           delete_dsn = TRUE, quiet = TRUE)
)
cat("GeoJSON written:", file.path(output_folder, paste0(display_name, ".geojson")), "\n")

############################################
# STANDALONE LEAFLET HTML MAP
############################################

cat("\n=== Writing standalone HTML map ===\n")

geojson_content <- readr::read_file(file.path(output_folder, paste0(display_name, ".geojson")))

leaflet_html <- paste0(
'<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>', display_name, ' \u2014 Ocean Monitoring Inventory</title>
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: "Segoe UI", Arial, sans-serif; }
    #map { width: 100vw; height: 100vh; }
    .leaflet-popup-content-wrapper { border-radius: 6px; padding: 0; overflow: hidden;
      box-shadow: 0 4px 18px rgba(0,0,0,0.22); min-width: 300px; max-width: 420px; }
    .leaflet-popup-content { margin: 0; width: 100% !important; }
    .pp-header   { background: #005e8c; color: #fff; padding: 14px 18px 12px; }
    .pp-title    { font-size: 16px; font-weight: 700; line-height: 1.3; margin-bottom: 3px; }
    .pp-subtitle { font-size: 11px; opacity: 0.8; font-style: italic; }
    .pp-narrative { background: #f0f6fb; padding: 9px 18px; font-size: 12px; color: #333;
                    border-bottom: 1px solid #d6e6f0; line-height: 1.55; }
    .pp-narrative b { color: #005e8c; }
    .pp-table { width: 100%; border-collapse: collapse; font-size: 12px; }
    .pp-table tr:nth-child(even) td { background: #f7f7f7; }
    .pp-table td { padding: 6px 18px; vertical-align: top;
                   border-bottom: 1px solid #ececec; line-height: 1.4; }
    .pp-table td:first-child { font-weight: 600; color: #555; white-space: nowrap; width: 38%; }
    .pp-section-head td { background: #e8f0f5 !important; font-size: 10px; font-weight: 700;
      text-transform: uppercase; letter-spacing: 0.06em; color: #005e8c; padding: 5px 18px; }
  </style>
</head>
<body>
<div id="map"></div>
<script>
var geojsonData = ', geojson_content, ';
var map = L.map("map").setView([36.2, -121.0], 7);
L.tileLayer(
  "https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}",
  { attribution: "Tiles &copy; Esri", maxZoom: 18 }
).addTo(map);

function val(v) {
  return (v===null||v===undefined||v==="NA"||v===""||v==="null")?null:v;
}
function row(label, value) {
  if (!val(value)) return "";
  return "<tr><td>" + label + "</td><td>" + value + "</td></tr>";
}
function sh(label) {
  return "<tr class=\'pp-section-head\'><td colspan=\'2\'>" + label + "</td></tr>";
}
function makePopup(p) {
  var narrative = val(p["Parameters"])
    ? "<div class=\'pp-narrative\'><b>Parameters:</b> " + p["Parameters"] + "</div>"
    : "";
  var depthRaw = val(p["Depth Range (m)"]);
  var depthDisplay = depthRaw;
  if (depthRaw) {
    var nums = depthRaw.match(/-?([\\d.]+)/g);
    if (nums && nums.length === 1)
      depthDisplay = "0\u2013" + Math.abs(parseFloat(nums[0])) + " m";
    else if (nums && nums.length >= 2)
      depthDisplay = Math.abs(parseFloat(nums[0])) + "\u2013" + Math.abs(parseFloat(nums[1])) + " m";
  }
  var tableRows =
    sh("Sampling") +
    row("Frequency",      p["Frequency"]) +
    row("Platform",       p["Platform"]) +
    row("First Year",     p["First Year"]) +
    row("Last Year",      p["Last Year"]) +
    row("Years Sampled",  p["Years Sampled"]) +
    row("Sample Locations", p["Sample Locations"]) +
    sh("Depth") +
    row("Overall Range",  depthDisplay) +
    row("Seafloor Depth", val(p["Gebco Mean Depth"])
      ? parseFloat(p["Gebco Mean Depth"]).toFixed(0) + " m" : null) +
    sh("Location") +
    row("Latitude",   val(p["Centroid Latitude"])  ? parseFloat(p["Centroid Latitude"]).toFixed(4)  : null) +
    row("Longitude",  val(p["Centroid Longitude"]) ? parseFloat(p["Centroid Longitude"]).toFixed(4) : null) +
    row("EOV Groups", p["EOV Groups"]);

  return "<div class=\'pp-header\'>" +
    "<div class=\'pp-title\'>"    + (val(p["Program Name"]) || "Station") + "</div>" +
    "<div class=\'pp-subtitle\'>" + (val(p["Full Program Name"]) || "")   + "</div>" +
    "</div>" + narrative + "<table class=\'pp-table\'>" + tableRows + "</table>";
}

L.geoJSON(geojsonData, {
  style: function() { return { color: "#005e8c", weight: 1, fillColor: "#0079c1", fillOpacity: 0.3 }; },
  onEachFeature: function(feature, layer) {
    layer.bindPopup(makePopup(feature.properties), { maxWidth: 420 });
    layer.on("mouseover", function() { this.setStyle({ fillOpacity: 0.55, weight: 2 }); });
    layer.on("mouseout",  function() { this.setStyle({ fillOpacity: 0.3,  weight: 1 }); });
  }
}).addTo(map);
</script>
</body>
</html>')

writeLines(leaflet_html,
           file.path(output_folder, paste0(display_name, "_map.html")),
           useBytes = FALSE)

cat("HTML map written:", file.path(output_folder, paste0(display_name, "_map.html")), "\n")

cat("\n===== WORKFLOW COMPLETE =====\n")
cat("Program:             ", program_name, "\n")
cat("CSV files read:      ", length(csv_files), "\n")
cat("Rows after year filter:", nrow(all_filtered_df), "\n")
cat("Rows after clip:     ", nrow(all_filtered_coords_df), "\n")
cat("Active hex rows:     ", nrow(program_hex_inventory), "\n")
cat("GeoJSON:             ", file.path(output_folder, paste0(display_name, ".geojson")), "\n")
cat("=================================\n")
