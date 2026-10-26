################################################################################
# SI Idiographic Prediction Project - Data Cleaning Script
#
# Companion to si_data_analysis.R, which consumes the output built here
# (all_features_cleaned.csv, the inclusion report, and three small report files).
# Run this script first.
#
# Folder layout expected under BASE_PATH (the working directory set below):
#   1. Code
#   2. Raw data
#   3. Cleaned data
#   4. Results
#
# The script runs top to bottom when sourced:
#   1. shared dynamic-index helpers
#   2. EMA
#   3. passive sensing data
#   4. feature assembly and analytic inclusion
#
# During the run, many intermediate files are written under "3. Cleaned data"
# and "4. Results"; a cleanup step at the very end removes them, keeping
# all_features_cleaned.csv, Table1_inclusion_participant_report.csv and three small report files
# (raw_passive_compliance_by_uid_sensor.csv, dropped_high_missing_columns.csv,
# dropped_no_variability_columns.csv) that the analysis script and the Methods use.
################################################################################

# ------------------------------------------------------------------------------
# 1) Set the working directory to the SI_Prediction project folder.
# ------------------------------------------------------------------------------
setwd("/Users/ahhyun/Library/CloudStorage/GoogleDrive-ahhyun.lee@dartmouth.edu/My Drive/SI_Prediction")       # your path from the previous file; replace or remove it before the code is shared
# Guard: the working directory must be the top-level SI_Prediction folder, not one of its subfolders
# (otherwise paths double up, e.g. ".../2. Raw data/2. Raw data/...").
if (basename(getwd()) %in% c("1. Code", "2. Raw data", "3. Cleaned data", "4. Results")) {
  setwd(dirname(getwd()))
  message("Working directory pointed at a subfolder; using its parent instead: ", getwd())
}
if (!dir.exists("2. Raw data"))
  stop("Could not find '2. Raw data' in the working directory: ", getwd(),
       "\nSet the working directory above to the SI_Prediction folder that contains '2. Raw data'.")
cat("Working directory:", getwd(), "\n")

# Design decisions that change the features (state them in the Methods).
ZERO_FILL_EVENT_STREAMS         <- TRUE   # SMS / calls / conversation: no event in a bin = 0, not missing
SUPPRESS_ENTROPY_EVENT_STREAMS  <- FALSE  # TRUE = do not compute sample entropy for those three streams
DEDUP_EVENT_LOGS                <- TRUE   # SMS and call logs: drop rows that are exact copies of another row before binning (no effect if there are none)
SLEEP_STRICT_PRIOR_NIGHTS       <- FALSE  # FALSE = the night before the assessment, as preregistered (OSF): day1 = that night, day7 = the 7 nights ending with it.
# TRUE  = only nights BEFORE the prompt date (stricter 4-hour separation, but not what the OSF registered)


# ------------------------------------------------------------------------------
# 2) Install any missing packages, then load them all.
# ------------------------------------------------------------------------------
packages <- c(
  "data.table", "stringr", "lubridate", "parallel", "pbmcapply", "geosphere",
  "tidyverse", "progress", "TSEntropies"
)
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# Create directories quietly, including missing parent folders.
make_dirs <- function(...) {
  invisible(lapply(c(...), dir.create, showWarnings = FALSE, recursive = TRUE))
}

# ------------------------------------------------------------------------------
# SHARED PROJECT PATHS
#
# All derived once from the working directory set above and reused by every
# section below.
# ------------------------------------------------------------------------------
BASE_PATH <- getwd()
CODE_DIR <- file.path(BASE_PATH, "1. Code")
RAW_DATA_DIR <- file.path(BASE_PATH, "2. Raw data")
CLEAN_DATA_DIR <- file.path(BASE_PATH, "3. Cleaned data")
CLEAN_SENSOR_DIR <- file.path(CLEAN_DATA_DIR, "per_sensors")
CLEAN_UID_DIR <- file.path(CLEAN_DATA_DIR, "per_uid")
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
results_dir <- file.path(RESULTS_ROOT_DIR, "report")
# Only the cleaned-data folders are created up front; the results folder is
# created later, right before the first results are written.
make_dirs(CLEAN_DATA_DIR, CLEAN_SENSOR_DIR, CLEAN_UID_DIR)

cat("\n===== SECTION: DATA CLEANING =====\n")
# ==============================================================================
# DATA CLEANING, 1 OF 4: SHARED DYNAMIC-INDEX HELPERS (Methods-aligned N thresholds)
# ==============================================================================
# Mean (N >= 1)
mean_thr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 1) return(NA_real_)
  base::mean(x)
}
# MAD (N >= 2)
mad_thr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  stats::mad(x)
}
# P90-P10 relative position: (current_x - P10) / (P90 - P10).
# P90 == P10 returns NA (undefined; the old code injected 0.5). Values are not clipped.
# NOTE: this is the RELATIVE POSITION of the current value, not a dispersion measure.
p90p10_relpos <- function(window_vals, current_x) {
  window_vals <- window_vals[!is.na(window_vals)]
  if (length(window_vals) < 20) return(NA_real_)
  if (is.na(current_x)) return(NA_real_)
  p10 <- unname(stats::quantile(window_vals, 0.10, type = 7))
  p90 <- unname(stats::quantile(window_vals, 0.90, type = 7))
  if (p90 == p10) return(NA_real_)
  (current_x - p10) / (p90 - p10)
}
# Sample entropy (N >= 250; m = 2, r = 0.2 * SD)
sample_entropy <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 250) return(NA_real_)
  r <- 0.2 * stats::sd(x)
  if (is.na(r) || r == 0) return(NA_real_)
  if (!requireNamespace("TSEntropies", quietly = TRUE)) return(NA_real_)
  v <- tryCatch(TSEntropies::FastSampEn(x, dim = 2, r = r), error = function(e) NA_real_)
  v <- suppressWarnings(as.numeric(v))
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v     # undefined when no template matches
}
################################################################################
# Make UIDs and dates consistent across data sources.
################################################################################
# Participant ID -> canonical 4-digit form ("11", "11@m.life", "pt-11" -> "0011").
normalize_uid <- function(x) {
  raw <- trimws(as.character(x))
  raw[raw %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  id <- sub("@.*$", "", raw)                  # drop anything after @
  digits <- gsub("\\D", "", id)               # keep only the digits
  out <- ifelse(!is.na(id) & nzchar(digits), digits, id)
  pad <- !is.na(out) & grepl("^\\d+$", out) & nchar(out) <= 4
  out[pad] <- stringr::str_pad(out[pad], 4, pad = "0")
  out
}
# Any common date format -> Date.
parse_date_any <- function(x) {
  raw <- trimws(as.character(x))
  raw[raw %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  out <- as.Date(suppressWarnings(lubridate::parse_date_time(
    raw,
    orders = c("Ymd", "mdY", "dmY", "Ymd HMS", "mdY HMS", "dmY HMS", "B d Y", "d B Y"),
    quiet = TRUE
  )))
  num <- suppressWarnings(as.numeric(raw))
  unix <- is.na(out) & !is.na(num) & num > 1e9          # Unix epoch seconds
  out[unix] <- as.Date(as.POSIXct(num[unix], origin = "1970-01-01", tz = "UTC"))
  excel <- is.na(out) & !is.na(num) & num >= 20000 & num <= 80000  # Excel serial days
  out[excel] <- as.Date(num[excel], origin = "1899-12-30")
  failed <- !is.na(raw) & is.na(out)
  if (any(failed)) {
    warning(sprintf("parse_date_any(): %d unrecognized date(s) set to NA (e.g., %s)",
                    sum(failed), paste(utils::head(unique(raw[failed]), 3), collapse = ", ")),
            call. = FALSE)
  }
  out
}
# ==============================================================================
# DATA CLEANING, 2 OF 4: EMA PROCESSING
# ==============================================================================
EMA_FILE          <- "3._TD_cohort_EMA_FINAL.csv"
ema_features_path <- file.path(CLEAN_DATA_DIR, "ema_features.csv")
derive_compliant_uids_from_raw <- function() {
  per_uid_dirs <- c("acceleration", "hr", "steps", "battery", "sleep", "raw_gps_clean")
  uid_hits <- character(0)
  for (d in per_uid_dirs) {
    dir_path <- file.path(RAW_DATA_DIR, d)
    if (!dir.exists(dir_path)) next
    fs <- list.files(dir_path, pattern = "\\.csv$", full.names = FALSE)
    uid_hits <- c(uid_hits, normalize_uid(tools::file_path_sans_ext(fs)))
  }
  convo_dir <- file.path(RAW_DATA_DIR, "convo")
  if (dir.exists(convo_dir)) {
    fs <- list.files(convo_dir, pattern = "\\.csv$", full.names = FALSE)
    uid_hits <- c(uid_hits, normalize_uid(stringr::str_extract(fs, "(?<=convo_)[^@._]+")))
  }
  for (combined in c("smslog.csv", "calllog.csv", "sms.csv", "calls.csv")) {
    fp <- file.path(RAW_DATA_DIR, combined)
    if (file.exists(fp)) {
      cols <- names(data.table::fread(fp, nrows = 0))
      id_col <- cols[grepl("^uid$|mlife_id|participant|subject|user", cols, ignore.case = TRUE)][1]
      if (!is.na(id_col)) {
        ids <- unique(data.table::fread(fp, select = id_col)[[1]])
        uid_hits <- c(uid_hits, normalize_uid(ids))
      }
    }
  }
  sort(unique(uid_hits[!is.na(uid_hits) & nzchar(uid_hits)]))
}
cat("\n===== Build EMA features =====\n")
compliant_uids <- derive_compliant_uids_from_raw()
cat(sprintf("Compliant UIDs from raw sensor data: %d\n", length(compliant_uids)))
write.csv(data.frame(uid = compliant_uids),
          file.path(CLEAN_DATA_DIR, "compliant_uids.csv"), row.names = FALSE)
ema_path <- file.path(RAW_DATA_DIR, EMA_FILE)
if (!file.exists(ema_path))
  stop("EMA file not found: ", ema_path,
       "\nCheck that the working directory (BASE_PATH) is the SI_Prediction folder and that EMA_FILE matches the file name.")
ema_raw <- read.csv(ema_path)
# DIAGNOSTIC: is `compliance` a participant-level constant or does it vary by row?
# The filter further down (compliance >= 0.75) is applied ROW-wise. If compliance varies within a
# participant, it drops individual days (and the later participant-level check becomes trivial).
{
  .cmp <- tapply(ema_raw$compliance, ema_raw$mlife_id, function(z) length(unique(z[!is.na(z)])))
  cat(sprintf("[CHECK] compliance takes >1 distinct value within %d of %d participants%s\n",
              sum(.cmp > 1, na.rm = TRUE), length(.cmp),
              if (any(.cmp > 1, na.rm = TRUE)) "  -> the 0.75 filter drops individual rows; describe it as such" else "  -> participant-level"))
  rm(.cmp)
}
ema_selected <- ema_raw %>%
  select(mlife_id, compliance, date_rep, ema_time, ema_timestamp,
         starts_with("PHQ."), starts_with("GAD."), starts_with("RespTime_")) %>%
  rename_with(~ str_replace(.x, "^PHQ\\.", "q"), starts_with("PHQ.")) %>%
  rename_with(~ str_replace(.x, "^GAD\\.", "g"), starts_with("GAD.")) %>%
  mutate(
    ema_day = select(., starts_with("RespTime_")) %>%
      apply(1, function(row) {
        ts <- unlist(str_extract_all(as.character(row), "(?<=;)[0-9.]+"))
        if (length(ts) == 0 || all(is.na(ts))) return(NA_real_)
        max(as.numeric(ts), na.rm = TRUE)
      }) %>% as.POSIXct(origin = "1970-01-01", tz = "UTC")
  ) %>%
  select(mlife_id, compliance, date_rep, ema_time, ema_timestamp, ema_day,
         all_of(c(paste0("q", 1:9), "g1", "g2"))) %>%
  mutate(compliance = pmin(compliance, 1)) %>%
  filter(compliance >= 0.75)
ema_complete <- ema_selected %>%
  filter(!is.na(date_rep), !is.na(ema_time)) %>%
  mutate(date_rep = as.Date(date_rep, format = "%m/%d/%Y"),
         ema_time = factor(ema_time, levels = c("morning","afternoon","evening"))) %>%
  group_by(mlife_id) %>%
  complete(date_rep = seq(min(date_rep, na.rm = TRUE),
                          max(date_rep, na.rm = TRUE), by = "day"),
           ema_time = c("morning","afternoon","evening")) %>%
  filter(!(date_rep == min(date_rep) & ema_time < first(ema_time[!is.na(compliance)])),
         !(date_rep == max(date_rep) & ema_time > last(ema_time[!is.na(compliance)]))) %>%
  ungroup()
cat(sprintf("EMA rows: %d | Participants: %d\n",
            nrow(ema_complete), uniqueN(ema_complete$mlife_id)))
vars <- c(paste0("q", 1:9), "g1", "g2")
dt <- as.data.table(ema_complete)
dt[, ema_num    := as.numeric(ema_day)]
dt[, anchor_end := ema_num - 3600 * 4]
dt[, row_id     := .I]
setkey(dt, mlife_id, ema_num)
n_rows <- nrow(dt)
d1_mat <- matrix(NA_real_, nrow = n_rows, ncol = length(vars),
                 dimnames = list(NULL, paste0(vars, "_mean_d1")))
d7_mat <- matrix(NA_real_, nrow = n_rows, ncol = length(vars) * 2,
                 dimnames = list(NULL, c(paste0(vars, "_mean_d7"),
                                         paste0(vars, "_mad_d7"))))
pb_ema <- progress::progress_bar$new(
  format = "  [:bar] :current/:total | :percent | eta: :eta",
  total = n_rows, clear = FALSE, width = 80)
for (idx in split(seq_len(n_rows), dt$mlife_id)) {
  v_all <- dt$ema_num[idx]
  ok    <- which(!is.na(v_all))
  v_ok  <- v_all[ok]
  vals  <- lapply(vars, function(v) dt[[v]][idx][ok]); names(vals) <- vars
  for (r in seq_along(idx)) {
    pb_ema$tick()
    i  <- idx[r]
    ae <- dt$anchor_end[i]
    if (is.na(ae) || !length(v_ok)) next
    hi  <- findInterval(ae,          v_ok, left.open = TRUE)       # ema_num <  ae
    lo1 <- findInterval(ae - 86400,  v_ok, left.open = TRUE) + 1L  # >= ae - 1d
    lo7 <- findInterval(ae - 604800, v_ok, left.open = TRUE) + 1L  # >= ae - 7d
    if (hi >= lo1) {
      w1 <- lo1:hi
      d1_mat[i, ] <- sapply(vars, function(v) mean_thr(vals[[v]][w1]))
    }
    if (hi >= lo7) {
      w7 <- lo7:hi
      d7_mat[i, seq_along(vars)]                <- sapply(vars, function(v) mean_thr(vals[[v]][w7]))
      d7_mat[i, seq_along(vars) + length(vars)] <- sapply(vars, function(v) mad_thr(vals[[v]][w7]))
    }
  }
}
dt <- cbind(dt, as.data.table(d1_mat), as.data.table(d7_mat))
dt[, c("ema_num","anchor_end","row_id") := NULL]
ema_features <- as.data.frame(dt) %>%
  mutate(uid = stringr::str_pad(str_extract(mlife_id, "^[^@]+"), width = 4, pad = "0"))
if (length(compliant_uids) > 0) {
  ema_features <- ema_features %>% filter(uid %in% compliant_uids)
}
data.table::fwrite(ema_features, ema_features_path)
cat(sprintf("Saved: ema_features.csv (%d rows | %d PTs)\n",
            nrow(ema_features), uniqueN(ema_features$uid)))
summary(ema_features)
# ==============================================================================
# DATA CLEANING, 3 OF 4: PASSIVE SENSING DATA
# ==============================================================================
# ------------------------------------------------------------------------------
# 3a. Shared passive-sensing setup
# ------------------------------------------------------------------------------
WORKERS           <- max(1, parallel::detectCores() - 1)   # auto: all cores minus 1
SEC_2020          <- 1577836800
ACC_EPOCH_S       <- 10
MAX_GPS_SPEED_MPS <- 50    # ~180 km/h; consecutive GPS fixes implying more are corrupt
WINDOWS           <- list(hr1 = 3600, hr4 = 14400, day1 = 86400, day7 = 604800)
data.table::setDTthreads(1) # 1 thread/worker -> avoid oversubscription across forks
uid_default <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  normalize_uid(b)
}
uid_convo <- function(fp) stringr::str_extract(basename(fp), "(?<=convo_)[^@]+")
# --- EMA anchors ----------------------------------------------------------------
build_ema_lookup <- function(path = file.path(CLEAN_DATA_DIR, "ema_features.csv")) {
  d <- data.table::fread(path)
  if (!inherits(d$ema_day, "POSIXct")) {
    raw <- as.character(d$ema_day)
    p   <- suppressWarnings(as.POSIXct(raw, tz = "UTC"))
    if (all(is.na(p))) p <- lubridate::ymd_hms(raw, tz = "UTC")
    data.table::set(d, j = "ema_day", value = p)
  }
  d <- d[!is.na(ema_day)]
  if (!"uid" %in% names(d)) d[, uid := mlife_id]
  d[, uid := normalize_uid(uid)]
  d[, anchor_end := ema_day - 3600 * 4]
  out <- d[, .(uid, ema_day, date_rep, ema_time, anchor_end)]
  data.table::setkey(out, uid); out[]
}
ema_lookup <- build_ema_lookup()
cat(sprintf("ema_lookup: %d rows, %d uids\n", nrow(ema_lookup), data.table::uniqueN(ema_lookup$uid)))
# --- uid alias map (raw sensor id variants -> canonical EMA uid) ----------------
build_uid_alias_lookup <- function(path = file.path(CLEAN_DATA_DIR, "ema_features.csv")) {
  d <- data.table::fread(path, colClasses = list(character = c("uid", "mlife_id")))
  if (!"uid" %in% names(d) && !"mlife_id" %in% names(d)) return(setNames(character(0), character(0)))
  canonical <- if ("uid" %in% names(d)) normalize_uid(d$uid) else normalize_uid(d$mlife_id)
  alias_cols <- intersect(c("uid", "mlife_id"), names(d))
  alias_rows <- data.table::rbindlist(lapply(alias_cols, function(col) {
    raw <- trimws(as.character(d[[col]]))
    data.table::data.table(
      canonical_uid = canonical,
      alias = c(raw, tolower(raw), normalize_uid(raw), tolower(normalize_uid(raw)),
                stringr::str_extract(raw, "^[^@]+"), tolower(stringr::str_extract(raw, "^[^@]+")))
    )
  }), use.names = TRUE, fill = TRUE)
  alias_rows <- alias_rows[!is.na(canonical_uid) & !is.na(alias) & nzchar(alias)]
  alias_rows <- unique(alias_rows)
  alias_rows <- alias_rows[!duplicated(alias)]
  stats::setNames(alias_rows$canonical_uid, alias_rows$alias)
}
UID_ALIAS <- build_uid_alias_lookup()
map_to_ema_uid <- function(x) {
  raw <- trimws(as.character(x))
  norm <- normalize_uid(raw)
  keys <- list(raw, tolower(raw), norm, tolower(norm),
               stringr::str_extract(raw, "^[^@]+"),
               tolower(stringr::str_extract(raw, "^[^@]+")))
  out <- rep(NA_character_, length(raw))
  for (key in keys) {
    hit <- !is.na(key) & key %in% names(UID_ALIAS) & is.na(out)
    out[hit] <- UID_ALIAS[key[hit]]
  }
  out[is.na(out)] <- norm[is.na(out)]
  out
}
# --- HOME LOCATIONS -------------------------------------------------------------
first_existing <- function(paths) { hit <- paths[file.exists(paths)]; if (length(hit)) hit[1] else NA_character_ }
DERIVED_HOME_FILE <- first_existing(c(file.path(CLEAN_DATA_DIR, "homes_derived.csv"),
                                      file.path(RAW_DATA_DIR,   "homes_derived.csv")))
HOME_FILE         <- first_existing(c(file.path(CLEAN_DATA_DIR, "time_at_home.csv"),
                                      file.path(RAW_DATA_DIR,   "time_at_home.csv")))
load_home_coords <- function() {
  if (!is.na(DERIVED_HOME_FILE)) {
    d <- data.table::fread(DERIVED_HOME_FILE, colClasses = list(character = "uid"))
    d[, uid := stringr::str_pad(uid, 4, pad = "0")]
    L <- lapply(split(d, d$uid), function(x) {
      m <- as.matrix(x[, .(home_lat, home_lon)]); colnames(m) <- c("lat", "lon"); m
    })
    cat(sprintf("homes loaded from %s for %d participants\n", DERIVED_HOME_FILE, length(L)))
    return(L)
  }
  if (!is.na(HOME_FILE)) {
    d <- unique(data.table::fread(HOME_FILE, select = c("uid", "home_coords")), by = "uid")
    L <- lapply(seq_len(nrow(d)), function(i) {
      nums <- suppressWarnings(as.numeric(unlist(
        regmatches(as.character(d$home_coords[i]),
                   gregexpr("-?[0-9]+\\.?[0-9]*(?:[eE]-?[0-9]+)?", as.character(d$home_coords[i]))))))
      nums <- nums[!is.na(nums)]
      if (length(nums) < 2) return(NULL)
      m <- matrix(nums[1:(2 * (length(nums) %/% 2))], ncol = 2, byrow = TRUE)
      colnames(m) <- c("lat", "lon"); m
    })
    names(L) <- stringr::str_pad(stringr::str_extract(d$uid, "^[^@]+"), 4, pad = "0")
    cat(sprintf("homes loaded from time_at_home.csv for %d participants\n",
                sum(!vapply(L, is.null, logical(1)))))
    return(L)
  }
  cat("WARN: no home file found (run SI_derive_homes.R first) -> gps_home will be EMPTY\n")
  list()
}
HOME_COORDS <- load_home_coords()
# ------------------------------------------------------------------------------
# UNIT SANITY CHECKS  (printed once per run; do NOT alter any output)
#   ACC   : resting vector-magnitude of x/y/z should be ~9.81 if the sensor logs
#           m/s^2 (ENMO = vm/9.80665 - 1 assumes this). A median near ~1.0 means
#           the data is already in g and the ENMO conversion is wrong.
#   STEPS : element-1 of each record should be a per-interval step COUNT (small
#           integers). A large, mostly-non-decreasing series means it is a
#           cumulative daily counter and summing it over windows is invalid.
# ------------------------------------------------------------------------------
sanity_check_units <- function() {
  # ACC
  af <- list.files(file.path(RAW_DATA_DIR, "acceleration"), pattern = "\\.csv$", full.names = TRUE)
  if (length(af)) {
    a <- tryCatch(data.table::fread(af[1], select = "data", nrows = 20000), error = function(e) NULL)
    if (!is.null(a) && "data" %in% names(a)) {
      parts <- data.table::tstrsplit(gsub("[][]", "", a$data), ",", fixed = TRUE, type.convert = TRUE)
      if (length(parts) >= 3) {
        vm <- sqrt(as.numeric(parts[[1]])^2 + as.numeric(parts[[2]])^2 + as.numeric(parts[[3]])^2)
        vm <- vm[is.finite(vm)]
        if (length(vm))
          cat(sprintf("[SANITY] ACC   %s: median |xyz| = %.2f  (expect ~9.81 for m/s^2; ~1.0 => data in g, ENMO conversion wrong)\n",
                      basename(af[1]), stats::median(vm)))
      }
    }
  }
  # STEPS
  sf <- list.files(file.path(RAW_DATA_DIR, "steps"), pattern = "\\.csv$", full.names = TRUE)
  if (length(sf)) {
    s <- tryCatch(data.table::fread(sf[1], nrows = 20000), error = function(e) NULL)
    if (!is.null(s) && "data" %in% names(s)) {
      m <- stringr::str_match(s$data,
                              "\\[\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*\\]")
      v <- suppressWarnings(as.numeric(m[, 2])); v <- v[is.finite(v)]
      if (length(v))
        cat(sprintf("[SANITY] STEPS %s: element-1 median = %.0f, max = %.0f  (per-interval counts stay small; a cumulative counter rises into the thousands)\n",
                    basename(sf[1]), stats::median(v), max(v)))
    }
  }
}
sanity_check_units()
# ------------------------------------------------------------------------------
# 3b. ACCELEROMETER
# ------------------------------------------------------------------------------
read_acc <- function(fp) {
  dt <- tryCatch(data.table::fread(fp, select = c("data", "day")),
                 error = function(e) NULL)
  if (is.null(dt) || !all(c("data", "day") %in% names(dt))) return(NULL)
  dt[, c("x", "y", "z", "acc") :=
       data.table::tstrsplit(gsub("[][]", "", data), ",", fixed = TRUE, type.convert = TRUE)]
  dt[, vm   := sqrt(x^2 + y^2 + z^2)]
  dt[, enmo := pmax(vm / 9.80665 - 1, 0)]
  dt[, day  := as.POSIXct(as.character(day), tz = "UTC")]
  dt <- dt[!is.na(enmo) & !is.na(day) & as.numeric(day) >= SEC_2020]
  if (!nrow(dt)) return(NULL)
  dt[, ep := floor(as.numeric(day) / 10) * 10]
  agg <- dt[, .(signal = mean(enmo)), by = ep]
  agg[, day := as.POSIXct(ep, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(agg, day); agg[, .(day, signal)]
}
output_dir <- file.path(CLEAN_UID_DIR, "acc_per_uid")
make_dirs(output_dir)
process_one_acc <- function(fp) {
  uid_this <- stringr::str_pad(as.character(uid_default(fp)), width = 4, pad = "0")
  out_path <- file.path(output_dir, sprintf("%s_dynamic_indices_acc.csv", uid_this))
  tryCatch({
    dt <- read_acc(fp)
    if (is.null(dt) || !nrow(dt)) return(paste("EMPTY", uid_this))
    ema_this <- ema_lookup[uid == uid_this]
    if (!nrow(ema_this)) return(paste("NOEMA", uid_this))
    max_t <- as.numeric(max(ema_this$anchor_end))
    min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
    dt <- dt[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
    if (!nrow(dt)) return(paste("NORANGE", uid_this))
    data.table::setkey(dt, day)
    day_vec <- as.numeric(dt$day); sig_vec <- dt$signal
    
    out <- vector("list", nrow(ema_this) * length(WINDOWS)); k <- 1L
    for (j in seq_len(nrow(ema_this))) {
      ae <- as.numeric(ema_this$anchor_end[j])
      ae_b    <- floor(ae / ACC_EPOCH_S) * ACC_EPOCH_S      # last complete 10-s bin boundary
      cur_len <- max(1800, 2 * ACC_EPOCH_S)
      cl <- findInterval(ae_b - cur_len, day_vec); ch <- findInterval(ae_b - ACC_EPOCH_S, day_vec)
      x_current <- mean_thr(if (ch > cl) sig_vec[(cl + 1L):ch] else numeric(0))
      for (wn in names(WINDOWS)) {
        ws <- WINDOWS[[wn]]
        wl <- findInterval(ae_b - ws, day_vec, left.open = TRUE); wh <- findInterval(ae_b - ACC_EPOCH_S, day_vec)
        xc <- { x <- if (wh > wl) sig_vec[(wl + 1L):wh] else numeric(0); x[!is.na(x)] }
        n_obs <- length(xc)
        cov   <- min(n_obs / (ws / ACC_EPOCH_S), 1)
        mean_v <- mean_thr(xc); sum_v <- if (n_obs) sum(xc) else NA_real_
        mad_v <- mad_thr(xc); p_v <- p90p10_relpos(xc, x_current); ent_v <- sample_entropy(xc)
        out[[k]] <- list(
          uid = ema_this$uid[j], ema_day = ema_this$ema_day[j],
          date_rep = ema_this$date_rep[j], ema_time = ema_this$ema_time[j],
          anchor_end = ema_this$anchor_end[j], window = wn,
          n_obs = n_obs, coverage = cov, mean_val = mean_v, sum_val = sum_v,
          mad_val = mad_v, p90p10 = p_v, entropy = ent_v)
        k <- k + 1L
      }
    }
    data.table::fwrite(data.table::rbindlist(out), out_path)
    paste("DONE", uid_this, nrow(ema_this))
  }, error = function(e) paste("ERROR", uid_this, conditionMessage(e)))
}
acc_final_path <- file.path(CLEAN_SENSOR_DIR, "dynamic_indices_acc.csv")
unlink(acc_final_path)
unlink(output_dir, recursive = TRUE)
make_dirs(output_dir)
acc_files <- list.files(file.path(RAW_DATA_DIR, "acceleration"),
                        pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("acc files: %d | workers: %d\n", length(acc_files), WORKERS))
t0 <- proc.time()
logs <- if (requireNamespace("pbmcapply", quietly = TRUE)) {
  pbmcapply::pbmclapply(acc_files, process_one_acc,
                        mc.cores = WORKERS, mc.preschedule = FALSE)
} else {
  parallel::mclapply(acc_files, process_one_acc,
                     mc.cores = WORKERS, mc.preschedule = FALSE)
}
cat(sprintf("Parallel acc done in %.1f min\n", (proc.time() - t0)["elapsed"] / 60))
status <- sub(" .*", "", unlist(logs))
print(table(status))
errs <- grep("^ERROR", unlist(logs), value = TRUE)
if (length(errs)) { cat("\nErrors:\n"); print(head(errs, 20)) }
done <- list.files(output_dir, pattern = "_dynamic_indices_acc\\.csv$", full.names = TRUE)
acc_res <- data.table::rbindlist(lapply(done, function(f)
  data.table::fread(f, colClasses = list(character = "uid"))))
acc_res[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
data.table::fwrite(acc_res, acc_final_path)
cat(sprintf("\nSaved dynamic_indices_acc.csv: %d rows, %d participants\n",
            nrow(acc_res), data.table::uniqueN(acc_res$uid)))
print(acc_res[, .(n = .N,
                  mean_coverage  = round(mean(coverage, na.rm = TRUE), 2),
                  pct_gated      = round(mean(is.na(mean_val)) * 100),
                  pct_entropy_na = round(mean(is.na(entropy)) * 100)), by = window])
data.table::setDTthreads(0)
# ------------------------------------------------------------------------------
# 3c. GENERIC PARALLEL SENSOR RUNNER
# ------------------------------------------------------------------------------
# ==============================================================================
# READERS  -> each returns data.table(day = POSIXct UTC, signal = numeric)
# ==============================================================================
read_hr <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !all(c("data","day") %in% names(dt))) return(NULL)
  dt[, day := as.POSIXct(as.character(day), tz = "UTC")]
  dt[, signal := suppressWarnings(as.numeric(data))]
  dt <- dt[!is.na(day) & as.numeric(day) >= SEC_2020 &
             !is.na(signal) & signal >= 30 & signal <= 220]
  if (!nrow(dt)) return(NULL)
  data.table::setkey(dt, day); dt[, .(day, signal)]
}
read_steps <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !"data" %in% names(dt)) return(NULL)
  m <- stringr::str_match(dt$data,
                          "\\[\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*,\\s*(-?\\d+)\\s*\\]")
  out <- data.table::data.table(
    signal = suppressWarnings(as.numeric(m[, 2])),     # step count (element 1)
    epoch  = suppressWarnings(as.numeric(m[, 5])))     # UTC epoch (element 4)
  out[, day := as.POSIXct(epoch, origin = "1970-01-01", tz = "UTC")]
  out <- out[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(signal) & signal >= 0]
  if (!nrow(out)) return(NULL)
  data.table::setkey(out, day); out[, .(day, signal)]
}
read_battery <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !all(c("data","day") %in% names(dt))) return(NULL)
  dt[, day := as.POSIXct(as.character(day), tz = "UTC")]
  dt[, signal := suppressWarnings(as.numeric(data))]
  dt <- dt[!is.na(day) & as.numeric(day) >= SEC_2020 &
             !is.na(signal) & signal >= 0 & signal <= 100]
  if (!nrow(dt)) return(NULL)
  data.table::setkey(dt, day); dt[, .(day, signal)]
}
read_convo <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !"event_data" %in% names(dt)) return(NULL)
  nums  <- stringr::str_extract_all(dt$event_data, "\\d+")
  start <- suppressWarnings(as.numeric(vapply(nums, function(z) z[1], character(1)))) / 1000
  end   <- suppressWarnings(as.numeric(vapply(nums, function(z) z[2], character(1)))) / 1000
  ev <- data.table::data.table(t = start, dur = end - start)
  ev <- ev[!is.na(t) & !is.na(dur) & dur >= 0 & dur <= 24*3600 & t >= SEC_2020]
  if (!nrow(ev)) return(NULL)
  ev[, bin := floor(t / 300) * 300]                         # 5-min epochs
  b <- ev[, .(signal = sum(dur)), by = bin]
  b[, day := as.POSIXct(bin, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(b, day); b[, .(day, signal)]
}
pick_col <- function(nms, patterns) {
  hit <- nms[Reduce(`|`, lapply(patterns, function(p) grepl(p, nms, ignore.case = TRUE)))]
  if (length(hit)) hit[1] else NA_character_
}
extract_uid_values <- function(dt) {
  uid_col <- pick_col(names(dt), c("^uid$", "mlife_id", "participant", "subject", "user"))
  if (is.na(uid_col)) return(NULL)
  map_to_ema_uid(dt[[uid_col]])
}
parse_event_time <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  xn <- suppressWarnings(as.numeric(x))
  if (any(is.finite(xn), na.rm = TRUE)) {
    if (stats::median(xn, na.rm = TRUE) > 1e12) xn <- xn / 1000
    return(as.POSIXct(xn, origin = "1970-01-01", tz = "UTC"))
  }
  suppressWarnings(as.POSIXct(as.character(x), tz = "UTC"))
}
read_sms <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  if (isTRUE(DEDUP_EVENT_LOGS)) {                    # repeated uploads of the same log would otherwise be summed again and again
    n0 <- nrow(dt); dt <- unique(dt)
    cat(sprintf("  [dedup] %s: removed %d exact duplicate rows of %d (%.1f%%)\n", basename(fp), n0 - nrow(dt), n0, 100 * (n0 - nrow(dt)) / n0))
  }
  nms <- names(dt)
  time_col <- pick_col(nms, c("^day$", "timestamp", "time", "date", "created"))
  if (is.na(time_col)) return(NULL)
  body_col <- pick_col(nms, c("body", "message", "text", "content", "^data$"))
  len_col  <- pick_col(nms, c("body_len", "message_len", "text_len", "length", "char", "chars"))
  uid <- extract_uid_values(dt)
  if (!is.null(uid)) dt[, uid := uid]
  dt[, day := parse_event_time(get(time_col))]
  if (!is.na(len_col) && suppressWarnings(any(is.finite(as.numeric(dt[[len_col]])), na.rm = TRUE))) {
    dt[, signal := suppressWarnings(as.numeric(get(len_col)))]
  } else if (!is.na(body_col)) {
    dt[, signal := nchar(as.character(get(body_col)), type = "chars", allowNA = TRUE)]
  } else {
    dt[, signal := 1]
  }
  dt <- dt[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(signal) & signal >= 0]
  if (!nrow(dt)) return(NULL)
  dt[, bin := floor(as.numeric(day) / 1800) * 1800]             # 30-min epochs
  by_cols <- if ("uid" %in% names(dt)) c("uid", "bin") else "bin"
  out <- dt[, .(signal = sum(signal, na.rm = TRUE)), by = by_cols]
  out[, day := as.POSIXct(bin, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day)
  if ("uid" %in% names(out)) out[, .(uid, day, signal)] else out[, .(day, signal)]
}
read_calls <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  if (isTRUE(DEDUP_EVENT_LOGS)) {                    # repeated uploads of the same log would otherwise be summed again and again
    n0 <- nrow(dt); dt <- unique(dt)
    cat(sprintf("  [dedup] %s: removed %d exact duplicate rows of %d (%.1f%%)\n", basename(fp), n0 - nrow(dt), n0, 100 * (n0 - nrow(dt)) / n0))
  }
  nms <- names(dt)
  time_col <- pick_col(nms, c("^day$", "timestamp", "start", "time", "date", "created"))
  dur_col  <- pick_col(nms, c("duration", "dur", "length", "seconds", "sec"))
  end_col  <- pick_col(nms, c("end", "stop"))
  data_col <- pick_col(nms, c("^data$", "event_data"))
  if (is.na(time_col)) return(NULL)
  uid <- extract_uid_values(dt)
  if (!is.null(uid)) dt[, uid := uid]
  dt[, day := parse_event_time(get(time_col))]
  if (!is.na(dur_col) && suppressWarnings(any(is.finite(as.numeric(dt[[dur_col]])), na.rm = TRUE))) {
    dt[, signal := suppressWarnings(as.numeric(get(dur_col)))]
  } else if (!is.na(end_col)) {
    dt[, signal := as.numeric(parse_event_time(get(end_col)) - day, units = "secs")]
  } else if (!is.na(data_col)) {
    nums <- stringr::str_extract_all(as.character(dt[[data_col]]), "\\d+")
    dt[, signal := suppressWarnings(as.numeric(vapply(nums, function(z) {
      if (length(z) >= 2) {
        a <- as.numeric(z[1]); b <- as.numeric(z[2])
        if (max(a, b, na.rm = TRUE) > 1e12) {
          abs(b - a) / 1000
        } else {
          abs(b - a)
        }
      } else if (length(z) == 1) {
        as.numeric(z[1])
      } else {
        NA_real_
      }
    }, numeric(1))))]
  } else {
    dt[, signal := 0]
  }
  dt <- dt[!is.na(day) & as.numeric(day) >= SEC_2020 &
             !is.na(signal) & signal >= 0 & signal <= 24*3600]
  if (!nrow(dt)) return(NULL)
  dt[, bin := floor(as.numeric(day) / 1800) * 1800]             # 30-min epochs
  by_cols <- if ("uid" %in% names(dt)) c("uid", "bin") else "bin"
  out <- dt[, .(signal = sum(signal, na.rm = TRUE)), by = by_cols]
  out[, day := as.POSIXct(bin, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day)
  if ("uid" %in% names(out)) out[, .(uid, day, signal)] else out[, .(day, signal)]
}
# --- GPS: distance travelled per 30-min epoch (speed-gated) --------------------
#   Sums haversine distance between consecutive fixes, but first drops
#   null-island (0,0) fixes and any segment implying an impossible speed
#   (> MAX_GPS_SPEED_MPS). Without these, a handful of corrupt fixes accounted
#   for ~all of each participant's measured distance. shift() replaces the old
#   [-n] indexing.
read_gps <- function(fp) {
  d <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  data.table::setnames(d, toupper(names(d)))
  if (!all(c("LATITUDE", "LONGITUDE", "DAY") %in% names(d))) return(NULL)
  d[, day := as.POSIXct(as.character(DAY), tz = "UTC")]
  d[, lat := suppressWarnings(as.numeric(LATITUDE))]
  d[, lon := suppressWarnings(as.numeric(LONGITUDE))]
  d <- d[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(lat) & !is.na(lon) &
           lat >= -90 & lat <= 90 & lon >= -180 & lon <= 180 &
           !(lat == 0 & lon == 0)]                     # drop null-island (0,0) fixes
  if (nrow(d) < 2) return(NULL)
  data.table::setkey(d, day)
  dist  <- geosphere::distHaversine(
    cbind(data.table::shift(d$lon), data.table::shift(d$lat)),
    cbind(d$lon, d$lat))                       # meters; first fix = NA
  dt_s  <- c(NA_real_, diff(as.numeric(d$day)))         # seconds between fixes
  speed <- dist / dt_s                                  # m/s
  dist[is.na(dist)] <- 0                                # first fix has no predecessor
  dist[is.finite(speed) & speed > MAX_GPS_SPEED_MPS] <- NA   # drop impossible jumps
  d[, dist_m := dist]
  d[, ep := floor(as.numeric(day) / 1800) * 1800]       # 30-min epoch
  out <- d[, .(signal = sum(dist_m, na.rm = TRUE)), by = ep]  # total m / epoch
  out[, day := as.POSIXct(ep, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day); out[, .(day, signal)]
}
# --- GPS TIME-AT-HOME: fraction of fixes within 40m of home x 1800 = sec/epoch --
#   Same coordinate sanity filter as read_gps (valid range + drop (0,0)); without
#   it, corrupt fixes were counted as "away" and deflated time-at-home.
read_gps_home <- function(fp, radius_m = 40) {
  uid_this <- stringr::str_pad(uid_default(fp), 4, pad = "0")
  homes <- HOME_COORDS[[uid_this]]
  if (is.null(homes) || !nrow(homes)) return(NULL)        # no home -> no TAH
  d <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  data.table::setnames(d, toupper(names(d)))
  if (!all(c("LATITUDE", "LONGITUDE", "DAY") %in% names(d))) return(NULL)
  d[, day := as.POSIXct(as.character(DAY), tz = "UTC")]
  d[, lat := suppressWarnings(as.numeric(LATITUDE))]
  d[, lon := suppressWarnings(as.numeric(LONGITUDE))]
  d <- d[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(lat) & !is.na(lon) &
           lat >= -90 & lat <= 90 & lon >= -180 & lon <= 180 &
           !(lat == 0 & lon == 0)]                     # drop null-island (0,0) fixes
  if (!nrow(d)) return(NULL)
  pts <- cbind(d$lon, d$lat)
  is_home <- rep(FALSE, nrow(d))
  for (h in seq_len(nrow(homes)))
    is_home <- is_home | (geosphere::distHaversine(pts, c(homes[h, "lon"], homes[h, "lat"])) <= radius_m)
  d[, ih := as.integer(is_home)]
  d[, ep := floor(as.numeric(day) / 1800) * 1800]
  out <- d[, .(signal = mean(ih) * 1800), by = ep]        # seconds at home / epoch
  out[, day := as.POSIXct(ep, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day); out[, .(day, signal)]
}
# ==============================================================================
# EVENT STREAMS (SMS, calls, conversation): a bin with no event is a ZERO, not missing.
# Fills every bin between the participant's first and last event with 0. Assumes the phone was
# logging throughout that span (state this in the Methods); set ZERO_FILL_EVENT_STREAMS <- FALSE
# to keep the old behaviour.
# ==============================================================================
zero_fill_events <- function(dt, bin_sec) {
  if (is.null(dt) || !nrow(dt)) return(dt)
  tt   <- as.numeric(dt$day)
  grid <- data.table::data.table(t = seq(min(tt), max(tt), by = bin_sec))
  obs  <- data.table::data.table(t = tt, signal = dt$signal)
  obs  <- obs[, .(signal = sum(signal, na.rm = TRUE)), by = t]
  out  <- merge(grid, obs, by = "t", all.x = TRUE)
  out[is.na(signal), signal := 0]
  out[, day := as.POSIXct(t, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day)
  out[, .(day, signal)]
}
# ==============================================================================
# CORE: compute dynamic indices for one participant's (day, signal) series
# ==============================================================================
# bin_sec: width of the bin that each `day` timestamp labels (0 for instantaneous readings).
# Bins are labelled by their START. The cutoff is snapped down to the last bin boundary
# (ae_b = floor(anchor_end / bin_sec) * bin_sec) and only bins that END at or before ae_b are used,
# so every bin is complete before t - 4 h and each window keeps its full number of bins.
compute_indices <- function(day_vec, sig_vec, ema_this, epoch_sec, bin_sec = 0,
                            entropy_ok = TRUE, windows = WINDOWS) {
  out <- vector("list", nrow(ema_this) * length(windows)); k <- 1L
  for (j in seq_len(nrow(ema_this))) {
    ae <- as.numeric(ema_this$anchor_end[j])
    # Snap the cutoff DOWN to the last complete bin boundary, so each window keeps its full
    # number of complete bins (e.g. 2 half-hour bins for hr1). bin_sec = 0 (instantaneous readings): no snapping.
    ae_b    <- if (bin_sec > 0) floor(ae / bin_sec) * bin_sec else ae
    cur_len <- max(1800, 2 * bin_sec)
    cl <- findInterval(ae_b - cur_len, day_vec); ch <- findInterval(ae_b - bin_sec, day_vec)
    x_current <- mean_thr(if (ch > cl) sig_vec[(cl + 1L):ch] else numeric(0))
    for (wn in names(windows)) {
      ws <- windows[[wn]]
      wl <- if (bin_sec > 0) findInterval(ae_b - ws, day_vec, left.open = TRUE) else findInterval(ae - ws, day_vec)
      wh <- findInterval(ae_b - bin_sec, day_vec)
      xc <- { x <- if (wh > wl) sig_vec[(wl + 1L):wh] else numeric(0); x[!is.na(x)] }
      n_obs <- length(xc)
      cov   <- if (!is.na(epoch_sec)) min(n_obs / (ws / epoch_sec), 1) else NA_real_
      mean_v <- mean_thr(xc); sum_v <- if (n_obs) sum(xc) else NA_real_
      mad_v <- mad_thr(xc); p_v <- p90p10_relpos(xc, x_current)
      ent_v <- if (isTRUE(entropy_ok)) sample_entropy(xc) else NA_real_
      out[[k]] <- list(
        uid = ema_this$uid[j], ema_day = ema_this$ema_day[j],
        date_rep = ema_this$date_rep[j], ema_time = ema_this$ema_time[j],
        anchor_end = ema_this$anchor_end[j], window = wn,
        n_obs = n_obs, coverage = cov, mean_val = mean_v, sum_val = sum_v,
        mad_val = mad_v, p90p10 = p_v, entropy = ent_v)
      k <- k + 1L
    }
  }
  data.table::rbindlist(out)
}
process_one <- function(fp, spec) {
  uid_this <- stringr::str_pad(as.character(spec$uid_fn(fp)), width = 4, pad = "0")
  out_dir  <- file.path(CLEAN_UID_DIR, paste0(spec$tag, "_per_uid"))
  out_path <- file.path(out_dir, sprintf("%s_dynamic_indices_%s.csv", uid_this, spec$tag))
  tryCatch({
    dt <- spec$read(fp)
    if (is.null(dt) || !nrow(dt)) return(paste("EMPTY", uid_this))
    if (isTRUE(spec$zero_fill)) dt <- zero_fill_events(dt, spec$bin_sec)
    ema_this <- ema_lookup[uid == uid_this]
    if (!nrow(ema_this)) return(paste("NOEMA", uid_this))
    min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
    max_t <- as.numeric(max(ema_this$anchor_end))
    dt <- dt[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
    if (!nrow(dt)) return(paste("NORANGE", uid_this))
    data.table::setkey(dt, day)
    res <- compute_indices(as.numeric(dt$day), dt$signal, ema_this, spec$epoch_sec,
                           bin_sec = spec$bin_sec, entropy_ok = !isFALSE(spec$entropy_ok))
    make_dirs(out_dir)
    data.table::fwrite(res, out_path)
    paste("DONE", uid_this, nrow(ema_this))
  }, error = function(e) paste("ERROR", uid_this, conditionMessage(e)))
}
sensor_search_roots <- function() {
  roots <- c(RAW_DATA_DIR, CLEAN_SENSOR_DIR, CLEAN_DATA_DIR, BASE_PATH)
  unique(roots[!is.na(roots) & nzchar(roots) & dir.exists(roots)])
}
sensor_path_candidates <- function(folder) {
  if (grepl("^/", folder) || grepl("^[A-Za-z]:", folder)) {
    direct <- folder
    with_ext <- if (grepl("\\.csv$", folder, ignore.case = TRUE)) direct else paste0(direct, ".csv")
    return(unique(c(direct, with_ext)))
  }
  unique(unlist(lapply(sensor_search_roots(), function(root) {
    direct <- file.path(root, folder)
    with_ext <- if (grepl("\\.csv$", folder, ignore.case = TRUE)) direct else paste0(direct, ".csv")
    c(direct, with_ext)
  }), use.names = FALSE))
}
case_insensitive_csv_matches <- function(folder) {
  stem <- tools::file_path_sans_ext(basename(folder))
  roots <- sensor_search_roots()
  unique(unlist(lapply(roots, function(root) {
    if (!dir.exists(root)) return(character(0))
    list.files(root, pattern = paste0("^", stem, "\\.csv$"),
               full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
}
list_sensor_files <- function(spec) {
  candidates <- unique(c(
    unlist(lapply(spec$folder, sensor_path_candidates), use.names = FALSE),
    unlist(lapply(spec$folder, case_insensitive_csv_matches), use.names = FALSE)
  ))
  files <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  dirs <- candidates[dir.exists(candidates)]
  if (length(dirs)) {
    files <- c(files, unlist(lapply(dirs, list.files, pattern = "\\.csv$", full.names = TRUE),
                             use.names = FALSE))
  }
  unique(files)
}
# ==============================================================================
# run_full(): parallel scale-up for one sensor
# ==============================================================================
run_full <- function(spec, workers = WORKERS) {
  final_path <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", spec$tag))
  unlink(final_path)
  unlink(file.path(CLEAN_UID_DIR, paste0(spec$tag, "_per_uid")), recursive = TRUE)
  files <- list_sensor_files(spec)
  if (!length(files)) {
    cat(sprintf("%s: no CSV files found in %s — skipping\n",
                spec$tag, paste(spec$folder, collapse = ", ")))
    return(invisible(NULL))
  }
  if (isTRUE(spec$combined_file)) {
    cat(sprintf("%s: %d combined CSV file(s)\n", spec$tag, length(files)))
    dt <- data.table::rbindlist(lapply(files, spec$read), use.names = TRUE, fill = TRUE)
    if (!nrow(dt)) {
      cat(sprintf("%s: reader returned no event rows\n", spec$tag))
      return(invisible(NULL))
    }
    if (!all(c("uid", "day", "signal") %in% names(dt))) {
      stop(sprintf("%s reader must return uid, day, and signal for combined files. Run inspect_sensor_input_files() and check the ID column name.",
                   spec$tag))
    }
    dt <- dt[!is.na(uid) & !is.na(day) & !is.na(signal)]
    dt[, uid := map_to_ema_uid(uid)]
    data.table::setkey(dt, uid, day)
    uids <- intersect(unique(dt$uid), unique(ema_lookup$uid))
    cat(sprintf("%s: %d PTs with events and EMA anchors\n", spec$tag, length(uids)))
    out_list <- vector("list", length(uids))
    for (i in seq_along(uids)) {
      uid_this <- uids[i]
      dt_this <- dt[uid == uid_this]
      if (isTRUE(spec$zero_fill)) dt_this <- zero_fill_events(dt_this, spec$bin_sec)
      ema_this <- ema_lookup[uid == uid_this]
      min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
      max_t <- as.numeric(max(ema_this$anchor_end))
      dt_this <- dt_this[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
      if (!nrow(dt_this)) next
      data.table::setkey(dt_this, day)
      out_list[[i]] <- compute_indices(as.numeric(dt_this$day), dt_this$signal,
                                       ema_this, spec$epoch_sec,
                                       bin_sec = spec$bin_sec, entropy_ok = !isFALSE(spec$entropy_ok))
    }
    final <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
    if (!nrow(final)) {
      cat(sprintf("%s: no rows overlapped EMA windows\n", spec$tag))
      return(invisible(NULL))
    }
    final[, uid := map_to_ema_uid(uid)]
    data.table::fwrite(final, final_path)
    cat(sprintf("Saved dynamic_indices_%s.csv: %d rows, %d PTs\n",
                spec$tag, nrow(final), data.table::uniqueN(final$uid)))
    return(invisible(final))
  }
  uids  <- stringr::str_pad(as.character(vapply(files, spec$uid_fn, character(1))), width = 4, pad = "0")
  if (anyDuplicated(uids))
    cat(sprintf("WARN %s: %d files, %d unique uids — dups would be skipped\n",
                spec$tag, length(files), length(unique(uids))))
  cat(sprintf("%s: %d files | %d workers\n", spec$tag, length(files), workers))
  data.table::setDTthreads(1)
  logs <- pbmcapply::pbmclapply(files, process_one, spec = spec,
                                mc.cores = workers, mc.preschedule = FALSE)
  data.table::setDTthreads(0)
  print(table(sub(" .*", "", unlist(logs))))
  errs <- grep("^ERROR", unlist(logs), value = TRUE); if (length(errs)) print(head(errs, 10))
  out_dir <- file.path(CLEAN_UID_DIR, paste0(spec$tag, "_per_uid"))
  done <- list.files(out_dir, pattern = sprintf("_dynamic_indices_%s\\.csv$", spec$tag), full.names = TRUE)
  final <- data.table::rbindlist(lapply(done, function(f)
    data.table::fread(f, colClasses = list(character = "uid"))))
  final[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  data.table::fwrite(final, final_path)
  cat(sprintf("Saved dynamic_indices_%s.csv: %d rows, %d PTs\n",
              spec$tag, nrow(final), data.table::uniqueN(final$uid)))
  invisible(final)
}
# ==============================================================================
# SLEEP
# ==============================================================================
NIGHT_OFFSET_DAYS <- 0L   # calendarDate = wake date (per Garmin documentation)
run_sleep <- function(folder = "sleep") {
  final_path <- file.path(CLEAN_SENSOR_DIR, "dynamic_indices_sleep.csv")
  unlink(final_path)
  sdir <- file.path(RAW_DATA_DIR, folder)
  files <- list.files(sdir, pattern = "\\.csv$", full.names = TRUE)
  cat(sprintf("sleep: %d files\n", length(files)))
  ema_for_sleep <- build_ema_lookup()
  if (!"uid" %in% names(ema_for_sleep)) {
    stop("EMA lookup for sleep is missing uid; rebuild ema_features.csv before run_sleep().")
  }
  ema_dates <- unique(ema_for_sleep[, .(uid, ema_day, date_rep, ema_time)])
  ad <- suppressWarnings(as.Date(as.character(ema_dates$date_rep)))            # ISO
  if (all(is.na(ad))) ad <- as.Date(as.character(ema_dates$date_rep), format = "%m/%d/%Y")
  ema_dates[, assess_date := ad]
  
  out_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    fp <- files[i]
    uid_this <- stringr::str_pad(as.character(uid_default(fp)), width = 4, pad = "0")
    s <- tryCatch(data.table::fread(fp), error = function(e) NULL)
    need <- c("calendarDate","durationInSeconds","remSleepInSeconds")
    if (is.null(s) || !all(need %in% names(s))) next
    sl <- data.table::data.table(
      uid = uid_this, sleep_date = as.Date(s$calendarDate),
      total_sleep = as.numeric(s$durationInSeconds) / 60,
      rem_sleep   = as.numeric(s$remSleepInSeconds) / 60)
    sl <- sl[!is.na(total_sleep)]
    sl <- sl[, .(total_sleep = mean(total_sleep),
                 rem_sleep   = mean(rem_sleep, na.rm = TRUE)),
             by = .(uid, sleep_date)]
    em <- ema_dates[uid == uid_this]; if (!nrow(em)) next
    em[, match_date := assess_date - NIGHT_OFFSET_DAYS]
    out_list[[i]] <- merge(em, sl, by.x = "match_date", by.y = "sleep_date",
                           all.x = TRUE, suffixes = c("", "_sleep"))[
                             , .(uid, ema_day, date_rep, ema_time, total_sleep, rem_sleep)]
  }
  final <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  data.table::fwrite(final, final_path)
  cat(sprintf("Saved dynamic_indices_sleep.csv: %d rows, %d PTs | %.0f%% nights matched\n",
              nrow(final), data.table::uniqueN(final$uid),
              100 * mean(!is.na(final$total_sleep))))
  invisible(final)
}
# ==============================================================================
# SENSOR SPECS
# ==============================================================================
SENSORS <- list(
  hr      = list(tag = "hr",      folder = "hr",      read = read_hr,      uid_fn = uid_default, epoch_sec = 10,   bin_sec = 0),
  steps   = list(tag = "steps",   folder = "steps",   read = read_steps,   uid_fn = uid_default, epoch_sec = 30,   bin_sec = 0),
  battery = list(tag = "battery", folder = "battery", read = read_battery, uid_fn = uid_default, epoch_sec = 3600, bin_sec = 0),
  convo   = list(tag = "convo",   folder = "convo",        read = read_convo,   uid_fn = uid_convo,   epoch_sec = NA_real_, bin_sec = 300,
                 zero_fill = ZERO_FILL_EVENT_STREAMS, entropy_ok = !SUPPRESS_ENTROPY_EVENT_STREAMS),
  gps      = list(tag = "gps",      folder = "raw_gps_clean", read = read_gps,      uid_fn = uid_default, epoch_sec = 1800, bin_sec = 1800),
  gps_home = list(tag = "gps_home", folder = "raw_gps_clean", read = read_gps_home, uid_fn = uid_default, epoch_sec = 1800, bin_sec = 1800),
  sms      = list(tag = "sms",      folder = c("smslog.csv", "smslog", "sms.csv", "sms"),
                  read = read_sms, uid_fn = uid_default, epoch_sec = NA_real_, bin_sec = 1800, combined_file = TRUE,
                  zero_fill = ZERO_FILL_EVENT_STREAMS, entropy_ok = !SUPPRESS_ENTROPY_EVENT_STREAMS),
  calls    = list(tag = "calls",    folder = c("calllog.csv", "calllog", "calls.csv", "calls"),
                  read = read_calls, uid_fn = uid_default, epoch_sec = NA_real_, bin_sec = 1800, combined_file = TRUE,
                  zero_fill = ZERO_FILL_EVENT_STREAMS, entropy_ok = !SUPPRESS_ENTROPY_EVENT_STREAMS)
)
# --- process every sensor -------------------------------------------------------
run_full(SENSORS$hr)
run_full(SENSORS$steps)
run_full(SENSORS$battery)
run_full(SENSORS$convo)
run_full(SENSORS$gps)
run_full(SENSORS$gps_home)
run_full(SENSORS$sms)
run_full(SENSORS$calls)
run_sleep()
# ==============================================================================
# DATA CLEANING, 4 OF 4: FEATURE ASSEMBLY + ANALYTIC INCLUSION
# ==============================================================================
make_dirs(results_dir)
reshape_dynamic_indices <- function(path, tag) {
  d <- data.table::fread(path, colClasses = list(character = "uid"))
  d[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  key_cols <- c("uid", "date_rep", "ema_time")
  if (!all(key_cols %in% names(d))) {
    warning(sprintf("%s is missing merge keys; skipping", basename(path)))
    return(NULL)
  }
  
  if (tag == "sleep") {
    # calendarDate = WAKE date. Features are built at the level of NIGHTS (one row per
    # uid x night), then joined back to prompts. Fixes in this version:
    #  * day7 mean is a true mean over the previous 7 nights (it used to equal the single-night value);
    #  * MAD is computed over distinct nights (it used to see 3 copies of each night, so one prior
    #    night gave N = 3 and MAD = 0 instead of NA);
    #  * SLEEP_STRICT_PRIOR_NIGHTS = FALSE (preregistered): the night before the assessment (wake date = prompt date) is
    #    included. TRUE: only nights BEFORE the prompt date, so every
    #    feature respects the t - 4 h cutoff (a same-day wake can fall after t - 4 h for early prompts).
    out <- d[, .(total_sleep = mean(total_sleep, na.rm = TRUE),
                 rem_sleep   = mean(rem_sleep,   na.rm = TRUE)), by = key_cols]
    out[, date_for_sleep := parse_date_any(date_rep)]
    nights <- unique(out[, .(uid, nd = date_for_sleep, ts = total_sleep, rs = rem_sleep)],
                     by = c("uid", "nd"))
    data.table::setorder(nights, uid, nd)
    upper_incl <- !SLEEP_STRICT_PRIOR_NIGHTS      # strict: window = nights D-7 .. D-1; preregistered (FALSE): nights D-6 .. D
    roll_feat <- function(dd, x, f) vapply(seq_along(dd), function(i) {
      hi <- if (upper_incl) dd[i] else dd[i] - 1L
      v  <- x[dd <= hi & dd >= hi - 6L]
      f(v[!is.na(v)])
    }, numeric(1))
    last_feat <- function(dd, x) vapply(seq_along(dd), function(i) {
      tgt <- if (upper_incl) dd[i] else dd[i] - 1L
      v <- x[dd == tgt]
      if (length(v)) v[1] else NA_real_
    }, numeric(1))
    nights[, `:=`(
      sleep_mean_total_sleep_day1 = last_feat(as.integer(nd), ts),
      sleep_mean_total_sleep_day7 = roll_feat(as.integer(nd), ts, mean_thr),
      sleep_mad_total_sleep_day7  = roll_feat(as.integer(nd), ts, mad_thr),
      sleep_mean_rem_sleep_day1   = last_feat(as.integer(nd), rs),
      sleep_mean_rem_sleep_day7   = roll_feat(as.integer(nd), rs, mean_thr),
      sleep_mad_rem_sleep_day7    = roll_feat(as.integer(nd), rs, mad_thr)
    ), by = uid]
    nights[, c("ts", "rs") := NULL]
    out[, c("total_sleep", "rem_sleep") := NULL]
    out <- merge(out, nights, by.x = c("uid", "date_for_sleep"), by.y = c("uid", "nd"), all.x = TRUE)
    out[, date_for_sleep := NULL]
    cat(sprintf("Sleep reshape columns: %s\n",
                paste(grep("^sleep_", names(out), value = TRUE), collapse = ", ")))
    chk <- out[, .(same = isTRUE(all.equal(sleep_mean_total_sleep_day1, sleep_mean_total_sleep_day7)))]
    if (isTRUE(chk$same)) warning("sleep day1 and day7 means are identical; check the sleep reshape.")
    return(unique(out, by = key_cols))
  }
  
  if (!"window" %in% names(d)) {
    warning(sprintf("%s is missing window; skipping", basename(path)))
    return(NULL)
  }
  measures <- intersect(c("coverage", "n_obs", "mean_val", "sum_val",
                          "mad_val", "p90p10", "entropy"), names(d))
  if (!length(measures)) return(NULL)
  d[, (measures) := lapply(.SD, as.numeric), .SDcols = measures]
  long <- data.table::melt(
    d,
    id.vars = c(key_cols, "window"),
    measure.vars = measures,
    variable.name = "measure",
    value.name = "value"
  )
  long[, feature := paste(tag, window, measure, sep = "_")]
  data.table::dcast(long, uid + date_rep + ema_time ~ feature,
                    value.var = "value", fun.aggregate = mean, na.rm = TRUE)
}
assemble_feature_matrix_local <- function(outfile = "all_features_combined.csv") {
  base <- data.table::fread(file.path(CLEAN_DATA_DIR, "ema_features.csv"), colClasses = list(character = "uid"))
  base[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  key_cols <- c("uid", "date_rep", "ema_time")
  if (!all(key_cols %in% names(base))) {
    stop("ema_features.csv must contain uid, date_rep, and ema_time.")
  }
  
  tags <- c("acc", "hr", "steps", "battery", "convo", "gps", "gps_home", "sleep", "sms", "calls")
  for (tag in tags) {
    fp <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", tag))
    if (!file.exists(fp)) {
      cat(sprintf("Missing dynamic_indices_%s.csv — skipping %s columns\n", tag, tag))
      next
    }
    wide <- reshape_dynamic_indices(fp, tag)
    if (is.null(wide)) next
    add_cols <- setdiff(setdiff(names(wide), key_cols), names(base))
    if (length(add_cols)) {
      base <- merge(base, wide[, c(key_cols, add_cols), with = FALSE],
                    by = key_cols, all.x = TRUE, sort = FALSE)
      cat(sprintf("Merged %s: %d columns\n", tag, length(add_cols)))
    } else {
      cat(sprintf("%s columns already present or empty\n", tag))
    }
  }
  data.table::fwrite(base, file.path(CLEAN_DATA_DIR, outfile))
  cat(sprintf("Saved %s: %d rows x %d cols (%d participants)\n",
              outfile, nrow(base), ncol(base), data.table::uniqueN(base$uid)))
  invisible(base)
}
assemble_feature_matrix_local()
validate_feature_matrix <- function(path = file.path(CLEAN_DATA_DIR, "all_features_combined.csv")) {
  d <- data.table::fread(path, nrows = 1)
  core_prefixes <- c("acc_", "hr_", "steps_", "battery_", "convo_", "gps_",
                     "gps_home_", "sleep_")
  optional_prefixes <- c("sms_", "calls_")
  core_present <- vapply(core_prefixes, function(p) any(grepl(paste0("^", p), names(d))), logical(1))
  optional_present <- vapply(optional_prefixes, function(p) any(grepl(paste0("^", p), names(d))), logical(1))
  if (any(!core_present)) {
    warning(sprintf("Feature matrix is missing expected core prefixes: %s",
                    paste(core_prefixes[!core_present], collapse = ", ")))
  }
  if (any(!optional_present)) {
    message(sprintf("Optional event-log prefixes not present: %s",
                    paste(optional_prefixes[!optional_present], collapse = ", ")))
  }
  invisible(c(core_present, optional_present))
}
validate_feature_matrix()
dat <- read.csv(file.path(CLEAN_DATA_DIR, "all_features_combined.csv"))
dynamic_sensor_days <- function(sensor) {
  fp <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", sensor))
  if (!file.exists(fp)) {
    return(data.table::data.table(uid = character(), sensor = sensor, date = as.Date(character())))
  }
  cols <- names(data.table::fread(fp, nrows = 0))
  select_cols <- intersect(c("uid", "date_rep", "n_obs", "total_sleep", "rem_sleep"), cols)
  d <- data.table::fread(fp, select = select_cols, colClasses = list(character = "uid"))
  if (!all(c("uid", "date_rep") %in% names(d))) {
    return(data.table::data.table(uid = character(), sensor = sensor, date = as.Date(character())))
  }
  d[, uid := normalize_uid(uid)]
  d[, date := parse_date_any(date_rep)]
  if ("n_obs" %in% names(d)) {
    d[, n_obs := suppressWarnings(as.numeric(n_obs))]
    d <- d[!is.na(n_obs) & n_obs > 0]
  } else if (sensor == "sleep") {
    d <- d[(!is.na(total_sleep)) | (!is.na(rem_sleep))]
  }
  d <- unique(d[!is.na(uid) & !is.na(date), .(uid, date)])
  if (!nrow(d)) return(data.table::data.table(uid = character(), sensor = sensor, date = as.Date(character())))
  d[, sensor := sensor]
  d[, .(uid, sensor, date)]
}
compute_raw_passive_compliance <- function(expected_days,
                                           passive_sensors = c("acc","steps","hr","battery","gps","gps_home"),
                                           outfile = file.path(results_dir,
                                                               "raw_passive_compliance_by_uid_sensor.csv")) {
  make_dirs(dirname(outfile))
  expected_days[, uid := normalize_uid(uid)]
  expected_days <- unique(expected_days[!is.na(uid) & !is.na(date), .(uid, date)])
  uid_expected <- expected_days[, .(expected_days = data.table::uniqueN(date)), by = uid]
  valid_days <- data.table::rbindlist(lapply(passive_sensors, function(s) {
    cat(sprintf("Computing passive compliance days: %s\n", s))
    dynamic_sensor_days(s)
  }), use.names = TRUE, fill = TRUE)
  valid_days <- unique(valid_days[uid %in% uid_expected$uid & !is.na(date), .(uid, sensor, date)])
  valid_counts <- valid_days[, .(valid_days = data.table::uniqueN(date)), by = .(uid, sensor)]
  
  grid <- data.table::CJ(uid = uid_expected$uid, sensor = passive_sensors, unique = TRUE)
  grid <- uid_expected[grid, on = "uid"]
  grid <- valid_counts[grid, on = c("uid", "sensor")]
  grid[is.na(valid_days), valid_days := 0L]
  grid[, passive_compliance := valid_days / expected_days]
  data.table::fwrite(grid, outfile)
  grid[]
}
################################################################################
# ANALYTIC INCLUSION + FEATURE-DROP CRITERIA
################################################################################
apply_inclusion_criteria <- function(
    infile = "all_features_combined.csv",
    outfile = "all_features_cleaned.csv",
    ema_compliance_threshold = 0.75,
    passive_data_threshold = 0.75,
    passive_sensors = c("acc","steps","hr","battery","gps","gps_home"),
    drop_sensors = c("stress"),
    min_si_pos = 10, min_si_neg = 20,
    min_folds = 10,
    train_days = 45, val_days = 10, test_days = 1, ema_per_day = 3,
    max_missing = 0.50) {
  
  d <- data.table::fread(file.path(CLEAN_DATA_DIR, infile))
  d[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  n_uid0 <- data.table::uniqueN(d$uid)
  cat(sprintf("Input: %d rows, %d cols, %d participants\n", nrow(d), ncol(d), n_uid0))
  
  for (s in drop_sensors) {
    cols <- grep(sprintf("^%s_", s), names(d), value = TRUE)
    if (length(cols)) d[, (cols) := NULL]
  }
  
  ema_comp <- d[, .(ema_compliance = max(suppressWarnings(as.numeric(compliance)), na.rm = TRUE)), by = uid]
  ema_comp[!is.finite(ema_compliance), ema_compliance := NA_real_]
  keep_ema <- ema_comp[ema_compliance >= ema_compliance_threshold, uid]
  
  expected_days <- unique(d[, .(uid, date = parse_date_any(date_rep))])
  uid_sensor <- compute_raw_passive_compliance(expected_days, passive_sensors)
  
  dead_sensors <- uid_sensor[, .(has_data = any(valid_days > 0)), by = sensor][
    has_data == FALSE, sensor]
  if (length(dead_sensors)) {
    cat(sprintf("NOTE: passive stream(s) with no data in the whole sample -> dropped from compliance and features: %s\n",
                paste(dead_sensors, collapse = ", ")))
    uid_sensor <- uid_sensor[!sensor %in% dead_sensors]
    for (s in dead_sensors) {
      cols <- grep(sprintf("^%s_(hr1|hr4|day1|day7)_", s), names(d), value = TRUE)
      if (length(cols)) d[, (cols) := NULL]
    }
  }
  
  uid_cov <- uid_sensor[, .(mean_cov = mean(passive_compliance, na.rm = TRUE)), by = uid]
  keep_cov <- uid_cov[is.finite(mean_cov) & mean_cov >= passive_data_threshold, uid]
  
  si <- d[!is.na(q9), .(n_pos = sum(q9 > 0), n_neg = sum(q9 == 0)), by = uid]
  keep_si <- si[n_pos >= min_si_pos & n_neg >= min_si_neg, uid]
  
  si_var <- d[!is.na(q9), .(q9_var = stats::var(q9)), by = uid]
  keep_var <- si_var[is.finite(q9_var) & q9_var > 0, uid]
  
  fold_n  <- (train_days + val_days + test_days) * ema_per_day   # 168
  slide_n <- test_days * ema_per_day                             # 3
  min_obs <- fold_n + (min_folds - 1) * slide_n                  # 195
  n_obs_pt <- d[!is.na(q9), .(n_obs = .N), by = uid]
  keep_obs <- n_obs_pt[n_obs >= min_obs, uid]
  
  report <- Reduce(function(a, b) merge(a, b, by = "uid", all = TRUE),
                   list(ema_comp, uid_cov, si, si_var, n_obs_pt))
  report <- data.table::as.data.table(report)
  report[, pass_ema     := uid %in% keep_ema]
  report[, pass_passive := uid %in% keep_cov]
  report[, pass_si      := uid %in% keep_si]
  report[, pass_var     := uid %in% keep_var]
  report[, pass_folds   := uid %in% keep_obs]
  make_dirs(results_dir)
  data.table::fwrite(report, file.path(results_dir, "Table1_inclusion_participant_report.csv"))
  
  cat(sprintf("EMA >= %.0f%%: %d | passive >= %.0f%%: %d | SI counts: %d | non-zero var: %d | >=%d folds (>=%d obs): %d  (of %d PTs)\n",
              ema_compliance_threshold * 100, length(keep_ema),
              passive_data_threshold * 100, length(keep_cov),
              length(keep_si), length(keep_var),
              min_folds, min_obs, length(keep_obs), n_uid0))
  
  keep <- Reduce(intersect, list(keep_ema, keep_cov, keep_si, keep_var, keep_obs))
  d <- d[uid %in% keep]
  cat(sprintf("==> %d participants retained\n", length(keep)))
  
  # 7. drop feature columns >max_missing missing -------------------------------
  protected <- c("uid", "date_rep", "ema_time", "mlife_id", "compliance",
                 "ema_timestamp", "ema_day", paste0("q", 1:9), "g1", "g2")
  feat <- setdiff(names(d), protected)
  pct_missing <- vapply(feat, function(cc) mean(is.na(d[[cc]])), numeric(1))
  drop_cols <- feat[pct_missing > max_missing]
  
  if (length(drop_cols)) {
    miss_report <- data.frame(column = drop_cols,
                              pct_missing = round(100 * pct_missing[drop_cols], 1),
                              row.names = NULL)
    miss_report <- miss_report[order(-miss_report$pct_missing), ]
    data.table::fwrite(miss_report,
                       file.path(results_dir, "dropped_high_missing_columns.csv"))
    d[, (drop_cols) := NULL]
  }
  
  # 7b. drop feature columns with no variability in the analytic sample -------
  feat_zv <- setdiff(names(d), protected)
  zv <- vapply(feat_zv, function(cc) {
    v <- suppressWarnings(stats::var(as.numeric(d[[cc]]), na.rm = TRUE))
    !is.finite(v) || v == 0
  }, logical(1))
  drop_zv <- feat_zv[zv]
  if (length(drop_zv)) {
    data.table::fwrite(data.frame(column = drop_zv),
                       file.path(results_dir, "dropped_no_variability_columns.csv"))
    d[, (drop_zv) := NULL]
    cat(sprintf("Dropped %d feature column(s) with zero variability in the analytic sample.\n",
                length(drop_zv)))
  }
  
  d[, .ord := match(ema_time, c("morning", "afternoon", "evening"))]
  data.table::setorderv(d, c("uid", "date_rep", ".ord"))
  d[, .ord := NULL]
  
  data.table::fwrite(d, file.path(CLEAN_DATA_DIR, outfile))
  cat(sprintf("Saved %s: %d rows x %d cols (%d participants) | dropped %d cols >%.0f%% missing\n",
              outfile, nrow(d), ncol(d), data.table::uniqueN(d$uid),
              length(drop_cols), max_missing * 100))
  invisible(list(data = d, kept_uids = keep, dropped_cols = drop_cols))
}
res <- apply_inclusion_criteria()
cat("\nData cleaning complete.\n")
dat <- res$data
if (interactive()) View(dat)
nrow(dplyr::distinct(dat, uid))
# ==============================================================================
# DATA CLEANING, 5 OF 5: FINAL CLEANUP
# ==============================================================================
final_output <- file.path(CLEAN_DATA_DIR, "all_features_cleaned.csv")
if (file.exists(final_output)) {
  intermediate_files <- file.path(
    CLEAN_DATA_DIR,
    c("compliant_uids.csv", "ema_features.csv", "all_features_combined.csv")
  )
  intermediate_dirs <- c(CLEAN_SENSOR_DIR, CLEAN_UID_DIR)
  # The three small report files in results_dir are KEPT on purpose: the per-sensor passive compliance
  # (Table 1 / coverage rule) and the lists of dropped feature columns are needed for the Methods.
  
  unlink(intermediate_files)
  unlink(intermediate_dirs, recursive = TRUE)
  
  for (dd in c(results_dir, RESULTS_ROOT_DIR)) {
    if (dir.exists(dd) &&
        length(list.files(dd, all.files = TRUE, no.. = TRUE)) == 0)
      unlink(dd, recursive = TRUE)
  }
  
  cat(sprintf("Cleanup done. Retained:\n  %s\n  %s  (manuscript Table 1)\n  plus the report files in %s\n",
              final_output,
              file.path(results_dir, "Table1_inclusion_participant_report.csv"), results_dir))
} else {
  warning("all_features_cleaned.csv not found — skipping cleanup so nothing is lost.")
}

# ------------------------------------------------------------------------------
# Record the environment and the design switches that produced this file
# ------------------------------------------------------------------------------
try(writeLines(c(
  sprintf("Run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("ZERO_FILL_EVENT_STREAMS = %s", ZERO_FILL_EVENT_STREAMS),
  sprintf("SUPPRESS_ENTROPY_EVENT_STREAMS = %s", SUPPRESS_ENTROPY_EVENT_STREAMS),
  sprintf("SLEEP_STRICT_PRIOR_NIGHTS = %s", SLEEP_STRICT_PRIOR_NIGHTS),
  sprintf("DEDUP_EVENT_LOGS = %s", DEDUP_EVENT_LOGS),
  "", capture.output(sessionInfo())),
  file.path(CLEAN_DATA_DIR, "cleaning_sessionInfo.txt")), silent = TRUE)