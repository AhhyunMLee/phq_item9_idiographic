################################################################################
# SI Idiographic Prediction Project - Final Single R Script
#
# Folder layout expected under BASE_PATH:
#   1. Code
#   2. Raw data
#   3. Cleaned data
#   4. Results
#
# To run from another R session, set flags before sourcing this file, for example:
#   Sys.setenv(SI_PROJECT_BASE = "/path/to/SI_Prediction")
#   RUN_MAIN_PIPELINE <- TRUE
#   source("si_idiographic_2026.R")
################################################################################

set_default_flag <- function(name, value) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv, inherits = FALSE))
  }
  assign(name, value, envir = .GlobalEnv)
  value
}

require_packages_or_stop <- function(packages, context = "this analysis") {
  missing <- packages[!(packages %in% rownames(installed.packages()))]
  if (length(missing)) {
    stop(
      sprintf(
        "Missing R packages required for %s: %s\nInstall them with:\ninstall.packages(c(%s))",
        context,
        paste(missing, collapse = ", "),
        paste(sprintf('"%s"', missing), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Create directories quietly, including missing parent folders.
make_dirs <- function(...) {
  invisible(lapply(c(...), dir.create, showWarnings = FALSE, recursive = TRUE))
}

# Public run flags. Existing values set before source() are preserved.
RUN_DATA_CLEANING <- set_default_flag("RUN_DATA_CLEANING", FALSE)
RUN_MAIN_PIPELINE <- set_default_flag("RUN_MAIN_PIPELINE", FALSE)
RUN_FINAL_MANUSCRIPT_OUTPUTS <- set_default_flag("RUN_FINAL_MANUSCRIPT_OUTPUTS", FALSE)
RUN_FEATURE_IMPORTANCE <- set_default_flag("RUN_FEATURE_IMPORTANCE", FALSE)
RUN_ALL_PT_CONTINUOUS_TIMESERIES <- set_default_flag("RUN_ALL_PT_CONTINUOUS_TIMESERIES", FALSE)
RUN_ALL_PT_HURDLE_TIMESERIES <- set_default_flag("RUN_ALL_PT_HURDLE_TIMESERIES", FALSE)
KEEP_FEATURE_IMPORTANCE_DETAILS <- set_default_flag("KEEP_FEATURE_IMPORTANCE_DETAILS", FALSE)
RUN_EXTRA_DIAGNOSTIC_FIGURES <- set_default_flag("RUN_EXTRA_DIAGNOSTIC_FIGURES", FALSE)
RUN_BASE_RATE_ANALYSIS <- set_default_flag("RUN_BASE_RATE_ANALYSIS", FALSE)

PROJECT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(PROJECT_DIR)
cat("\nSI Prediction Project - Final Single R Script\n")
cat(sprintf("Code directory: %s\n", PROJECT_DIR))

resolve_project_base <- function() {
  env_base <- Sys.getenv("SI_PROJECT_BASE", unset = "")
  cloud_roots <- Sys.glob(file.path(path.expand("~"), "Library/CloudStorage", "GoogleDrive*"))
  google_drive_bases <- unique(c(
    file.path(cloud_roots, "My Drive", "SI_Prediction"),
    file.path(cloud_roots, "Shared drives", "SI_Prediction"),
    file.path(cloud_roots, "SI_Prediction")
  ))
  nearby_bases <- unique(c(
    if (basename(PROJECT_DIR) == "1. Code") dirname(PROJECT_DIR) else character(0),
    dirname(PROJECT_DIR),
    PROJECT_DIR
  ))
  candidates <- unique(c(
    if (nzchar(env_base)) env_base else character(0),
    google_drive_bases,
    nearby_bases
  ))
  candidates <- normalizePath(candidates, mustWork = FALSE)
  is_project_root <- vapply(candidates, function(x) {
    dir.exists(file.path(x, "2. Raw data")) ||
      dir.exists(file.path(x, "3. Cleaned data")) ||
      dir.exists(file.path(x, "4. Results"))
  }, logical(1))
  if (any(is_project_root)) return(candidates[which(is_project_root)[1]])
  deep_hits <- unlist(lapply(cloud_roots, function(root) {
    tryCatch(
      list.files(root, pattern = "^SI_Prediction$", recursive = TRUE,
                 full.names = TRUE, include.dirs = TRUE),
      error = function(e) character()
    )
  }), use.names = FALSE)
  deep_hits <- unique(normalizePath(deep_hits, mustWork = FALSE))
  deep_is_project_root <- vapply(deep_hits, function(x) {
    dir.exists(file.path(x, "2. Raw data")) ||
      dir.exists(file.path(x, "3. Cleaned data")) ||
      dir.exists(file.path(x, "4. Results"))
  }, logical(1))
  if (any(deep_is_project_root)) return(deep_hits[which(deep_is_project_root)[1]])
  PROJECT_DIR
}

discover_project_roots <- function() {
  env_base <- Sys.getenv("SI_PROJECT_BASE", unset = "")
  cloud_roots <- Sys.glob(file.path(path.expand("~"), "Library/CloudStorage", "GoogleDrive*"))
  roots <- unique(c(
    if (nzchar(env_base)) env_base else character(0),
    resolve_project_base(),
    file.path(cloud_roots, "My Drive", "SI_Prediction"),
    file.path(cloud_roots, "Shared drives", "SI_Prediction"),
    file.path(cloud_roots, "SI_Prediction")
  ))
  deep_hits <- unlist(lapply(cloud_roots, function(root) {
    tryCatch(
      list.files(root, pattern = "^SI_Prediction$", recursive = TRUE,
                 full.names = TRUE, include.dirs = TRUE),
      error = function(e) character()
    )
  }), use.names = FALSE)
  roots <- unique(normalizePath(c(roots, deep_hits), mustWork = FALSE))
  roots[vapply(roots, dir.exists, logical(1))]
}

# ------------------------------------------------------------------------------
# Shared dynamic-index helpers (Methods-aligned N thresholds)
# ------------------------------------------------------------------------------
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
# P90 == P10 returns 0.5; values are not clipped.
p90p10_relpos <- function(window_vals, current_x) {
  window_vals <- window_vals[!is.na(window_vals)]
  if (length(window_vals) < 20) return(NA_real_)
  if (is.na(current_x)) return(NA_real_)
  p10 <- unname(stats::quantile(window_vals, 0.10, type = 7))
  p90 <- unname(stats::quantile(window_vals, 0.90, type = 7))
  if (p90 == p10) return(0.5)
  (current_x - p10) / (p90 - p10)
}

# Sample entropy (N >= 250; m = 2, r = 0.2 * SD)
sample_entropy <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 250) return(NA_real_)
  r <- 0.2 * stats::sd(x)
  if (is.na(r) || r == 0) return(NA_real_)
  if (!requireNamespace("TSEntropies", quietly = TRUE)) return(NA_real_)
  tryCatch(TSEntropies::FastSampEn(x, dim = 2, r = r), error = function(e) NA_real_)
}

if (RUN_DATA_CLEANING) {
cat("\n===== OPTIONAL SECTION: DATA CLEANING =====\n")
################################################################################
# RUN ACCELERATION PROCESSING IN PARALLEL
#
# Computes ENMO, 10-second epochs, and dynamic indices with the 4-hour buffer.
# Participant files are processed in parallel.
# Resumable: each participant writes its own per-uid file; finished ones skip.
################################################################################

require_packages_or_stop("pbmcapply", "data cleaning")
suppressPackageStartupMessages({
  library(data.table); library(stringr); library(lubridate); library(parallel)
  library(pbmcapply)
})

DEFAULT_DATA_PATH <- resolve_project_base()

resolve_base_path <- function() {
  existing_base <- if (exists("BASE_PATH", inherits = TRUE)) get("BASE_PATH", inherits = TRUE) else NA_character_
  candidates <- unique(c(
    Sys.getenv("SI_DATA_PATH", unset = NA_character_),
    getwd(),
    existing_base,
    DEFAULT_DATA_PATH
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & dir.exists(candidates)]
  has_ema <- candidates[
    file.exists(file.path(candidates, "3. Cleaned data", "ema_features.csv")) |
      file.exists(file.path(candidates, "ema_features.csv"))
  ]
  if (length(has_ema)) return(has_ema[1])
  has_sensor_logs <- candidates[
    file.exists(file.path(candidates, "2. Raw data", "smslog.csv")) |
      file.exists(file.path(candidates, "2. Raw data", "calllog.csv")) |
      file.exists(file.path(candidates, "smslog.csv")) |
      file.exists(file.path(candidates, "calllog.csv"))
  ]
  if (length(has_sensor_logs)) return(has_sensor_logs[1])
  if (!is.na(existing_base) && nzchar(existing_base)) return(existing_base)
  getwd()
}

BASE_PATH <- resolve_base_path()
CODE_DIR <- file.path(BASE_PATH, "1. Code")
RAW_DATA_DIR <- file.path(BASE_PATH, "2. Raw data")
CLEAN_DATA_DIR <- file.path(BASE_PATH, "3. Cleaned data")
CLEAN_SENSOR_DIR <- file.path(CLEAN_DATA_DIR, "per_sensors")
CLEAN_UID_DIR <- file.path(CLEAN_DATA_DIR, "per_uid")
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
if (!exists("results_dir", inherits = TRUE)) {
  results_dir <- file.path(RESULTS_ROOT_DIR, "report")
}
make_dirs(CLEAN_DATA_DIR, CLEAN_SENSOR_DIR, CLEAN_UID_DIR, results_dir)
cat(sprintf("Using BASE_PATH: %s\n", BASE_PATH))

normalize_uid <- function(x) {
  raw <- trimws(as.character(x))
  raw[raw %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  before_at <- sub("@.*$", "", raw)
  digits <- ifelse(grepl("^\\d{1,4}$", before_at), before_at, NA_character_)
  digits[is.na(digits)] <- stringr::str_extract(raw[is.na(digits)], "(?<!\\d)\\d{1,4}(?!\\d)")
  compact_digits <- gsub("\\D", "", raw)
  has_compact <- is.na(digits) & !is.na(compact_digits) & nzchar(compact_digits)
  digits[has_compact] <- compact_digits[has_compact]
  fallback <- stringr::str_extract(raw, "^[^@_ -]+")
  out <- ifelse(!is.na(digits), digits, fallback)
  numeric_id <- !is.na(out) & grepl("^\\d+$", out)
  out[numeric_id] <- ifelse(
    nchar(out[numeric_id]) <= 4,
    stringr::str_pad(out[numeric_id], width = 4, pad = "0"),
    out[numeric_id]
  )
  out
}

parse_date_any <- function(x) {
  raw <- as.character(x)
  out <- suppressWarnings(as.Date(raw, format = "%Y-%m-%d"))
  miss <- is.na(out)
  out[miss] <- suppressWarnings(as.Date(raw[miss], format = "%m/%d/%Y"))
  miss <- is.na(out)
  out[miss] <- suppressWarnings(lubridate::ymd(raw[miss]))
  miss <- is.na(out)
  out[miss] <- suppressWarnings(lubridate::mdy(raw[miss]))
  out
}

WORKERS      <- 10
SEC_2020     <- 1577836800
ACC_EPOCH_S  <- 10
# Coverage is stored descriptively only. Window validity is governed by the
# minimum-observation thresholds for each dynamic index and by participant-level
# passive compliance criteria.
data.table::setDTthreads(1) # 1 thread/worker -> avoid oversubscription across forks

uid_default <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  normalize_uid(b)
}

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
WINDOWS <- list(hr1 = 3600, hr4 = 14400, day1 = 86400, day7 = 604800)

# --- per-participant worker -----------------------------------------------------
process_one_acc <- function(fp) {
  uid_this <- stringr::str_pad(as.character(uid_default(fp)), width = 4, pad = "0")
  out_path <- file.path(output_dir, sprintf("%s_dynamic_indices_acc.csv", uid_this))
  if (file.exists(out_path)) return(paste("SKIP", uid_this))
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
      cl <- findInterval(ae - 1800, day_vec); ch <- findInterval(ae, day_vec)
      x_current <- mean_thr(if (ch > cl) sig_vec[(cl + 1L):ch] else numeric(0))
      for (wn in names(WINDOWS)) {
        ws <- WINDOWS[[wn]]
        wl <- findInterval(ae - ws, day_vec); wh <- findInterval(ae, day_vec)
        xc <- { x <- if (wh > wl) sig_vec[(wl + 1L):wh] else numeric(0); x[!is.na(x)] }
        n_obs <- length(xc)
        cov   <- min(n_obs / (ws / ACC_EPOCH_S), 1)
        mean_v <- mean_thr(xc); sum_v <- if (n_obs) sum(xc) else NA_real_
        mad_v <- mad_thr(xc); p_v <- p90p10_relpos(xc, x_current); ent_v <- sample_entropy(xc)
        out[[k]] <- data.frame(
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

# --- run in parallel ------------------------------------------------------------
acc_final_path <- file.path(CLEAN_SENSOR_DIR, "dynamic_indices_acc.csv")
if (file.exists(acc_final_path)) {
  cat("dynamic_indices_acc.csv already exists — skipping accelerometer processing\n")
} else {
  acc_files <- list.files(file.path(RAW_DATA_DIR, "acceleration"),
                          pattern = "\\.csv$", full.names = TRUE)
  cat(sprintf("acc files: %d | workers: %d\n", length(acc_files), WORKERS))
  t0 <- proc.time()
  # pbmclapply = mclapply + live progress bar/ETA; fall back to mclapply if needed
  logs <- if (requireNamespace("pbmcapply", quietly = TRUE)) {
    pbmcapply::pbmclapply(acc_files, process_one_acc,
                          mc.cores = WORKERS, mc.preschedule = FALSE)
  } else {
    parallel::mclapply(acc_files, process_one_acc,
                       mc.cores = WORKERS, mc.preschedule = FALSE)
  }
  cat(sprintf("Parallel acc done in %.1f min\n", (proc.time() - t0)["elapsed"] / 60))
  
  # status tally (DONE / SKIP / EMPTY / NOEMA / NORANGE / ERROR)
  status <- sub(" .*", "", unlist(logs))
  print(table(status))
  errs <- grep("^ERROR", unlist(logs), value = TRUE)
  if (length(errs)) { cat("\nErrors:\n"); print(head(errs, 20)) }
  
  # --- combine to dynamic_indices_acc.csv ----------------------------------------
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
}

data.table::setDTthreads(0)   # restore default threading
################################################################################


################################################################################
# GENERIC PARALLEL SENSOR RUNNER  (hr, steps, battery, convo, gps, gps_home)
#
# Sensor processing is parameterized by sensor specification and writes resumable
# per-participant output files.
# Coverage is stored descriptively only. Window validity is governed by the
# minimum-observation thresholds for each dynamic index and by participant-level
# passive compliance criteria.
################################################################################

require_packages_or_stop("pbmcapply", "data cleaning")
suppressPackageStartupMessages({
  library(data.table); library(stringr); library(lubridate)
  library(parallel);   library(pbmcapply); library(geosphere)
})

SEC_2020     <- 1577836800
WINDOWS      <- list(hr1 = 3600, hr4 = 14400, day1 = 86400, day7 = 604800)
uid_default <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  normalize_uid(b)
}
uid_convo   <- function(fp) stringr::str_extract(basename(fp), "(?<=convo_)[^@]+")

parse_pydict <- function(s) {
  if (is.na(s) || s == "") return(NULL)
  s <- gsub("u'", "'", s, fixed = TRUE); s <- gsub("'", '"', s, fixed = TRUE)
  tryCatch(jsonlite::fromJSON(s, simplifyVector = TRUE), error = function(e) NULL)
}

# --- anchors --------------------------------------------------------------------
build_ema_lookup <- function(path = file.path(CLEAN_DATA_DIR, "ema_features.csv")) {
  d <- data.table::fread(path)
  if (!inherits(d$ema_day, "POSIXct")) {
    raw <- as.character(d$ema_day); p <- suppressWarnings(as.POSIXct(raw, tz = "UTC"))
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

# --- HOME LOCATIONS (reuse Dawson's validated homes; do NOT re-run DBSCAN) ------
# Reads per-uid home_coords from time_at_home.csv (column "home_coords" holds a
# python list of (lat, lon) tuples, incl. the hand-assigned ones). Builds a
# named list HOME_COORDS[[uid]] = matrix(cols = lat, lon).
DERIVED_HOME_FILE <- file.path(CLEAN_DATA_DIR, "homes_derived.csv")   # from SI_derive_homes.R
HOME_FILE         <- file.path(CLEAN_DATA_DIR, "time_at_home.csv")    # Dawson's output (fallback)

load_home_coords <- function() {
  # 1) preferred: homes_derived.csv (long: uid, home_lat, home_lon)
  if (file.exists(DERIVED_HOME_FILE)) {
    d <- data.table::fread(DERIVED_HOME_FILE, colClasses = list(character = "uid"))
    d[, uid := stringr::str_pad(uid, 4, pad = "0")]
    L <- lapply(split(d, d$uid), function(x) {
      m <- as.matrix(x[, .(home_lat, home_lon)]); colnames(m) <- c("lat", "lon"); m
    })
    cat(sprintf("homes loaded from homes_derived.csv for %d participants\n", length(L)))
    return(L)
  }
  # 2) fallback: parse home_coords tuples from time_at_home.csv
  if (file.exists(HOME_FILE)) {
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
  out <- dt[, .(signal = sum(signal, na.rm = TRUE)), by = by_cols]  # body length/count per epoch
  out[, day := as.POSIXct(bin, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day)
  if ("uid" %in% names(out)) out[, .(uid, day, signal)] else out[, .(day, signal)]
}

read_calls <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
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
  out <- dt[, .(signal = sum(signal, na.rm = TRUE)), by = by_cols]  # total call duration per epoch
  out[, day := as.POSIXct(bin, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day)
  if ("uid" %in% names(out)) out[, .(uid, day, signal)] else out[, .(day, signal)]
}

read_stress <- function(fp) {
  dt <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  need <- c("calendarDate","startTimeInSeconds","durationInSeconds","timeOffsetStressLevelValues")
  if (is.null(dt) || !all(need %in% names(dt))) return(NULL)
  data.table::setorder(dt, calendarDate, -durationInSeconds)
  dt <- dt[!duplicated(calendarDate)]                       # dedup per day
  rows <- lapply(seq_len(nrow(dt)), function(i) {
    v <- parse_pydict(dt$timeOffsetStressLevelValues[i])
    if (is.null(v) || !length(v)) return(NULL)
    off <- suppressWarnings(as.numeric(names(v))); lvl <- suppressWarnings(as.numeric(v))
    keep <- !is.na(off) & !is.na(lvl) & lvl >= 0 & lvl <= 100   # drop -1/-2 no-data
    if (!any(keep)) return(NULL)
    data.table::data.table(
      day = as.POSIXct(dt$startTimeInSeconds[i] + off[keep], origin = "1970-01-01", tz = "UTC"),
      signal = lvl[keep])
  })
  out <- data.table::rbindlist(rows)
  if (!nrow(out)) return(NULL)
  data.table::setkey(out, day); out[]
}

# --- GPS: LATITUDE/LONGITUDE/ACCURACY/DAY -> distance traveled / 30-min epoch ----
#   distance = haversine between consecutive fixes (meters), summed per 30-min bin.
#   Filters low-accuracy fixes and caps implausible-speed jumps (GPS error).
read_gps <- function(fp, max_accuracy_m = 100, max_speed_ms = 70) {
  d <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  data.table::setnames(d, toupper(names(d)))
  if (!all(c("LATITUDE", "LONGITUDE", "DAY") %in% names(d))) return(NULL)
  d[, day := as.POSIXct(as.character(DAY), tz = "UTC")]
  d[, lat := suppressWarnings(as.numeric(LATITUDE))]
  d[, lon := suppressWarnings(as.numeric(LONGITUDE))]
  d[, acc := if ("ACCURACY" %in% names(d)) suppressWarnings(as.numeric(ACCURACY)) else NA_real_]
  d <- d[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(lat) & !is.na(lon) &
           lat >= -90 & lat <= 90 & lon >= -180 & lon <= 180 &
           (is.na(acc) | acc <= max_accuracy_m)]
  if (nrow(d) < 2) return(NULL)
  data.table::setkey(d, day); n <- nrow(d)
  dist <- c(0, geosphere::distHaversine(cbind(d$lon[-n], d$lat[-n]),
                                        cbind(d$lon[-1], d$lat[-1])))   # meters
  tgap <- c(NA, diff(as.numeric(d$day)))                               # seconds
  spd  <- dist / tgap
  dist[!is.na(spd) & spd > max_speed_ms] <- NA                         # drop GPS jumps
  d[, dist_m := dist]
  d[, ep := floor(as.numeric(day) / 1800) * 1800]                      # 30-min epoch
  out <- d[, .(signal = sum(dist_m, na.rm = TRUE)), by = ep]           # total m/epoch
  out[, day := as.POSIXct(ep, origin = "1970-01-01", tz = "UTC")]
  data.table::setkey(out, day); out[, .(day, signal)]
}

# --- GPS TIME-AT-HOME: fraction of fixes within 40m of home x 1800 = sec/epoch --
#   Reuses Dawson's per-uid home_coords (HOME_COORDS), applies the 40m flag to
#   raw_gps_clean fixes (haversine ~ geodesic at this scale), aggregates per
#   30-min epoch. signal = seconds at home per epoch (0..1800).
read_gps_home <- function(fp, radius_m = 40, max_accuracy_m = 100) {
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
  d[, acc := if ("ACCURACY" %in% names(d)) suppressWarnings(as.numeric(ACCURACY)) else NA_real_]
  d <- d[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(lat) & !is.na(lon) &
           (is.na(acc) | acc <= max_accuracy_m)]
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
# CORE: compute dynamic indices for one participant's (day, signal) series
# ==============================================================================
compute_indices <- function(day_vec, sig_vec, ema_this, epoch_sec,
                            windows = WINDOWS) {
  out <- vector("list", nrow(ema_this) * length(windows)); k <- 1L
  for (j in seq_len(nrow(ema_this))) {
    ae <- as.numeric(ema_this$anchor_end[j])
    cl <- findInterval(ae - 1800, day_vec); ch <- findInterval(ae, day_vec)
    x_current <- mean_thr(if (ch > cl) sig_vec[(cl + 1L):ch] else numeric(0))
    for (wn in names(windows)) {
      ws <- windows[[wn]]
      wl <- findInterval(ae - ws, day_vec); wh <- findInterval(ae, day_vec)
      xc <- { x <- if (wh > wl) sig_vec[(wl + 1L):wh] else numeric(0); x[!is.na(x)] }
      n_obs <- length(xc)
      cov   <- if (!is.na(epoch_sec)) min(n_obs / (ws / epoch_sec), 1) else NA_real_
      mean_v <- mean_thr(xc); sum_v <- if (n_obs) sum(xc) else NA_real_
      mad_v <- mad_thr(xc); p_v <- p90p10_relpos(xc, x_current); ent_v <- sample_entropy(xc)
      out[[k]] <- data.frame(
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

# --- one participant -> writes per-uid csv, returns status ----------------------
process_one <- function(fp, spec) {
  uid_this <- stringr::str_pad(as.character(spec$uid_fn(fp)), width = 4, pad = "0")
  out_dir  <- file.path(CLEAN_UID_DIR, paste0(spec$tag, "_per_uid"))
  out_path <- file.path(out_dir, sprintf("%s_dynamic_indices_%s.csv", uid_this, spec$tag))
  if (file.exists(out_path)) return(paste("SKIP", uid_this))
  tryCatch({
    dt <- spec$read(fp)
    if (is.null(dt) || !nrow(dt)) return(paste("EMPTY", uid_this))
    ema_this <- ema_lookup[uid == uid_this]
    if (!nrow(ema_this)) return(paste("NOEMA", uid_this))
    min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
    max_t <- as.numeric(max(ema_this$anchor_end))
    dt <- dt[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
    if (!nrow(dt)) return(paste("NORANGE", uid_this))
    data.table::setkey(dt, day)
    res <- compute_indices(as.numeric(dt$day), dt$signal, ema_this, spec$epoch_sec)
    make_dirs(out_dir)
    data.table::fwrite(res, out_path)
    paste("DONE", uid_this, nrow(ema_this))
  }, error = function(e) paste("ERROR", uid_this, conditionMessage(e)))
}

sensor_search_roots <- function() {
  roots <- c(
    RAW_DATA_DIR,
    CLEAN_SENSOR_DIR,
    CLEAN_DATA_DIR,
    BASE_PATH,
    Sys.getenv("SI_DATA_PATH", unset = NA_character_),
    resolve_project_base()
  )
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

inspect_sensor_input_files <- function(folders = c("sleep", "smslog", "calllog"), n = 3) {
  cat(sprintf("Search roots:\n  %s\n", paste(sensor_search_roots(), collapse = "\n  ")))
  out <- lapply(folders, function(folder) {
    paths <- unique(c(sensor_path_candidates(folder), case_insensitive_csv_matches(folder)))
    files <- paths[file.exists(paths) & !dir.exists(paths)]
    dirs <- paths[dir.exists(paths)]
    if (length(dirs)) {
      files <- c(files, unlist(lapply(dirs, list.files, pattern = "\\.csv$", full.names = TRUE),
                               use.names = FALSE))
    }
    files <- unique(files)
    cat(sprintf("\n[%s] %d CSV files\n", folder, length(files)))
    if (!length(files)) return(data.table::data.table(folder = folder, files = 0L))
    rows <- lapply(utils::head(files, n), function(fp) {
      d <- tryCatch(data.table::fread(fp, nrows = 0), error = function(e) NULL)
      cols <- if (is.null(d)) "<cannot read>" else paste(names(d), collapse = ", ")
      cat(sprintf("  %s\n    uid guess: %s\n    columns: %s\n",
                  basename(fp), uid_default(fp), cols))
      data.table::data.table(
        folder = folder,
        file = basename(fp),
        uid_guess = uid_default(fp),
        columns = cols
      )
    })
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  })
  res <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  out_csv <- file.path(results_dir, "sensor_input_file_inspection.csv")
  make_dirs(dirname(out_csv))
  data.table::fwrite(res, out_csv)
  invisible(res)
}

diagnose_event_uid_overlap <- function(spec) {
  files <- list_sensor_files(spec)
  cat(sprintf("\n=== UID overlap diagnostic: %s ===\n", spec$tag))
  cat(sprintf("Files found: %d\n", length(files)))
  if (!length(files)) return(invisible(NULL))
  dt <- data.table::rbindlist(lapply(files, spec$read), use.names = TRUE, fill = TRUE)
  if (!nrow(dt)) {
    cat("Reader returned zero usable rows after parsing time, uid, and signal.\n")
    return(invisible(NULL))
  }
  if (!"uid" %in% names(dt)) {
    cat("Reader output does not include uid. Check the source ID column name.\n")
    return(invisible(NULL))
  }
  dt[, uid := map_to_ema_uid(uid)]
  event_uids <- sort(unique(dt[!is.na(uid), uid]))
  ema_uids <- sort(unique(ema_lookup[!is.na(uid), uid]))
  overlap <- intersect(event_uids, ema_uids)
  cat(sprintf("Event UIDs: %d | EMA UIDs: %d | overlap: %d\n",
              length(event_uids), length(ema_uids), length(overlap)))
  cat(sprintf("Event time range: %s .. %s\n", min(dt$day, na.rm = TRUE), max(dt$day, na.rm = TRUE)))
  cat(sprintf("EMA anchor range: %s .. %s\n", min(ema_lookup$anchor_end, na.rm = TRUE),
              max(ema_lookup$anchor_end, na.rm = TRUE)))
  if (!length(overlap)) {
    cat("No UID overlap after EMA alias mapping. Compare these locally:\n")
    cat(sprintf("  event uid examples: %s\n", paste(utils::head(event_uids, 10), collapse = ", ")))
    cat(sprintf("  EMA uid examples:   %s\n", paste(utils::head(ema_uids, 10), collapse = ", ")))
  }
  invisible(list(event_uids = event_uids, ema_uids = ema_uids, overlap = overlap))
}

# ==============================================================================
# test_one(): verify a reader on ONE participant (prints, does not write)
# ==============================================================================
test_one <- function(spec, uid_sel = "0011") {
  files <- list_sensor_files(spec)
  if (isTRUE(spec$combined_file)) {
    cat(sprintf("\n=== TEST %s (uid %s) ===\n", spec$tag, uid_sel))
    if (!length(files)) { cat("  no combined file found\n"); return(invisible()) }
    dt <- data.table::rbindlist(lapply(files, spec$read), use.names = TRUE, fill = TRUE)
    if (is.null(dt) || !nrow(dt) || !"uid" %in% names(dt)) {
      cat("  reader returned EMPTY or no uid column\n")
      return(invisible())
    }
    dt[, uid := map_to_ema_uid(uid)]
    dt <- dt[uid == uid_sel]
    if (!nrow(dt)) { cat("  no events for this uid in combined file\n"); return(invisible()) }
    cat(sprintf("  combined files: %d\n  signal rows: %d | range: %s .. %s\n",
                length(files), nrow(dt), min(dt$day), max(dt$day)))
    print(summary(dt$signal))
    ema_this <- ema_lookup[uid == uid_sel]
    if (!nrow(ema_this)) { cat("  no EMA for this uid\n"); return(invisible()) }
    min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
    max_t <- as.numeric(max(ema_this$anchor_end))
    dt <- dt[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
    data.table::setkey(dt, day)
    res <- compute_indices(as.numeric(dt$day), dt$signal, ema_this, spec$epoch_sec)
    print(utils::head(res))
    print(res[, .(n = .N, mean_cov = round(mean(coverage, na.rm = TRUE), 2),
                  pct_gated = round(mean(is.na(mean_val)) * 100),
                  pct_ent_na = round(mean(is.na(entropy)) * 100)), by = window])
    return(invisible(res))
  }
  hit <- files[stringr::str_pad(as.character(vapply(files, spec$uid_fn, character(1))), width = 4, pad = "0") == uid_sel]
  cat(sprintf("\n=== TEST %s (uid %s) ===\n", spec$tag, uid_sel))
  if (!length(hit)) { cat("  no file for this uid\n"); return(invisible()) }
  dt <- spec$read(hit[1])
  if (is.null(dt) || !nrow(dt)) { cat("  reader returned EMPTY\n"); return(invisible()) }
  cat(sprintf("  file: %s\n  signal rows: %d | range: %s .. %s\n",
              basename(hit[1]), nrow(dt), min(dt$day), max(dt$day)))
  print(summary(dt$signal))
  ema_this <- ema_lookup[uid == uid_sel]
  if (!nrow(ema_this)) { cat("  no EMA for this uid\n"); return(invisible()) }
  min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
  max_t <- as.numeric(max(ema_this$anchor_end))
  dt <- dt[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
  data.table::setkey(dt, day)
  res <- compute_indices(as.numeric(dt$day), dt$signal, ema_this, spec$epoch_sec)
  print(utils::head(res))
  print(res[, .(n = .N, mean_cov = round(mean(coverage, na.rm = TRUE), 2),
                pct_gated = round(mean(is.na(mean_val)) * 100),
                pct_ent_na = round(mean(is.na(entropy)) * 100)), by = window])
  invisible(res)
}

# ==============================================================================
# run_full(): parallel scale-up for one sensor
# ==============================================================================
run_full <- function(spec, workers = 10, force = FALSE) {
  final_path <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", spec$tag))
  if (file.exists(final_path) && !force) {
    cat(sprintf("dynamic_indices_%s.csv already exists — skipping %s processing\n",
                spec$tag, spec$tag))
    return(invisible(data.table::fread(final_path, colClasses = list(character = "uid"))))
  }
  if (force && file.exists(final_path)) {
    cat(sprintf("force=TRUE — recalculating %s and replacing %s\n", spec$tag, basename(final_path)))
    unlink(final_path)
    unlink(file.path(CLEAN_UID_DIR, paste0(spec$tag, "_per_uid")), recursive = TRUE)
  }
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
      ema_this <- ema_lookup[uid == uid_this]
      min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
      max_t <- as.numeric(max(ema_this$anchor_end))
      dt_this <- dt_this[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
      if (!nrow(dt_this)) next
      data.table::setkey(dt_this, day)
      out_list[[i]] <- compute_indices(as.numeric(dt_this$day), dt_this$signal,
                                       ema_this, spec$epoch_sec)
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
# SLEEP  — NOT windowed. One record per night, joined to each assessment by date.
#   total_sleep + rem_sleep (minutes), deduped per night, night-before of the
#   assessment date. No sliding windows / coverage / entropy apply here.
#   Fast & single file per uid, so no parallelism needed.
# ==============================================================================
NIGHT_OFFSET_DAYS <- 0L   # calendarDate assumed = wake date; set 1 if onset date
.val_rank <- function(v) ifelse(grepl("ENHANCED", v), 2L, 1L) * 10L +
  ifelse(grepl("FINAL", v), 2L, 1L)

run_sleep <- function(folder = "sleep", force = FALSE) {
  final_path <- file.path(CLEAN_SENSOR_DIR, "dynamic_indices_sleep.csv")
  if (file.exists(final_path) && !force) {
    cat("dynamic_indices_sleep.csv already exists — skipping sleep processing\n")
    return(invisible(data.table::fread(final_path, colClasses = list(character = "uid"))))
  }
  if (force && file.exists(final_path)) {
    cat("force=TRUE — recalculating sleep and replacing dynamic_indices_sleep.csv\n")
    unlink(final_path)
  }
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
    need <- c("calendarDate","durationInSeconds","remSleepInSeconds","startTimeInSeconds","validation")
    if (is.null(s) || !all(need %in% names(s))) next
    s[, .vr := .val_rank(validation)]
    data.table::setorder(s, calendarDate, -.vr, -durationInSeconds, -startTimeInSeconds)
    s <- s[!duplicated(calendarDate)]                                          # 1 record/night
    sl <- data.table::data.table(
      uid = uid_this, sleep_date = as.Date(s$calendarDate),
      total_sleep = as.numeric(s$durationInSeconds) / 60,
      rem_sleep   = as.numeric(s$remSleepInSeconds) / 60)
    sl <- sl[!is.na(total_sleep) & total_sleep >= 60 & total_sleep <= 14*60 & rem_sleep <= total_sleep]
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
#   acc is processed in the accelerometer block above because it uses a separate
#   per-uid cache directory. Sleep is processed by run_sleep() because it is a
#   night-before join rather than a rolling-window sensor.
# ==============================================================================
SENSORS <- list(
  hr      = list(tag = "hr",      folder = "hr",      read = read_hr,      uid_fn = uid_default, epoch_sec = 10),
  steps   = list(tag = "steps",   folder = "steps",   read = read_steps,   uid_fn = uid_default, epoch_sec = 30),
  battery = list(tag = "battery", folder = "battery", read = read_battery, uid_fn = uid_default, epoch_sec = 3600),
  convo   = list(tag = "convo",   folder = "convo",        read = read_convo,   uid_fn = uid_convo,   epoch_sec = NA_real_),
  gps      = list(tag = "gps",      folder = "raw_gps_clean", read = read_gps,      uid_fn = uid_default, epoch_sec = 1800),
  gps_home = list(tag = "gps_home", folder = "raw_gps_clean", read = read_gps_home, uid_fn = uid_default, epoch_sec = 1800),
  sms      = list(tag = "sms",      folder = c("smslog.csv", "smslog", "sms.csv", "sms"),
                  read = read_sms, uid_fn = uid_default, epoch_sec = NA_real_, combined_file = TRUE),
  calls    = list(tag = "calls",    folder = c("calllog.csv", "calllog", "calls.csv", "calls"),
                  read = read_calls, uid_fn = uid_default, epoch_sec = NA_real_, combined_file = TRUE)
)

# ------------------------------------------------------------------------------
# ON SOURCE: test ONE participant for each sensor (verify readers, no scaling)
# ------------------------------------------------------------------------------
for (nm in names(SENSORS)) {
  tryCatch(test_one(SENSORS[[nm]], "0011"),
           error = function(e) cat(sprintf("TEST ERROR %s: %s\n", nm, conditionMessage(e))))
}

cat("\n--- to scale a sensor up once it looks right: ---\n")
cat("# acc is processed by the accelerometer block above\n")
cat("# run_full(SENSORS$hr)      # ~minutes; hr files are far smaller than acc\n")
cat("# run_full(SENSORS$steps)\n# run_full(SENSORS$battery)\n")
cat("# run_full(SENSORS$convo)\n# run_full(SENSORS$sms)\n# run_full(SENSORS$calls)\n")
cat("# run_sleep()               # sleep is separate (night-before join, not windowed)\n")
################################################################################

run_full(SENSORS$hr)        # biggest after acc (~610k rows/uid) — a few min
run_full(SENSORS$steps)
run_full(SENSORS$battery)
run_full(SENSORS$convo)
run_full(SENSORS$gps)
run_full(SENSORS$gps_home)
run_full(SENSORS$sms, force = FALSE)
run_full(SENSORS$calls, force = FALSE)
run_sleep()

reshape_dynamic_indices <- function(path, tag) {
  d <- data.table::fread(path, colClasses = list(character = "uid"))
  d[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  key_cols <- c("uid", "date_rep", "ema_time")
  if (!all(key_cols %in% names(d))) {
    warning(sprintf("%s is missing merge keys; skipping", basename(path)))
    return(NULL)
  }
  
  if (tag == "sleep") {
    out <- d[, .(
      total_sleep = mean(total_sleep, na.rm = TRUE),
      rem_sleep = mean(rem_sleep, na.rm = TRUE)
    ), by = key_cols]
    out[, date_for_sleep := parse_date_any(date_rep)]
    data.table::setorder(out, uid, date_for_sleep, ema_time)
    out[, sleep_mad_total_sleep_day7 := {
      vapply(seq_len(.N), function(i) {
        x <- total_sleep[date_for_sleep < date_for_sleep[i] &
                           date_for_sleep >= date_for_sleep[i] - 7]
        mad_thr(x)
      }, numeric(1))
    }, by = uid]
    out[, sleep_mad_rem_sleep_day7 := {
      vapply(seq_len(.N), function(i) {
        x <- rem_sleep[date_for_sleep < date_for_sleep[i] &
                         date_for_sleep >= date_for_sleep[i] - 7]
        mad_thr(x)
      }, numeric(1))
    }, by = uid]
    for (w in names(WINDOWS)) {
      if ("total_sleep" %in% names(out)) out[[paste0("sleep_mean_total_sleep_", w)]] <- out$total_sleep
      if ("rem_sleep" %in% names(out)) out[[paste0("sleep_mean_rem_sleep_", w)]] <- out$rem_sleep
    }
    out[, intersect(c("total_sleep", "rem_sleep", "date_for_sleep"), names(out)) := NULL]
    cat(sprintf("Sleep reshape columns: %s\n",
                paste(grep("^sleep_", names(out), value = TRUE), collapse = ", ")))
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

verify_sleep_mad_features <- function(path = file.path(CLEAN_DATA_DIR, "all_features_combined.csv")) {
  d <- data.table::fread(path, nrows = 1)
  cols <- grep("^sleep_mad_", names(d), value = TRUE)
  if (!length(cols)) {
    warning("No sleep_mad_ columns found. Re-run the data-cleaning section, then run assemble_feature_matrix_local().")
  } else {
    cat(sprintf("Found sleep MAD columns: %s\n", paste(cols, collapse = ", ")))
  }
  invisible(cols)
}

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

read_sensor_day_file <- function(sensor, fp) {
  tryCatch({
    if (sensor == "sleep") {
      d <- data.table::fread(fp, select = c("calendarDate", "durationInSeconds",
                                            "remSleepInSeconds", "validation"))
      if (!"calendarDate" %in% names(d)) return(NULL)
      d[, uid := normalize_uid(uid_default(fp))]
      d[, date := parse_date_any(calendarDate)]
      d <- d[!is.na(date)]
      if ("durationInSeconds" %in% names(d)) {
        d[, durationInSeconds := suppressWarnings(as.numeric(durationInSeconds))]
        d <- d[!is.na(durationInSeconds) & durationInSeconds >= 60 * 60 &
                 durationInSeconds <= 14 * 60 * 60]
      }
      return(unique(d[, .(uid, date)]))
    }
    uid_this <- if (sensor == "acc") normalize_uid(uid_default(fp)) else normalize_uid(SENSORS[[sensor]]$uid_fn(fp))
    d <- if (sensor == "acc") read_acc(fp) else SENSORS[[sensor]]$read(fp)
    if (is.null(d) || !nrow(d) || !"day" %in% names(d)) return(NULL)
    data.table::data.table(uid = uid_this, date = as.Date(d$day))
  }, error = function(e) NULL)
}

raw_sensor_days <- function(sensor, workers = 10) {
  if (sensor == "sleep") {
    sdir <- file.path(RAW_DATA_DIR, "sleep")
    files <- list.files(sdir, pattern = "\\.csv$", full.names = TRUE)
  } else if (sensor == "acc") {
    files <- list.files(file.path(RAW_DATA_DIR, "acceleration"), pattern = "\\.csv$", full.names = TRUE)
  } else if (!is.null(SENSORS[[sensor]])) {
    files <- list_sensor_files(SENSORS[[sensor]])
  } else {
    files <- character(0)
  }
  if (!length(files)) {
    return(data.table::data.table(uid = character(), sensor = sensor, date = as.Date(character())))
  }
  
  data.table::setDTthreads(1)
  cores <- max(1L, min(workers, parallel::detectCores() - 1L, length(files)))
  rows <- if (requireNamespace("pbmcapply", quietly = TRUE) && cores > 1L) {
    pbmcapply::pbmclapply(files, read_sensor_day_file, sensor = sensor,
                          mc.cores = cores, mc.preschedule = FALSE)
  } else if (cores > 1L) {
    parallel::mclapply(files, read_sensor_day_file, sensor = sensor,
                       mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(files, read_sensor_day_file, sensor = sensor)
  }
  data.table::setDTthreads(0)
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(data.table::data.table(uid = character(), sensor = sensor, date = as.Date(character())))
  out <- out[!is.na(uid) & !is.na(date)]
  out[, sensor := sensor]
  unique(out[, .(uid, sensor, date)])
}

compute_raw_passive_compliance <- function(expected_days,
                                           passive_sensors = c("acc","steps","hr","battery","gps","gps_home"),
                                           force = FALSE,
                                           method = c("dynamic_n_obs", "raw_files"),
                                           workers = 10,
                                           outfile = file.path(results_dir,
                                                               "raw_passive_compliance_by_uid_sensor.csv")) {
  method <- match.arg(method)
  if (file.exists(outfile) && !force) {
    return(data.table::fread(outfile, colClasses = list(character = "uid")))
  }
  make_dirs(dirname(outfile))
  expected_days[, uid := normalize_uid(uid)]
  expected_days <- unique(expected_days[!is.na(uid) & !is.na(date), .(uid, date)])
  uid_expected <- expected_days[, .(expected_days = data.table::uniqueN(date)), by = uid]
  valid_days <- data.table::rbindlist(lapply(passive_sensors, function(s) {
    cat(sprintf("Computing passive compliance days: %s (%s)\n", s, method))
    if (method == "dynamic_n_obs") dynamic_sensor_days(s) else raw_sensor_days(s, workers = workers)
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
# ANALYTIC INCLUSION + FEATURE-DROP CRITERIA  (Methods 2.4.1 / 2.4.3)
#
# Applied to all_features_combined.csv -> all_features_included.csv
#
# Drop criteria implemented:
#   1. DROP stress entirely (sparse, unreliable).
#   2. PARTICIPANT passive sensing data filter: mean data availability across
#      passive streams, averaged equally across streams. SMS and calls are
#      excluded from this threshold because they are sparse event logs.
#      Default = 0.75, aligned with the Methods section.
#   3. DROP feature columns whose values are ALL NA (e.g. indices never estimable
#      — like event-based coverage, or sensors with no data).
#   4. (optional) SI-count criterion: keep PTs with >= 10 SI-positive and
#      >= 20 SI-negative q9 observations (Methods 2.4.1; Peduzzi/Riley).
#
################################################################################
suppressPackageStartupMessages({ library(data.table); library(stringr) })

apply_inclusion_criteria <- function(
    infile             = "all_features_combined.csv",
    outfile            = "all_features_included.csv",
    ema_compliance_threshold = 0.75,
    passive_data_threshold = 0.75,
    passive_sensors    = c("acc","steps","hr","battery","gps","gps_home"),
    drop_sensors       = c("stress"),          # dropped entirely
    apply_ema_compliance = TRUE,
    apply_passive_compliance = TRUE,
    passive_compliance_method = c("dynamic_n_obs", "raw_files"),
    passive_compliance_workers = WORKERS,
    force_recompute_passive_compliance = FALSE,
    apply_si_counts    = TRUE, min_si_pos = 10, min_si_neg = 20) {
  passive_compliance_method <- match.arg(passive_compliance_method)
  
  input_path <- file.path(CLEAN_DATA_DIR, infile)
  cat(sprintf("Reading inclusion input: %s\n", input_path))
  d <- data.table::fread(input_path)
  d[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  n_uid0 <- data.table::uniqueN(d$uid); n_col0 <- ncol(d)
  cat(sprintf("Input: %d rows, %d cols, %d participants\n", nrow(d), n_col0, n_uid0))
  
  ## 1. DROP stress (and any other listed sensors) entirely -------------------
  for (s in drop_sensors) {
    cols <- grep(sprintf("^%s_", s), names(d), value = TRUE)
    if (length(cols)) { d[, (cols) := NULL]; cat(sprintf("Dropped %s: %d cols\n", s, length(cols))) }
  }
  
  ## 2. EMA compliance from participant-level EMA completion rate --------------
  keep_ema <- unique(d$uid)
  ema_comp <- data.table::data.table(uid = unique(d$uid), ema_compliance = NA_real_)
  if (apply_ema_compliance && "compliance" %in% names(d)) {
    ema_comp <- d[, .(ema_compliance = max(suppressWarnings(as.numeric(compliance)), na.rm = TRUE)), by = uid]
    ema_comp[!is.finite(ema_compliance), ema_compliance := NA_real_]
    keep_ema <- ema_comp[!is.na(ema_compliance) & ema_compliance >= ema_compliance_threshold, uid]
    cat(sprintf("EMA compliance >= %.0f%%: %d / %d participants kept (median compliance: %.2f)\n",
                ema_compliance_threshold * 100, length(keep_ema), n_uid0,
                stats::median(ema_comp$ema_compliance, na.rm = TRUE)))
  }
  
  ## 3. Raw/minimally processed passive sensing compliance --------------------
  expected_days <- unique(d[, .(uid, date = parse_date_any(date_rep))])
  keep_cov <- unique(d$uid)
  uid_cov <- data.table::data.table(uid = unique(d$uid), mean_cov = NA_real_)
  uid_sensor <- data.table::data.table(uid = character(), sensor = character(),
                                       expected_days = integer(), valid_days = integer(),
                                       passive_compliance = numeric())
  if (apply_passive_compliance) {
    uid_sensor <- compute_raw_passive_compliance(
      expected_days = expected_days,
      passive_sensors = passive_sensors,
      force = force_recompute_passive_compliance,
      method = passive_compliance_method,
      workers = passive_compliance_workers
    )
    uid_cov <- uid_sensor[, .(mean_cov = mean(passive_compliance, na.rm = TRUE)), by = uid]
    keep_cov <- uid_cov[is.finite(mean_cov) & mean_cov >= passive_data_threshold, uid]
    cat(sprintf("Raw passive sensing compliance >= %.0f%%: %d / %d participants kept (median compliance: %.2f)\n",
                passive_data_threshold * 100, length(keep_cov), n_uid0,
                stats::median(uid_cov$mean_cov, na.rm = TRUE)))
  }
  sensor_availability <- uid_sensor[, .(
    n_participants = .N,
    median_availability = stats::median(passive_compliance, na.rm = TRUE),
    mean_availability = mean(passive_compliance, na.rm = TRUE),
    n_below_75 = sum(passive_compliance < passive_data_threshold, na.rm = TRUE),
    pct_below_75 = mean(passive_compliance < passive_data_threshold, na.rm = TRUE) * 100
  ), by = sensor][order(-n_below_75, median_availability)]
  cat("Raw passive sensor compliance summary (largest drop first):\n")
  print(sensor_availability)
  
  ## 4. SI-count criterion (q9 raw 0-100; SI-positive = q9 > 0) ----------------
  keep_si <- unique(d$uid)
  if (apply_si_counts && "q9" %in% names(d)) {
    si <- d[!is.na(q9), .(n_pos = sum(q9 > 0), n_neg = sum(q9 == 0)), by = uid]
    keep_si <- si[n_pos >= min_si_pos & n_neg >= min_si_neg, uid]
    cat(sprintf("SI counts (>=%d pos, >=%d neg): %d participants kept\n",
                min_si_pos, min_si_neg, length(keep_si)))
  }
  
  ## Participant-level inclusion breakdown ------------------------------------
  all_uids <- sort(unique(d$uid))
  si_counts <- if ("q9" %in% names(d)) {
    d[!is.na(q9), .(n_si_pos = sum(q9 > 0), n_si_neg = sum(q9 == 0)), by = uid]
  } else {
    data.table::data.table(uid = all_uids, n_si_pos = NA_integer_, n_si_neg = NA_integer_)
  }
  inclusion_report <- data.table::data.table(uid = all_uids)[
    uid_cov, on = "uid"
  ][
    ema_comp, on = "uid"
  ][
    si_counts, on = "uid"
  ]
  inclusion_report[, pass_ema := uid %in% keep_ema]
  inclusion_report[, pass_passive := uid %in% keep_cov]
  inclusion_report[, pass_si := uid %in% keep_si]
  inclusion_report[, fail_si_positive := !is.na(n_si_pos) & n_si_pos < min_si_pos]
  inclusion_report[, fail_si_negative := !is.na(n_si_neg) & n_si_neg < min_si_neg]
  inclusion_report[, si_failure_reason := data.table::fcase(
    pass_si, "pass_si_count",
    fail_si_positive & fail_si_negative, "fail_both_si_positive_and_non_si",
    fail_si_positive, "fail_si_positive_observations",
    fail_si_negative, "fail_non_si_observations",
    default = "unknown_or_missing_q9"
  )]
  inclusion_report[, inclusion_group := data.table::fcase(
    pass_ema & pass_passive & pass_si, "pass_all",
    !pass_ema & pass_passive & pass_si, "fail_ema_only",
    pass_ema & !pass_passive & pass_si, "fail_passive_only",
    pass_ema & pass_passive & !pass_si, "fail_si_only",
    default = "pass_neither"
  )]
  breakdown <- inclusion_report[, .N, by = inclusion_group][order(inclusion_group)]
  cat("Inclusion overlap breakdown:\n")
  print(breakdown)
  si_keystone <- inclusion_report[pass_ema & pass_passive, .(
    denominator_passed_ema_passive = .N,
    retained_after_si_count = sum(pass_si, na.rm = TRUE),
    excluded_by_si_count = sum(!pass_si, na.rm = TRUE),
    excluded_due_to_low_si_positive = sum(!pass_si & fail_si_positive, na.rm = TRUE),
    excluded_due_to_low_si_positive_only = sum(!pass_si & fail_si_positive & !fail_si_negative, na.rm = TRUE),
    excluded_due_to_low_non_si_only = sum(!pass_si & !fail_si_positive & fail_si_negative, na.rm = TRUE),
    excluded_due_to_both_low_si_and_low_non_si = sum(!pass_si & fail_si_positive & fail_si_negative, na.rm = TRUE)
  )]
  si_reason_breakdown <- inclusion_report[pass_ema & pass_passive & !pass_si,
                                          .N, by = si_failure_reason][order(-N)]
  cat("\nKEYSTONE SI-count exclusion summary among participants who passed EMA + passive criteria:\n")
  print(si_keystone)
  cat("SI-count failure reason breakdown:\n")
  print(si_reason_breakdown)
  si_keystone_path <- file.path(results_dir, "keystone_si_count_exclusion_summary.csv")
  data.table::fwrite(si_keystone, si_keystone_path)
  si_reason_path <- file.path(results_dir, "keystone_si_count_failure_reasons.csv")
  data.table::fwrite(si_reason_breakdown, si_reason_path)
  cat(sprintf("Saved keystone SI-count summary: %s\n", si_keystone_path))
  cat(sprintf("Saved keystone SI-count reason breakdown: %s\n", si_reason_path))
  report_path <- file.path(results_dir, "inclusion_participant_report.csv")
  make_dirs(dirname(report_path))
  data.table::fwrite(inclusion_report, report_path)
  cat(sprintf("Saved inclusion participant report: %s\n", report_path))
  sensor_report_path <- file.path(results_dir, "passive_sensor_availability_report.csv")
  data.table::fwrite(sensor_availability, sensor_report_path)
  cat(sprintf("Saved passive sensor availability report: %s\n", sensor_report_path))
  
  ## intersect all participant criteria ---------------------------------------
  keep <- Reduce(intersect, list(keep_ema, keep_cov, keep_si))
  d <- d[uid %in% keep]
  cat(sprintf("==> %d participants retained after all criteria\n", length(keep)))
  
  ## 5. DROP all-NA feature columns (computed AFTER row filtering) -------------
  id_cols  <- intersect(c("uid","ema_day","date_rep","ema_time"), names(d))
  feat_cols <- setdiff(names(d), id_cols)
  all_na <- feat_cols[vapply(feat_cols, function(c) all(is.na(d[[c]])), logical(1))]
  if (length(all_na)) {
    d[, (all_na) := NULL]
    cat(sprintf("Dropped %d all-NA feature columns\n", length(all_na)))
  }
  
  data.table::fwrite(d, file.path(CLEAN_DATA_DIR, outfile))
  cat(sprintf("\nSaved %s: %d rows x %d cols (%d participants)\n",
              outfile, nrow(d), ncol(d), data.table::uniqueN(d$uid)))
  invisible(list(data = d, kept_uids = keep, coverage = uid_cov,
                 sensor_coverage = uid_sensor,
                 sensor_availability = sensor_availability,
                 dropped_cols = all_na))
}

res <- apply_inclusion_criteria(
  outfile = "all_features_included.csv"
)

make_sensor_quality_distributions <- function(
    infile = "all_features_combined.csv",
    out_dir = results_dir,
    tag = tools::file_path_sans_ext(basename(infile))) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping sensor quality distribution plots.")
    return(invisible(NULL))
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    warning("scales is not installed; skipping sensor quality distribution plots.")
    return(invisible(NULL))
  }
  make_dirs(out_dir)
  
  d <- data.table::fread(file.path(CLEAN_DATA_DIR, infile), colClasses = list(character = "uid"))
  d[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  sensor_labels <- c(
    acc = "Accelerometer", hr = "Heart rate", steps = "Steps",
    battery = "Battery", convo = "Conversation", gps = "GPS / location",
    gps_home = "Time at home", sleep = "Sleep", sms = "SMS", calls = "Calls"
  )
  sensor_tags <- names(sensor_labels)
  
  sensor_cols <- function(sensor, pattern) {
    grep(sprintf("^%s_.*%s", sensor, pattern), names(d), value = TRUE)
  }
  
  value_long <- data.table::rbindlist(lapply(sensor_tags, function(sensor) {
    cols <- if (sensor == "sleep") {
      grep("^sleep_mean_", names(d), value = TRUE)
    } else {
      sensor_cols(sensor, "_mean_val$")
    }
    if (!length(cols)) return(NULL)
    d[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
    x <- data.table::melt(
      d[, c("uid", cols), with = FALSE],
      id.vars = "uid",
      variable.name = "feature",
      value.name = "value"
    )
    x[, sensor := sensor_labels[[sensor]]]
    x
  }), use.names = TRUE, fill = TRUE)
  
  availability_long <- data.table::rbindlist(lapply(sensor_tags, function(sensor) {
    cov_cols <- sensor_cols(sensor, "_coverage$")
    mean_cols <- if (sensor == "sleep") grep("^sleep_mean_", names(d), value = TRUE) else sensor_cols(sensor, "_mean_val$")
    if (length(cov_cols)) {
      data.table::data.table(
        uid = d$uid,
        sensor = sensor_labels[[sensor]],
        availability = rowMeans(as.matrix(d[, ..cov_cols]), na.rm = TRUE)
      )
    } else if (length(mean_cols)) {
      data.table::data.table(
        uid = d$uid,
        sensor = sensor_labels[[sensor]],
        availability = rowMeans(!is.na(as.data.frame(d[, ..mean_cols])))
      )
    } else {
      NULL
    }
  }), use.names = TRUE, fill = TRUE)
  
  summary_rows <- data.table::rbindlist(lapply(sensor_tags, function(sensor) {
    mean_cols <- if (sensor == "sleep") grep("^sleep_mean_", names(d), value = TRUE) else sensor_cols(sensor, "_mean_val$")
    cov_cols <- sensor_cols(sensor, "_coverage$")
    data.table::data.table(
      sensor = sensor_labels[[sensor]],
      n_features = length(mean_cols),
      n_coverage_features = length(cov_cols),
      pct_rows_any_value = if (length(mean_cols)) {
        mean(rowSums(!is.na(as.data.frame(d[, ..mean_cols]))) > 0) * 100
      } else NA_real_,
      median_availability = if (length(cov_cols)) {
        median(rowMeans(as.matrix(d[, ..cov_cols]), na.rm = TRUE), na.rm = TRUE)
      } else if (length(mean_cols)) {
        median(rowMeans(!is.na(as.data.frame(d[, ..mean_cols]))), na.rm = TRUE)
      } else NA_real_
    )
  }), use.names = TRUE)
  data.table::fwrite(summary_rows, file.path(out_dir, sprintf("sensor_quality_summary_%s.csv", tag)))
  
  if (nrow(value_long)) {
    value_long <- value_long[is.finite(value)]
    if (nrow(value_long)) {
      p_values <- ggplot2::ggplot(value_long, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 40, fill = "#7fb5a3", color = "white") +
        ggplot2::facet_wrap(~ sensor, scales = "free", ncol = 2) +
        ggplot2::labs(
          title = "Sensor feature value distributions",
          subtitle = sprintf("%s | mean-level features across EMA rows", infile),
          x = "Feature value",
          y = "Rows"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
      ggplot2::ggsave(file.path(out_dir, sprintf("sensor_value_distributions_%s.png", tag)),
                      p_values, width = 12, height = 14, dpi = 200, bg = "white")
    }
  }
  
  if (nrow(availability_long)) {
    availability_long <- availability_long[is.finite(availability)]
    p_avail <- ggplot2::ggplot(availability_long, ggplot2::aes(x = availability)) +
      ggplot2::geom_histogram(bins = 30, fill = "#5c7a8a", color = "white") +
      ggplot2::facet_wrap(~ sensor, scales = "free_y", ncol = 2) +
      ggplot2::scale_x_continuous(labels = scales::label_percent()) +
      ggplot2::labs(
        title = "Sensor data availability distributions",
        subtitle = sprintf("%s | coverage when available, otherwise non-missing feature proportion", infile),
        x = "Availability",
        y = "Rows"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    ggplot2::ggsave(file.path(out_dir, sprintf("sensor_availability_distributions_%s.png", tag)),
                    p_avail, width = 12, height = 14, dpi = 200, bg = "white")
  }
  
  cat(sprintf("Saved sensor QA plots and summary to: %s\n", out_dir))
  invisible(list(values = value_long, availability = availability_long, summary = summary_rows))
}

make_sensor_quality_distributions("all_features_combined.csv", tag = "combined")
if (file.exists(file.path(CLEAN_DATA_DIR, "all_features_included.csv"))) {
  make_sensor_quality_distributions("all_features_included.csv", tag = "included")
}

make_dynamic_index_quality_distributions <- function(
    infile = "all_features_included.csv",
    out_dir = results_dir,
    tag = tools::file_path_sans_ext(basename(infile))) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping dynamic-index QA plots.")
    return(invisible(NULL))
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    warning("scales is not installed; skipping dynamic-index QA plots.")
    return(invisible(NULL))
  }
  make_dirs(out_dir)
  d <- data.table::fread(file.path(CLEAN_DATA_DIR, infile), colClasses = list(character = "uid"))
  d[, uid := stringr::str_pad(as.character(uid), 4, pad = "0")]
  
  passive_tags <- c("gps_home", "battery", "steps", "convo", "calls", "sleep",
                    "acc", "gps", "sms", "hr")
  feature_meta <- function(cols) {
    rows <- lapply(cols, function(col) {
      sensor <- passive_tags[vapply(passive_tags, function(s) {
        grepl(paste0("^", s, "_"), col)
      }, logical(1))][1]
      if (is.na(sensor)) return(NULL)
      rest <- sub(paste0("^", sensor, "_"), "", col)
      parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
      if (sensor == "sleep") {
        idx <- "mean"
        win <- parts[length(parts)]
      } else {
        win <- parts[1]
        idx <- paste(parts[-1], collapse = "_")
      }
      idx <- sub("_val$", "", idx)
      data.table::data.table(feature = col, sensor = sensor, window = win, index = idx)
    })
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  }
  
  dynamic_cols <- grep("^(acc|hr|steps|battery|convo|gps|gps_home|sleep|sms|calls)_",
                       names(d), value = TRUE)
  dynamic_cols <- dynamic_cols[!grepl("_(coverage|n_obs)$", dynamic_cols)]
  meta <- feature_meta(dynamic_cols)
  if (!nrow(meta)) {
    warning(sprintf("No dynamic-index feature columns found in %s", infile))
    return(invisible(NULL))
  }
  
  d[, (meta$feature) := lapply(.SD, as.numeric), .SDcols = meta$feature]
  long <- data.table::melt(
    d[, c("uid", meta$feature), with = FALSE],
    id.vars = "uid",
    variable.name = "feature",
    value.name = "value"
  )
  long <- meta[long, on = "feature"]
  
  summary <- long[, .(
    n_rows = .N,
    n_participants = data.table::uniqueN(uid),
    pct_missing = mean(is.na(value)) * 100,
    pct_observed = mean(!is.na(value)) * 100,
    median = stats::median(value, na.rm = TRUE),
    iqr = stats::IQR(value, na.rm = TRUE),
    p01 = stats::quantile(value, 0.01, na.rm = TRUE, names = FALSE),
    p99 = stats::quantile(value, 0.99, na.rm = TRUE, names = FALSE)
  ), by = .(sensor, index, window, feature)][order(-pct_missing, sensor, index, window)]
  data.table::fwrite(summary, file.path(out_dir, sprintf("dynamic_index_quality_summary_%s.csv", tag)))
  
  missing_by_group <- summary[, .(
    n_features = .N,
    median_pct_missing = stats::median(pct_missing, na.rm = TRUE),
    max_pct_missing = max(pct_missing, na.rm = TRUE)
  ), by = .(sensor, index, window)][order(-median_pct_missing, sensor, index, window)]
  data.table::fwrite(missing_by_group, file.path(out_dir, sprintf("dynamic_index_missingness_by_group_%s.csv", tag)))
  
  p_missing <- ggplot2::ggplot(missing_by_group,
                               ggplot2::aes(x = window, y = median_pct_missing, fill = index)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ sensor, ncol = 2) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(scale = 1), limits = c(0, 100)) +
    ggplot2::labs(
      title = "Dynamic-index missingness by sensor, index, and window",
      subtitle = infile,
      x = "Window",
      y = "Median missingness across features",
      fill = "Index"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(out_dir, sprintf("dynamic_index_missingness_%s.png", tag)),
                  p_missing, width = 12, height = 14, dpi = 200, bg = "white")
  
  value_long <- long[!is.na(value) & is.finite(value)]
  if (nrow(value_long)) {
    p_values <- ggplot2::ggplot(value_long, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = 40, fill = "#7fb5a3", color = "white") +
      ggplot2::facet_grid(sensor ~ index, scales = "free") +
      ggplot2::labs(
        title = "Dynamic-index value distributions",
        subtitle = sprintf("%s | faceted by sensor and index", infile),
        x = "Feature value",
        y = "Rows"
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    ggplot2::ggsave(file.path(out_dir, sprintf("dynamic_index_value_distributions_%s.png", tag)),
                    p_values, width = 16, height = 18, dpi = 200, bg = "white")
  }
  
  cat(sprintf("Saved dynamic-index QA plots and summaries to: %s\n", out_dir))
  invisible(list(long = long, summary = summary, missing_by_group = missing_by_group))
}

inspect_dynamic_index_merge <- function(
    feature_file = "all_features_combined.csv",
    out_dir = results_dir) {
  
  make_dirs(out_dir)
  feature_path <- file.path(CLEAN_DATA_DIR, feature_file)
  if (!file.exists(feature_path)) stop(sprintf("Missing feature file: %s", feature_path))
  
  fm <- data.table::fread(feature_path, nrows = 1)
  combined_path <- file.path(CLEAN_DATA_DIR, "all_features_combined.csv")
  combined_names <- if (file.exists(combined_path)) names(data.table::fread(combined_path, nrows = 1)) else character(0)
  tags <- c("acc", "hr", "steps", "battery", "convo", "gps", "gps_home", "sleep", "sms", "calls")
  tag_feature_count <- function(tag, nms) {
    if (!length(nms)) return(NA_integer_)
    pattern <- if (tag == "gps") "^gps_(?!home_)" else sprintf("^%s_", tag)
    sum(grepl(pattern, nms, perl = TRUE))
  }
  expected <- data.table::rbindlist(lapply(tags, function(tag) {
    fp <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", tag))
    file_exists <- file.exists(fp)
    n_rows <- NA_integer_
    n_pts <- NA_integer_
    source_cols <- 0L
    
    if (file_exists) {
      src <- data.table::fread(fp, colClasses = list(character = "uid"))
      n_rows <- nrow(src)
      if ("uid" %in% names(src)) n_pts <- data.table::uniqueN(src$uid)
      reshaped <- reshape_dynamic_indices(fp, tag)
      if (!is.null(reshaped)) {
        source_cols <- length(setdiff(names(reshaped), c("uid", "date_rep", "ema_time")))
      }
    }
    
    merged_cols <- tag_feature_count(tag, names(fm))
    combined_cols <- tag_feature_count(tag, combined_names)
    merge_status <- if (!file_exists) {
      "missing_dynamic_file"
    } else if (source_cols == 0L) {
      "no_reshaped_columns"
    } else if (merged_cols >= source_cols) {
      "ok"
    } else if (!identical(feature_file, "all_features_combined.csv") &&
               !is.na(combined_cols) && combined_cols >= source_cols && merged_cols > 0L) {
      "ok_after_filter_drop"
    } else if (merged_cols > 0L) {
      "partial"
    } else {
      "not_merged"
    }
    data.table::data.table(
      sensor = tag,
      dynamic_file_exists = file_exists,
      dynamic_rows = n_rows,
      dynamic_participants = n_pts,
      expected_feature_columns = source_cols,
      combined_feature_columns = combined_cols,
      merged_feature_columns = merged_cols,
      merge_status = merge_status
    )
  }), use.names = TRUE)
  
  out_csv <- file.path(out_dir, sprintf("dynamic_index_merge_inspection_%s.csv",
                                        tools::file_path_sans_ext(basename(feature_file))))
  data.table::fwrite(expected, out_csv)
  print(expected)
  if (any(!expected$merge_status %in% c("ok", "ok_after_filter_drop"))) {
    warning(sprintf("Some dynamic index sources are missing or incomplete. See: %s", out_csv))
  } else {
    cat(sprintf("All dynamic index sources are represented in %s\n", feature_file))
  }
  invisible(expected)
}

inspect_dynamic_index_merge("all_features_combined.csv")
if (file.exists(file.path(CLEAN_DATA_DIR, "all_features_included.csv"))) {
  inspect_dynamic_index_merge("all_features_included.csv")
}
}

if (RUN_MAIN_PIPELINE || RUN_FINAL_MANUSCRIPT_OUTPUTS ||
    RUN_ALL_PT_CONTINUOUS_TIMESERIES || RUN_ALL_PT_HURDLE_TIMESERIES) {
cat("\n===== SECTION: MAIN PREDICTION PIPELINE =====\n")
################################################################################
# SI Prediction Pipeline
################################################################################
# ==============================================================================
# CONFIG
# ==============================================================================
BASE_PATH  <- resolve_project_base()
cat(sprintf("Project base: %s\n", BASE_PATH))
CODE_DIR   <- file.path(BASE_PATH, "1. Code")
RAW_DATA_DIR <- file.path(BASE_PATH, "2. Raw data")
CLEAN_DATA_DIR <- file.path(BASE_PATH, "3. Cleaned data")
CLEAN_SENSOR_DIR <- file.path(CLEAN_DATA_DIR, "per_sensors")
CLEAN_UID_DIR <- file.path(CLEAN_DATA_DIR, "per_uid")
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
EMA_FILE   <- "3._TD_cohort_EMA_FINAL.csv"
MODEL_FILE <- "all_features_included.csv"
DEMO_FILE  <- "9._FinalDemographics.csv"
N_CORES    <- 12
NTREE      <- 700
MIN_FOLDS  <- 10
REBUILD_EMA_FEATURES <- FALSE
REBUILD_MODEL_CACHE <- FALSE
RERUN_MODEL_RESULTS <- FALSE
RERUN_F2_HURDLE <- TRUE
RERUN_F2_REFERENCE_STYLE <- FALSE
RUN_CONTINUOUS_MODELS <- TRUE
RUN_HURDLE_MODELS <- TRUE
RERUN_COMPLETED_CONFIGS <- FALSE
RUN_MODEL_DIAGNOSTICS <- FALSE
SENSOR_FOLDERS <- list(
  acc="acceleration", hr="hr", steps="steps", convo="convo", gps="raw_gps_clean",
  gps_home="raw_gps_clean", sms="sms", calls="calls", battery="battery",
  sleep="sleep"
)
CV_VARIANTS <- list(
  main  = list(train_days=45, val_days=10, test_days=1, ema_per_day=3)
)
results_dir       <- file.path(RESULTS_ROOT_DIR, "report")
plots_dir         <- results_dir
report_dir        <- results_dir
report_all_pt_dir <- file.path(report_dir, "all_PT")
report_all_pt_continuous_dir <- file.path(report_all_pt_dir, "continuous")
report_all_pt_hurdle_dir <- file.path(report_all_pt_dir, "hurdle")
model_config_cache_dir <- file.path(RESULTS_ROOT_DIR, "model_config_cache")
ema_features_path <- file.path(CLEAN_DATA_DIR, "ema_features.csv")
make_dirs(
  CLEAN_DATA_DIR,
  CLEAN_SENSOR_DIR,
  CLEAN_UID_DIR,
  results_dir,
  plots_dir,
  report_dir,
  report_all_pt_dir,
  report_all_pt_continuous_dir,
  report_all_pt_hurdle_dir,
  model_config_cache_dir
)
# ==============================================================================
# SECTION 0: PACKAGES + SHARED HELPERS
# ==============================================================================
required_packages <- c(
  "tidyverse", "data.table", "lubridate", "stringr", "tidyr",
  "TSEntropies", "bit64", "parallel", "future.apply", "furrr", "progress",
  "RANN", "geosphere", "caret", "randomForest", "foreach", "doParallel",
  "effsize", "patchwork", "gbm", "pROC",
  "ggplot2", "ggridges", "scales"
)
require_packages_or_stop(required_packages, "the main prediction pipeline")
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(lubridate); library(stringr)
  library(tidyr); library(TSEntropies); library(bit64); library(progress)
  library(future.apply); library(furrr); library(RANN); library(geosphere)
  library(caret); library(randomForest); library(foreach); library(doParallel)
  library(effsize); library(patchwork); library(gbm); library(pROC)
  library(ggplot2); library(ggridges); library(scales)
})
setwd(BASE_PATH)
# Per-person min-max scaling
minmax_scale <- function(x) {
  mn <- min(x, na.rm = TRUE); mx <- max(x, na.rm = TRUE)
  if (is.na(mn) || is.na(mx) || is.infinite(mn) || is.infinite(mx) || mx == mn) return(x)
  (x - mn) / (mx - mn)
}
# Logger factory
make_logger <- function(log_path) {
  function(msg) {
    out <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
    cat(out, "\n"); cat(out, "\n", file = log_path, append = TRUE)
  }
}
# EMA lookup builder (used by run_sensor_pipeline)
build_ema_lookup <- function(ema_path) {
  read.csv(ema_path) %>%
    filter(!is.na(ema_day)) %>%
    mutate(uid        = str_extract(mlife_id, "^[^@]+"),
           ema_day    = as.POSIXct(ema_day, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
           anchor_end = ema_day - hours(4)) %>%
    select(uid, ema_day, date_rep, ema_time, anchor_end) %>%
    as.data.table() %>% { setkey(., uid); . }
}
# Cache short-circuit
run_if_missing <- function(sensor_tag, run_fn) {
  out_csv <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", sensor_tag))
  if (file.exists(out_csv)) {
    cat(sprintf("%-10s already done — skipping\n", sensor_tag)); return(invisible(NULL))
  }
  run_fn()
}
# within-fold per-feature min-max (learned on TRAIN rows only, applied to all splits)
apply_minmax_fold <- function(ref_df, target_df, cols) {
  for (cc in cols) {
    if (!cc %in% names(ref_df) || !cc %in% names(target_df)) next
    if (!is.numeric(ref_df[[cc]]) || !is.numeric(target_df[[cc]])) next
    mn <- suppressWarnings(min(ref_df[[cc]], na.rm = TRUE))
    mx <- suppressWarnings(max(ref_df[[cc]], na.rm = TRUE))
    if (is.finite(mn) && is.finite(mx) && mx > mn)
      target_df[[cc]] <- (target_df[[cc]] - mn) / (mx - mn)
  }
  target_df
}
candidate_predictors <- function(pt) {
  meta <- c("uid", "day", outcome, "mlife_id", "ema_day", "ema_timestamp",
            "date_rep", "ema_time")
  feats <- setdiff(names(pt), meta)
  feats <- feats[vapply(feats, function(f) is.numeric(pt[[f]]) && !all(is.na(pt[[f]])), logical(1))]
  feats
}
# ------------------------------------------------------------------------------
# run_sensor_pipeline — time-based windows + P90-P10 relative position
# ------------------------------------------------------------------------------
run_sensor_pipeline <- function(file_list, ema_lookup, sensor_tag,
                                read_fn, signal_col, output_dir, log_fn,
                                windows = list(hr1  = 3600,
                                               hr4  = 14400,
                                               day1 = 86400,
                                               day7 = 604800)) {
  make_dirs(output_dir)
  pb <- progress::progress_bar$new(
    format = "  [:bar] :current/:total | :uid | eta: :eta",
    total = length(file_list), clear = FALSE, width = 90)
  log_fn(sprintf("Pipeline start: %s", sensor_tag))
  t0 <- proc.time()
  
  for (i in seq_along(file_list)) {
    fp <- file_list[i]
    uid_this <- str_extract(basename(fp), "^[^@]+")
    out_path <- file.path(output_dir,
                          sprintf("%s_dynamic_indices_%s.csv", uid_this, sensor_tag))
    pb$tick(tokens = list(uid = uid_this))
    if (file.exists(out_path)) {
      log_fn(sprintf("SKIP %s — already done", uid_this)); next
    }
    
    tryCatch({
      dt <- read_fn(fp)
      if (is.null(dt) || nrow(dt) == 0L) {
        log_fn(sprintf("SKIP %s — empty", uid_this)); return(NULL)
      }
      ema_this <- ema_lookup[uid == uid_this]
      if (nrow(ema_this) == 0L) {
        log_fn(sprintf("SKIP %s — no EMA", uid_this)); return(NULL)
      }
      
      max_t <- as.numeric(max(ema_this$anchor_end))
      min_t <- as.numeric(min(ema_this$anchor_end)) - 604800
      dt <- dt[as.numeric(day) >= min_t & as.numeric(day) <= max_t]
      if (nrow(dt) == 0L) {
        log_fn(sprintf("SKIP %s — no range", uid_this)); return(NULL)
      }
      
      data.table::setkey(dt, day)
      day_vec <- as.numeric(dt$day)
      sig_vec <- dt[[signal_col]]
      rm(dt); gc()
      
      n_ema <- nrow(ema_this)
      results <- vector("list", n_ema * length(windows))
      k <- 1L
      
      for (j in seq_len(n_ema)) {
        ae <- as.numeric(ema_this$anchor_end[j])
        
        # current_x = mean of obs in [ae - 1800, ae)
        cur_lo <- findInterval(ae - 1800, day_vec)
        cur_hi <- findInterval(ae,        day_vec)
        x_cur_window <- if (cur_hi > cur_lo) sig_vec[(cur_lo + 1L):cur_hi] else numeric(0)
        x_current <- mean_thr(x_cur_window)
        
        # Each window selected by TIME range
        for (wn in names(windows)) {
          win_sec <- windows[[wn]]
          win_lo  <- findInterval(ae - win_sec, day_vec)
          win_hi  <- findInterval(ae,           day_vec)
          x <- if (win_hi > win_lo) sig_vec[(win_lo + 1L):win_hi] else numeric(0)
          x_clean <- x[!is.na(x)]
          
          results[[k]] <- data.frame(
            uid        = ema_this$uid[j],
            ema_day    = ema_this$ema_day[j],
            date_rep   = ema_this$date_rep[j],
            ema_time   = ema_this$ema_time[j],
            anchor_end = ema_this$anchor_end[j],
            window     = wn,
            n_obs      = length(x_clean),
            mean_val   = mean_thr(x_clean),
            mad_val    = mad_thr(x_clean),
            p90p10     = p90p10_relpos(x_clean, x_current),
            entropy    = sample_entropy(x_clean)
          )
          k <- k + 1L
        }
        if (j %% 20 == 0) gc()
      }
      
      out <- data.table::rbindlist(results)
      data.table::fwrite(out, out_path)
      rm(day_vec, sig_vec, results, out); gc()
      log_fn(sprintf("DONE %s — %d EMA rows", uid_this, n_ema))
    }, error = function(e) log_fn(sprintf("ERROR %s: %s", uid_this, e$message)))
  }
  
  log_fn(sprintf("Pipeline done: %s — %.1f min",
                 sensor_tag, (proc.time() - t0)["elapsed"] / 60))
  
  done_files <- list.files(output_dir,
                           pattern = sprintf("_dynamic_indices_%s\\.csv$", sensor_tag),
                           full.names = TRUE)
  final <- data.table::rbindlist(lapply(done_files, data.table::fread))
  out_csv <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", sensor_tag))
  data.table::fwrite(final, out_csv)
  log_fn(sprintf("Saved: dynamic_indices_%s.csv (%d rows, %d participants)",
                 sensor_tag, nrow(final), uniqueN(final$uid)))
  invisible(final)
}
# ==============================================================================
# SECTION 1: EMA PROCESSING
# ==============================================================================
if (file.exists(ema_features_path) && !REBUILD_EMA_FEATURES) {
  cat(sprintf("\n===== EMA features cached — using %s =====\n", basename(ema_features_path)))
} else {
  cat("\n===== Rebuild EMA features =====\n")
  # 1. Derive compliant_uids from existing sensor outputs (restored)
  sensor_csvs <- list.files(CLEAN_SENSOR_DIR, pattern = "^dynamic_indices_.*\\.csv$", full.names = TRUE)
  cat(sprintf("Found %d sensor CSVs\n", length(sensor_csvs)))
  if (length(sensor_csvs) > 0) {
    all_sensor_uids <- unique(unlist(lapply(sensor_csvs, function(fp) {
      uids <- data.table::fread(fp, select = "uid", nThread = 2)$uid
      as.character(uids)
    })))
    all_sensor_uids <- stringr::str_pad(all_sensor_uids, width = 4, pad = "0")
    compliant_uids <- sort(unique(all_sensor_uids))
  } else {
    cat("WARNING: no sensor CSVs found yet — compliant_uids set from EMA file only.\n")
    compliant_uids <- character(0)
  }
  cat(sprintf("Compliant UIDs from sensor data: %d\n", length(compliant_uids)))
  write.csv(data.frame(uid = compliant_uids),
            file.path(CLEAN_DATA_DIR, "compliant_uids.csv"), row.names = FALSE)
  # 2. Load and clean EMA
  ema_raw <- read.csv(file.path(RAW_DATA_DIR, EMA_FILE))
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
    select(-starts_with("RespTime_")) %>%
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
  # 3. Compute dynamic indices (uses Section 0 helpers: mean_thr, mad_thr, p90p10_relpos)
  vars <- c(paste0("q", 1:9), "g1", "g2")
  dt <- as.data.table(ema_complete)
  dt[, ema_num    := as.numeric(ema_day)]
  dt[, anchor_end := ema_num - 3600 * 4]
  dt[, row_id     := .I]
  setkey(dt, mlife_id, ema_num)
  src <- dt[, c("row_id","mlife_id","ema_num", vars), with = FALSE]
  setkey(src, mlife_id, ema_num)
  n_rows <- nrow(dt)
  d1_mat <- matrix(NA_real_, nrow = n_rows, ncol = length(vars),
                   dimnames = list(NULL, paste0(vars, "_mean_d1")))
  d7_mat <- matrix(NA_real_, nrow = n_rows, ncol = length(vars) * 2,
                   dimnames = list(NULL, c(paste0(vars, "_mean_d7"),
                                           paste0(vars, "_mad_d7"))))
  pb_ema <- progress::progress_bar$new(
    format = "  [:bar] :current/:total | :percent | eta: :eta",
    total = n_rows, clear = FALSE, width = 80)
  for (i in seq_len(n_rows)) {
    pb_ema$tick()
    uid_i <- dt$mlife_id[i]
    ae    <- dt$anchor_end[i]
    
    # 1-day window: mean only
    p1 <- src[mlife_id == uid_i & ema_num >= (ae - 86400) & ema_num < ae]
    if (nrow(p1) > 0) {
      d1_mat[i, ] <- sapply(vars, function(v) mean_thr(p1[[v]]))
    }
    
    # 7-day window: mean + MAD (no P90-P10 for EMA)
    p7 <- src[mlife_id == uid_i & ema_num >= (ae - 604800) & ema_num < ae]
    if (nrow(p7) > 0) {
      d7_mat[i, seq_along(vars)]                  <- sapply(vars, function(v) mean_thr(p7[[v]]))
      d7_mat[i, seq_along(vars) + length(vars)]   <- sapply(vars, function(v) mad_thr(p7[[v]]))
    }
    
    if (i %% 1000 == 0) gc()
  }
  dt <- cbind(dt, as.data.table(d1_mat), as.data.table(d7_mat))
  dt[, c("ema_num","anchor_end","row_id") := NULL]
  # 4. Filter to compliant uids (if any) and save
  ema_features <- as.data.frame(dt) %>%
    mutate(uid = stringr::str_pad(str_extract(mlife_id, "^[^@]+"), width = 4, pad = "0"))
  if (length(compliant_uids) > 0) {
    ema_features <- ema_features %>% filter(uid %in% compliant_uids)
  }
  data.table::fwrite(ema_features, ema_features_path)
  cat(sprintf("Saved: ema_features.csv (%d rows | %d PTs)\n",
              nrow(ema_features), uniqueN(ema_features$uid)))
  summary(ema_features)
}
# ==============================================================================
# SECTION 2: DIGITAL PHENOTYPE PROCESSING
# ==============================================================================
# Digital phenotype cleaning, EMA feature construction, passive dynamic-index
# processing, and analytic inclusion criteria are implemented above in this same
# single script. Stress is dropped because it was not retained for analysis.
# Participants are retained if they meet EMA compliance, passive compliance
# across continuously sampled streams (acc, steps, HR, battery, GPS distance,
# GPS time at home), and SI-count criteria. The exact retained N is printed and
# saved when apply_inclusion_criteria() is run.
# ==============================================================================
# SECTION 3: RANDOM FOREST PREDICTION
# ==============================================================================
# cor (Pearson r), r2 (coefficient of determination), mse, mae, rmse, auc.
# Figures and R²>0 subsets use the TRUE r2 column (not sign(cor)*cor^2).
# ==============================================================================
cat("\n===== SECTION 3: Random Forest Prediction =====\n")
set.seed(12361488)
model_source_path <- file.path(CLEAN_DATA_DIR, MODEL_FILE)
data_model_path <- file.path(
  CLEAN_DATA_DIR,
  sprintf("%s_filtered.csv", tools::file_path_sans_ext(basename(MODEL_FILE)))
)
cat(sprintf("Model source file: %s\n", basename(model_source_path)))
cat(sprintf("Model cache file: %s\n", basename(data_model_path)))
parse_date_rep <- function(x) {
  raw <- as.character(x)
  out <- suppressWarnings(as.Date(raw, format = "%Y-%m-%d"))
  miss <- is.na(out)
  out[miss] <- suppressWarnings(as.Date(raw[miss], format = "%m/%d/%Y"))
  miss <- is.na(out)
  out[miss] <- suppressWarnings(lubridate::ymd(raw[miss]))
  miss <- is.na(out)
  out[miss] <- suppressWarnings(lubridate::mdy(raw[miss]))
  out
}
n_pt <- function(d) dplyr::n_distinct(d$uid[!is.na(d$uid)])
log_model_step <- function(label, d) {
  cat(sprintf("  %-38s %7d rows | %4d participants\n", label, nrow(d), n_pt(d)))
}
lag_previous_observed <- function(x) {
  out <- rep(NA_real_, length(x))
  last <- NA_real_
  for (i in seq_along(x)) {
    out[i] <- last
    if (!is.na(x[i])) last <- x[i]
  }
  out
}
model_cache_stale <- file.exists(data_model_path) &&
  file.exists(model_source_path) &&
  file.info(model_source_path)$mtime > file.info(data_model_path)$mtime
if (REBUILD_MODEL_CACHE || model_cache_stale || !file.exists(data_model_path)) {
  if (file.exists(data_model_path)) {
    cat(sprintf("Rebuilding %s from %s\n", basename(data_model_path), basename(model_source_path)))
  }
  data_raw <- data.table::fread(model_source_path, colClasses = list(character = "uid")) %>%
    as.data.frame()
  log_model_step("input included feature matrix", data_raw)
  
  data_filtered <- data_raw %>%
    mutate(uid = stringr::str_pad(as.character(uid), width = 4, pad = "0"),
           date_rep = parse_date_rep(date_rep),
           ema_time = factor(ema_time, levels = c("morning", "afternoon", "evening"))) %>%
    filter(!is.na(date_rep))
  log_model_step("after date parsing", data_filtered)
  
  data_filtered <- data_filtered %>%
    group_by(uid) %>% filter(!all(is.na(date_rep))) %>%
    complete(date_rep = seq(min(date_rep, na.rm = TRUE),
                            max(date_rep, na.rm = TRUE), by = "day"),
             ema_time = c("morning", "afternoon", "evening")) %>%
    arrange(uid, date_rep, ema_time) %>% ungroup()
  log_model_step("after complete daily EMA grid", data_filtered)
  
  # Lag all 11 EMA items by one prompt within each participant.
  ema_items <- c("q1","q2","q3","q4","q5","q6","q7","q8","q9","g1","g2")
  
  data_after_burnin <- data_filtered %>%
    group_by(uid) %>%
    mutate(first_date = min(date_rep, na.rm = TRUE)) %>%
    filter(date_rep >= first_date + 7) %>%
    select(-first_date) %>%
    ungroup()
  log_model_step("after 7-day burn-in", data_after_burnin)
  
  data_after_lag <- data_after_burnin %>%
    group_by(uid) %>%
    arrange(date_rep, ema_time, .by_group = TRUE) %>%
    mutate(across(all_of(ema_items), lag_previous_observed, .names = "{.col}_lag")) %>%
    ungroup() %>%
    # Forward-fill q3 and q3_lag from morning to afternoon/evening of same day
    group_by(uid, date_rep) %>%
    mutate(q3     = dplyr::first(q3[ema_time     == "morning"], default = NA_real_),
           q3_lag = dplyr::first(q3_lag[ema_time == "morning"], default = NA_real_)) %>%
    ungroup()
  log_model_step("after EMA lag construction", data_after_lag)
  
  data_model_full <- data_after_lag %>%
    filter(!is.na(q9))
  log_model_step("after current q9 filter", data_model_full)
  retention_report <- data_after_lag %>%
    group_by(uid) %>%
    summarise(
      rows_after_lag = n(),
      observed_q9 = sum(!is.na(q9)),
      observed_q9_with_lag = sum(!is.na(q9) & !is.na(q9_lag)),
      model_rows = sum(!is.na(q9)),
      .groups = "drop"
    )
  retention_path <- file.path(results_dir, "model_filter_retention_by_uid.csv")
  write.csv(retention_report, retention_path, row.names = FALSE)
  cat(sprintf("Saved model retention report: %s\n", retention_path))
  
  data_model_full <- data_model_full %>%
    rename(day = ema_time) %>%
    select(-any_of(c("date_rep", "ema_time"))) %>%
    # Outcome rescaling only (constant /100 -> no leakage). Per-feature min-max is
    # applied WITHIN each training split inside the fit functions (apply_minmax_fold).
    mutate(q9     = q9     / 100,
           q9_lag = ifelse(is.na(q9_lag), NA_real_, q9_lag / 100)) %>%
    ungroup()
  
  write.csv(data_model_full, data_model_path, row.names = FALSE)
  cat(sprintf("Saved: %s (%d rows | %d participants)\n",
              basename(data_model_path),
              nrow(data_model_full), uniqueN(data_model_full$uid)))
}
data_model_base <- data.table::fread(data_model_path, colClasses = list(character = "uid")) %>%
  as.data.frame()
# Target = current q9 (PATCHED: was q9_lag)
outcome <- "q9"
# Sliding-window CV factory (unchanged)
make_sliding_cv <- function(train_days, val_days, test_days, ema_per_day = 3) {
  function(data) {
    n <- nrow(data)
    train_n <- train_days * ema_per_day
    val_n   <- val_days   * ema_per_day
    test_n  <- test_days  * ema_per_day
    slide_n <- test_n
    if (n < train_n + val_n + test_n) return(NULL)
    trains <- list(); vals <- list(); tests <- list(); k <- 1L; start <- 1L
    while ((start + train_n + val_n + test_n - 1) <= n) {
      trains[[k]] <- start:(start + train_n - 1)
      vals[[k]]   <- (start + train_n):(start + train_n + val_n - 1)
      tests[[k]]  <- (start + train_n + val_n):(start + train_n + val_n + test_n - 1)
      k <- k + 1L; start <- start + slide_n
    }
    if (length(trains) == 0) return(NULL)
    list(trains = trains, vals = vals, tests = tests)
  }
}
summarize_si_non_si_counts <- function(data, variant_name = "main") {
  cv_params <- CV_VARIANTS[[variant_name]]
  cv_fn <- make_sliding_cv(cv_params$train_days, cv_params$val_days,
                           cv_params$test_days, cv_params$ema_per_day)
  test_y <- unlist(lapply(split(data, data$uid), function(pt) {
    sl <- cv_fn(pt)
    if (is.null(sl) || !length(sl$tests)) return(numeric(0))
    pt[[outcome]][unlist(sl$tests)]
  }))
  out <- data.frame(
    variant = variant_name,
    analytic_participants = dplyr::n_distinct(data$uid),
    analytic_rows = nrow(data),
    analytic_SI_obs = sum(data[[outcome]] > 0, na.rm = TRUE),
    analytic_nonSI_obs = sum(data[[outcome]] == 0, na.rm = TRUE),
    test_SI_obs = sum(test_y > 0, na.rm = TRUE),
    test_nonSI_obs = sum(test_y == 0, na.rm = TRUE),
    test_total_obs = length(test_y),
    stringsAsFactors = FALSE
  )
  out_path <- file.path(results_dir, sprintf("si_non_si_counts_%s.csv", variant_name))
  write.csv(out, out_path, row.names = FALSE)
  cat(sprintf(
    "\nSI/non-SI counts [%s]: analytic SI=%s, analytic non-SI=%s | test SI=%s, test non-SI=%s\n",
    variant_name,
    format(out$analytic_SI_obs, big.mark = ","),
    format(out$analytic_nonSI_obs, big.mark = ","),
    format(out$test_SI_obs, big.mark = ","),
    format(out$test_nonSI_obs, big.mark = ",")
  ))
  cat(sprintf("Saved SI/non-SI count summary: %s\n", out_path))
  invisible(out)
}
si_non_si_counts_main <- summarize_si_non_si_counts(data_model_base, "main")
get_passive_mean_cols <- function(sensor, window, data) {
  if (sensor == "sleep") grep(sprintf("^sleep_mean_.*_%s$", window),    names(data), value = TRUE)
  else                   grep(sprintf("^%s_%s_(mean_val|sum_val)$", sensor, window), names(data), value = TRUE)
}
get_passive_dynamic_cols <- function(sensor, window, data) {
  if (sensor == "sleep") grep(sprintf("^sleep_mad_.*_%s$", window),                   names(data), value = TRUE)
  else                   grep(sprintf("^%s_%s_(mad_val|p90p10|entropy)$", sensor, window), names(data), value = TRUE)
}
# All 11 EMA predictors are LAGGED (PATCHED).
ema_mean <- c("q1_lag","q2_lag","q3_lag","q4_lag","q5_lag",
              "q6_lag","q7_lag","q8_lag","q9_lag","g1_lag","g2_lag")
# Window features are unchanged (already pre-anchor)
ema_dyn  <- grep("^(q|g)[0-9]+_(mean_d|mad_d|p90p10_d)", names(data_model_base), value = TRUE)
sensors  <- c("acc","hr","steps","convo","gps","gps_home","sms","calls","battery","sleep")
pmean    <- function(w, d) unlist(lapply(sensors, get_passive_mean_cols,    window = w, data = d))
pdyn     <- function(w, d) unlist(lapply(sensors, get_passive_dynamic_cols, window = w, data = d))
configs <- list(
  EMA_mean                     = function(w, d) ema_mean,
  EMA_mean_dynamic             = function(w, d) c(ema_mean, ema_dyn),
  Opt_Passive_mean             = function(w, d) pmean(w, d),
  Opt_Passive_mean_dynamic     = function(w, d) c(pmean(w, d), pdyn(w, d)),
  EMA_Opt_Passive_mean         = function(w, d) c(ema_mean, pmean(w, d)),
  EMA_Opt_Passive_mean_dynamic = function(w, d) c(ema_mean, ema_dyn, pmean(w, d), pdyn(w, d))
)
MODEL_CONFIGS <- names(configs)  # The six feature-set models; baseline is a benchmark, not counted as a model.
WINDOWS    <- c("hr1","hr4","day1","day7")
MTRY_GRID  <- c(2, 5, 10)
NTREE_GRID <- c(200, 500)
GBM_CLS    <- data.frame(n.trees=100, interaction.depth=3, shrinkage=0.1, n.minobsinnode=10)
GBM_REG    <- data.frame(n.trees=100, interaction.depth=3, shrinkage=0.1, n.minobsinnode=5)
is_passive_feature <- function(x) {
  grepl("^(acc|hr|steps|convo|gps_home|gps|sms|calls|battery|sleep)_", x)
}
drop_near_zero_variance_features <- function(cols, training_data,
                                             freqCut = 19,
                                             uniqueCut = 10) {
  cols <- intersect(cols, names(training_data))
  if (length(cols) < 1) return(cols)
  numeric_cols <- cols[vapply(cols, function(cc) is.numeric(training_data[[cc]]), logical(1))]
  non_numeric_cols <- setdiff(cols, numeric_cols)
  if (length(numeric_cols) < 1) return(cols)
  nzv <- tryCatch(
    caret::nearZeroVar(training_data[, numeric_cols, drop = FALSE],
                       freqCut = freqCut, uniqueCut = uniqueCut),
    error = function(e) integer(0)
  )
  if (length(nzv)) numeric_cols <- numeric_cols[-nzv]
  unique(c(non_numeric_cols, numeric_cols))
}
# ------------------------------------------------------------------------------
reg_metrics <- function(pred, true) {
  ok   <- !is.na(pred) & !is.na(true)
  pred <- pred[ok]; true <- true[ok]
  n <- length(true)
  if (n < 1) return(list(cor=NA_real_, r2=NA_real_, mse=NA_real_,
                         mae=NA_real_, rmse=NA_real_))
  resid  <- pred - true
  mse    <- mean(resid^2)
  mae    <- mean(abs(resid))
  rmse   <- sqrt(mse)
  cor_v  <- if (n >= 2) suppressWarnings(stats::cor(pred, true)) else NA_real_
  ss_tot <- sum((true - mean(true))^2)
  r2     <- if (ss_tot > 0) 1 - sum(resid^2) / ss_tot else NA_real_
  list(cor=cor_v, r2=r2, mse=mse, mae=mae, rmse=rmse)
}
# [UPDATED] Row schema — carries r2, mse, mae alongside cor, rmse, auc.
empty_row <- function(uid, status) data.frame(
  Participant=stringr::str_pad(as.character(uid), width = 4, pad = "0"),
  cor=NA_real_, r2=NA_real_, mse=NA_real_, mae=NA_real_,
  rmse=NA_real_, auc=NA_real_, sensitivity=NA_real_, specificity=NA_real_,
  balanced_acc=NA_real_, n_obs=0, n_SI=NA_integer_, n_nonSI=NA_integer_,
  n_folds=0, status=status)
ok_row <- function(uid, n_folds, n_obs,
                   cor_v=NA, r2_v=NA, mse_v=NA, mae_v=NA, rmse_v=NA, auc_v=NA,
                   sensitivity_v=NA, specificity_v=NA, balanced_acc_v=NA,
                   n_SI_v=NA_integer_, n_nonSI_v=NA_integer_)
  data.frame(
    Participant=stringr::str_pad(as.character(uid), width = 4, pad = "0"),
    cor=cor_v, r2=r2_v, mse=mse_v, mae=mae_v,
    rmse=rmse_v, auc=auc_v, sensitivity=sensitivity_v, specificity=specificity_v,
    balanced_acc=balanced_acc_v, n_obs=n_obs, n_SI=n_SI_v, n_nonSI=n_nonSI_v,
    n_folds=n_folds, status="ok")
binary_metrics <- function(true_bin, p_si, threshold = 0.5) {
  ok <- !is.na(true_bin) & !is.na(p_si)
  true_bin <- as.integer(true_bin[ok])
  pred_bin <- as.integer(p_si[ok] >= threshold)
  if (!length(true_bin) || length(unique(true_bin)) < 2) {
    return(list(sensitivity=NA_real_, specificity=NA_real_, balanced_acc=NA_real_))
  }
  tp <- sum(pred_bin == 1 & true_bin == 1)
  fn <- sum(pred_bin == 0 & true_bin == 1)
  tn <- sum(pred_bin == 0 & true_bin == 0)
  fp <- sum(pred_bin == 1 & true_bin == 0)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  list(sensitivity=sens, specificity=spec, balanced_acc=mean(c(sens, spec), na.rm=TRUE))
}
ensure_result_metric_cols <- function(d) {
  if ("Participant" %in% names(d)) {
    d$Participant <- stringr::str_pad(as.character(d$Participant), width = 4, pad = "0")
  }
  for (cc in c("sensitivity", "specificity", "balanced_acc")) {
    if (!cc %in% names(d)) d[[cc]] <- NA_real_
  }
  for (cc in c("n_SI", "n_nonSI")) {
    if (!cc %in% names(d)) d[[cc]] <- NA_integer_
  }
  d
}
# --- CONTINUOUS ---
fit_continuous <- function(pt, fb, cv_fn) {
  feats <- candidate_predictors(pt)
  if (length(feats)<1 || var(pt[[outcome]],na.rm=TRUE) %in% c(0,NA))
    return(empty_row(pt$uid[1], "skipped"))
  sl <- cv_fn(pt); if (is.null(sl) || length(sl$trains)<1) return(empty_row(pt$uid[1], "no_slices"))
  preds <- list()
  for (i in seq_along(sl$trains)) {
    tr <- pt[sl$trains[[i]],]; vl <- pt[sl$vals[[i]],]; te <- pt[sl$tests[[i]],]
    scale_cols <- feats   # within-fold min-max
    tr_ref <- tr                                               # original train values
    tr <- apply_minmax_fold(tr_ref, tr, scale_cols)
    vl <- apply_minmax_fold(tr_ref, vl, scale_cols)
    te <- apply_minmax_fold(tr_ref, te, scale_cols)
    best <- list(score=Inf)
    for (w in WINDOWS) {
      pv <- intersect(fb(w, data_model_base), feats); if (length(pv)<2) next
      pv <- drop_near_zero_variance_features(pv, tr, freqCut = 19, uniqueCut = 10); if (length(pv)<2) next
      pp <- tryCatch(preProcess(tr[,pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
      if (is.null(pp)) next
      xtr <- tryCatch(predict(pp, tr[,pv,drop=FALSE]), error=function(e) NULL)
      xvl <- tryCatch(predict(pp, vl[,pv,drop=FALSE]), error=function(e) NULL)
      if (is.null(xtr) || is.null(xvl)) next
      d_tr <- cbind(xtr, .y=tr[[outcome]]); names(d_tr)[ncol(d_tr)] <- outcome
      for (mt in MTRY_GRID) for (nt in NTREE_GRID) {
        s <- tryCatch({
          m <- randomForest(as.formula(paste(outcome,"~ .")), data=d_tr,
                            mtry=min(mt,ncol(xtr)), ntree=nt)
          mean(abs(predict(m, xvl) - vl[[outcome]]), na.rm=TRUE)
        }, error=function(e) NA_real_)
        if (!is.na(s) && s < best$score) best <- list(score=s, w=w, pv=pv, mt=mt, nt=nt)
      }
    }
    if (best$score==Inf) next
    tv <- rbind(tr, vl)
    pp2 <- tryCatch(preProcess(tv[,best$pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
    if (is.null(pp2)) next
    x_tv <- tryCatch(predict(pp2, tv[,best$pv,drop=FALSE]), error=function(e) NULL)
    x_te <- tryCatch(predict(pp2, te[,best$pv,drop=FALSE]), error=function(e) NULL)
    if (is.null(x_tv) || is.null(x_te)) next
    d_tv <- cbind(x_tv, .y=tv[[outcome]]); names(d_tv)[ncol(d_tv)] <- outcome
    mf <- tryCatch(randomForest(as.formula(paste(outcome,"~ .")), data=d_tv,
                                mtry=min(best$mt,ncol(x_tv)), ntree=best$nt), error=function(e) NULL)
    if (is.null(mf)) next
    p <- tryCatch(predict(mf, x_te), error=function(e) NULL); if (is.null(p)) next
    preds[[i]] <- data.frame(true=te[[outcome]], pred=as.numeric(p))
  }
  if (length(preds)==0) return(empty_row(pt$uid[1], "no_successful_folds"))
  pdf_out <- bind_rows(preds)
  met <- reg_metrics(pdf_out$pred, pdf_out$true)
  ok_row(pt$uid[1], length(preds), nrow(pdf_out),
         cor_v=met$cor, r2_v=met$r2, mse_v=met$mse, mae_v=met$mae, rmse_v=met$rmse)
}
fit_continuous_baseline <- function(pt, cv_fn) {
  if (var(pt[[outcome]],na.rm=TRUE) %in% c(0,NA)) return(empty_row(pt$uid[1], "skipped"))
  sl <- cv_fn(pt); if (is.null(sl)) return(empty_row(pt$uid[1], "no_slices"))
  preds <- lapply(seq_along(sl$trains), function(i) {
    tm <- mean(pt[[outcome]][sl$trains[[i]]], na.rm=TRUE); if (is.na(tm)) return(NULL)
    data.frame(true=pt[[outcome]][sl$tests[[i]]], pred=rep(tm, length(sl$tests[[i]])))
  })
  preds <- Filter(Negate(is.null), preds)
  if (length(preds)==0) return(empty_row(pt$uid[1], "no_successful_folds"))
  pdf_out <- bind_rows(preds)
  met <- reg_metrics(pdf_out$pred, pdf_out$true)
  ok_row(pt$uid[1], length(preds), nrow(pdf_out),
         cor_v=met$cor, r2_v=met$r2, mse_v=met$mse, mae_v=met$mae, rmse_v=met$rmse)
}
# --- HURDLE ---
hurdle_one_fit <- function(tr_full, x_eval, ff) {
  cls <- tryCatch(suppressWarnings(train(
    as.formula(paste(".bin ~", paste(ff, collapse="+"))), data=tr_full, method="gbm",
    trControl=trainControl(method="none", classProbs=TRUE),
    tuneGrid=GBM_CLS, verbose=FALSE)), error=function(e) NULL)
  if (is.null(cls)) return(NULL)
  p_si <- tryCatch(predict(cls, x_eval, type="prob")[,"SI"], error=function(e) NULL)
  if (is.null(p_si)) return(NULL)
  sev <- rep(0, nrow(x_eval))
  tr_si <- tr_full[tr_full[[outcome]] > 0,]
  if (nrow(tr_si)>=10 && var(tr_si[[outcome]],na.rm=TRUE) > 0) {
    reg <- tryCatch(suppressWarnings(train(
      as.formula(paste(outcome,"~", paste(ff, collapse="+"))),
      data=tr_si[,c(outcome,ff)], method="gbm",
      trControl=trainControl(method="none"), tuneGrid=GBM_REG, verbose=FALSE)),
      error=function(e) NULL)
    if (!is.null(reg)) sev <- tryCatch(predict(reg, x_eval), error=function(e) rep(0, nrow(x_eval)))
  }
  list(pred=pmax(0, pmin(1, p_si * sev)), p_si=p_si)
}
fit_hurdle <- function(pt, fb, cv_fn) {
  feats <- candidate_predictors(pt)
  if (length(feats)<1 || var(pt[[outcome]],na.rm=TRUE) %in% c(0,NA))
    return(empty_row(pt$uid[1], "skipped"))
  n_si_pos <- sum(pt[[outcome]] > 0, na.rm = TRUE)
  n_si_neg <- sum(pt[[outcome]] == 0, na.rm = TRUE)
  if (n_si_pos < 10 || n_si_neg < 20) return(empty_row(pt$uid[1], "too_few_SI_counts"))
  sl <- cv_fn(pt); if (is.null(sl) || length(sl$trains)<1)
    return(empty_row(pt$uid[1], "no_slices"))
  preds <- list()
  for (i in seq_along(sl$trains)) {
    tr <- pt[sl$trains[[i]],]; vl <- pt[sl$vals[[i]],]; te <- pt[sl$tests[[i]],]
    scale_cols <- feats   # within-fold min-max
    tr_ref <- tr
    tr <- apply_minmax_fold(tr_ref, tr, scale_cols)
    vl <- apply_minmax_fold(tr_ref, vl, scale_cols)
    te <- apply_minmax_fold(tr_ref, te, scale_cols)
    tv <- rbind(tr, vl); best <- list(score=Inf)
    for (w in WINDOWS) {
      pv <- intersect(fb(w, data_model_base), feats); if (length(pv)<2) next
      pv <- drop_near_zero_variance_features(pv, tr, freqCut = 19, uniqueCut = 10); if (length(pv)<2) next
      pp <- tryCatch(preProcess(tr[,pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
      if (is.null(pp)) next
      xtr <- tryCatch(predict(pp, tr[,pv,drop=FALSE]), error=function(e) NULL)
      xvl <- tryCatch(predict(pp, vl[,pv,drop=FALSE]), error=function(e) NULL)
      if (is.null(xtr) || is.null(xvl)) next
      tr_full <- cbind(xtr, .y=tr[[outcome]]); names(tr_full)[ncol(tr_full)] <- outcome
      tr_full$.bin <- factor(ifelse(tr_full[[outcome]] > 0, "SI","noSI"), levels=c("noSI","SI"))
      if (length(unique(tr_full$.bin)) < 2) next
      ff <- setdiff(names(tr_full), c(outcome, ".bin"))
      h <- hurdle_one_fit(tr_full, xvl, ff); if (is.null(h)) next
      s <- mean(abs(h$pred - vl[[outcome]]), na.rm=TRUE)
      if (!is.na(s) && s < best$score) best <- list(score=s, w=w, pv=pv)
    }
    if (best$score==Inf) next
    pp2 <- tryCatch(preProcess(tv[,best$pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
    if (is.null(pp2)) next
    x_tv <- tryCatch(predict(pp2, tv[,best$pv,drop=FALSE]), error=function(e) NULL)
    x_te <- tryCatch(predict(pp2, te[,best$pv,drop=FALSE]), error=function(e) NULL)
    if (is.null(x_tv) || is.null(x_te)) next
    tv_full <- cbind(x_tv, .y=tv[[outcome]]); names(tv_full)[ncol(tv_full)] <- outcome
    tv_full$.bin <- factor(ifelse(tv_full[[outcome]] > 0, "SI","noSI"), levels=c("noSI","SI"))
    if (length(unique(tv_full$.bin)) < 2) next
    ff <- setdiff(names(tv_full), c(outcome, ".bin"))
    h <- hurdle_one_fit(tv_full, x_te, ff); if (is.null(h)) next
    preds[[i]] <- data.frame(true=te[[outcome]], pred=h$pred,
                             p_si=h$p_si, true_bin=as.integer(te[[outcome]]>0))
  }
  if (length(preds)==0) return(empty_row(pt$uid[1], "no_successful_folds"))
  pdf_out <- bind_rows(preds)
  auc_v <- tryCatch({
    if (length(unique(pdf_out$true_bin))<2) NA_real_
    else as.numeric(pROC::auc(pROC::roc(pdf_out$true_bin, pdf_out$p_si, quiet=TRUE)))
  }, error=function(e) NA_real_)
  met <- reg_metrics(pdf_out$pred, pdf_out$true)
  bin <- binary_metrics(pdf_out$true_bin, pdf_out$p_si, threshold = 0.5)
  ok_row(pt$uid[1], length(preds), nrow(pdf_out),
         cor_v=met$cor, r2_v=met$r2, mse_v=met$mse, mae_v=met$mae, rmse_v=met$rmse,
         auc_v=auc_v, sensitivity_v=bin$sensitivity, specificity_v=bin$specificity,
         balanced_acc_v=bin$balanced_acc,
         n_SI_v=sum(pdf_out$true_bin == 1, na.rm = TRUE),
         n_nonSI_v=sum(pdf_out$true_bin == 0, na.rm = TRUE))
}
fit_hurdle_baseline <- function(pt, cv_fn) {
  if (var(pt[[outcome]],na.rm=TRUE) %in% c(0,NA)) return(empty_row(pt$uid[1], "skipped"))
  sl <- cv_fn(pt); if (is.null(sl)) return(empty_row(pt$uid[1], "no_slices"))
  preds <- lapply(seq_along(sl$trains), function(i) {
    tr_y <- pt[[outcome]][sl$trains[[i]]]
    p_si <- mean(tr_y > 0, na.rm=TRUE)
    sev  <- mean(tr_y[tr_y > 0], na.rm=TRUE); if (is.na(sev)) sev <- 0
    final <- max(0, min(1, p_si * sev))
    te_y <- pt[[outcome]][sl$tests[[i]]]
    data.frame(true=te_y, pred=rep(final, length(te_y)),
               p_si=rep(p_si, length(te_y)), true_bin=as.integer(te_y>0))
  })
  preds <- Filter(Negate(is.null), preds)
  if (length(preds)==0) return(empty_row(pt$uid[1], "no_successful_folds"))
  pdf_out <- bind_rows(preds)
  auc_v <- tryCatch({
    if (length(unique(pdf_out$true_bin))<2) NA_real_
    else as.numeric(pROC::auc(pROC::roc(pdf_out$true_bin, pdf_out$p_si, quiet=TRUE)))
  }, error=function(e) NA_real_)
  met <- reg_metrics(pdf_out$pred, pdf_out$true)
  bin <- binary_metrics(pdf_out$true_bin, pdf_out$p_si, threshold = 0.5)
  ok_row(pt$uid[1], length(preds), nrow(pdf_out),
         cor_v=met$cor, r2_v=met$r2, mse_v=met$mse, mae_v=met$mae, rmse_v=met$rmse,
         auc_v=auc_v, sensitivity_v=bin$sensitivity, specificity_v=bin$specificity,
         balanced_acc_v=bin$balanced_acc,
         n_SI_v=sum(pdf_out$true_bin == 1, na.rm = TRUE),
         n_nonSI_v=sum(pdf_out$true_bin == 0, na.rm = TRUE))
}
# Generic parallel runner — parameterized by CV fn
run_one <- function(model_type, config_name, fb=NULL, cv_fn) {
  cat(sprintf("===== %s | %s =====\n", model_type, config_name)); t0 <- Sys.time()
  uids <- unique(data_model_base$uid)
  cl <- makeCluster(min(N_CORES, detectCores()-1, length(uids)))
  registerDoParallel(cl); on.exit(stopCluster(cl))
  out <- foreach(u=uids, .combine=bind_rows, .errorhandling="remove",
                 .packages=c("caret","randomForest","gbm","pROC","dplyr"),
                 .export=c("data_model_base","outcome","WINDOWS","MTRY_GRID","NTREE_GRID",
                           "GBM_CLS","GBM_REG","apply_minmax_fold",
                           "candidate_predictors",
                           "drop_near_zero_variance_features","is_passive_feature",
                           "reg_metrics","binary_metrics",
                           "fit_continuous","fit_continuous_baseline",
                           "fit_hurdle","fit_hurdle_baseline","hurdle_one_fit",
                           "empty_row","ok_row",
                           "ema_mean","ema_dyn","sensors","pmean","pdyn",
                           "get_passive_mean_cols","get_passive_dynamic_cols",
                           "configs","cv_fn")) %dopar% {
                             pt <- data_model_base[data_model_base$uid==u, , drop=FALSE]
                             tryCatch({
                               if (model_type == "Continuous") {
                                 if (config_name == "Baseline") fit_continuous_baseline(pt, cv_fn) else fit_continuous(pt, fb, cv_fn)
                               } else if (model_type == "Hurdle") {
                                 if (config_name == "Baseline") fit_hurdle_baseline(pt, cv_fn) else fit_hurdle(pt, fb, cv_fn)
                               }
                             }, error = function(e) empty_row(u, paste0("error: ", conditionMessage(e))))
                           }
  out$model_type <- model_type
  out$config     <- config_name
  out <- ensure_result_metric_cols(out)
  cat(sprintf("  %.1f min, n_ok=%d/%d\n",
              as.numeric(difftime(Sys.time(),t0,units="mins")),
              sum(out$status=="ok"), nrow(out)))
  out
}
# Run all configs for a given CV variant
run_variant <- function(variant_name, cv_params) {
  results_path <- file.path(model_config_cache_dir, sprintf("All_models_results_%s.csv", variant_name))
  results_stale <- file.exists(results_path) &&
    file.exists(data_model_path) &&
    file.info(data_model_path)$mtime > file.info(results_path)$mtime
  if (file.exists(results_path) && !RERUN_MODEL_RESULTS && !results_stale) {
    cat(sprintf("\n--- Variant '%s' cached, skipping ---\n", variant_name))
    return(ensure_result_metric_cols(read.csv(results_path, stringsAsFactors=FALSE)))
  }
  cat(sprintf("\n========== Variant: %s (train=%d, val=%d, test=%d) ==========\n",
              variant_name, cv_params$train_days, cv_params$val_days, cv_params$test_days))
  cv_fn <- make_sliding_cv(cv_params$train_days, cv_params$val_days,
                           cv_params$test_days, cv_params$ema_per_day)
  backfill_hurdle_class_counts <- function(out) {
    out <- ensure_result_metric_cols(out)
    if (!nrow(out) || !("model_type" %in% names(out)) || !any(out$model_type == "Hurdle")) return(out)
    out$Participant <- stringr::str_pad(as.character(out$Participant), width = 4, pad = "0")
    needs_counts <- out$model_type == "Hurdle" &
      (is.na(out$n_SI) | is.na(out$n_nonSI)) &
      out$status == "ok"
    if (!any(needs_counts)) return(out)
    for (rr in which(needs_counts)) {
      u <- stringr::str_pad(as.character(out$Participant[rr]), width = 4, pad = "0")
      pt <- data_model_base[data_model_base$uid == u, , drop = FALSE]
      if (!nrow(pt)) next
      sl <- cv_fn(pt)
      if (is.null(sl) || !length(sl$tests)) next
      y <- pt[[outcome]][unlist(sl$tests)]
      if (is.finite(out$n_obs[rr]) && out$n_obs[rr] > 0 && out$n_obs[rr] < length(y)) {
        y <- y[seq_len(out$n_obs[rr])]
      }
      out$n_SI[rr] <- sum(y > 0, na.rm = TRUE)
      out$n_nonSI[rr] <- sum(y == 0, na.rm = TRUE)
    }
    out
  }
  run_one_cached <- function(model_type, config_name, fb = NULL) {
    cache_file <- file.path(
      model_config_cache_dir,
      sprintf("model_results_%s_%s_%s.csv",
              variant_name,
              gsub("[^A-Za-z0-9]+", "_", model_type),
              gsub("[^A-Za-z0-9]+", "_", config_name))
    )
    cache_stale <- file.exists(cache_file) &&
      file.exists(data_model_path) &&
      file.info(data_model_path)$mtime > file.info(cache_file)$mtime
    if (file.exists(cache_file) && !RERUN_COMPLETED_CONFIGS && !cache_stale) {
      cat(sprintf("===== %s | %s =====\n", model_type, config_name))
      cat(sprintf("  cached, skipping: %s\n", cache_file))
      out <- ensure_result_metric_cols(read.csv(cache_file, stringsAsFactors = FALSE))
      out$model_type <- model_type
      out$config <- config_name
      out <- backfill_hurdle_class_counts(out)
      write.csv(out, cache_file, row.names = FALSE)
      return(out)
    }
    out <- run_one(model_type, config_name, fb = fb, cv_fn = cv_fn)
    write.csv(out, cache_file, row.names = FALSE)
    cat(sprintf("  cached result saved: %s\n", cache_file))
    ensure_result_metric_cols(out)
  }
  all_res <- list()
  # Run one baseline benchmark plus the six planned feature-set models.
  # Baseline is retained for comparison statistics but is not counted as a feature-set model.
  if (RUN_CONTINUOUS_MODELS) {
    all_res[["Cont_Baseline"]] <- run_one_cached("Continuous", "Baseline")
    for (cn in MODEL_CONFIGS)
      all_res[[paste0("Cont_", cn)]] <- run_one_cached("Continuous", cn, fb=configs[[cn]])
  }
  if (RUN_HURDLE_MODELS) {
    all_res[["Hurdle_Baseline"]] <- run_one_cached("Hurdle", "Baseline")
    for (cn in MODEL_CONFIGS)
      all_res[[paste0("Hurdle_", cn)]] <- run_one_cached("Hurdle", cn, fb=configs[[cn]])
  }
  combined <- bind_rows(all_res) %>%
    select(model_type, config, Participant,
           cor, r2, mse, mae, rmse, auc, sensitivity, specificity, balanced_acc,
           n_obs, n_SI, n_nonSI, n_folds, status) %>%
    mutate(variant=variant_name)
  combined <- backfill_hurdle_class_counts(combined)
  write.csv(combined, results_path, row.names=FALSE)
  cat(sprintf("✅ Cached combined model results: %s (%d rows)\n", results_path, nrow(combined)))
  ensure_result_metric_cols(combined)
}
results_main  <- run_variant("main",  CV_VARIANTS$main)
print_summary <- function(combined, label) {
  cat(sprintf("\n===== Summary [%s] (n_folds >= %d, status==ok) =====\n", label, MIN_FOLDS))
  s <- combined %>% filter(status=="ok", n_folds>=MIN_FOLDS) %>%
    group_by(model_type, config) %>%
    summarise(n=n(),
              median_cor =round(median(cor,  na.rm=TRUE),3),
              mean_r2    =round(mean(r2,    na.rm=TRUE),3),
              median_r2  =round(median(r2,   na.rm=TRUE),3),
              median_mse =round(median(mse,  na.rm=TRUE),3),
              median_mae =round(median(mae,  na.rm=TRUE),3),
              median_rmse=round(median(rmse, na.rm=TRUE),3),
              median_auc =round(median(auc,  na.rm=TRUE),3),
              median_sensitivity =round(median(sensitivity,  na.rm=TRUE),3),
              median_specificity =round(median(specificity,  na.rm=TRUE),3),
              median_balanced_acc=round(median(balanced_acc, na.rm=TRUE),3),
              .groups="drop") %>%
    arrange(model_type, desc(coalesce(median_auc, median_cor)))
  print(s, n=Inf); s
}
summary_main  <- print_summary(results_main,  "main")

diagnose_continuous_model_failures <- function(combined, variant_name = "main",
                                               min_folds = MIN_FOLDS) {
  analytic_uids <- data_model_base %>%
    distinct(uid) %>%
    mutate(uid = as.character(uid)) %>%
    pull(uid)
  cont <- combined %>%
    filter(variant == variant_name,
           model_type == "Continuous",
           config %in% MODEL_CONFIGS) %>%
    mutate(
      evaluable = status == "ok" & n_folds >= min_folds,
      Participant = as.character(Participant)
    )
  outcome_diag <- data_model_base %>%
    group_by(uid) %>%
    summarise(
      model_rows = n(),
      observed_q9 = sum(!is.na(q9)),
      q9_positive = sum(q9 > 0, na.rm = TRUE),
      q9_zero = sum(q9 == 0, na.rm = TRUE),
      q9_variance = stats::var(q9, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(uid = as.character(uid))
  predictor_diag <- bind_rows(lapply(analytic_uids, function(u) {
    pt <- data_model_base[data_model_base$uid == u, , drop = FALSE]
    data.frame(
      uid = u,
      candidate_predictors = length(candidate_predictors(pt)),
      ema_mean_predictors = length(intersect(ema_mean, names(pt))),
      ema_dynamic_predictors = length(intersect(ema_dyn, names(pt))),
      passive_mean_predictors = length(intersect(unique(unlist(lapply(WINDOWS, pmean, d = data_model_base))), names(pt))),
      passive_dynamic_predictors = length(intersect(unique(unlist(lapply(WINDOWS, pdyn, d = data_model_base))), names(pt))),
      stringsAsFactors = FALSE
    )
  }))
  by_config <- cont %>%
    left_join(outcome_diag, by = c("Participant" = "uid")) %>%
    left_join(predictor_diag, by = c("Participant" = "uid")) %>%
    arrange(config, desc(evaluable), status, Participant)
  by_participant <- by_config %>%
    group_by(Participant) %>%
    summarise(
      n_configs_ok = sum(evaluable, na.rm = TRUE),
      configs_failed = paste(config[!evaluable], collapse = ", "),
      statuses = paste(unique(status[!evaluable]), collapse = "; "),
      min_folds_observed = min(n_folds, na.rm = TRUE),
      max_folds_observed = max(n_folds, na.rm = TRUE),
      model_rows = dplyr::first(model_rows),
      observed_q9 = dplyr::first(observed_q9),
      q9_positive = dplyr::first(q9_positive),
      q9_zero = dplyr::first(q9_zero),
      q9_variance = dplyr::first(q9_variance),
      candidate_predictors = dplyr::first(candidate_predictors),
      .groups = "drop"
    ) %>%
    mutate(all_six_configs_ok = n_configs_ok == length(MODEL_CONFIGS)) %>%
    arrange(n_configs_ok, Participant)
  out_config <- file.path(report_dir, sprintf("continuous_model_failure_diagnostic_by_config_%s.csv", variant_name))
  out_pt <- file.path(report_dir, sprintf("continuous_model_failure_diagnostic_by_participant_%s.csv", variant_name))
  write.csv(by_config, out_config, row.names = FALSE)
  write.csv(by_participant, out_pt, row.names = FALSE)
  cat(sprintf("Saved continuous model diagnostics: %s\n", out_config))
  cat(sprintf("Saved continuous model diagnostics: %s\n", out_pt))
  cat("\nContinuous model evaluability by config:\n")
  print(cont %>%
          group_by(config, status) %>%
          summarise(n_pt = n_distinct(Participant),
                    median_folds = median(n_folds, na.rm = TRUE),
                    .groups = "drop") %>%
          arrange(config, status), n = Inf)
  invisible(list(by_config = by_config, by_participant = by_participant))
}

if (isTRUE(RUN_MODEL_DIAGNOSTICS)) {
  continuous_failure_diagnostics_main <- diagnose_continuous_model_failures(results_main, "main", MIN_FOLDS)
} else {
  continuous_failure_diagnostics_main <- NULL
  cat("Skipping continuous model diagnostics (RUN_MODEL_DIAGNOSTICS=FALSE)\n")
}

compare_models_to_baseline <- function(combined, variant_name = "main",
                                       min_folds = MIN_FOLDS) {
  metrics <- tibble::tribble(
    ~metric, ~higher_is_better,
    "r2",    TRUE,
    "cor",   TRUE,
    "mae",   FALSE,
    "rmse",  FALSE,
    "mse",   FALSE,
    "auc",   TRUE,
    "sensitivity", TRUE,
    "specificity", TRUE,
    "balanced_acc", TRUE
  )
  valid <- combined %>%
    filter(status == "ok", n_folds >= min_folds, variant == variant_name)
  out <- list()
  k <- 1L
  for (mt in unique(valid$model_type)) {
    base <- valid %>%
      filter(model_type == mt, config == "Baseline") %>%
      select(Participant, all_of(metrics$metric))
    if (!nrow(base)) next
    for (cfg in setdiff(unique(valid$config[valid$model_type == mt]), "Baseline")) {
      cfg_dat <- valid %>%
        filter(model_type == mt, config == cfg) %>%
        select(Participant, all_of(metrics$metric))
      paired <- inner_join(cfg_dat, base, by = "Participant", suffix = c("_model", "_baseline"))
      for (ii in seq_len(nrow(metrics))) {
        mm <- metrics$metric[ii]
        hi <- metrics$higher_is_better[ii]
        model_col <- paste0(mm, "_model")
        base_col  <- paste0(mm, "_baseline")
        if (!model_col %in% names(paired) || !base_col %in% names(paired)) next
        d <- paired %>%
          transmute(
            diff_raw = .data[[model_col]] - .data[[base_col]],
            diff_better = if (hi) diff_raw else -diff_raw
          ) %>%
          filter(is.finite(diff_raw), is.finite(diff_better))
        if (nrow(d) < 5 || length(unique(d$diff_better)) < 2) next
        wt <- tryCatch(
          suppressWarnings(stats::wilcox.test(d$diff_better, mu = 0, alternative = "greater",
                                              exact = FALSE, paired = FALSE)),
          error = function(e) NULL
        )
        out[[k]] <- data.frame(
          variant = variant_name,
          model_type = mt,
          config = cfg,
          metric = mm,
          higher_is_better = hi,
          n_paired = nrow(d),
          median_model_minus_baseline = median(d$diff_raw, na.rm = TRUE),
          median_improvement = median(d$diff_better, na.rm = TRUE),
          mean_improvement = mean(d$diff_better, na.rm = TRUE),
          pct_better_than_baseline = mean(d$diff_better > 0, na.rm = TRUE) * 100,
          p_value = if (is.null(wt)) NA_real_ else wt$p.value,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  res <- bind_rows(out)
  if (nrow(res)) {
    res <- res %>%
      group_by(model_type, metric) %>%
      mutate(p_fdr = p.adjust(p_value, method = "BH"),
             significant_fdr_05 = !is.na(p_fdr) & p_fdr < 0.05) %>%
      ungroup() %>%
      arrange(model_type, metric, p_fdr, p_value)
    cat("\n===== Models Significantly Better Than Baseline (FDR < .05) =====\n")
    print(res %>% filter(significant_fdr_05) %>%
            arrange(model_type, metric, p_fdr), n = Inf)
  } else {
    cat("No valid paired baseline comparisons available.\n")
  }
  invisible(res)
}

baseline_tests_main <- compare_models_to_baseline(results_main, "main", MIN_FOLDS)

make_predictive_performance_tables <- function(combined, baseline_tests = NULL,
                                               variant_name = "main",
                                               min_folds = MIN_FOLDS) {
  combined <- ensure_result_metric_cols(combined)
  config_label_map <- c(
    Baseline = "Baseline",
    EMA_mean = "EMA",
    EMA_mean_dynamic = "EMA (dynamic)",
    Opt_Passive_mean = "Passive",
    Opt_Passive_mean_dynamic = "Passive (dynamic)",
    EMA_Opt_Passive_mean = "EMA + Passive",
    EMA_Opt_Passive_mean_dynamic = "EMA + Passive (dynamic)"
  )
  config_order <- c("Baseline", MODEL_CONFIGS)
  fmt_num <- function(x, digits = 3, plus = FALSE) {
    if (is.na(x) || !is.finite(x)) return("")
    out <- sprintf(paste0("%.", digits, "f"), x)
    if (plus && x > 0) out <- paste0("+", out)
    out
  }
  fmt_iqr <- function(x, digits = 3, plus = FALSE) {
    x <- x[is.finite(x)]
    if (!length(x)) return("")
    qs <- stats::quantile(x, probs = c(.25, .5, .75), na.rm = TRUE)
    sprintf("%s [%s, %s]",
            fmt_num(qs[2], digits, plus),
            fmt_num(qs[1], digits, plus),
            fmt_num(qs[3], digits, plus))
  }
  pct_gt <- function(x, threshold) {
    x <- x[is.finite(x)]
    if (!length(x)) return("")
    sprintf("%.0f%%", mean(x > threshold, na.rm = TRUE) * 100)
  }
  summarize_config <- function(model_type_i, config_i) {
    combined %>%
      filter(status == "ok", n_folds >= min_folds, variant == variant_name,
             model_type == model_type_i, config == config_i)
  }
  cont_raw <- bind_rows(lapply(config_order, function(cfg) {
    d <- summarize_config("Continuous", cfg)
    data.frame(
      `Model type` = "Continuous",
      Configuration = unname(config_label_map[cfg]),
      n = nrow(d),
      cor_median = median(d$cor, na.rm = TRUE),
      r2_median = median(d$r2, na.rm = TRUE),
      mae_median = median(d$mae, na.rm = TRUE),
      rmse_median = median(d$rmse, na.rm = TRUE),
      `% r>0` = pct_gt(d$cor, 0),
      `r, median [IQR]` = fmt_iqr(d$cor, 3, plus = TRUE),
      `R², median [IQR]` = fmt_iqr(d$r2, 3, plus = TRUE),
      `MAE, median [IQR]` = fmt_iqr(d$mae, 4),
      `RMSE, median [IQR]` = fmt_iqr(d$rmse, 4),
      check.names = FALSE
    )
  }))
  hurdle_raw <- bind_rows(lapply(config_order, function(cfg) {
    d <- summarize_config("Hurdle", cfg)
    data.frame(
      `Model type` = "Hurdle",
      Configuration = unname(config_label_map[cfg]),
      n = nrow(d),
      auc_median = median(d$auc, na.rm = TRUE),
      sensitivity_median = median(d$sensitivity, na.rm = TRUE),
      specificity_median = median(d$specificity, na.rm = TRUE),
      balanced_acc_median = median(d$balanced_acc, na.rm = TRUE),
      mae_median = median(d$mae, na.rm = TRUE),
      `AUC, median [IQR]` = fmt_iqr(d$auc, 2),
      `% AUC>0.5` = pct_gt(d$auc, 0.5),
      `Sensitivity, median [IQR]` = fmt_iqr(d$sensitivity, 2),
      `Specificity, median [IQR]` = fmt_iqr(d$specificity, 2),
      `Balanced acc., median [IQR]` = fmt_iqr(d$balanced_acc, 2),
      `MAE, median [IQR]` = fmt_iqr(d$mae, 4),
      check.names = FALSE
    )
  }))
  cont_table <- cont_raw %>%
    select(`Model type`, Configuration, n, `r, median [IQR]`, `% r>0`,
           `R², median [IQR]`, `MAE, median [IQR]`, `RMSE, median [IQR]`)
  hurdle_table <- hurdle_raw %>%
    select(`Model type`, Configuration, n, `AUC, median [IQR]`, `% AUC>0.5`,
           `Sensitivity, median [IQR]`, `Specificity, median [IQR]`,
           `Balanced acc., median [IQR]`, `MAE, median [IQR]`)
  perf_table <- bind_rows(cont_table, hurdle_table)

  cont_report_table <- cont_table %>% select(-`Model type`)
  hurdle_report_table <- hurdle_table %>% select(-`Model type`)
  cont_csv <- file.path(report_dir, sprintf("Table2_continuous_model_results_%s.csv", variant_name))
  hurdle_csv <- file.path(report_dir, sprintf("Table3_hurdle_model_results_%s.csv", variant_name))
  write.csv(cont_report_table, cont_csv, row.names = FALSE)
  write.csv(hurdle_report_table, hurdle_csv, row.names = FALSE)
  cat(sprintf("Saved report continuous table CSV: %s\n", cont_csv))
  cat(sprintf("Saved report hurdle table CSV: %s\n", hurdle_csv))
  invisible(perf_table)
}

table2_predictive_performance_main <- make_predictive_performance_tables(
  results_main, baseline_tests_main, "main", MIN_FOLDS
)
# Auto-detect best configs from MAIN variant
# Best continuous config = highest MEAN R² (per request)
auto_best_cont <- summary_main %>% filter(model_type=="Continuous", config != "Baseline") %>%
  slice_max(mean_r2, n=1, with_ties=FALSE) %>% pull(config)
if (length(auto_best_cont)>0)   BEST_CONT_CFG   <- auto_best_cont
if (RUN_HURDLE_MODELS) {
  auto_best_hurdle <- summary_main %>% filter(model_type=="Hurdle", config != "Baseline") %>%
    slice_max(median_auc, n=1, with_ties=FALSE) %>% pull(config)
  if (length(auto_best_hurdle)>0) BEST_HURDLE_CFG <- auto_best_hurdle
} else {
  BEST_HURDLE_CFG <- NA_character_
}
cat(sprintf("\nBest Continuous: %s | Best Hurdle: %s\n",
            ifelse(exists("BEST_CONT_CFG"), BEST_CONT_CFG, "not run"),
            ifelse(!is.na(BEST_HURDLE_CFG), BEST_HURDLE_CFG, "not run")))
# ==============================================================================
# SECTION 4: REPORTING & VISUALIZATION
# ==============================================================================
cat("\n===== SECTION 5: Reporting & Visualization =====\n")
config_levels <- c("Baseline","EMA_mean","EMA_mean_dynamic",
                   "Opt_Passive_mean","Opt_Passive_mean_dynamic",
                   "EMA_Opt_Passive_mean","EMA_Opt_Passive_mean_dynamic")
config_labels <- c("Baseline","EMA mean","EMA mean + dyn",
                   "Passive mean","Passive mean + dyn",
                   "EMA + Passive mean","EMA + Passive\nmean + dyn")
config_family <- c("Baseline","EMA","EMA","Passive","Passive","Combined","Combined")
family_palette <- c(Baseline="#bfbfbf", EMA="#8aae8b", Passive="#e9c578", Combined="#5c7a8a")
sensor_labels <- c(acc="Accelerometer", hr="Heart rate", steps="Steps",
                   convo="Conversation", gps="GPS / location", gps_home="Time at home",
                   sms="SMS", calls="Calls", battery="Battery", sleep="Sleep")
ema_item_labels <- c(q1="q1: anhedonia", q2="q2: depression", q3="q3: sleep",
                     q4="q4: fatigue", q5="q5: appetite", q6="q6: worthlessness",
                     q7="q7: concentration", q8="q8: psychomotor", q9="q9: SI severity",
                     g1="g1: anxiety", g2="g2: worry")
window_labels   <- c(hr1="1 hr", hr4="4 hr", day1="1 day", day7="7 day")
window_palette  <- c("1 hr"="#e6b800", "4 hr"="#e07a3c", "1 day"="#7fb5d8", "7 day"="#2d5e7e")
ema_time_palette<- c("prior prompt"="#9b6fa3","prior 1d"="#7fb5d8","prior 7d"="#2d5e7e")
# --- F1: model performance ---
make_F1 <- function(combined, variant_name) {
  perf_df <- combined %>%
    filter(status=="ok", n_folds>=MIN_FOLDS, model_type %in% c("Continuous","Hurdle")) %>%
    mutate(config=factor(config, levels=config_levels, labels=config_labels),
           family=factor(config_family[match(config, config_labels)],
                         levels=c("Baseline","EMA","Passive","Combined")),
           r2_disp=r2)   # true coefficient of determination (can be negative)
  
  get_winner <- function(df, metric, higher=TRUE) {
    base_med <- df %>% filter(config=="Baseline") %>%
      summarise(b=median(.data[[metric]], na.rm=TRUE)) %>% pull(b)
    s <- df %>% filter(config != "Baseline") %>% group_by(config) %>%
      summarise(med=median(.data[[metric]], na.rm=TRUE), .groups="drop") %>%
      mutate(delta=med - base_med)
    if (higher) s <- arrange(s, desc(delta)) else s <- arrange(s, delta)
    list(winner=as.character(s$config[1]), delta=s$delta[1])
  }
  
  plot_panel <- function(df, metric, ylab, mt, higher=TRUE, ref_line=NULL) {
    d <- df %>% filter(model_type==mt, !is.na(.data[[metric]]))
    if (nrow(d)==0) return(patchwork::plot_spacer())
    d <- d %>% mutate(config=droplevels(config))   # drop empty configs so the winner box aligns
    w <- get_winner(d, metric, higher)
    meds <- d %>% group_by(config) %>%
      summarise(m=median(.data[[metric]], na.rm=TRUE), .groups="drop")
    p <- ggplot(d, aes(x=config, y=.data[[metric]], fill=family))
    if (!is.null(ref_line))
      p <- p + geom_hline(yintercept=ref_line, linetype="dashed", color="grey50")
    p <- p +
      geom_jitter(aes(color=family), width=0.18, alpha=0.25, size=0.6) +
      geom_boxplot(outlier.shape=NA, alpha=0.85, width=0.65, color="grey25") +
      geom_point(data=meds, aes(x=config, y=m), inherit.aes=FALSE,
                 shape=18, size=2.6, color="black") +
      geom_text(data=meds, aes(x=config, y=m, label=sprintf("%.3f", m)),
                inherit.aes=FALSE, vjust=-1.1, size=2.6,
                color=ifelse(meds$config==w$winner, "red", "black")) +
      scale_fill_manual(values=family_palette, guide="none") +
      scale_color_manual(values=family_palette, guide="none") +
      labs(title=sprintf("%s — %s", mt, ylab), y=ylab, x=NULL) +
      theme_minimal(base_size=10) +
      theme(axis.text.x=element_text(angle=30, hjust=1, size=8),
            plot.title=element_text(face="bold", size=11),
            panel.grid.minor=element_blank())
    winner_idx <- which(levels(d$config) == w$winner)
    if (length(winner_idx))
      p <- p + geom_rect(xmin=winner_idx-0.4, xmax=winner_idx+0.4,
                         ymin=-Inf, ymax=Inf, fill=NA, color="red",
                         linewidth=0.6, inherit.aes=FALSE)
    p + annotate("label", x=Inf, y=-Inf, hjust=1.05, vjust=-0.3,
                 label=sprintf("winner: %s\nΔ vs baseline: %+.3f", w$winner, w$delta),
                 size=2.5, color="red", fill="white",
                 label.padding=unit(0.15,"lines"))
  }
  
  fig1 <- (plot_panel(perf_df, "r2_disp",   "R²",         "Continuous", TRUE) |
             plot_panel(perf_df, "cor",       "Pearson r",  "Continuous", TRUE) |
             plot_panel(perf_df, "mae",       "MAE",        "Continuous", FALSE)) /
    (plot_panel(perf_df, "auc",       "AUC",        "Hurdle",     TRUE, ref_line=0.5) |
       plot_panel(perf_df, "cor",       "Pearson r",  "Hurdle",     TRUE) |
       plot_panel(perf_df, "mae",       "MAE",        "Hurdle",     FALSE)) +
    plot_annotation(
      title=sprintf("Model performance across configurations [%s variant]", variant_name),
      subtitle=sprintf("CV: train=%d days, val=%d days, test=%d day | Boxes = PT distribution; diamond = median",
                       CV_VARIANTS[[variant_name]]$train_days,
                       CV_VARIANTS[[variant_name]]$val_days,
                       CV_VARIANTS[[variant_name]]$test_days),
      theme=theme(plot.title=element_text(face="bold", size=14)))
  ggsave(file.path(plots_dir, sprintf("F1_model_performance_%s.png", variant_name)),
         fig1, width=14, height=9, dpi=200, bg="white")
  cat(sprintf("Saved F1_model_performance_%s.png\n", variant_name))
}
if (isTRUE(RUN_EXTRA_DIAGNOSTIC_FIGURES)) make_F1(results_main, "main")
# --- F2: within-person performance for best continuous (MAIN variant only) ---
make_F2 <- function(combined, variant_name) {
  best_df <- combined %>%
    filter(status=="ok", n_folds>=MIN_FOLDS,
           model_type=="Continuous", config==BEST_CONT_CFG)
  if (nrow(best_df)==0) { cat("⚠ No rows for best continuous — skipping F2\n"); return() }
  pct_pos <- mean(best_df$r2 > 0, na.rm=TRUE) * 100
  dot_df <- best_df %>% arrange(r2) %>%
    mutate(rank=row_number()) %>%
    select(rank, r2, cor, mae) %>%
    pivot_longer(c(r2, cor, mae), names_to="metric", values_to="value") %>%
    mutate(metric=recode(metric, r2="R²", cor="Pearson r", mae="MAE (norm)"),
           metric=factor(metric, levels=c("R²","Pearson r","MAE (norm)")))
  metric_palette <- c("R²"="grey20","Pearson r"="#7fb5a3","MAE (norm)"="#e09a3c")
  
  p_dot <- ggplot(dot_df, aes(x=value, y=rank, color=metric)) +
    geom_vline(xintercept=0, color="grey50") +
    geom_point(alpha=0.75, size=1.2) +
    scale_color_manual(values=metric_palette) +
    annotate("text", x=-0.45, y=max(dot_df$rank)*0.05,
             label=sprintf("%.0f%%", pct_pos), size=4, fontface="bold") +
    labs(x="Value", y="Participants", color=NULL,
         title=sprintf("Within-person performance — %s [%s variant]", BEST_CONT_CFG, variant_name),
         subtitle=sprintf("Continuous regression | %.0f%% R²>0 | Mean R²=%.3f | Median MAE=%.3f | Test obs: %s",
                          pct_pos, mean(best_df$r2, na.rm=TRUE),
                          median(best_df$mae, na.rm=TRUE),
                          format(sum(best_df$n_obs), big.mark=","))) +
    theme_minimal(base_size=10) +
    theme(legend.position="top", plot.title=element_text(face="bold"),
          panel.grid.major.y=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  obs_df <- best_df %>% arrange(r2) %>% mutate(rank=row_number())
  p_obs <- ggplot(obs_df, aes(x=n_obs, y=rank)) +
    geom_col(fill="grey60", width=0.8) +
    labs(x="Test observations", y=NULL) +
    theme_minimal(base_size=10) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.grid.major.y=element_blank())
  p_ridge <- ggplot(dot_df, aes(x=value, y=metric, fill=metric)) +
    geom_density_ridges(alpha=0.85, scale=1.3, rel_min_height=0.005) +
    geom_jitter(height=0.05, alpha=0.4, size=0.4, color="black") +
    geom_vline(xintercept=0, color="grey50") +
    scale_fill_manual(values=metric_palette, guide="none") +
    labs(x="Value", y=NULL,
         caption=sprintf("N = %d participants", nrow(best_df))) +
    theme_minimal(base_size=10) +
    theme(plot.caption=element_text(hjust=1, size=8, color="grey40"))
  fig2 <- p_dot + p_obs + p_ridge + plot_layout(widths=c(2.3, 1, 1.5))
  ggsave(file.path(plots_dir, sprintf("F2_withinperson_continuous_%s.png", variant_name)),
         fig2, width=14, height=9, dpi=200, bg="white")
  cat(sprintf("Saved F2_withinperson_continuous_%s.png\n", variant_name))
}
if (isTRUE(RUN_EXTRA_DIAGNOSTIC_FIGURES)) make_F2(results_main, "main")
# ==============================================================================
# SECTION 6: DEMOGRAPHICS TABLE 1
# ==============================================================================
cat("\n===== SECTION 6: Demographics Table 1 =====\n")
demo_path <- file.path(RAW_DATA_DIR, DEMO_FILE)
if (!file.exists(demo_path)) {
  cat(sprintf("⚠ Demographics file not found: %s — skipping Table 1\n", demo_path))
} else {
  demo_raw <- read.csv(demo_path, stringsAsFactors=FALSE)
  cat(sprintf("Loaded demographics: %d rows × %d cols\n", nrow(demo_raw), ncol(demo_raw)))
  cat("Available columns:\n"); print(names(demo_raw))
  
  finalized_uids <- data_model_base %>%
    distinct(uid) %>%
    pull(uid) %>%
    as.character() %>%
    str_pad(width=4, pad="0")
  cat(sprintf("Analytic feature-matrix sample for demographics: %d PTs\n", length(finalized_uids)))

  normalize_demo_uid <- function(x) {
    out <- stringr::str_extract(as.character(x), "[0-9]+")
    stringr::str_pad(out, width = 4, pad = "0")
  }
  uid_col_candidates <- intersect(
    c("uid","UID","mlife_id","Participant","participant","id","ID","subject"),
    names(demo_raw)
  )
  if (!length(uid_col_candidates)) uid_col_candidates <- names(demo_raw)
  uid_match_counts <- vapply(uid_col_candidates, function(cc) {
    sum(unique(normalize_demo_uid(demo_raw[[cc]])) %in% finalized_uids, na.rm = TRUE)
  }, integer(1))
  uid_col <- names(sort(uid_match_counts, decreasing = TRUE))[1]
  cat(sprintf("Using '%s' as participant ID column (%d UID matches).\n",
              uid_col, uid_match_counts[uid_col]))

  demo_raw$uid <- normalize_demo_uid(demo_raw[[uid_col]])

  # Descriptive Table 1 over the analytic feature-matrix sample. Deduplicate first:
  # one participant should contribute one row to demographics.
  demo_sample <- demo_raw %>%
    filter(uid %in% finalized_uids) %>%
    arrange(uid) %>%
    distinct(uid, .keep_all = TRUE)
  demo_n <- dplyr::n_distinct(demo_sample$uid)
  cat(sprintf("Demographic participants matched: %d / %d finalized PTs\n",
              demo_n, length(finalized_uids)))
  missing_demo_uids <- setdiff(finalized_uids, unique(demo_sample$uid))
  if (length(missing_demo_uids)) {
    cat(sprintf("Missing demographics for %d finalized PTs: %s\n",
                length(missing_demo_uids), paste(missing_demo_uids, collapse=", ")))
  }
  
  # Auto-detect demographic columns by name pattern
  detect_col <- function(df, patterns) {
    cols <- names(df)
    for (p in patterns) {
      hit <- cols[grepl(p, cols, ignore.case=TRUE)]
      if (length(hit)>0) return(hit[1])
    }
    return(NA_character_)
  }
  col_age    <- detect_col(demo_sample, c("^age$", "age_at"))
  col_sex    <- detect_col(demo_sample, c("^sex$", "^gender$", "gender_identity"))
  col_race   <- detect_col(demo_sample, c("^race$", "race_eth", "ethnicity"))
  col_edu    <- detect_col(demo_sample, c("^education", "edu_level", "highest_education"))
  col_income <- detect_col(demo_sample, c("^income", "household_income"))
  
  detected <- c(age=col_age, sex=col_sex, race=col_race, education=col_edu, income=col_income)
  detected <- detected[!is.na(detected)]
  cat("Detected demographic columns:\n"); print(detected)
  
  # Build publication-style Table 1 summary:
  # one characteristic column, one overall column, category levels indented.
  summarize_var <- function(df, var, label) {
    if (is.na(var) || !var %in% names(df)) return(NULL)
    df <- df %>% distinct(uid, .keep_all = TRUE)
    x <- df[[var]]
    if (is.numeric(x)) {
      v <- x[!is.na(x)]
      ov <- if (length(v)==0) "" else sprintf("%.1f (%.1f)", mean(v), stats::sd(v))
      data.frame(Characteristic=sprintf("%s, mean (SD)", label),
                 Overall=ov, check.names=FALSE)
    } else {
      x <- as.character(x); x[is.na(x) | x==""] <- "Unknown"
      lvls <- sort(unique(x))
      tab <- table(factor(x, levels=lvls))
      pct <- round(prop.table(tab)*100, 1)
      bind_rows(
        data.frame(Characteristic=label, Overall="", check.names=FALSE),
        data.frame(
          Characteristic=paste0("  ", lvls),
          Overall=sprintf("%d (%.1f%%)", as.integer(tab), pct),
          check.names=FALSE)
      )
    }
  }
  
  table1_parts <- list()
  # N row
  table1_parts[["N"]] <- data.frame(
    Characteristic="N",
    Overall=as.character(demo_n),
    check.names=FALSE)
  for (nm in names(detected))
    table1_parts[[nm]] <- summarize_var(demo_sample, detected[[nm]],
                                        tools::toTitleCase(nm))
  table1 <- bind_rows(table1_parts)
  cat("\n── Table 1: Demographics ──\n"); print(table1)
  
  # Save one demographic table CSV only.
  table1_csv  <- file.path(report_dir, "Table1_demographics.csv")
  write.csv(table1, table1_csv, row.names=FALSE)
  cat(sprintf("Saved: %s\n", table1_csv))
}
cat("\n✅ Pipeline complete.\n")
cat(sprintf("   • Final report folder: %s\n", report_dir))
################################################################################
# EXTRA FIGURES — within-person performance and time-series figures
#
# Source AFTER your Section 3-6 has run (needs these objects/functions in scope):
#   data_model_base, outcome, WINDOWS, configs, hurdle_one_fit, make_sliding_cv,
#   BEST_HURDLE_CFG, BEST_CONT_CFG, MIN_FOLDS, results_main, CV_VARIANTS,
#   plots_dir, results_dir, BASE_PATH, DEMO_FILE, N_CORES,
#   and the naming patch applied (so passive features resolve).
# Packages: caret, gbm, pROC, foreach, doParallel, dplyr, tidyr, ggplot2,
#           ggridges, patchwork, stringr.
#
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggridges)
  library(patchwork); library(foreach); library(doParallel); library(stringr)
})
# ==============================================================================
# A. WITHIN-PERSON HURDLE FIGURE  (AUC / Pearson r / MAE + test observations)
#    Refits ONLY the best hurdle config per participant to recover the per-fold
#    held-out predictions and SI probabilities.
# ==============================================================================
hurdle_pt_clsmetrics <- function(pt, fb, cv_fn, thr = 0.5) {
  feats <- candidate_predictors(pt)
  if (length(feats) < 1 || var(pt[[outcome]], na.rm = TRUE) %in% c(0, NA)) return(NULL)
  n_si_pos <- sum(pt[[outcome]] > 0, na.rm = TRUE)
  n_si_neg <- sum(pt[[outcome]] == 0, na.rm = TRUE)
  if (n_si_pos < 10 || n_si_neg < 20) return(NULL)
  sl <- cv_fn(pt); if (is.null(sl) || length(sl$trains) < 1) return(NULL)
  preds <- list()
  for (i in seq_along(sl$trains)) {
    tr <- pt[sl$trains[[i]],]; vl <- pt[sl$vals[[i]],]; te <- pt[sl$tests[[i]],]
    scale_cols <- feats   # within-fold min-max
    tr_ref <- tr
    tr <- apply_minmax_fold(tr_ref, tr, scale_cols)
    vl <- apply_minmax_fold(tr_ref, vl, scale_cols)
    te <- apply_minmax_fold(tr_ref, te, scale_cols)
    tv <- rbind(tr, vl); best <- list(score = Inf)
    for (w in WINDOWS) {
      pv <- intersect(fb(w, data_model_base), feats); if (length(pv) < 2) next
      pv <- drop_near_zero_variance_features(pv, tr, freqCut = 19, uniqueCut = 10); if (length(pv) < 2) next
      pp <- tryCatch(caret::preProcess(tr[,pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
      if (is.null(pp)) next
      xtr <- tryCatch(predict(pp, tr[,pv,drop=FALSE]), error=function(e) NULL)
      xvl <- tryCatch(predict(pp, vl[,pv,drop=FALSE]), error=function(e) NULL)
      if (is.null(xtr) || is.null(xvl)) next
      trf <- cbind(xtr, .y=tr[[outcome]]); names(trf)[ncol(trf)] <- outcome
      trf$.bin <- factor(ifelse(trf[[outcome]]>0,"SI","noSI"), levels=c("noSI","SI"))
      if (length(unique(trf$.bin)) < 2) next
      ff <- setdiff(names(trf), c(outcome, ".bin"))
      h <- hurdle_one_fit(trf, xvl, ff); if (is.null(h)) next
      s <- mean(abs(h$pred - vl[[outcome]]), na.rm = TRUE)
      if (!is.na(s) && s < best$score) best <- list(score=s, w=w, pv=pv)
    }
    if (best$score == Inf) next
    pp2 <- tryCatch(caret::preProcess(tv[,best$pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
    if (is.null(pp2)) next
    x_tv <- tryCatch(predict(pp2, tv[,best$pv,drop=FALSE]), error=function(e) NULL)
    x_te <- tryCatch(predict(pp2, te[,best$pv,drop=FALSE]), error=function(e) NULL)
    if (is.null(x_tv) || is.null(x_te)) next
    tvf <- cbind(x_tv, .y=tv[[outcome]]); names(tvf)[ncol(tvf)] <- outcome
    tvf$.bin <- factor(ifelse(tvf[[outcome]]>0,"SI","noSI"), levels=c("noSI","SI"))
    if (length(unique(tvf$.bin)) < 2) next
    ff <- setdiff(names(tvf), c(outcome, ".bin"))
    h <- hurdle_one_fit(tvf, x_te, ff); if (is.null(h)) next
    preds[[i]] <- data.frame(true=te[[outcome]], pred=h$pred,
                             p_si=h$p_si, true_bin=as.integer(te[[outcome]]>0))
  }
  if (length(preds) == 0) return(NULL)
  pd <- do.call(rbind, preds)
  if (length(unique(pd$true_bin)) < 2) return(NULL)
  pb <- as.integer(pd$p_si >= thr)
  TP <- sum(pb==1 & pd$true_bin==1); FN <- sum(pb==0 & pd$true_bin==1)
  TN <- sum(pb==0 & pd$true_bin==0); FP <- sum(pb==1 & pd$true_bin==0)
  sens <- if ((TP+FN)>0) TP/(TP+FN) else NA_real_
  spec <- if ((TN+FP)>0) TN/(TN+FP) else NA_real_
  auc  <- tryCatch(as.numeric(pROC::auc(pROC::roc(pd$true_bin, pd$p_si, quiet=TRUE))), error=function(e) NA_real_)
  met <- reg_metrics(pd$pred, pd$true)
  data.frame(Participant=pt$uid[1], auc=auc, sensitivity=sens, specificity=spec,
             balanced_acc=mean(c(sens,spec), na.rm=TRUE),
             cor=met$cor, r2=met$r2, mae=met$mae,
             n_obs=nrow(pd), n_SI=sum(pd$true_bin==1), n_nonSI=sum(pd$true_bin==0),
             mean_pred_SI=mean(pd$pred[pd$true_bin==1], na.rm=TRUE),
             mean_pred_nonSI=mean(pd$pred[pd$true_bin==0], na.rm=TRUE))
}
hurdle_pt_predictions_one_config <- function(pt, fb, config_name, cv_fn) {
  feats <- candidate_predictors(pt)
  if (length(feats) < 1 || var(pt[[outcome]], na.rm = TRUE) %in% c(0, NA)) return(NULL)
  n_si_pos <- sum(pt[[outcome]] > 0, na.rm = TRUE)
  n_si_neg <- sum(pt[[outcome]] == 0, na.rm = TRUE)
  if (n_si_pos < 10 || n_si_neg < 20) return(NULL)
  sl <- cv_fn(pt); if (is.null(sl) || length(sl$trains) < 1) return(NULL)
  preds <- list()
  for (i in seq_along(sl$trains)) {
    tr <- pt[sl$trains[[i]],]; vl <- pt[sl$vals[[i]],]; te <- pt[sl$tests[[i]],]
    scale_cols <- feats
    tr_ref <- tr
    tr <- apply_minmax_fold(tr_ref, tr, scale_cols)
    vl <- apply_minmax_fold(tr_ref, vl, scale_cols)
    te <- apply_minmax_fold(tr_ref, te, scale_cols)
    tv <- rbind(tr, vl); best <- list(score = Inf)
    for (w in WINDOWS) {
      pv <- intersect(fb(w, data_model_base), feats); if (length(pv) < 2) next
      pv <- drop_near_zero_variance_features(pv, tr, freqCut = 19, uniqueCut = 10); if (length(pv) < 2) next
      pp <- tryCatch(caret::preProcess(tr[,pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
      if (is.null(pp)) next
      xtr <- tryCatch(predict(pp, tr[,pv,drop=FALSE]), error=function(e) NULL)
      xvl <- tryCatch(predict(pp, vl[,pv,drop=FALSE]), error=function(e) NULL)
      if (is.null(xtr) || is.null(xvl)) next
      trf <- cbind(xtr, .y=tr[[outcome]]); names(trf)[ncol(trf)] <- outcome
      trf$.bin <- factor(ifelse(trf[[outcome]] > 0, "SI", "noSI"), levels=c("noSI","SI"))
      if (length(unique(trf$.bin)) < 2) next
      ff <- setdiff(names(trf), c(outcome, ".bin"))
      h <- hurdle_one_fit(trf, xvl, ff); if (is.null(h)) next
      s <- mean(abs(h$pred - vl[[outcome]]), na.rm = TRUE)
      if (!is.na(s) && s < best$score) best <- list(score=s, w=w, pv=pv)
    }
    if (best$score == Inf) next
    pp2 <- tryCatch(caret::preProcess(tv[,best$pv,drop=FALSE], method=c("zv","knnImpute"), k=5), error=function(e) NULL)
    if (is.null(pp2)) next
    x_tv <- tryCatch(predict(pp2, tv[,best$pv,drop=FALSE]), error=function(e) NULL)
    x_te <- tryCatch(predict(pp2, te[,best$pv,drop=FALSE]), error=function(e) NULL)
    if (is.null(x_tv) || is.null(x_te)) next
    tvf <- cbind(x_tv, .y=tv[[outcome]]); names(tvf)[ncol(tvf)] <- outcome
    tvf$.bin <- factor(ifelse(tvf[[outcome]] > 0, "SI", "noSI"), levels=c("noSI","SI"))
    if (length(unique(tvf$.bin)) < 2) next
    ff <- setdiff(names(tvf), c(outcome, ".bin"))
    h <- hurdle_one_fit(tvf, x_te, ff); if (is.null(h)) next
    preds[[i]] <- data.frame(
      Participant = pt$uid[1],
      config = config_name,
      fold = i,
      row_index = sl$tests[[i]],
      day_in_study = ceiling(sl$tests[[i]] / 3),
      true = te[[outcome]],
      pred = as.numeric(h$pred),
      p_si = as.numeric(h$p_si),
      true_bin = as.integer(te[[outcome]] > 0),
      stringsAsFactors = FALSE
    )
  }
  preds <- Filter(Negate(is.null), preds)
  if (!length(preds)) return(NULL)
  bind_rows(preds)
}

summarize_hurdle_prediction_panel <- function(pd) {
  if (is.null(pd) || !nrow(pd)) {
    return(data.frame(cor=NA_real_, rmse=NA_real_, auc=NA_real_))
  }
  met <- reg_metrics(pd$pred, pd$true)
  auc <- tryCatch({
    if (length(unique(pd$true_bin)) < 2) NA_real_
    else as.numeric(pROC::auc(pROC::roc(pd$true_bin, pd$p_si, quiet=TRUE)))
  }, error = function(e) NA_real_)
  data.frame(cor = met$cor, rmse = met$rmse * 100, auc = auc)
}

continuous_pt_predictions_one_config <- function(pt, fb, config_name, cv_fn) {
  feats <- candidate_predictors(pt)
  if (length(feats) < 1 || var(pt[[outcome]], na.rm = TRUE) %in% c(0, NA)) return(NULL)
  sl <- cv_fn(pt); if (is.null(sl) || length(sl$trains) < 1) return(NULL)
  preds <- list()
  for (i in seq_along(sl$trains)) {
    tr <- pt[sl$trains[[i]],]; vl <- pt[sl$vals[[i]],]; te <- pt[sl$tests[[i]],]
    scale_cols <- feats
    tr_ref <- tr
    tr <- apply_minmax_fold(tr_ref, tr, scale_cols)
    vl <- apply_minmax_fold(tr_ref, vl, scale_cols)
    te <- apply_minmax_fold(tr_ref, te, scale_cols)
    tv <- rbind(tr, vl); best <- list(score = Inf)
    for (w in WINDOWS) {
      pv <- intersect(fb(w, data_model_base), feats); if (length(pv) < 2) next
      pv <- drop_near_zero_variance_features(pv, tr, freqCut = 19, uniqueCut = 10); if (length(pv) < 2) next
      pp <- tryCatch(caret::preProcess(tr[, pv, drop = FALSE], method = c("zv", "knnImpute"), k = 5),
                     error = function(e) NULL)
      if (is.null(pp)) next
      xtr <- tryCatch(predict(pp, tr[, pv, drop = FALSE]), error = function(e) NULL)
      xvl <- tryCatch(predict(pp, vl[, pv, drop = FALSE]), error = function(e) NULL)
      if (is.null(xtr) || is.null(xvl)) next
      trf <- cbind(xtr, .y = tr[[outcome]]); names(trf)[ncol(trf)] <- outcome
      for (mt in MTRY_GRID) for (nt in NTREE_GRID) {
        s <- tryCatch({
          m <- randomForest::randomForest(
            as.formula(paste(outcome, "~ .")),
            data = trf,
            mtry = min(mt, ncol(xtr)),
            ntree = nt
          )
          mean(abs(predict(m, xvl) - vl[[outcome]]), na.rm = TRUE)
        }, error = function(e) NA_real_)
        if (!is.na(s) && s < best$score) best <- list(score = s, w = w, pv = pv, mt = mt, nt = nt)
      }
    }
    if (best$score == Inf) next
    pp2 <- tryCatch(caret::preProcess(tv[, best$pv, drop = FALSE], method = c("zv", "knnImpute"), k = 5),
                    error = function(e) NULL)
    if (is.null(pp2)) next
    x_tv <- tryCatch(predict(pp2, tv[, best$pv, drop = FALSE]), error = function(e) NULL)
    x_te <- tryCatch(predict(pp2, te[, best$pv, drop = FALSE]), error = function(e) NULL)
    if (is.null(x_tv) || is.null(x_te)) next
    tvf <- cbind(x_tv, .y = tv[[outcome]]); names(tvf)[ncol(tvf)] <- outcome
    mf <- tryCatch(randomForest::randomForest(
      as.formula(paste(outcome, "~ .")),
      data = tvf,
      mtry = min(best$mt, ncol(x_tv)),
      ntree = best$nt
    ), error = function(e) NULL)
    if (is.null(mf)) next
    p <- tryCatch(predict(mf, x_te), error = function(e) NULL)
    if (is.null(p)) next
    preds[[i]] <- data.frame(
      Participant = pt$uid[1],
      config = config_name,
      fold = i,
      row_index = sl$tests[[i]],
      day_in_study = ceiling(sl$tests[[i]] / 3),
      true = te[[outcome]],
      pred = pmax(0, pmin(1, as.numeric(p))),
      stringsAsFactors = FALSE
    )
  }
  preds <- Filter(Negate(is.null), preds)
  if (!length(preds)) return(NULL)
  bind_rows(preds)
}

summarize_continuous_prediction_panel <- function(pd) {
  if (is.null(pd) || !nrow(pd)) {
    return(data.frame(cor=NA_real_, rmse=NA_real_, r2=NA_real_))
  }
  met <- reg_metrics(pd$pred, pd$true)
  data.frame(cor = met$cor, rmse = met$rmse * 100, r2 = met$r2)
}

make_continuous_timeseries_for_participant <- function(uid, cv_params = CV_VARIANTS$main,
                                                       out_dir = report_all_pt_continuous_dir,
                                                       overwrite = TRUE) {
  make_dirs(out_dir)
  uid <- stringr::str_pad(as.character(uid), width = 4, pad = "0")
  out <- file.path(out_dir, sprintf("F3_continuous_timeseries_PT%s.png", uid))
  if (file.exists(out) && !isTRUE(overwrite)) {
    cat(sprintf("  exists, skipping PT %s: %s\n", uid, out))
    return(invisible(out))
  }
  pt <- data_model_base[data_model_base$uid == uid, , drop = FALSE]
  if (!nrow(pt)) return(invisible(NULL))
  cv_fn <- make_sliding_cv(cv_params$train_days, cv_params$val_days,
                           cv_params$test_days, cv_params$ema_per_day)
  config_order_local <- c("EMA_mean", "EMA_mean_dynamic",
                          "Opt_Passive_mean", "Opt_Passive_mean_dynamic",
                          "EMA_Opt_Passive_mean", "EMA_Opt_Passive_mean_dynamic")
  panel_preds <- bind_rows(lapply(config_order_local, function(cfg) {
    tryCatch(
      continuous_pt_predictions_one_config(pt, configs[[cfg]], cfg, cv_fn),
      error = function(e) {
        cat(sprintf("  skipped PT %s config %s: %s\n", uid, cfg, conditionMessage(e)))
        NULL
      }
    )
  }))
  if (!nrow(panel_preds)) return(invisible(NULL))
  panel_preds <- panel_preds %>%
    mutate(
      observed_100 = true * 100,
      predicted_100 = pred * 100
    )
  panel_stats <- panel_preds %>%
    group_by(config) %>%
    group_modify(~summarize_continuous_prediction_panel(.x)) %>%
    ungroup() %>%
    mutate(
      config_label = config_labels[match(config, config_levels)],
      facet_label = sprintf("%s\nr = %.2f | RMSE = %.2f | R² = %.2f",
                            config_label, cor, rmse, r2)
    )
  panel_preds <- panel_preds %>%
    left_join(panel_stats %>% select(config, facet_label), by = "config") %>%
    mutate(facet_label = factor(facet_label, levels = panel_stats$facet_label[match(config_order_local, panel_stats$config)]))
  best_row <- results_main %>%
    filter(status == "ok", n_folds >= MIN_FOLDS, model_type == "Continuous",
           Participant == uid, !is.na(r2)) %>%
    arrange(desc(r2)) %>%
    slice_head(n = 1)
  best_txt <- if (nrow(best_row)) {
    sprintf("Best config: %s (R² = %.2f)",
            config_labels[match(best_row$config[1], config_levels)], best_row$r2[1])
  } else {
    "Best config: not available"
  }
  day_range <- range(panel_preds$day_in_study, na.rm = TRUE)
  plot_df <- panel_preds %>%
    select(facet_label, day_in_study, observed_100, predicted_100) %>%
    pivot_longer(c(observed_100, predicted_100), names_to = "series", values_to = "q9_100") %>%
    mutate(series = recode(series, observed_100 = "Observed", predicted_100 = "Predicted"),
           series = factor(series, levels = c("Observed", "Predicted")))
  p <- ggplot(plot_df, aes(x = day_in_study, y = q9_100, color = series, group = series)) +
    geom_line(alpha = 0.48, linewidth = 0.45) +
    geom_point(size = 1.1, alpha = 0.9) +
    facet_wrap(~facet_label, ncol = 2) +
    scale_color_manual(values = c("Observed" = "#c4514a", "Predicted" = "#2f6473"), name = NULL) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(
      title = sprintf("Within-person time series — Participant %s (Continuous model)", uid),
      subtitle = sprintf("Predicted vs observed q9 over test days %.0f-%.0f (sliding-window CV requires %d prior days) | %s",
                         day_range[1], day_range[2], cv_params$train_days + cv_params$val_days, best_txt),
      x = "Day in study",
      y = "q9 (SI severity, 0-100)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "#fbfaf4", color = NA),
      panel.background = element_rect(fill = "#fbfaf4", color = NA),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank()
    )
  ggsave(out, p, width = 14, height = 8, dpi = 220, bg = "#fbfaf4")
  invisible(out)
}

make_continuous_timeseries_all_participants <- function(combined, variant_name = "main",
                                                        cv_params = CV_VARIANTS[[variant_name]],
                                                        participant_ids = NULL,
                                                        overwrite = TRUE) {
  uids <- combined %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Continuous", config != "Baseline") %>%
    pull(Participant) %>%
    unique() %>%
    as.character() %>%
    stringr::str_pad(width = 4, pad = "0")
  if (!is.null(participant_ids)) {
    participant_ids <- stringr::str_pad(as.character(participant_ids), width = 4, pad = "0")
    uids <- intersect(uids, participant_ids)
  }
  out_dir <- report_all_pt_continuous_dir
  make_dirs(out_dir)
  cat(sprintf("Saving continuous time-series plots for %d participants to: %s\n", length(uids), out_dir))
  outs <- lapply(uids, function(u) {
    tryCatch(make_continuous_timeseries_for_participant(u, cv_params, out_dir, overwrite = overwrite),
             error = function(e) {
               cat(sprintf("  skipped PT %s: %s\n", u, conditionMessage(e)))
               NULL
             })
  })
  outs <- unlist(outs)
  cat(sprintf("Saved %d continuous participant time-series plots.\n", length(outs)))
  invisible(outs)
}

make_hurdle_timeseries_for_participant <- function(uid, cv_params = CV_VARIANTS$main,
                                                   out_dir = report_all_pt_hurdle_dir,
                                                   overwrite = TRUE) {
  make_dirs(out_dir)
  uid <- stringr::str_pad(as.character(uid), width = 4, pad = "0")
  out <- file.path(out_dir, sprintf("F3_hurdle_timeseries_PT%s.png", uid))
  if (file.exists(out) && !isTRUE(overwrite)) {
    cat(sprintf("  exists, skipping PT %s: %s\n", uid, out))
    return(invisible(out))
  }
  pt <- data_model_base[data_model_base$uid == uid, , drop = FALSE]
  if (!nrow(pt)) return(invisible(NULL))
  cv_fn <- make_sliding_cv(cv_params$train_days, cv_params$val_days,
                           cv_params$test_days, cv_params$ema_per_day)
  config_order_local <- c("EMA_mean", "EMA_mean_dynamic",
                          "Opt_Passive_mean", "Opt_Passive_mean_dynamic",
                          "EMA_Opt_Passive_mean", "EMA_Opt_Passive_mean_dynamic")
  panel_preds <- bind_rows(lapply(config_order_local, function(cfg) {
    tryCatch(
      hurdle_pt_predictions_one_config(pt, configs[[cfg]], cfg, cv_fn),
      error = function(e) {
        cat(sprintf("  skipped PT %s config %s: %s\n", uid, cfg, conditionMessage(e)))
        NULL
      }
    )
  }))
  if (!nrow(panel_preds)) return(invisible(NULL))
  panel_preds <- panel_preds %>%
    mutate(
      observed_100 = true * 100,
      predicted_100 = pred * 100
    )
  panel_stats <- panel_preds %>%
    group_by(config) %>%
    group_modify(~summarize_hurdle_prediction_panel(.x)) %>%
    ungroup() %>%
    mutate(
      config_label = config_labels[match(config, config_levels)],
      facet_label = sprintf("%s\nr = %.2f | RMSE = %.2f | AUC = %.2f",
                            config_label, cor, rmse, auc)
    )
  panel_preds <- panel_preds %>%
    left_join(panel_stats %>% select(config, facet_label), by = "config") %>%
    mutate(facet_label = factor(facet_label, levels = panel_stats$facet_label[match(config_order_local, panel_stats$config)]))
  best_row <- results_main %>%
    filter(status == "ok", n_folds >= MIN_FOLDS, model_type == "Hurdle",
           Participant == uid, !is.na(auc)) %>%
    arrange(desc(auc)) %>%
    slice_head(n = 1)
  best_txt <- if (nrow(best_row)) {
    sprintf("Best config: %s (AUC = %.2f)",
            config_labels[match(best_row$config[1], config_levels)], best_row$auc[1])
  } else {
    "Best config: not available"
  }
  day_range <- range(panel_preds$day_in_study, na.rm = TRUE)
  plot_df <- panel_preds %>%
    select(facet_label, day_in_study, observed_100, predicted_100) %>%
    pivot_longer(c(observed_100, predicted_100), names_to = "series", values_to = "q9_100") %>%
    mutate(series = recode(series, observed_100 = "Observed", predicted_100 = "Predicted"),
           series = factor(series, levels = c("Observed", "Predicted")))
  p <- ggplot(plot_df, aes(x = day_in_study, y = q9_100, color = series, group = series)) +
    geom_line(alpha = 0.48, linewidth = 0.45) +
    geom_point(size = 1.1, alpha = 0.9) +
    facet_wrap(~facet_label, ncol = 2) +
    scale_color_manual(values = c("Observed" = "#c4514a", "Predicted" = "#2f6473"), name = NULL) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(
      title = sprintf("Within-person time series — Participant %s (Hurdle model)", uid),
      subtitle = sprintf("Predicted vs observed q9 over test days %.0f-%.0f (sliding-window CV requires %d prior days) | %s",
                         day_range[1], day_range[2], cv_params$train_days + cv_params$val_days, best_txt),
      x = "Day in study",
      y = "q9 (SI severity, 0-100)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "#fbfaf4", color = NA),
      panel.background = element_rect(fill = "#fbfaf4", color = NA),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank()
    )
  ggsave(out, p, width = 14, height = 8, dpi = 220, bg = "#fbfaf4")
  invisible(out)
}

make_hurdle_timeseries_all_participants <- function(combined, variant_name = "main",
                                                    cv_params = CV_VARIANTS[[variant_name]],
                                                    participant_ids = NULL,
                                                    overwrite = TRUE) {
  uids <- combined %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Hurdle", config != "Baseline") %>%
    pull(Participant) %>%
    unique() %>%
    as.character() %>%
    stringr::str_pad(width = 4, pad = "0")
  if (!is.null(participant_ids)) {
    participant_ids <- stringr::str_pad(as.character(participant_ids), width = 4, pad = "0")
    uids <- intersect(uids, participant_ids)
  }
  out_dir <- report_all_pt_hurdle_dir
  make_dirs(out_dir)
  cat(sprintf("Saving hurdle time-series plots for %d participants to: %s\n", length(uids), out_dir))
  outs <- lapply(uids, function(u) {
    tryCatch(make_hurdle_timeseries_for_participant(u, cv_params, out_dir, overwrite = overwrite),
             error = function(e) {
               cat(sprintf("  skipped PT %s: %s\n", u, conditionMessage(e)))
               NULL
             })
  })
  outs <- unlist(outs)
  cat(sprintf("Saved %d participant time-series plots.\n", length(outs)))
  invisible(outs)
}

make_F2_hurdle <- function(combined, variant_name, cv_params = NULL) {
  combined <- ensure_result_metric_cols(combined)
  if (is.null(cv_params)) cv_params <- CV_VARIANTS[[variant_name]]
  config_label <- dplyr::recode(
    BEST_HURDLE_CFG,
    "EMA_mean" = "EMA — Mean",
    "EMA_mean_dynamic" = "EMA — Mean + Dynamic",
    "Opt_Passive_mean" = "Passive — Mean",
    "Opt_Passive_mean_dynamic" = "Passive — Mean + Dynamic",
    "EMA_Opt_Passive_mean" = "EMA + Passive — Mean",
    "EMA_Opt_Passive_mean_dynamic" = "EMA + Passive — Mean + Dynamic",
    .default = gsub("_", " + ", BEST_HURDLE_CFG)
  )
  best_df <- combined %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Hurdle", config == BEST_HURDLE_CFG,
           !is.na(auc), !is.na(sensitivity), !is.na(specificity), !is.na(balanced_acc))
  if (nrow(best_df) == 0) {
    cat("⚠ No rows for best hurdle — skipping F2_hurdle\n")
    return(invisible(NULL))
  }

  cv_fn <- make_sliding_cv(cv_params$train_days, cv_params$val_days,
                           cv_params$test_days, cv_params$ema_per_day)
  obs_counts <- bind_rows(lapply(best_df$Participant, function(u) {
    pt <- data_model_base[data_model_base$uid == u, , drop = FALSE]
    sl <- cv_fn(pt)
    if (is.null(sl) || !length(sl$tests)) {
      return(data.frame(Participant = u, n_SI = NA_integer_, n_nonSI = NA_integer_))
    }
    test_idx <- unique(unlist(sl$tests))
    y <- pt[[outcome]][test_idx]
    data.frame(
      Participant = u,
      n_SI = sum(y > 0, na.rm = TRUE),
      n_nonSI = sum(y == 0, na.rm = TRUE)
    )
  }))
  best_df <- best_df %>%
    select(-any_of(c("n_SI_plot", "n_nonSI_plot"))) %>%
    left_join(
      obs_counts %>% rename(n_SI_plot = n_SI, n_nonSI_plot = n_nonSI),
      by = "Participant"
    ) %>%
    mutate(
      n_SI = coalesce(.data$n_SI, .data$n_SI_plot),
      n_nonSI = coalesce(.data$n_nonSI, .data$n_nonSI_plot),
      n_obs_plot = coalesce(.data$n_SI, 0L) + coalesce(.data$n_nonSI, 0L)
    ) %>%
    select(-any_of(c("n_SI_plot", "n_nonSI_plot")))

  pct_auc <- mean(best_df$auc > 0.5, na.rm = TRUE) * 100
  metric_map <- c(auc = "AUC", sensitivity = "Sensitivity",
                  specificity = "Specificity", balanced_acc = "Balanced acc.")
  metric_palette <- c("AUC" = "grey20", "Sensitivity" = "#7fb5a3",
                      "Specificity" = "#e7b64b", "Balanced acc." = "#e66a7a")
  dot_df <- best_df %>% arrange(auc) %>%
    mutate(rank = row_number()) %>%
    select(rank, all_of(names(metric_map))) %>%
    pivot_longer(all_of(names(metric_map)), names_to = "metric_raw", values_to = "value") %>%
    mutate(metric = factor(metric_map[metric_raw], levels = unname(metric_map)))
  best_df <- best_df %>% arrange(auc) %>% mutate(rank = row_number())
  total_si <- sum(best_df$n_SI, na.rm = TRUE)
  total_nonsi <- sum(best_df$n_nonSI, na.rm = TRUE)

  p_dot <- ggplot(dot_df, aes(x = value, y = rank, color = metric)) +
    geom_hline(yintercept = seq_len(nrow(best_df)), color = "#e7e5dc", linewidth = 0.25) +
    geom_vline(xintercept = 0.5, color = "grey55", linetype = "dashed") +
    geom_point(alpha = 0.75, size = 1.2) +
    scale_color_manual(values = metric_palette) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    annotate("text", x = 0.04, y = max(dot_df$rank) * 0.06,
             label = sprintf("%.0f%%", pct_auc), size = 4, fontface = "bold") +
    labs(x = "Value", y = "Participants", color = NULL,
         title = sprintf("Within-person performance — %s (Hurdle model)", config_label),
         subtitle = sprintf("Binary classification  |  threshold: any SI (q9 > 0)\n%.0f%% of participants with an AUC > 0.5\nMean AUC (within-person) = %.2f\nTotal observations: SI = %s | Non-SI = %s  |  N = %d participants",
                            pct_auc, mean(best_df$auc, na.rm = TRUE),
                            format(total_si, big.mark = ","),
                            format(total_nonsi, big.mark = ","),
                            nrow(best_df))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())

  obs_df <- best_df %>%
    transmute(
      rank,
      n_nonSI = coalesce(n_nonSI, 0L),
      n_SI = coalesce(n_SI, 0L),
      n_total = n_nonSI + n_SI
    )
  p_obs <- ggplot(obs_df, aes(y = rank)) +
    geom_col(aes(x = n_total, fill = "Non-SI"), width = 0.85) +
    geom_col(aes(x = n_SI, fill = "SI"), width = 0.85) +
    scale_fill_manual(values = c("Non-SI" = "#bdbbb3", "SI" = "#e34a3a"), name = NULL) +
    labs(x = "Observations", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  density_order <- c("Specificity", "AUC", "Sensitivity", "Balanced acc.")
  density_panels <- lapply(density_order, function(mm) {
    dd <- dot_df %>% filter(metric == mm, is.finite(value))
    med <- median(dd$value, na.rm = TRUE)
    ggplot(dd, aes(x = value, y = 0, fill = metric)) +
      geom_density(alpha = 0.85, color = "#555555", linewidth = 0.25) +
      geom_jitter(aes(y = -0.015), height = 0.012, alpha = 0.35, size = 0.45,
                  color = metric_palette[[mm]]) +
      geom_vline(xintercept = med, color = "#444444", linewidth = 0.35) +
      annotate("text", x = med, y = 0.78, label = sprintf("%.2f", med),
               size = 2.8, hjust = -0.05) +
      annotate("text", x = -Inf, y = 0.45, label = mm, hjust = -0.05, size = 3.3) +
      scale_fill_manual(values = metric_palette, guide = "none") +
      coord_cartesian(xlim = c(0, 1.2), ylim = c(-0.04, 1), clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_void(base_size = 9) +
      theme(plot.margin = ggplot2::margin(2, 8, 2, 8))
  })
  p_density <- wrap_plots(density_panels, ncol = 1) +
    plot_annotation(caption = sprintf("N = %d participants.", nrow(best_df)),
                    theme = theme(plot.caption = element_text(hjust = 1, size = 8, color = "grey35")))

  fig <- p_dot + p_obs + p_density + plot_layout(widths = c(2.2, 1.15, 1.35)) +
    plot_annotation(theme = theme(
      plot.background = element_rect(fill = "#fbfaf4", color = "#eadb9a", linewidth = 1.1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ))
  report_out <- file.path(report_dir, sprintf("Figure_hurdle_model_performance_%s.png", variant_name))
  ggsave(report_out, fig, width = 10, height = 8.5, dpi = 220, bg = "#fbfaf4")
  cat(sprintf("Saved report hurdle figure: %s\n", report_out))
  invisible(report_out)
}

if (RUN_HURDLE_MODELS && RERUN_F2_HURDLE) make_F2_hurdle(results_main,  "main",  CV_VARIANTS$main)

make_F2_reference_style <- function(combined, variant_name = "main") {
  combined <- ensure_result_metric_cols(combined)
  build_row <- function(df, model_label, title, subtitle, metrics, palette,
                        order_metric, pct_metric, pct_cut, pct_label,
                        out_caption = TRUE) {
    df <- df %>% arrange(.data[[order_metric]]) %>% mutate(rank = row_number())
    long <- df %>%
      select(rank, n_obs, all_of(names(metrics))) %>%
      pivot_longer(all_of(names(metrics)), names_to = "metric_raw", values_to = "value") %>%
      mutate(metric = factor(metrics[metric_raw], levels = unname(metrics)))
    pct <- mean(df[[pct_metric]] > pct_cut, na.rm = TRUE) * 100
    stat_lines <- sprintf(
      "%.0f%% of participants with %s > %s\nMean %s = %.3f\nMedian MAE = %.3f\nTotal test observations: %s",
      pct, pct_label, pct_cut,
      ifelse(pct_metric == "auc", "AUC", "R²"),
      mean(df[[pct_metric]], na.rm = TRUE),
      median(df$mae, na.rm = TRUE),
      format(sum(df$n_obs), big.mark = ",")
    )
    x_rng <- range(long$value, na.rm = TRUE)
    if (model_label == "continuous") {
      x_left <- max(-0.5, stats::quantile(long$value, 0.01, na.rm = TRUE))
      x_right <- min(1, max(1, stats::quantile(long$value, 0.99, na.rm = TRUE)))
    } else {
      x_left <- 0
      x_right <- 1
    }
    p_dot <- ggplot(long, aes(x = value, y = rank, color = metric)) +
      geom_hline(yintercept = seq_len(nrow(df)), color = "#e7e5dc", linewidth = 0.25) +
      geom_vline(xintercept = ifelse(model_label == "hurdle", 0.5, 0),
                 color = "#777777", linewidth = 0.45) +
      geom_point(alpha = 0.86, size = 1.15) +
      annotate("text", x = ifelse(model_label == "hurdle", 0.08, x_left + 0.05),
               y = max(long$rank, na.rm = TRUE) * 0.08,
               label = sprintf("%.0f%%", pct), size = 4, fontface = "bold") +
      scale_color_manual(values = palette, name = NULL) +
      coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
      labs(x = "Value", y = "Participants", title = title, subtitle = paste0(subtitle, "\n\n", stat_lines)) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.5, lineheight = 1.05),
        legend.position = "top",
        legend.justification = "left",
        legend.key.width = unit(0.35, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = ggplot2::margin(8, 8, 8, 8)
      )

    p_obs <- ggplot(df, aes(x = n_obs, y = rank)) +
      geom_col(fill = "#a8a69e", width = 0.85) +
      labs(x = "Test observations", y = NULL) +
      theme_minimal(base_size = 9) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = ggplot2::margin(62, 8, 8, 8)
      )

    dens <- lapply(unname(metrics), function(mm) {
      dd <- long %>% filter(metric == mm, is.finite(value))
      med <- median(dd$value, na.rm = TRUE)
      ggplot(dd, aes(x = value, y = 0, fill = metric)) +
        geom_density(alpha = 0.85, color = "#555555", linewidth = 0.25) +
        geom_jitter(aes(y = -0.015), height = 0.012, alpha = 0.35, size = 0.35, color = "#333333") +
        geom_vline(xintercept = med, color = "#444444", linewidth = 0.35) +
        annotate("text", x = med, y = 0.78, label = sprintf("%.3f", med), size = 2.5, hjust = -0.05) +
        annotate("text", x = -Inf, y = 0.45, label = mm, hjust = -0.05, size = 3) +
        scale_fill_manual(values = palette, guide = "none") +
        coord_cartesian(xlim = c(x_left, x_right), ylim = c(-0.04, 1), clip = "off") +
        labs(x = NULL, y = NULL) +
        theme_void(base_size = 8) +
        theme(plot.margin = ggplot2::margin(2, 8, 2, 8))
    })
    p_density <- wrap_plots(dens, ncol = 1) +
      plot_annotation(caption = if (out_caption) sprintf("N = %d participants.", nrow(df)) else NULL,
                      theme = theme(plot.caption = element_text(hjust = 1, size = 8, color = "grey35")))

    p_dot + p_obs + p_density + plot_layout(widths = c(2.2, 1.15, 1.35))
  }

  cont_df <- combined %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Continuous", config == BEST_CONT_CFG,
           !is.na(r2), !is.na(cor), !is.na(mae))
  hurdle_df <- combined %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Hurdle", config == BEST_HURDLE_CFG,
           !is.na(auc))
  has_hurdle_class_metrics <- any(is.finite(hurdle_df$sensitivity)) &&
    any(is.finite(hurdle_df$specificity)) &&
    any(is.finite(hurdle_df$balanced_acc))
  if (has_hurdle_class_metrics) {
    hurdle_df <- hurdle_df %>%
      filter(!is.na(sensitivity), !is.na(specificity), !is.na(balanced_acc))
    hurdle_metrics <- c(auc = "AUC", sensitivity = "Sensitivity",
                        specificity = "Specificity", balanced_acc = "Balanced acc.")
    hurdle_palette <- c("AUC" = "#252525", "Sensitivity" = "#7fb5a3",
                        "Specificity" = "#e7b64b", "Balanced acc." = "#e66a7a")
  } else {
    hurdle_df <- hurdle_df %>% filter(!is.na(cor), !is.na(mae))
    hurdle_metrics <- c(auc = "AUC", cor = "Pearson r", mae = "MAE (norm)")
    hurdle_palette <- c("AUC" = "#252525", "Pearson r" = "#7fb5a3", "MAE (norm)" = "#e7b64b")
  }
  if (!nrow(cont_df) || !nrow(hurdle_df)) {
    cat("⚠ Missing continuous or hurdle rows — skipping reference-style F2\n")
    return(invisible(NULL))
  }
  continuous_row <- build_row(
    cont_df,
    model_label = "continuous",
    title = sprintf("Within-person performance — %s", gsub("_", " + ", BEST_CONT_CFG)),
    subtitle = "Continuous regression | outcome: SI severity (q9, scaled 0-1)",
    metrics = c(r2 = "R² (signed)", cor = "Pearson r", mae = "MAE (norm)"),
    palette = c("R² (signed)" = "#252525", "Pearson r" = "#7fb5a3", "MAE (norm)" = "#e7b64b"),
    order_metric = "r2",
    pct_metric = "r2",
    pct_cut = 0,
    pct_label = "R²",
    out_caption = FALSE
  )
  hurdle_row <- build_row(
    hurdle_df,
    model_label = "hurdle",
    title = sprintf("Within-person performance — %s", gsub("_", " + ", BEST_HURDLE_CFG)),
    subtitle = "Hurdle model | binary SI occurrence plus conditional severity",
    metrics = hurdle_metrics,
    palette = hurdle_palette,
    order_metric = "auc",
    pct_metric = "auc",
    pct_cut = 0.5,
    pct_label = "AUC",
    out_caption = TRUE
  )
  fig <- continuous_row / hurdle_row +
    plot_annotation(theme = theme(
      plot.background = element_rect(fill = "#fbfaf4", color = "#eadb9a", linewidth = 1.1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ))
  out <- file.path(plots_dir, sprintf("F2_withinperson_reference_style_%s.png", variant_name))
  ggsave(out, fig, width = 10, height = 15, dpi = 220, bg = "#fbfaf4")
  cat(sprintf("Saved reference-style within-person figure: %s\n", out))
  invisible(out)
}

if (isTRUE(RERUN_F2_REFERENCE_STYLE)) {
  make_F2_reference_style(results_main, "main")
}
################################################################################
cat(sprintf("   • Table 1:       %s\n", file.path(results_dir, "Table1_demographics.csv")))

build_final_report_folder <- function() {
  make_dirs(
    report_dir,
    report_all_pt_dir,
    report_all_pt_continuous_dir,
    report_all_pt_hurdle_dir
  )

  copy_if_exists <- function(from, to) {
    if (file.exists(from)) {
      if (normalizePath(from, mustWork = FALSE) != normalizePath(to, mustWork = FALSE)) {
        file.copy(from, to, overwrite = TRUE)
      }
      TRUE
    } else {
      warning(sprintf("Missing expected report source file: %s", from))
      FALSE
    }
  }

  stale_report_files <- c(
    file.path(report_dir, "Figure_demographic_description.png"),
    file.path(report_dir, "Table1_demographics.docx"),
    file.path(report_dir, "Table1_demographics_continuous_model_sample.csv"),
    file.path(report_dir, "Table1_demographics_continuous_model_sample.docx"),
    file.path(report_dir, "Table1_demographics_continuous_model_sample.png"),
    file.path(report_dir, "README_report_folder.txt"),
    file.path(report_dir, "report_manifest.csv"),
    file.path(report_dir, "final_manuscript_output_manifest.csv"),
    file.path(report_dir, "All_models_results_main.csv"),
    file.path(report_dir, "baseline_comparison_tests_main.csv"),
    file.path(report_dir, "continuous_model_failure_diagnostic_by_config_main.csv"),
    file.path(report_dir, "continuous_model_failure_diagnostic_by_participant_main.csv"),
    file.path(report_dir, "Table1_demographics_PT143.png"),
    file.path(report_dir, "Table2_predictive_performance_main.csv"),
    file.path(report_dir, "Table2_predictive_performance_main.docx"),
    file.path(report_dir, "Table2_continuous_model_results_main.docx"),
    file.path(report_dir, "Table3_hurdle_model_results_main.docx"),
    file.path(plots_dir, "F1_model_performance_main.png"),
    file.path(plots_dir, "F2_withinperson_continuous_main.png"),
    file.path(plots_dir, "F2_withinperson_hurdle_main.png"),
    file.path(plots_dir, "F2_reference_style_main.png"),
    file.path(plots_dir, "F3_si_obs_counts_main.png"),
    file.path(plots_dir, "F8_demographics.png")
  )
  stale_report_files <- stale_report_files[file.exists(stale_report_files)]
  if (length(stale_report_files)) {
    unlink(stale_report_files)
    cat(sprintf("Removed %d stale report files that are no longer part of final outputs.\n",
                length(stale_report_files)))
  }
  # 1. Two results tables: continuous and hurdle models.
  if (exists("table2_predictive_performance_main")) {
    perf_table <- table2_predictive_performance_main
    write.csv(perf_table %>% filter(`Model type` == "Continuous") %>% select(-`Model type`),
              file.path(report_dir, "Table2_continuous_model_results_main.csv"),
              row.names = FALSE)
    write.csv(perf_table %>% filter(`Model type` == "Hurdle") %>% select(-`Model type`),
              file.path(report_dir, "Table3_hurdle_model_results_main.csv"),
              row.names = FALSE)
  } else {
    warning("table2_predictive_performance_main does not exist; result CSV tables were not rebuilt.")
  }
  if (!file.exists(file.path(report_dir, "Figure_hurdle_model_performance_main.png"))) {
    warning("Missing expected hurdle model figure: Figure_hurdle_model_performance_main.png")
  }

  # 3. all_PT folders for participant-level prediction plots.
  all_pt_continuous_files <- list.files(
    report_all_pt_continuous_dir,
    pattern = "^F3_continuous_timeseries_PT.*\\.png$",
    full.names = TRUE
  )
  all_pt_hurdle_files <- list.files(
    report_all_pt_hurdle_dir,
    pattern = "^F3_hurdle_timeseries_PT.*\\.png$",
    full.names = TRUE
  )
  if (!length(all_pt_continuous_files)) {
    message(sprintf(
      "No continuous all_PT plots found yet in %s. Set RUN_ALL_PT_CONTINUOUS_TIMESERIES <- TRUE and source this script to populate it.",
      report_all_pt_continuous_dir
    ))
  }
  if (!length(all_pt_hurdle_files)) {
    message(sprintf(
      "No hurdle all_PT plots found yet in %s. Set RUN_ALL_PT_HURDLE_TIMESERIES <- TRUE and source this script to populate it.",
      report_all_pt_hurdle_dir
    ))
  }

  cont_n <- results_main %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Continuous", config %in% MODEL_CONFIGS) %>%
    summarise(n = n_distinct(Participant), .groups = "drop") %>%
    pull(n)
  hurdle_n <- results_main %>%
    filter(status == "ok", n_folds >= MIN_FOLDS,
           model_type == "Hurdle", config == BEST_HURDLE_CFG) %>%
    summarise(n = n_distinct(Participant), .groups = "drop") %>%
    pull(n)

  manifest <- data.frame(
    item = c(
      "Demographic table CSV",
      "Continuous model results CSV",
      "Hurdle model results CSV",
      "Hurdle model figure",
      "All participant continuous time-series folder",
      "All participant hurdle time-series folder"
    ),
    path = c(
      file.path(report_dir, "Table1_demographics.csv"),
      file.path(report_dir, "Table2_continuous_model_results_main.csv"),
      file.path(report_dir, "Table3_hurdle_model_results_main.csv"),
      file.path(report_dir, "Figure_hurdle_model_performance_main.png"),
      report_all_pt_continuous_dir,
      report_all_pt_hurdle_dir
    ),
    exists = c(
      file.exists(file.path(report_dir, "Table1_demographics.csv")),
      file.exists(file.path(report_dir, "Table2_continuous_model_results_main.csv")),
      file.exists(file.path(report_dir, "Table3_hurdle_model_results_main.csv")),
      file.exists(file.path(report_dir, "Figure_hurdle_model_performance_main.png")),
      dir.exists(report_all_pt_continuous_dir),
      dir.exists(report_all_pt_hurdle_dir)
    ),
    stringsAsFactors = FALSE
  )
  cat(sprintf("Final report folder assembled: %s\n", report_dir))
  print(manifest)
  invisible(manifest)
}

final_report_manifest <- build_final_report_folder()
}

if (RUN_FEATURE_IMPORTANCE) {
cat("\n===== SECTION: EXPLORATORY SHAP FEATURE IMPORTANCE =====\n")
################################################################################
# Exploratory feature importance among successful idiographic models
#
# This script does not rerun the full prediction pipeline.
# It reads:
#   - all_features_included_filtered.csv
#   - Results/All_models_results_main.csv
#
# It estimates participant-level approximate SHAP/Shapley feature importance for:
#   1. Continuous models with R² > 0
#   2. Hurdle models with AUC > 0.5
#
# Configs:
#   - EMA_mean_dynamic
#   - Opt_Passive_mean_dynamic
#   - EMA_Opt_Passive_mean_dynamic
#
# Interpretation:
#   Exploratory/subgroup interpretability only. This identifies which feature
#   families were important among participants whose model already performed
#   above the selected criterion. Importance is estimated as mean absolute
#   approximate SHAP value from participant-level random forests. Grouped
#   winners use per-feature normalized importance to avoid favoring
#   sources/sensors with more columns.
################################################################################

BASE_PATH <- resolve_project_base()
CLEAN_DATA_DIR <- file.path(BASE_PATH, "3. Cleaned data")
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
RESULTS_DIR <- file.path(RESULTS_ROOT_DIR, "report")
MODEL_READY_FILE <- file.path(CLEAN_DATA_DIR, "all_features_included_filtered.csv")
RESULTS_FILE <- file.path(RESULTS_ROOT_DIR, "model_config_cache", "All_models_results_main.csv")
if (!file.exists(RESULTS_FILE)) {
  RESULTS_FILE <- file.path(RESULTS_DIR, "All_models_results_main.csv")
}
OUT_DIR <- if (isTRUE(KEEP_FEATURE_IMPORTANCE_DETAILS)) {
  file.path(RESULTS_DIR, "feature_selection_analysis", "exploratory_successful_models")
} else {
  file.path(tempdir(), "si_prediction_feature_importance_work")
}
FEATURE_IMPORTANCE_WORK_DIR <- OUT_DIR
STATUS_FILE <- file.path(OUT_DIR, "exploratory_feature_importance_status.csv")
make_dirs(OUT_DIR)

required_packages <- c("data.table", "dplyr", "stringr", "caret", "randomForest", "ggplot2")
require_packages_or_stop(required_packages, "feature importance")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(caret)
  library(randomForest)
  library(ggplot2)
})

outcome <- "q9"
set.seed(12361488)
SHAP_NSIM <- 10
SHAP_MAX_ROWS <- 40
importance_failures <- list()

record_importance_failure <- function(participant, model_type, config, reason) {
  key <- sprintf("%s_%s_%s_%s", model_type, config, participant, length(importance_failures) + 1L)
  importance_failures[[key]] <<- data.frame(
    Participant = as.character(participant),
    model_type = model_type,
    config = config,
    reason = as.character(reason),
    stringsAsFactors = FALSE
  )
  NULL
}

sensor_labels <- c(
  acc = "Accelerometer",
  hr = "Heart rate",
  steps = "Steps",
  battery = "Battery",
  convo = "Conversation",
  gps = "GPS / location",
  gps_home = "Time at home",
  sleep = "Sleep",
  sms = "SMS",
  calls = "Calls"
)

ema_labels <- c(
  q1 = "Anhedonia",
  q2 = "Depressed mood",
  q3 = "Sleep",
  q4 = "Low energy",
  q5 = "Appetite",
  q6 = "Guilt/worthlessness",
  q7 = "Concentration",
  q8 = "Psychomotor",
  q9 = "SI/self-harm",
  g1 = "Anxiety 1",
  g2 = "Anxiety 2"
)

drop_near_zero_variance_features <- function(cols, training_data,
                                             freqCut = 19,
                                             uniqueCut = 10) {
  cols <- intersect(cols, names(training_data))
  if (!length(cols)) return(character())
  numeric_cols <- cols[vapply(cols, function(cc) is.numeric(training_data[[cc]]), logical(1))]
  if (!length(numeric_cols)) return(character())
  nzv <- tryCatch(
    caret::nearZeroVar(training_data[, numeric_cols, drop = FALSE],
                       freqCut = freqCut, uniqueCut = uniqueCut),
    error = function(e) integer(0)
  )
  if (length(nzv)) numeric_cols <- numeric_cols[-nzv]
  numeric_cols
}

approx_mean_abs_shap <- function(object, X, pred_wrapper, nsim = SHAP_NSIM) {
  X <- as.data.frame(X)
  p <- ncol(X)
  n <- nrow(X)
  if (!p || !n) return(numeric())
  importance_sum <- setNames(rep(0, p), names(X))
  importance_n <- setNames(rep(0, p), names(X))

  for (ii in seq_len(n)) {
    x_interest <- X[ii, , drop = FALSE]
    for (ss in seq_len(nsim)) {
      background <- X[sample.int(n, 1L), , drop = FALSE]
      feature_order <- sample.int(p)
      previous <- background
      previous_pred <- pred_wrapper(object, previous)

      for (jj in feature_order) {
        current <- previous
        current[[jj]] <- x_interest[[jj]]
        current_pred <- pred_wrapper(object, current)
        contribution <- as.numeric(current_pred - previous_pred)
        importance_sum[[jj]] <- importance_sum[[jj]] + abs(contribution)
        importance_n[[jj]] <- importance_n[[jj]] + 1
        previous <- current
        previous_pred <- current_pred
      }
    }
  }
  importance_sum / pmax(importance_n, 1)
}

# Match the main prediction section above:
# EMA_mean = prior-prompt lagged EMA items, including q9_lag.
# EMA_dynamic = rolling EMA means/MADs such as q*_mean_d1, q*_mean_d7, q*_mad_d7.
ema_mean <- c(paste0("q", 1:9, "_lag"), "g1_lag", "g2_lag")

ema_dynamic_cols <- function(data_names) {
  grep("^(q|g)[0-9]+_(mean_d|mad_d|p90p10_d)", data_names, value = TRUE)
}

passive_cols <- function(data_names) {
  grep("^(acc|hr|steps|battery|convo|gps_home|gps|sleep|sms|calls)_", data_names, value = TRUE)
}

passive_mean_cols <- function(data_names) {
  pc <- passive_cols(data_names)
  pc[grepl("(_mean_val$|_sum_val$|^sleep_mean_)", pc)]
}

passive_dynamic_cols <- function(data_names) {
  pc <- passive_cols(data_names)
  pc[grepl("(_mad|_p90p10|_entropy)", pc)]
}

get_config_predictors <- function(config, data_names) {
  ema_all <- intersect(c(ema_mean, ema_dynamic_cols(data_names)), data_names)
  passive_all <- intersect(c(passive_mean_cols(data_names), passive_dynamic_cols(data_names)), data_names)
  if (config == "EMA_mean_dynamic") return(ema_all)
  if (config == "Opt_Passive_mean_dynamic") return(passive_all)
  if (config == "EMA_Opt_Passive_mean_dynamic") return(unique(c(ema_all, passive_all)))
  stop(sprintf("Unsupported config for this exploratory script: %s", config))
}

classify_source <- function(feature) {
  if (!is.na(classify_ema_item(feature))) return("EMA")
  if (!is.na(classify_sensor(feature))) return("Passive sensor")
  "Other"
}

classify_sensor <- function(feature) {
  if (grepl("^gps_home_", feature)) return("Time at home")
  for (s in names(sensor_labels)[order(nchar(names(sensor_labels)), decreasing = TRUE)]) {
    if (grepl(paste0("^", s, "_"), feature)) return(sensor_labels[[s]])
  }
  NA_character_
}

classify_ema_item <- function(feature) {
  hit <- stringr::str_extract(feature, "^(q[1-9]|g[12])")
  if (is.na(hit)) return(NA_character_)
  ema_labels[[hit]]
}

classify_window <- function(feature) {
  if (grepl("(^|_)hr1($|_)", feature)) return("1 hour")
  if (grepl("(^|_)hr4($|_)", feature)) return("4 hours")
  if (grepl("(^|_)day1($|_)", feature) || grepl("_mean_d1$|_d1$", feature)) return("1 day")
  if (grepl("(^|_)day7($|_)", feature) || grepl("_mean_d7$|_mad_d7$|_d7$", feature)) return("7 days")
  if (grepl("_lag$", feature)) return("Prior prompt")
  NA_character_
}

classify_index <- function(feature) {
  if (grepl("entropy", feature)) return("Entropy")
  if (grepl("p90p10", feature)) return("P90-P10")
  if (grepl("mad", feature)) return("MAD")
  if (grepl("mean", feature)) return("Mean")
  if (grepl("sum", feature)) return("Sum")
  if (grepl("_lag$", feature)) return("Lagged EMA")
  NA_character_
}

fit_participant_importance <- function(pt, predictors, model_type, config) {
  uid_this <- if ("uid" %in% names(pt) && nrow(pt)) as.character(pt$uid[1]) else NA_character_
  pt <- pt[!is.na(pt[[outcome]]), , drop = FALSE]
  if (nrow(pt) < 10) {
    record_importance_failure(uid_this, model_type, config, "fewer than 10 observed outcome rows")
    return(NULL)
  }
  if (model_type == "Continuous") {
    if (is.na(stats::var(pt[[outcome]], na.rm = TRUE)) ||
        stats::var(pt[[outcome]], na.rm = TRUE) == 0) {
      record_importance_failure(uid_this, model_type, config, "zero continuous outcome variance")
      return(NULL)
    }
  } else {
    if (length(unique(as.integer(pt[[outcome]] > 0))) < 2) {
      record_importance_failure(uid_this, model_type, config, "only one binary outcome class")
      return(NULL)
    }
  }

  pv <- drop_near_zero_variance_features(predictors, pt, freqCut = 19, uniqueCut = 10)
  if (length(pv) < 2) {
    record_importance_failure(uid_this, model_type, config, "fewer than 2 predictors after nearZeroVar")
    return(NULL)
  }

  k_use <- max(1, min(5, nrow(pt) - 1))
  pp <- tryCatch(
    caret::preProcess(pt[, pv, drop = FALSE], method = c("zv", "knnImpute"), k = k_use),
    error = function(e) {
      record_importance_failure(uid_this, model_type, config, paste("preProcess failed:", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(pp)) return(NULL)

  x <- tryCatch(
    predict(pp, pt[, pv, drop = FALSE]),
    error = function(e) {
      record_importance_failure(uid_this, model_type, config, paste("imputation/predict failed:", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(x) || ncol(x) < 2) {
    record_importance_failure(uid_this, model_type, config, "fewer than 2 predictors after imputation")
    return(NULL)
  }

  if (model_type == "Continuous") {
    fit_df <- cbind(x, .y = pt[[outcome]])
    names(fit_df)[ncol(fit_df)] <- outcome
    rf <- tryCatch(
      randomForest::randomForest(
        as.formula(paste(outcome, "~ .")),
        data = fit_df,
        mtry = min(5, ncol(x)),
        ntree = 500,
        importance = TRUE
      ),
      error = function(e) {
        record_importance_failure(uid_this, model_type, config, paste("random forest failed:", conditionMessage(e)))
        NULL
      }
    )
  } else {
    fit_df <- cbind(x, .bin = factor(ifelse(pt[[outcome]] > 0, "SI", "noSI"),
                                     levels = c("noSI", "SI")))
    rf <- tryCatch(
      randomForest::randomForest(
        .bin ~ .,
        data = fit_df,
        mtry = min(5, ncol(x)),
        ntree = 500,
        importance = TRUE
      ),
      error = function(e) {
        record_importance_failure(uid_this, model_type, config, paste("random forest failed:", conditionMessage(e)))
        NULL
      }
    )
  }
  if (is.null(rf)) return(NULL)

  shap_x <- x
  if (nrow(shap_x) > SHAP_MAX_ROWS) {
    shap_x <- shap_x[sample(seq_len(nrow(shap_x)), SHAP_MAX_ROWS), , drop = FALSE]
  }
  pred_wrapper <- if (model_type == "Continuous") {
    function(object, newdata) as.numeric(predict(object, newdata = newdata))
  } else {
    function(object, newdata) as.numeric(predict(object, newdata = newdata, type = "prob")[, "SI"])
  }
  shap_importance <- tryCatch(
    approx_mean_abs_shap(
      object = rf,
      X = as.data.frame(shap_x),
      pred_wrapper = pred_wrapper,
      nsim = SHAP_NSIM
    ),
    error = function(e) {
      record_importance_failure(uid_this, model_type, config, paste("approximate SHAP failed:", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(shap_importance) || !length(shap_importance)) return(NULL)

  data.frame(
    feature = names(shap_importance),
    importance = as.numeric(shap_importance),
    importance_positive = as.numeric(shap_importance),
    importance_column = "mean_abs_approx_shap",
    modal_window = vapply(names(shap_importance), classify_window, character(1)),
    source = vapply(names(shap_importance), classify_source, character(1)),
    sensor = vapply(names(shap_importance), classify_sensor, character(1)),
    ema_item = vapply(names(shap_importance), classify_ema_item, character(1)),
    index = vapply(names(shap_importance), classify_index, character(1)),
    model_type = model_type,
    config = config,
    importance_method = sprintf("participant_rf_monte_carlo_approx_shap_nsim%s_nrows%s", SHAP_NSIM, nrow(shap_x)),
    row.names = NULL
  )
}

summarize_importance <- function(importance, model_type, config, criterion_label, out_stub) {
  if (!nrow(importance)) {
    empty <- data.frame()
    write.csv(empty, file.path(OUT_DIR, sprintf("summary_feature_%s.csv", out_stub)), row.names = FALSE)
    return(invisible(NULL))
  }

  group_summary <- function(group_col) {
    importance %>%
      filter(!is.na(.data[[group_col]]), .data[[group_col]] != "") %>%
      group_by(Participant, level = .data[[group_col]], feature) %>%
      summarise(feature_importance = mean(importance_positive, na.rm = TRUE), .groups = "drop") %>%
      group_by(Participant, level) %>%
      summarise(
        total_importance = sum(feature_importance, na.rm = TRUE),
        mean_importance_per_feature = mean(feature_importance, na.rm = TRUE),
        n_features = n_distinct(feature),
        .groups = "drop"
      ) %>%
      group_by(level) %>%
      summarise(
        n_participants = n_distinct(Participant),
        median_total_importance = median(total_importance, na.rm = TRUE),
        mean_total_importance = mean(total_importance, na.rm = TRUE),
        median_importance_per_feature = median(mean_importance_per_feature, na.rm = TRUE),
        mean_importance_per_feature = mean(mean_importance_per_feature, na.rm = TRUE),
        median_n_features = median(n_features, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(median_importance_per_feature), desc(mean_importance_per_feature))
  }

  summaries <- list(
    feature = importance %>%
      group_by(feature) %>%
      summarise(
        n_participants = n_distinct(Participant),
        median_importance = median(importance_positive, na.rm = TRUE),
        mean_importance = mean(importance_positive, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        source = vapply(feature, classify_source, character(1)),
        sensor = vapply(feature, classify_sensor, character(1)),
        ema_item = vapply(feature, classify_ema_item, character(1)),
        window = vapply(feature, classify_window, character(1)),
        index = vapply(feature, classify_index, character(1))
      ) %>%
      arrange(desc(median_importance), desc(mean_importance)),
    source = group_summary("source"),
    sensor = group_summary("sensor"),
    ema_item = group_summary("ema_item"),
    window = group_summary("modal_window"),
    index = group_summary("index")
  )

  for (nm in names(summaries)) {
    write.csv(
      summaries[[nm]],
      file.path(OUT_DIR, sprintf("summary_%s_%s.csv", nm, out_stub)),
      row.names = FALSE
    )
  }

  winners <- bind_rows(lapply(names(summaries), function(nm) {
    x <- summaries[[nm]]
    if (!nrow(x)) return(NULL)
    data.frame(
      model_type = model_type,
      config = config,
      criterion = criterion_label,
      domain = nm,
      winner = if ("feature" %in% names(x)) x$feature[1] else x$level[1],
      winner_metric = if ("feature" %in% names(x)) {
        "median_feature_importance"
      } else {
        "median_importance_per_feature"
      },
      median_importance = if ("feature" %in% names(x)) x$median_importance[1] else x$median_importance_per_feature[1],
      mean_importance = if ("feature" %in% names(x)) x$mean_importance[1] else x$mean_importance_per_feature[1],
      median_total_importance = if ("median_total_importance" %in% names(x)) x$median_total_importance[1] else NA_real_,
      median_n_features = if ("median_n_features" %in% names(x)) x$median_n_features[1] else NA_real_,
      n_participants = x$n_participants[1],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(winners, file.path(OUT_DIR, sprintf("winners_%s.csv", out_stub)), row.names = FALSE)

  invisible(winners)
}

data_model <- data.table::fread(MODEL_READY_FILE, colClasses = list(character = "uid")) %>%
  as.data.frame() %>%
  mutate(uid = stringr::str_pad(as.character(uid), width = 4, pad = "0"))

results <- read.csv(RESULTS_FILE, stringsAsFactors = FALSE) %>%
  mutate(Participant = stringr::str_pad(as.character(Participant), width = 4, pad = "0"))

analyses <- data.frame(
  model_type = c(rep("Continuous", 3), rep("Hurdle", 3)),
  config = rep(c("EMA_mean_dynamic", "Opt_Passive_mean_dynamic", "EMA_Opt_Passive_mean_dynamic"), 2),
  metric = c(rep("r2", 3), rep("auc", 3)),
  threshold = c(rep(0, 3), rep(0.5, 3)),
  criterion = c(rep("R² > 0", 3), rep("AUC > 0.5", 3)),
  stringsAsFactors = FALSE
)

status_rows <- list()
winner_rows <- list()
predictor_audit_rows <- list()

for (ii in seq_len(nrow(analyses))) {
  a <- analyses[ii, ]
  out_stub <- sprintf(
    "%s_%s_%s",
    tolower(a$model_type),
    a$config,
    ifelse(a$model_type == "Continuous", "R2gt0", "AUCgt05")
  )
  cat(sprintf("\n[%s] %s | %s\n", out_stub, a$model_type, a$criterion))

  eligible <- results %>%
    filter(
      model_type == a$model_type,
      config == a$config,
      status == "ok",
      !is.na(.data[[a$metric]]),
      .data[[a$metric]] > a$threshold
    ) %>%
    select(Participant, config, r2, cor, mse, mae, auc, sensitivity, specificity,
           balanced_acc, n_obs, n_SI, n_nonSI, n_folds) %>%
    arrange(desc(.data[[a$metric]]))

  write.csv(
    eligible,
    file.path(OUT_DIR, sprintf("participants_%s.csv", out_stub)),
    row.names = FALSE
  )

  predictors <- get_config_predictors(a$config, names(data_model))
  predictor_audit_rows[[out_stub]] <- data.frame(
    model_type = a$model_type,
    config = a$config,
    feature = predictors,
    source = vapply(predictors, classify_source, character(1)),
    sensor = vapply(predictors, classify_sensor, character(1)),
    ema_item = vapply(predictors, classify_ema_item, character(1)),
    window = vapply(predictors, classify_window, character(1)),
    index = vapply(predictors, classify_index, character(1)),
    stringsAsFactors = FALSE
  )
  importance_rows <- list()
  for (u in eligible$Participant) {
    pt <- data_model[data_model$uid == u, , drop = FALSE]
    imp <- fit_participant_importance(pt, predictors, a$model_type, a$config)
    if (is.null(imp) || !nrow(imp)) next
    imp$Participant <- u
    importance_rows[[u]] <- imp
  }
  importance <- bind_rows(importance_rows)
  if (!nrow(importance)) {
    importance <- data.frame(
      Participant = character(),
      feature = character(),
      importance = numeric(),
      importance_positive = numeric(),
      importance_column = character(),
      modal_window = character(),
      source = character(),
      sensor = character(),
      ema_item = character(),
      index = character(),
      model_type = character(),
      config = character(),
      importance_method = character()
    )
  }

  importance_file <- file.path(OUT_DIR, sprintf("importance_%s.csv", out_stub))
  write.csv(importance, importance_file, row.names = FALSE)

  winners <- summarize_importance(importance, a$model_type, a$config, a$criterion, out_stub)
  if (!is.null(winners)) winner_rows[[out_stub]] <- winners

  status_rows[[out_stub]] <- data.frame(
    model_type = a$model_type,
    config = a$config,
    criterion = a$criterion,
    eligible_participants = nrow(eligible),
    importance_rows = nrow(importance),
    importance_participants = dplyr::n_distinct(importance$Participant),
    predictor_count_before_nzv = length(predictors),
    importance_file = importance_file,
    stringsAsFactors = FALSE
  )
  cat(sprintf(
    "  eligible=%d | importance rows=%d | importance PTs=%d | predictors=%d\n",
    nrow(eligible), nrow(importance), dplyr::n_distinct(importance$Participant),
    length(predictors)
  ))
}

status <- bind_rows(status_rows)
winners_all <- bind_rows(winner_rows)
predictor_audit <- bind_rows(predictor_audit_rows)
predictor_audit_summary <- predictor_audit %>%
  group_by(model_type, config, source, sensor, ema_item, window, index) %>%
  summarise(n_predictors = n_distinct(feature), .groups = "drop") %>%
  arrange(model_type, config, source, desc(n_predictors))

write.csv(status, STATUS_FILE, row.names = FALSE)
write.csv(winners_all, file.path(OUT_DIR, "exploratory_feature_importance_winners.csv"), row.names = FALSE)
write.csv(predictor_audit, file.path(OUT_DIR, "exploratory_predictor_audit.csv"), row.names = FALSE)
write.csv(predictor_audit_summary, file.path(OUT_DIR, "exploratory_predictor_audit_summary.csv"), row.names = FALSE)
failure_log <- dplyr::bind_rows(importance_failures)
if (!nrow(failure_log)) {
  failure_log <- data.frame(
    Participant = character(),
    model_type = character(),
    config = character(),
    reason = character()
  )
}
write.csv(failure_log, file.path(OUT_DIR, "exploratory_feature_importance_failures.csv"), row.names = FALSE)

cat(sprintf("\nSaved exploratory feature-importance status: %s\n", STATUS_FILE))
cat(sprintf("Saved exploratory feature-importance outputs to: %s\n", OUT_DIR))
if (nrow(failure_log)) {
  cat("\nFeature-importance failure reasons:\n")
  failure_summary <- failure_log %>% count(model_type, config, reason, sort = TRUE)
  utils::write.table(
    utils::head(as.data.frame(failure_summary), 30),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    na = ""
  )
}
cat("\nFeature-importance run status:\n")
utils::write.table(
  as.data.frame(status),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  na = ""
)
}

if (RUN_FINAL_MANUSCRIPT_OUTPUTS) {
cat("\n===== SECTION: FINAL MANUSCRIPT OUTPUTS =====\n")
################################################################################
# Final manuscript outputs only
#
# This is the clean runner for the final project.
# It produces the outputs needed for the Results section:
#
#   1. Participant flow / Table 1 demographic outputs
#   2. Continuous-model results table
#   3. Hurdle-model results table
#   4. Hurdle model within-person performance figure
#   5. Participant 0300 hurdle time-series figure, if enabled/generated
#   6. Feature-importance summaries using mean absolute approximate SHAP values
#
# Assumes the model-ready CSV and model result CSV already exist unless the
# main pipeline flags are changed manually.
################################################################################

message("Running final manuscript output pipeline...")


if (exists("make_F2_hurdle") && exists("results_main")) {
  message("Creating final hurdle within-person performance figure...")
  make_F2_hurdle(results_main, "main", CV_VARIANTS$main)
}

if (exists("make_hurdle_timeseries_for_participant")) {
  message("Creating final Participant 0300 hurdle time-series figure...")
  make_hurdle_timeseries_for_participant("0300", CV_VARIANTS$main, report_all_pt_hurdle_dir)
}

if (exists("make_continuous_timeseries_for_participant")) {
  message("Creating final Participant 0300 continuous time-series figure...")
  make_continuous_timeseries_for_participant("0300", CV_VARIANTS$main, report_all_pt_continuous_dir)
}

message("Building final SHAP feature-importance table from existing feature-importance outputs...")


BASE_PATH <- resolve_project_base()
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
RESULTS_DIR <- file.path(RESULTS_ROOT_DIR, "report")
REPORT_DIR <- RESULTS_DIR
PLOTS_DIR <- RESULTS_DIR
FEATURE_DIR <- if (exists("FEATURE_IMPORTANCE_WORK_DIR", inherits = FALSE)) {
  FEATURE_IMPORTANCE_WORK_DIR
} else {
  file.path(RESULTS_DIR, "feature_selection_analysis", "exploratory_successful_models")
}

table4_path <- file.path(REPORT_DIR, "Table4_feature_importance_shap.csv")
winners_path <- file.path(FEATURE_DIR, "exploratory_feature_importance_winners.csv")
status_path <- file.path(FEATURE_DIR, "exploratory_feature_importance_status.csv")
if (file.exists(winners_path)) {
  winners <- read.csv(winners_path, stringsAsFactors = FALSE)
  status <- if (file.exists(status_path)) read.csv(status_path, stringsAsFactors = FALSE) else data.frame()
  feature_value <- function(model_type, config, domain, value_col = "winner") {
    x <- winners[
      winners$model_type == model_type &
        winners$config == config &
        winners$domain == domain,
      ,
      drop = FALSE
    ]
    if (!nrow(x)) return(NA_character_)
    if (value_col == "winner") return(as.character(x$winner[1]))
    sprintf("%0.6f", x[[value_col]][1])
  }
  status_value <- function(model_type, config, value_col) {
    if (!nrow(status)) return(NA_integer_)
    x <- status[
      status$model_type == model_type &
        status$config == config,
      ,
      drop = FALSE
    ]
    if (!nrow(x) || !value_col %in% names(x)) return(NA_integer_)
    x[[value_col]][1]
  }
  table4_specs <- data.frame(
    Model = c("Continuous", "Continuous", "Continuous", "Hurdle", "Hurdle", "Hurdle"),
    Configuration = c(
      "EMA (dynamic)",
      "Passive (dynamic)",
      "EMA + Passive (dynamic)",
      "EMA (dynamic)",
      "Passive (dynamic)",
      "EMA + Passive (dynamic)"
    ),
    config_id = c(
      "EMA_mean_dynamic",
      "Opt_Passive_mean_dynamic",
      "EMA_Opt_Passive_mean_dynamic",
      "EMA_mean_dynamic",
      "Opt_Passive_mean_dynamic",
      "EMA_Opt_Passive_mean_dynamic"
    ),
    stringsAsFactors = FALSE
  )
  table4 <- table4_specs
  table4$Criterion <- ifelse(table4$Model == "Continuous", "R² > 0", "AUC > 0.5")
  table4$`Eligible N` <- mapply(status_value, table4$Model, table4$config_id, MoreArgs = list(value_col = "eligible_participants"))
  table4$`SHAP N` <- mapply(status_value, table4$Model, table4$config_id, MoreArgs = list(value_col = "importance_participants"))
  table4$`Top individual predictor` <- mapply(feature_value, table4$Model, table4$config_id, "feature")
  table4$`Feature median mean absolute approximate SHAP` <- mapply(feature_value, table4$Model, table4$config_id, "feature", MoreArgs = list(value_col = "median_importance"))
  table4$`Top grouped source` <- mapply(feature_value, table4$Model, table4$config_id, "source")
  table4$`Top grouped EMA item` <- mapply(feature_value, table4$Model, table4$config_id, "ema_item")
  table4$`Top grouped passive sensor` <- mapply(feature_value, table4$Model, table4$config_id, "sensor")
  table4$`Top grouped time window` <- mapply(feature_value, table4$Model, table4$config_id, "window")
  table4$`Top grouped index` <- mapply(feature_value, table4$Model, table4$config_id, "index")
  table4$config_id <- NULL
  write.csv(table4, table4_path, row.names = FALSE)
} else if (file.exists(table4_path)) {
  message(sprintf("Using existing SHAP Table 4: %s", table4_path))
} else {
  warning(sprintf("Missing SHAP winner file, so Table 4 was not created: %s", winners_path))
}

}

if (RUN_ALL_PT_CONTINUOUS_TIMESERIES || RUN_ALL_PT_HURDLE_TIMESERIES) {
cat("\n===== OPTIONAL SECTION: ALL-PT TIME-SERIES FIGURES =====\n")
################################################################################
# Generate time-series plots for PT0300 and all eligible participants.
#
# This section uses cached model results from the main prediction section above.
# It does not rerun the full continuous/hurdle model pipeline unless the cache
# rerun/cache flags above are changed.
################################################################################

if (RUN_ALL_PT_CONTINUOUS_TIMESERIES) {
  out_dir <- report_all_pt_continuous_dir
  make_dirs(out_dir)

  cat("\nCreating PT0300 continuous time-series figure...\n")
  make_continuous_timeseries_for_participant(
    uid = "0300",
    cv_params = CV_VARIANTS$main,
    out_dir = out_dir,
    overwrite = TRUE
  )

  cat("\nCreating continuous time-series figures for all remaining eligible participants...\n")
  make_continuous_timeseries_all_participants(
    combined = results_main,
    variant_name = "main",
    cv_params = CV_VARIANTS$main,
    participant_ids = NULL,
    overwrite = FALSE
  )

  cat(sprintf("\nContinuous figures are saved in:\n%s\n", out_dir))
}

if (RUN_ALL_PT_HURDLE_TIMESERIES) {
  out_dir <- report_all_pt_hurdle_dir
  make_dirs(out_dir)

  cat("\nCreating PT0300 hurdle time-series figure...\n")
  make_hurdle_timeseries_for_participant(
    uid = "0300",
    cv_params = CV_VARIANTS$main,
    out_dir = out_dir,
    overwrite = TRUE
  )

  cat("\nCreating hurdle time-series figures for all remaining eligible participants...\n")
  make_hurdle_timeseries_all_participants(
    combined = results_main,
    variant_name = "main",
    cv_params = CV_VARIANTS$main,
    participant_ids = NULL,
    overwrite = FALSE
  )

  cat(sprintf("\nHurdle figures are saved in:\n%s\n", out_dir))
}
}

if (RUN_FEATURE_IMPORTANCE || RUN_FINAL_MANUSCRIPT_OUTPUTS) {
cat("\n===== SECTION: FEATURE IMPORTANCE ONE-CSV EXPORT =====\n")
################################################################################
# Combine feature-importance outputs into one aggregate CSV
#
# This creates one clean, manuscript-friendly CSV from the feature-importance
# output folder. It includes aggregate winners, status, failure summaries, and
# top domain summaries. It intentionally does NOT merge participant-level
# importance_*.csv rows into the public/report CSV because those contain
# participant identifiers.
################################################################################

BASE_PATH <- resolve_project_base()
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
REPORT_DIR <- file.path(RESULTS_ROOT_DIR, "report")
FEATURE_DIR <- file.path(
  REPORT_DIR,
  "feature_selection_analysis",
  "exploratory_successful_models"
)
if (exists("FEATURE_IMPORTANCE_WORK_DIR", inherits = FALSE)) {
  FEATURE_DIR <- FEATURE_IMPORTANCE_WORK_DIR
}
OUT_FILE <- file.path(REPORT_DIR, "feature_importance_all_in_one_aggregate.csv")

make_dirs(REPORT_DIR)

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

as_chr <- function(x) {
  if (is.null(x)) return(NA_character_)
  as.character(x)
}

num_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

rows <- list()

# 1. Compact winners by model/config/domain.
winners <- read_if_exists(file.path(FEATURE_DIR, "exploratory_feature_importance_winners.csv"))
if (!is.null(winners) && nrow(winners)) {
  rows[["winners"]] <- data.frame(
    section = "winner",
    model_type = winners$model_type,
    config = winners$config,
    criterion = winners$criterion,
    domain = winners$domain,
    item = winners$winner,
    metric = winners$winner_metric,
    value = num_or_na(winners$median_importance),
    mean_value = num_or_na(winners$mean_importance),
    n_participants = num_or_na(winners$n_participants),
    eligible_n = NA_real_,
    shap_n = NA_real_,
    note = "Top aggregate winner by domain",
    source_file = "exploratory_feature_importance_winners.csv",
    stringsAsFactors = FALSE
  )
}

# 2. Eligible N / SHAP N status.
status <- read_if_exists(file.path(FEATURE_DIR, "exploratory_feature_importance_status.csv"))
if (!is.null(status) && nrow(status)) {
  rows[["status"]] <- data.frame(
    section = "status",
    model_type = status$model_type,
    config = status$config,
    criterion = status$criterion,
    domain = "run_status",
    item = "eligible_and_shap_n",
    metric = "counts",
    value = num_or_na(status$importance_rows),
    mean_value = NA_real_,
    n_participants = num_or_na(status$importance_participants),
    eligible_n = num_or_na(status$eligible_participants),
    shap_n = num_or_na(status$importance_participants),
    note = paste0("Predictors before near-zero-variance filtering: ", status$predictor_count_before_nzv),
    source_file = "exploratory_feature_importance_status.csv",
    stringsAsFactors = FALSE
  )
}

# 3. Failure reasons summarized without participant IDs.
failures <- read_if_exists(file.path(FEATURE_DIR, "exploratory_feature_importance_failures.csv"))
if (!is.null(failures) && nrow(failures)) {
  failure_summary <- aggregate(
    Participant ~ model_type + config + reason,
    data = failures,
    FUN = function(x) length(unique(x))
  )
  names(failure_summary)[names(failure_summary) == "Participant"] <- "n_failed_participants"
  rows[["failure_summary"]] <- data.frame(
    section = "failure_summary",
    model_type = failure_summary$model_type,
    config = failure_summary$config,
    criterion = NA_character_,
    domain = "failure_reason",
    item = failure_summary$reason,
    metric = "n_failed_participants",
    value = num_or_na(failure_summary$n_failed_participants),
    mean_value = NA_real_,
    n_participants = num_or_na(failure_summary$n_failed_participants),
    eligible_n = NA_real_,
    shap_n = NA_real_,
    note = "Summarized without participant IDs",
    source_file = "exploratory_feature_importance_failures.csv",
    stringsAsFactors = FALSE
  )
}

# 4. Top rows from domain-level summary files.
summary_files <- list.files(
  FEATURE_DIR,
  pattern = "^summary_(feature|source|ema_item|sensor|window|index)_.*\\.csv$",
  full.names = TRUE
)
summary_rows <- list()
for (fp in summary_files) {
  x <- read_if_exists(fp)
  if (is.null(x) || !nrow(x)) next
  base <- basename(fp)
  domain <- sub("^summary_([^_]+(?:_[^_]+)?)_.*$", "\\1", base)
  if (grepl("^summary_ema_item_", base)) domain <- "ema_item"
  if (grepl("^summary_feature_", base)) domain <- "feature"
  if (grepl("^summary_source_", base)) domain <- "source"
  if (grepl("^summary_sensor_", base)) domain <- "sensor"
  if (grepl("^summary_window_", base)) domain <- "window"
  if (grepl("^summary_index_", base)) domain <- "index"

  # Keep top 10 per summary file to keep the all-in-one CSV readable.
  x <- x[seq_len(min(10, nrow(x))), , drop = FALSE]
  item <- if ("feature" %in% names(x)) x$feature else if ("level" %in% names(x)) x$level else NA_character_
  val <- if ("median_importance" %in% names(x)) {
    x$median_importance
  } else if ("median_importance_per_feature" %in% names(x)) {
    x$median_importance_per_feature
  } else {
    NA_real_
  }
  mean_val <- if ("mean_importance" %in% names(x)) {
    x$mean_importance
  } else if ("mean_importance_per_feature" %in% names(x)) {
    x$mean_importance_per_feature
  } else {
    NA_real_
  }
  summary_rows[[base]] <- data.frame(
    section = "summary_top10",
    model_type = if ("model_type" %in% names(x)) x$model_type else NA_character_,
    config = if ("config" %in% names(x)) x$config else NA_character_,
    criterion = NA_character_,
    domain = domain,
    item = item,
    metric = "median_importance",
    value = num_or_na(val),
    mean_value = num_or_na(mean_val),
    n_participants = if ("n_participants" %in% names(x)) num_or_na(x$n_participants) else NA_real_,
    eligible_n = NA_real_,
    shap_n = NA_real_,
    note = "Top 10 rows from detailed summary file",
    source_file = base,
    stringsAsFactors = FALSE
  )
}
if (length(summary_rows)) rows[["summary_top10"]] <- do.call(rbind, summary_rows)

all_in_one <- if (length(rows)) {
  do.call(rbind, rows)
} else {
  data.frame(
    section = character(),
    model_type = character(),
    config = character(),
    criterion = character(),
    domain = character(),
    item = character(),
    metric = character(),
    value = numeric(),
    mean_value = numeric(),
    n_participants = numeric(),
    eligible_n = numeric(),
    shap_n = numeric(),
    note = character(),
    source_file = character(),
    stringsAsFactors = FALSE
  )
}

write.csv(all_in_one, OUT_FILE, row.names = FALSE)
cat(sprintf("Saved all-in-one aggregate feature-importance CSV:\n%s\n", OUT_FILE))
cat(sprintf("Rows: %d\n", nrow(all_in_one)))

if (!isTRUE(KEEP_FEATURE_IMPORTANCE_DETAILS)) {
  cleanup_feature_dirs <- unique(c(
    FEATURE_DIR,
    file.path(REPORT_DIR, "feature_selection_analysis", "exploratory_successful_models"),
    file.path(REPORT_DIR, "feature_selection_analysis")
  ))
  cleanup_feature_dirs <- cleanup_feature_dirs[dir.exists(cleanup_feature_dirs)]
  if (length(cleanup_feature_dirs)) {
    unlink(cleanup_feature_dirs, recursive = TRUE, force = TRUE)
    cat(sprintf(
      "Removed detailed feature-importance working files. Set KEEP_FEATURE_IMPORTANCE_DETAILS <- TRUE before sourcing if you want to retain them.\n"
    ))
  }
}
}

if (RUN_BASE_RATE_ANALYSIS) {
cat("\n===== OPTIONAL SECTION: OUTCOME BASE-RATE / ZERO-INFLATION ANALYSIS =====\n")
################################################################################
# Quantify zero-inflation/restricted range in item 9 and test whether
# within-person predictability scales with each participant's item 9 base rate.
#
# Default input files are intentionally simple and editable:
#   raw_outcome.csv             one row per participant-day or EMA occasion
#   perperson_performance.csv   one row per participant x feature configuration
#
# Output:
#   base_rate_analysis.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

BASE_PATH <- resolve_project_base()
PROJECT_ROOTS <- discover_project_roots()
if (!length(PROJECT_ROOTS)) PROJECT_ROOTS <- BASE_PATH
cat(sprintf("Searching project roots:\n%s\n", paste(PROJECT_ROOTS, collapse = "\n")))
CLEAN_DATA_DIR <- file.path(BASE_PATH, "3. Cleaned data")
RESULTS_ROOT_DIR <- file.path(BASE_PATH, "4. Results")
REPORT_DIR <- file.path(RESULTS_ROOT_DIR, "report")
make_dirs(REPORT_DIR)

# ---- User-tunable inputs and thresholds ---------------------------------------
DEFAULT_RAW_OUTCOME_FILE <- file.path(
  BASE_PATH,
  "3. Cleaned data",
  "all_features_included_filtered.csv"
)
DEFAULT_PERFORMANCE_FILE <- file.path(
  BASE_PATH,
  "4. Results",
  "report",
  "All_models_results_main.csv"
)

find_first_existing_file <- function(candidates, search_roots = PROJECT_ROOTS, pattern = NULL) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[1])
  if (!is.null(pattern)) {
    found <- unlist(lapply(search_roots[dir.exists(search_roots)], function(root) {
      tryCatch(
        list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE),
        error = function(e) character()
      )
    }), use.names = FALSE)
    found <- found[file.exists(found)]
    if (length(found)) return(found[1])
  }
  NA_character_
}

outcome_candidates <- c(
  DEFAULT_RAW_OUTCOME_FILE,
  file.path(REPORT_DIR, "raw_outcome.csv"),
  file.path(CLEAN_DATA_DIR, "all_features_included_filtered.csv"),
  file.path(CLEAN_DATA_DIR, "all_features_included.csv"),
  file.path(BASE_PATH, "all_features_included_filtered.csv"),
  file.path(BASE_PATH, "all_features_included.csv")
)
performance_candidates <- c(
  DEFAULT_PERFORMANCE_FILE,
  file.path(REPORT_DIR, "perperson_performance.csv"),
  file.path(RESULTS_ROOT_DIR, "model_config_cache", "All_models_results_main.csv"),
  file.path(REPORT_DIR, "All_models_results_main.csv"),
  file.path(BASE_PATH, "All_models_results_main.csv")
)

RAW_OUTCOME_FILE <- if (exists("RAW_OUTCOME_FILE", inherits = FALSE)) {
  RAW_OUTCOME_FILE
} else {
  find_first_existing_file(
    outcome_candidates,
    search_roots = PROJECT_ROOTS,
    pattern = "^(raw_outcome|all_features_included_filtered|all_features_included)\\.csv$"
  )
}
PERFORMANCE_FILE <- if (exists("PERFORMANCE_FILE", inherits = FALSE)) {
  PERFORMANCE_FILE
} else {
  find_first_existing_file(
    performance_candidates,
    search_roots = PROJECT_ROOTS,
    pattern = "^(perperson_performance|All_models_results_main)\\.csv$"
  )
}

outcome_source <- "file"
performance_source <- "file"
raw_outcome <- NULL
performance <- NULL
if (exists("data_model_base", inherits = TRUE)) {
  raw_outcome <- as.data.frame(get("data_model_base", inherits = TRUE))
  outcome_source <- "existing object: data_model_base"
}
if (exists("results_main", inherits = TRUE)) {
  performance <- as.data.frame(get("results_main", inherits = TRUE))
  performance_source <- "existing object: results_main"
}

OUTCOME_ID_COL <- if (exists("OUTCOME_ID_COL", inherits = FALSE)) OUTCOME_ID_COL else "participant_id"
OUTCOME_DATE_COL <- if (exists("OUTCOME_DATE_COL", inherits = FALSE)) OUTCOME_DATE_COL else "date"
OUTCOME_ITEM9_COL <- if (exists("OUTCOME_ITEM9_COL", inherits = FALSE)) OUTCOME_ITEM9_COL else "item9"

PERF_ID_COL <- if (exists("PERF_ID_COL", inherits = FALSE)) PERF_ID_COL else "participant_id"
PERF_CONFIG_COL <- if (exists("PERF_CONFIG_COL", inherits = FALSE)) PERF_CONFIG_COL else "config"

NEAR_INVARIANT_MAX_IDEATION_DAYS <- if (exists("NEAR_INVARIANT_MAX_IDEATION_DAYS", inherits = FALSE)) {
  NEAR_INVARIANT_MAX_IDEATION_DAYS
} else {
  2
}
LOW_BASE_RATE_1 <- if (exists("LOW_BASE_RATE_1", inherits = FALSE)) LOW_BASE_RATE_1 else 0.05
LOW_BASE_RATE_2 <- if (exists("LOW_BASE_RATE_2", inherits = FALSE)) LOW_BASE_RATE_2 else 0.10
MIN_SI_POSITIVE_FOR_MODEL <- if (exists("MIN_SI_POSITIVE_FOR_MODEL", inherits = FALSE)) {
  MIN_SI_POSITIVE_FOR_MODEL
} else {
  10
}
MIN_NON_SI_FOR_MODEL <- if (exists("MIN_NON_SI_FOR_MODEL", inherits = FALSE)) {
  MIN_NON_SI_FOR_MODEL
} else {
  20
}

normalize_analysis_uid <- function(x) {
  x <- as.character(x)
  x <- sub("@.*$", "", x)
  digits <- stringr::str_extract(x, "\\d+")
  out <- ifelse(!is.na(digits), digits, x)
  out <- ifelse(grepl("^\\d+$", out) & nchar(out) <= 4,
                stringr::str_pad(out, 4, pad = "0"), out)
  out
}

pick_col <- function(data, requested, candidates, label) {
  if (!is.null(requested) && requested %in% names(data)) return(requested)
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) return(hit[1])
  stop(sprintf("Could not find %s column. Tried: %s",
               label, paste(unique(c(requested, candidates)), collapse = ", ")))
}

if (is.null(raw_outcome) && !file.exists(RAW_OUTCOME_FILE)) {
  stop(sprintf(
    "Outcome file not found by R and no data_model_base object is available.\nTried:\n%s\n\nThis add-on used to work when the main pipeline had loaded data_model_base. Either run the main pipeline once, set RAW_OUTCOME_FILE to the exact local path, or download/make the file available offline.",
    paste(outcome_candidates, collapse = "\n")
  ))
}
if (is.null(performance) && !file.exists(PERFORMANCE_FILE)) {
  stop(sprintf(
    "Performance file not found by R and no results_main object is available.\nTried:\n%s\n\nThis add-on used to work when the main pipeline had loaded results_main. Either run the main pipeline once, set PERFORMANCE_FILE to the exact local path, or download/make the file available offline.",
    paste(performance_candidates, collapse = "\n")
  ))
}

if (is.null(raw_outcome)) {
  raw_outcome <- read.csv(RAW_OUTCOME_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  outcome_source <- RAW_OUTCOME_FILE
}
if (is.null(performance)) {
  performance <- read.csv(PERFORMANCE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  performance_source <- PERFORMANCE_FILE
}
cat(sprintf("Using outcome source: %s\n", outcome_source))
cat(sprintf("Using performance source: %s\n", performance_source))

out_id_col <- pick_col(raw_outcome, OUTCOME_ID_COL,
                       c("participant_id", "Participant", "uid", "mlife_id", "id"),
                       "raw outcome participant ID")
out_date_col <- pick_col(raw_outcome, OUTCOME_DATE_COL,
                         c("date", "date_rep", "ema_day", "ema_date", "day"),
                         "raw outcome date")
out_item_col <- pick_col(raw_outcome, OUTCOME_ITEM9_COL,
                         c("item9", "q9", "phq9_item9", "PHQ9_item9"),
                         "item 9 outcome")

perf_id_col <- pick_col(performance, PERF_ID_COL,
                        c("participant_id", "Participant", "uid", "mlife_id", "id"),
                        "performance participant ID")
perf_config_col <- pick_col(performance, PERF_CONFIG_COL,
                            c("config", "Configuration", "feature_configuration", "model_config"),
                            "performance configuration")

daily_outcome <- raw_outcome %>%
  transmute(
    participant_id = normalize_analysis_uid(.data[[out_id_col]]),
    date = as.Date(.data[[out_date_col]]),
    item9 = suppressWarnings(as.numeric(.data[[out_item_col]]))
  ) %>%
  filter(!is.na(participant_id), !is.na(date), !is.na(item9)) %>%
  group_by(participant_id, date) %>%
  summarise(
    daily_occurrence = as.integer(any(item9 > 0, na.rm = TRUE)),
    daily_intensity = max(item9, na.rm = TRUE),
    .groups = "drop"
  )

perperson_baserate <- daily_outcome %>%
  group_by(participant_id) %>%
  summarise(
    n_days = n(),
    n_ideation_days = sum(daily_occurrence == 1, na.rm = TRUE),
    n_non_si_days = sum(daily_occurrence == 0, na.rm = TRUE),
    base_rate = n_ideation_days / n_days,
    sd_item9 = stats::sd(daily_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_item9 = ifelse(is.na(sd_item9), 0, sd_item9),
    near_invariant = sd_item9 == 0 | n_ideation_days <= NEAR_INVARIANT_MAX_IDEATION_DAYS,
    fail_low_si_positive = n_ideation_days < MIN_SI_POSITIVE_FOR_MODEL,
    fail_low_non_si = n_non_si_days < MIN_NON_SI_FOR_MODEL,
    pass_si = !fail_low_si_positive & !fail_low_non_si,
    base_rate_below_5pct = base_rate < LOW_BASE_RATE_1,
    base_rate_below_10pct = base_rate < LOW_BASE_RATE_2
  )

q_base <- stats::quantile(perperson_baserate$base_rate, probs = c(.25, .50, .75), na.rm = TRUE)
base_min <- min(perperson_baserate$base_rate, na.rm = TRUE)
base_max <- max(perperson_baserate$base_rate, na.rm = TRUE)
pct_near_invariant <- mean(perperson_baserate$near_invariant, na.rm = TRUE) * 100
pct_below_5 <- mean(perperson_baserate$base_rate < LOW_BASE_RATE_1, na.rm = TRUE) * 100
pct_below_10 <- mean(perperson_baserate$base_rate < LOW_BASE_RATE_2, na.rm = TRUE) * 100
n_zero_ideation <- sum(perperson_baserate$n_ideation_days == 0, na.rm = TRUE)
n_low_si_positive <- sum(perperson_baserate$fail_low_si_positive, na.rm = TRUE)
n_low_si_positive_only <- sum(perperson_baserate$fail_low_si_positive & !perperson_baserate$fail_low_non_si, na.rm = TRUE)
n_low_non_si_only <- sum(!perperson_baserate$fail_low_si_positive & perperson_baserate$fail_low_non_si, na.rm = TRUE)
n_low_both <- sum(perperson_baserate$fail_low_si_positive & perperson_baserate$fail_low_non_si, na.rm = TRUE)
keystone_baserate_summary <- data.frame(
  unit = "participant_day_after_daily_aggregation",
  denominator_participants = nrow(perperson_baserate),
  min_si_positive_required = MIN_SI_POSITIVE_FOR_MODEL,
  min_non_si_required = MIN_NON_SI_FOR_MODEL,
  excluded_due_to_low_si_positive = n_low_si_positive,
  excluded_due_to_low_si_positive_only = n_low_si_positive_only,
  excluded_due_to_low_non_si_only = n_low_non_si_only,
  excluded_due_to_both_low_si_and_low_non_si = n_low_both,
  stringsAsFactors = FALSE
)

cat("\nOutcome base-rate summary:\n")
cat(sprintf("  Participants: %d\n", nrow(perperson_baserate)))
cat(sprintf("  Median base rate: %.3f [IQR %.3f, %.3f]\n", q_base[[2]], q_base[[1]], q_base[[3]]))
cat(sprintf("  Base-rate range: %.3f to %.3f\n", base_min, base_max))
cat(sprintf("  Near-invariant participants: %.1f%%\n", pct_near_invariant))
cat(sprintf("  Base rate < %.0f%%: %.1f%%\n", LOW_BASE_RATE_1 * 100, pct_below_5))
cat(sprintf("  Base rate < %.0f%%: %.1f%%\n", LOW_BASE_RATE_2 * 100, pct_below_10))
cat(sprintf("  Participants with zero ideation days: %d\n", n_zero_ideation))
cat(sprintf("  KEYSTONE: participants with fewer than %d SI-positive days/observations: %d / %d\n",
            MIN_SI_POSITIVE_FOR_MODEL, n_low_si_positive, nrow(perperson_baserate)))
cat(sprintf("  KEYSTONE only low SI-positive, not low non-SI: %d\n", n_low_si_positive_only))

metric_candidates <- c("corr", "cor", "r", "r2", "auc", "balanced_accuracy", "balanced_acc")
metric_cols <- intersect(metric_candidates, names(performance))
metric_map <- c(corr = "corr", cor = "corr", r = "corr",
                r2 = "r2", auc = "auc",
                balanced_accuracy = "balanced_accuracy",
                balanced_acc = "balanced_accuracy")

performance_clean <- performance %>%
  mutate(
    participant_id = normalize_analysis_uid(.data[[perf_id_col]]),
    config = as.character(.data[[perf_config_col]])
  )
for (mc in metric_cols) {
  performance_clean[[metric_map[[mc]]]] <- suppressWarnings(as.numeric(performance_clean[[mc]]))
}
metric_cols_standard <- intersect(c("corr", "r2", "auc", "balanced_accuracy"), names(performance_clean))

merged_perf <- performance_clean %>%
  select(participant_id, config, all_of(metric_cols_standard)) %>%
  left_join(perperson_baserate, by = "participant_id")

zero_ideation_with_corr_missing <- 0L
if ("corr" %in% names(merged_perf)) {
  zero_ideation_with_corr_missing <- merged_perf %>%
    filter(n_ideation_days == 0, !is.finite(corr)) %>%
    summarise(n = n_distinct(participant_id), .groups = "drop") %>%
    pull(n)
  if (!length(zero_ideation_with_corr_missing)) zero_ideation_with_corr_missing <- 0L
}

analysis_long <- merged_perf %>%
  pivot_longer(all_of(metric_cols_standard), names_to = "metric", values_to = "metric_value") %>%
  filter(!is.na(base_rate), is.finite(metric_value))

run_assoc <- function(df, predictor) {
  df <- df %>%
    filter(is.finite(.data[[predictor]]), is.finite(metric_value), is.finite(n_days))
  if (nrow(df) < 5 || length(unique(df[[predictor]])) < 2 || length(unique(df$metric_value)) < 2) {
    return(data.frame(
      spearman_rho = NA_real_,
      p = NA_real_,
      n = nrow(df),
      ols_beta = NA_real_,
      ols_p = NA_real_
    ))
  }
  sp <- suppressWarnings(stats::cor.test(df[[predictor]], df$metric_value,
                                         method = "spearman", exact = FALSE))
  df$.predictor <- df[[predictor]]
  fit <- tryCatch(
    stats::lm(metric_value ~ .predictor + n_days, data = df),
    error = function(e) NULL
  )
  beta <- NA_real_
  pval <- NA_real_
  if (!is.null(fit)) {
    sm <- summary(fit)$coefficients
    pred_row <- grep("^\\.predictor$", rownames(sm))
    if (length(pred_row)) {
      beta <- sm[pred_row[1], "Estimate"]
      pval <- sm[pred_row[1], "Pr(>|t|)"]
    }
  }
  data.frame(
    spearman_rho = unname(sp$estimate),
    p = sp$p.value,
    n = nrow(df),
    ols_beta = beta,
    ols_p = pval
  )
}

assoc_results <- bind_rows(lapply(split(analysis_long, list(analysis_long$config, analysis_long$metric), drop = TRUE), function(df) {
  key <- df[1, c("config", "metric"), drop = FALSE]
  bind_rows(
    cbind(key, predictor = "base_rate", run_assoc(df, "base_rate")),
    cbind(key, predictor = "sd_item9", run_assoc(df, "sd_item9"))
  )
})) %>%
  mutate(
    ols_beta_baserate = ifelse(predictor == "base_rate", ols_beta, NA_real_),
    ols_beta_sd_item9 = ifelse(predictor == "sd_item9", ols_beta, NA_real_)
  ) %>%
  arrange(config, metric, predictor)

median_perf <- analysis_long %>%
  group_by(config, metric) %>%
  summarise(
    n_all = n_distinct(participant_id),
    median_all = median(metric_value, na.rm = TRUE),
    n_not_near_invariant = n_distinct(participant_id[!near_invariant]),
    median_not_near_invariant = median(metric_value[!near_invariant], na.rm = TRUE),
    n_base_rate_ge_10pct = n_distinct(participant_id[base_rate >= LOW_BASE_RATE_2]),
    median_base_rate_ge_10pct = median(metric_value[base_rate >= LOW_BASE_RATE_2], na.rm = TRUE),
    rises_after_dropping_near_invariant = median_not_near_invariant > median_all,
    rises_after_base_rate_ge_10pct = median_base_rate_ge_10pct > median_all,
    .groups = "drop"
  )

cat(sprintf("  Zero-ideation participants with undefined corr in performance file: %d\n",
            zero_ideation_with_corr_missing))
cat("\nSensitivity check: median performance after dropping near-invariant / base_rate >= 10%:\n")
print(median_perf, n = Inf)

cat("\nHeadline numbers for abstract/results:\n")
cat(sprintf("  Median base rate = %.1f%% [IQR %.1f%%, %.1f%%]\n",
            100 * q_base[[2]], 100 * q_base[[1]], 100 * q_base[[3]]))
cat(sprintf("  Near-invariant = %.1f%%; base rate < %.0f%% = %.1f%%; base rate < %.0f%% = %.1f%%\n",
            pct_near_invariant, LOW_BASE_RATE_1 * 100, pct_below_5,
            LOW_BASE_RATE_2 * 100, pct_below_10))

base_rate_output_path <- file.path(REPORT_DIR, "base_rate_analysis.csv")

empty_base_rate_row <- function(n = 0L) {
  data.frame(
    section = rep(NA_character_, n),
    participant_id = rep(NA_character_, n),
    config = rep(NA_character_, n),
    metric = rep(NA_character_, n),
    predictor = rep(NA_character_, n),
    value = rep(NA_real_, n),
    p = rep(NA_real_, n),
    n = rep(NA_real_, n),
    n_days = rep(NA_real_, n),
    n_si_pos = rep(NA_real_, n),
    n_non_si = rep(NA_real_, n),
    base_rate = rep(NA_real_, n),
    sd_item9 = rep(NA_real_, n),
    near_invariant = rep(NA, n),
    fail_low_si_positive = rep(NA, n),
    fail_low_non_si = rep(NA, n),
    pass_si = rep(NA, n),
    spearman_rho = rep(NA_real_, n),
    ols_beta = rep(NA_real_, n),
    ols_p = rep(NA_real_, n),
    median_all = rep(NA_real_, n),
    median_not_near_invariant = rep(NA_real_, n),
    median_base_rate_ge_10pct = rep(NA_real_, n),
    note = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

participant_rows <- empty_base_rate_row(nrow(perperson_baserate))
participant_rows$section <- "participant_base_rate"
participant_rows$participant_id <- perperson_baserate$participant_id
participant_rows$value <- perperson_baserate$base_rate
participant_rows$n <- perperson_baserate$n_days
participant_rows$n_days <- perperson_baserate$n_days
participant_rows$n_si_pos <- perperson_baserate$n_ideation_days
participant_rows$n_non_si <- perperson_baserate$n_non_si_days
participant_rows$base_rate <- perperson_baserate$base_rate
participant_rows$sd_item9 <- perperson_baserate$sd_item9
participant_rows$near_invariant <- perperson_baserate$near_invariant
participant_rows$fail_low_si_positive <- perperson_baserate$fail_low_si_positive
participant_rows$fail_low_non_si <- perperson_baserate$fail_low_non_si
participant_rows$pass_si <- perperson_baserate$pass_si

summary_metrics <- c(
  participants = nrow(perperson_baserate),
  median_base_rate = unname(q_base[[2]]),
  base_rate_iqr_25 = unname(q_base[[1]]),
  base_rate_iqr_75 = unname(q_base[[3]]),
  min_base_rate = base_min,
  max_base_rate = base_max,
  pct_near_invariant = pct_near_invariant,
  pct_base_rate_below_5 = pct_below_5,
  pct_base_rate_below_10 = pct_below_10,
  zero_ideation_participants = n_zero_ideation,
  zero_ideation_with_corr_missing = zero_ideation_with_corr_missing
)
summary_rows <- empty_base_rate_row(length(summary_metrics))
summary_rows$section <- "overall_summary"
summary_rows$metric <- names(summary_metrics)
summary_rows$value <- as.numeric(summary_metrics)

exclusion_metrics <- unlist(keystone_baserate_summary[1, ], use.names = TRUE)
exclusion_rows <- empty_base_rate_row(length(exclusion_metrics))
exclusion_rows$section <- "si_count_exclusion_summary"
exclusion_rows$metric <- names(exclusion_metrics)
exclusion_rows$value <- suppressWarnings(as.numeric(exclusion_metrics))
exclusion_rows$note <- as.character(exclusion_metrics)

association_rows <- empty_base_rate_row(nrow(assoc_results))
association_rows$section <- "predictability_association"
association_rows$config <- assoc_results$config
association_rows$metric <- assoc_results$metric
association_rows$predictor <- assoc_results$predictor
association_rows$value <- assoc_results$spearman_rho
association_rows$p <- assoc_results$p
association_rows$n <- assoc_results$n
association_rows$spearman_rho <- assoc_results$spearman_rho
association_rows$ols_beta <- assoc_results$ols_beta
association_rows$ols_p <- assoc_results$ols_p

sensitivity_rows <- empty_base_rate_row(nrow(median_perf))
sensitivity_rows$section <- "sensitivity_medians"
sensitivity_rows$config <- median_perf$config
sensitivity_rows$metric <- median_perf$metric
sensitivity_rows$n <- median_perf$n_all
sensitivity_rows$value <- median_perf$median_all
sensitivity_rows$median_all <- median_perf$median_all
sensitivity_rows$median_not_near_invariant <- median_perf$median_not_near_invariant
sensitivity_rows$median_base_rate_ge_10pct <- median_perf$median_base_rate_ge_10pct
sensitivity_rows$note <- sprintf(
  "n_not_near_invariant=%s; n_base_rate_ge_10pct=%s; rises_after_dropping_near_invariant=%s; rises_after_base_rate_ge_10pct=%s",
  median_perf$n_not_near_invariant,
  median_perf$n_base_rate_ge_10pct,
  median_perf$rises_after_dropping_near_invariant,
  median_perf$rises_after_base_rate_ge_10pct
)

base_rate_analysis <- dplyr::bind_rows(
  participant_rows,
  summary_rows,
  exclusion_rows,
  association_rows,
  sensitivity_rows
)
write.csv(base_rate_analysis, base_rate_output_path, row.names = FALSE)
cat(sprintf("Saved consolidated base-rate analysis CSV: %s\n", base_rate_output_path))
}

cat("\nSingle-file workflow complete.\n")
