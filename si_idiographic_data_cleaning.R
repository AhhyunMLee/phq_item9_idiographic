################################################################################
# SI Idiographic Prediction Project - Data Cleaning Script
#
# Companion to si_idiographic_data_analysis.R, which consumes the outputs
# built here (ema_features.csv, per-sensor dynamic indices, and
# all_features_included.csv). Run this script first.
#
# Folder layout expected under BASE_PATH (the working directory set below):
#   1. Code
#   2. Raw data
#   3. Cleaned data
#   4. Results
#
# The script runs top to bottom when sourced:
#   1. shared dynamic-index helpers
#   2. EMA processing
#   3. passive sensing data
#   4. data inspection / visualization
################################################################################

# ------------------------------------------------------------------------------
# 1) Set the working directory to the SI_Prediction project folder.
# ------------------------------------------------------------------------------
setwd("/Users/ahhyun/Library/CloudStorage/GoogleDrive-ahhyun.lee@dartmouth.edu/My Drive/SI_Prediction")   # <-- EDIT: put the full path to your SI_Prediction folder here

# ------------------------------------------------------------------------------
# 2) Install any missing packages, then load them all.
# ------------------------------------------------------------------------------
packages <- c(
  "data.table", "stringr", "lubridate", "parallel", "pbmcapply", "geosphere",
  "jsonlite", "dplyr", "tidyr", "progress", "tidyverse", "TSEntropies",
  "bit64", "future.apply", "furrr", "RANN", "caret", "randomForest",
  "foreach", "doParallel", "effsize", "patchwork", "gbm", "pROC",
  "ggplot2", "ggridges", "scales"
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

################################################################################
# Make UIDs and dates consistent across data sources.
# normalize_uid(): every participant-ID variant -> canonical 4-digit "0011" form.
# parse_date_any(): every date format -> Date; warns if a format is unrecognized.
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

# Any common date format -> Date. NA only as a last resort (with a warning):
# tries text formats first, then numeric timestamps (Unix seconds, Excel days).
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
#
# Builds ema_features.csv. This runs before passive-sensing processing (below)
# because build_ema_lookup() there reads ema_features.csv to anchor each
# sensor's rolling windows to EMA prompt times. Processed per-sensor
# dynamic_indices_*.csv files do not exist yet at this point in the pipeline,
# so participant compliance here is derived from RAW sensor folders instead.
# ==============================================================================
EMA_FILE <- "3._TD_cohort_EMA_FINAL.csv"
ema_features_path <- file.path(CLEAN_DATA_DIR, "ema_features.csv")

# Scan RAW per-participant sensor folders (not yet-processed dynamic_indices_*.csv
# files, since passive-sensing processing hasn't run yet) for which UIDs have any
# raw passive data at all. Used only as a light pre-filter on EMA rows; the
# authoritative passive-compliance and SI-count criteria are re-applied later in
# apply_inclusion_criteria() during the main pipeline.
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
      if (!is.na(id_col) && !is.na(id_col)) {
        ids <- unique(data.table::fread(fp, select = id_col)[[1]])
        uid_hits <- c(uid_hits, normalize_uid(ids))
      }
    }
  }
  sort(unique(uid_hits[!is.na(uid_hits) & nzchar(uid_hits)]))
}


cat("\n===== Build EMA features =====\n")
  # 1. Derive compliant_uids from RAW sensor folders (passive-sensing processing
  #    has not produced per-sensor dynamic_indices_*.csv files yet at this point).
  compliant_uids <- derive_compliant_uids_from_raw()
  cat(sprintf("Compliant UIDs from raw sensor data: %d\n", length(compliant_uids)))
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
    # Keep only the 11 expected items — this drops stray raw columns such as
    # the character "q9." follow-up field swept in by starts_with("PHQ.").
    select(mlife_id, compliance, date_rep, ema_time, ema_timestamp, ema_day,
           all_of(c(paste0("q", 1:9), "g1", "g2"))) %>%
    # Compliance is a proportion; a few raw values slightly exceed 1 — cap at 1.
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
  # 3. Compute dynamic indices (uses shared helpers: mean_thr, mad_thr, p90p10_relpos)
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
# ==============================================================================
# END EMA PROCESSING
# ==============================================================================

# ==============================================================================
# DATA CLEANING, 3 OF 4: PASSIVE SENSING DATA
# ==============================================================================
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
# Always rebuild from scratch: clear any previous output first.
unlink(acc_final_path)
unlink(output_dir, recursive = TRUE)
make_dirs(output_dir)
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

# status tally (DONE / EMPTY / NOEMA / NORANGE / ERROR)
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

data.table::setDTthreads(0)   # restore default threading
################################################################################


################################################################################
# GENERIC PARALLEL SENSOR RUNNER  (hr, steps, battery, convo, gps, gps_home)
#
# Sensor processing is parameterized by sensor specification and writes
# per-participant output files (rebuilt from scratch on every run).
# Coverage is stored descriptively only. Window validity is governed by the
# minimum-observation thresholds for each dynamic index and by participant-level
# passive compliance criteria.
################################################################################

SEC_2020     <- 1577836800
WINDOWS      <- list(hr1 = 3600, hr4 = 14400, day1 = 86400, day7 = 604800)
uid_default <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  normalize_uid(b)
}
uid_convo   <- function(fp) stringr::str_extract(basename(fp), "(?<=convo_)[^@]+")


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


# --- GPS: LATITUDE/LONGITUDE/DAY -> distance traveled / 30-min epoch ----
#   distance = haversine between consecutive fixes (meters), summed per 30-min bin.
read_gps <- function(fp) {
  d <- tryCatch(data.table::fread(fp), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  data.table::setnames(d, toupper(names(d)))
  if (!all(c("LATITUDE", "LONGITUDE", "DAY") %in% names(d))) return(NULL)
  d[, day := as.POSIXct(as.character(DAY), tz = "UTC")]
  d[, lat := suppressWarnings(as.numeric(LATITUDE))]
  d[, lon := suppressWarnings(as.numeric(LONGITUDE))]
  d <- d[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(lat) & !is.na(lon) &
           lat >= -90 & lat <= 90 & lon >= -180 & lon <= 180]
  if (nrow(d) < 2) return(NULL)
  data.table::setkey(d, day); n <- nrow(d)
  dist <- c(0, geosphere::distHaversine(cbind(d$lon[-n], d$lat[-n]),
                                        cbind(d$lon[-1], d$lat[-1])))   # meters
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
  d <- d[!is.na(day) & as.numeric(day) >= SEC_2020 & !is.na(lat) & !is.na(lon)]
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
    BASE_PATH
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

# ==============================================================================
# run_full(): parallel scale-up for one sensor
# ==============================================================================
run_full <- function(spec, workers = 10) {
  final_path <- file.path(CLEAN_SENSOR_DIR, sprintf("dynamic_indices_%s.csv", spec$tag))
  # Always rebuild from scratch: clear any previous output first.
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
# SLEEP
#   run_sleep() joins each assessment to the previous night's total_sleep and
#   rem_sleep (minutes). Garmin's calendarDate is the wake date, so the record
#   dated the same day as the assessment IS the night before. Multiple records
#   for one night are averaged downstream (reshape_dynamic_indices).
#   Sleep time-scales are built downstream as well: the night-before value
#   becomes the day1/day7 mean features, and a 7-day rolling MAD of total/REM
#   sleep captures past-week variability. Coverage/entropy do not apply.
# ==============================================================================
NIGHT_OFFSET_DAYS <- 0L   # calendarDate = wake date (per Garmin documentation)

run_sleep <- function(folder = "sleep") {
  final_path <- file.path(CLEAN_SENSOR_DIR, "dynamic_indices_sleep.csv")
  # Always rebuild from scratch: clear any previous output first.
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
    # Average multiple records for the same night (e.g., re-synced Garmin
    # summaries) so each night contributes exactly one row to the join.
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

# --- process every sensor -------------------------------------------------------
run_full(SENSORS$hr)        # biggest after acc (~610k rows/uid) — a few min
run_full(SENSORS$steps)
run_full(SENSORS$battery)
run_full(SENSORS$convo)
run_full(SENSORS$gps)
run_full(SENSORS$gps_home)
run_full(SENSORS$sms)
run_full(SENSORS$calls)
run_sleep()

# ==============================================================================
# DATA CLEANING, 4 OF 4: DATA INSPECTION / VISUALIZATION
# ==============================================================================
make_dirs(results_dir)   # first results are written from here on
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
    # Sleep is a once-per-night value, not a rolling window like the other
    # sensors, so only day1 (the night before) and day7 (rolling week) are
    # meaningful here -- unlike the other passive sensors, hr1/hr4 don't apply.
    for (w in c("day1", "day7")) {
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
    apply_si_counts    = TRUE, min_si_pos = 10, min_si_neg = 20) {
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
      passive_sensors = passive_sensors
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

cat("\nData cleaning complete.\n")

