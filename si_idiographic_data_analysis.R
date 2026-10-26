################################################################################
# SI Idiographic Prediction Project - Data Analysis Script
#
# Continuous RF model and Hurdle GBM model, fit on the same participants and the same test days
#   Baseline                  (predict the training mean)
#   EMA                       (EMA mean features)
#   EMA (dynamic)             (EMA mean + MAD features)
#   Passive                   (all passive sensors pooled, mean features)
#   Passive (dynamic)         (all passive sensors pooled, mean + dynamic)
#   EMA + Passive             (EMA mean + pooled passive mean)
#   EMA + Passive (dynamic)   (all of the above)
#   Reference predictors      (persistence, previous prompt, previous-day and 7-day means, always zero)
#
# Feature-set selection is over TIME WINDOW only; sensors are pooled within a
# window. Window and mtry are chosen on the validation split of each fold
# (lowest validation MAE), then refit on train+val and evaluated once on the
# held-out test day.
#
# Input : all_features_cleaned.csv   (already the analytic sample: participant
#         inclusion criteria are applied in the data-cleaning script,
#         si_idiographic_data_cleaning_1.R, which is run first)
# Output: tables/  one CSV per manuscript table (Tables 1, 3-8; Supplementary Tables S2, S3, S6-S11), Figure 3 and
#         reporting_answers.txt. These are the files you need.
#         Working files: model_results_matched/ and posthoc/working_files/ (can be deleted after the run; set KEEP_WORKING_FILES),
#         checkpoints/ (keep: they make a rerun fast).
#
# STAGE 1  Matched model fitting (writes matched_predictions.csv and logs)
#   * Calendar-day folds: train 45 d -> validate 10 d -> test 1 d, sliding by 1 day, after a real
#     6-day run-in (OSF schedule: first training day = study day 7, test days 62-90). Fold ids are test-day offsets,
#     so they are comparable across models.
#   * Each framework can be fit with either learner (RF, GBM) through COMBOS. The reported analysis uses Continuous-RF and
#     Hurdle-GBM, the two models named in the OSF registration. Every model selects its feature window (and RF mtry)
#     with the SAME criterion: validation MAE of the point prediction.
#   * Hurdle never silently skips a fold (fallbacks: prevalence for p(SI); mean positive value for severity).
#   * Reference predictors: Persistence (last q9 >= 4 h before the prompt), previous prompt, q9 window means.
#     Persistence and the window means respect the 4 h gap, so they are the fair naive comparators. PrevPrompt may
#     use a response from LESS than 4 h earlier, so it is an optimistic reference, not a like-for-like forecast.
#   * Per-participant seeding; every skipped participant/fold is logged with a reason; resumable checkpoints.
#   * Optional single-window runs (RUN_FIXED_WINDOWS) for the adaptive-vs-fixed test of H4.
#
# STAGE 2  Post-hoc statistics (reads Stage 1 output; no refitting)
#   PART A  size of the ROC-direction bias in the OLD outputs (optional)
#   PART B  common-support person-level metrics (AUC with fixed direction, AP, MAE overall/zero/positive,
#           RMSE, R2, r, Brier), summaries, paired contrasts (H1, H2, vs baseline, vs persistence,
#           H3 with matched learner, learner effect, H4 adaptive vs fixed), per-person block-permutation
#           tests, window stability, variance partition
#   PART C  participant flow, included vs excluded, outcome descriptives, skipped participants/folds
#
# STAGE 3  Reporting items for the manuscript (reads the Stage 1 and 2 files; nothing is refit; every block runs on its own)
#   D1 participant accounting (who is in every analysis and why not)   D2 prompt times and the 4-hour gap
#   D3 per-person comparison with the training-mean baseline            D4 models vs the training-mean and previous-day-mean benchmarks (S7)
#   D5 correlation and balanced accuracy (S6)                            D6 bootstrap CIs for the Table 4 medians
#   D7 higher SI-positive cutoffs, predictions re-scored                 D8 higher-cutoff runs that were refit (SI_CUT and RUN_TAG)
#   D9 source for the feature table (S2)                                 D10 Figure 3
#
# STAGE 4  One CSV per manuscript table, written to 4. Results/tables (Tables 1, 3 to 8 and Supplementary Tables S2, S3, S6 to S11)
#   The tables folder holds the files you need. Everything else (posthoc/working_files, model_results_matched) is a working file;
#   set KEEP_WORKING_FILES <- FALSE to delete it afterwards (the checkpoints stay).
#
# Run order: si_data_cleaning.R -> this script. Set the switches below; stages are independent once the
# Stage 1 output exists. mclapply uses forking (macOS/Linux); set N_CORES <- 1 on Windows.
# Settings of the reported analysis (section 0): RUNIN_DAYS = 6, STUDY_DAYS = 90, 45/10/1-day folds, COMBOS = Continuous-RF +
# Hurdle-GBM, SI_CUT = 0, MIN_NONMISSING = 0.5, MAX_PREDICTORS = NA, RUN_FIXED_WINDOWS = TRUE, SEED = 12361488. Each run prints a
# settings fingerprint ("Settings fingerprint for checkpoints"); the reported run printed 788D0354.
################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(caret)
library(randomForest)
library(gbm)
library(parallel)
library(data.table)

# caret::preProcess(method = "knnImpute") needs RANN. Without it, fit_preprocess() would SILENTLY
# fall back to median imputation, which changes the method.
if (!requireNamespace("RANN", quietly = TRUE))
  stop("Package 'RANN' is required (knnImpute). install.packages('RANN') and rerun.")

# ------------------------------------------------------------------------------
# 1. SET THE WORKING DIRECTORY to the project folder (the one that holds the cleaned data and "4. Results")
# ------------------------------------------------------------------------------
setwd("/Users/ahhyun/Library/CloudStorage/GoogleDrive-ahhyun.lee-2.hse@dartmouth.edu/My Drive/SI_Prediction/3. Cleaned data")

# ------------------------------------------------------------------------------
# 0. SETTINGS (edit)
# ------------------------------------------------------------------------------
# The project folder is the working directory set above. (On a cluster, set the environment variable SI_PROJECT_DIR to override it.)
PROJECT_DIR <- Sys.getenv("SI_PROJECT_DIR", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR)
if (basename(PROJECT_DIR) %in% c("1. Code", "2. Raw data", "3. Cleaned data", "4. Results"))
  PROJECT_DIR <- dirname(PROJECT_DIR)           # must be the top-level project folder
# Files are looked for in the usual subfolder first, then directly in the project folder (as in the previous version of this script).
locate <- function(name, subdirs) { cands <- file.path(PROJECT_DIR, subdirs, name); hit <- cands[file.exists(cands)]; if (length(hit)) hit[1] else cands[1] }

# TEST MODE: set to e.g. 6 to run everything on the first 6 participants only (minutes, not hours).
# Test outputs go to separate *_debug folders, so they can never be mistaken for (or reused as) real results.
DEBUG_N_UIDS <- NA_integer_
# SENSITIVITY RUNS: give each variant its own tag so its outputs and checkpoints never mix with the main run,
# e.g. RUN_TAG <- "_cut10" (with SI_CUT <- 10), "_nm80" (MIN_NONMISSING <- 0.8) or "_top20" (MAX_PREDICTORS <- 20).
RUN_TAG <- ""

RUN_STAGE1        <- TRUE    # fit the models (slow; checkpointed, safe to stop and resume)
RUN_STAGE2        <- TRUE    # post-hoc statistics (needs Stage 1 output)
RUN_STAGE3        <- TRUE    # reporting items D1-D10 and Figure 3 (needs the Stage 1 and 2 output; a few minutes)
RUN_STAGE4        <- TRUE    # one CSV per manuscript table, written to 4. Results/tables (needs Stages 1 to 3)
KEEP_WORKING_FILES <- TRUE   # FALSE = delete the working files (posthoc/working_files and the Stage 1 CSVs) after the tables are written; the checkpoints stay, so a rerun is fast
RUN_PART_A        <- TRUE    # ROC-bias check on the old prediction files (skipped if they are absent)
RUN_FIXED_WINDOWS <- TRUE    # TRUE adds the single-window fits for H4 (reported in the manuscript; much slower; checkpointed, so it can also be run after the main fits)

DATA_FILE <- locate("all_features_cleaned.csv", c("3. Cleaned data", "."))
OLD_DIR   <- file.path(PROJECT_DIR, "4. Results", "model_results")           # old continuous/hurdle prediction CSVs
.sfx      <- paste0(if (is.na(DEBUG_N_UIDS)) "" else "_debug", RUN_TAG)
OUT_DIR   <- file.path(PROJECT_DIR, "4. Results", paste0("model_results_matched", .sfx))   # Stage 1 output
NEW_DIR   <- OUT_DIR
OUT       <- file.path(PROJECT_DIR, "4. Results", paste0("posthoc", .sfx), "working_files")   # Stages 2 and 3: working files (intermediate CSVs)
OUT_TABLES <- file.path(PROJECT_DIR, "4. Results", paste0("tables", .sfx))                   # Stage 4: the final outputs (one CSV per table, Figure 3)
INCL_FILE <- local({ cands <- file.path(PROJECT_DIR, "4. Results", "report", c("Table1_inclusion_participant_report.csv", "inclusion_participant_report.csv"))
if (any(file.exists(cands))) cands[file.exists(cands)][1] else cands[1] })
# Participant-characteristics variables (Part C5; Table 3 in the manuscript). Paths assume the files are in "2. Raw data"; edit them if they live elsewhere.
# IDs: MLife_ID in the demographics file -> uid. The C-SSRS file is keyed by src_subject_id, which is linked to uid
# through Participant_ID in the demographics file.
DEMO_FILE    <- locate("9._FinalDemographics.csv", c("2. Raw data", "."))
DEMO_UID_COL <- "MLife_ID"             # participant ID in the demographics file (becomes uid)
DEMO_PID_COL <- "Participant_ID"       # links the demographics file to the C-SSRS file
DEMO_VARS    <- c("age_years", "Gender", "Transgender", "Sexual_Orientation", "Race_For_Reporting", "Hispanic_Latino",
                  "Highest_Level_Education", "Employment", "Income")   # characteristics-table demographics (age_years is built from "Age at Study Start"). Empty = detect automatically
TABLE2_ONLY  <- FALSE                  # TRUE (with RUN_STAGE1 <- FALSE): rerun ONLY Part C5 (characteristics table), which takes seconds instead of the full Stage 2
ONE_ROW_PER_PERSON <- TRUE            # TRUE: the C-SSRS and PHQ-9 files have one row per participant, so no baseline-row selection is done
CSSRS_FILE   <- locate("CSSRS_modified.csv", c("2. Raw data", "."))
CSSRS_TIMEFRAME <- "screen"            # "screen" = the screening C-SSRS in your file (cssrs_add_1 to cssrs_add_7: ideation coded 0-3, behaviour lifetime and past 3 months).
# "3mo" / "ever" = the older item sets (si1p3 ..., si1l ...) if a file has them
TREATMENT_COL   <- "treatmentgroup"    # study treatment condition; looked for in the demographics, then PHQ-9, then C-SSRS file
CSSRS_ID_COL <- "src_subject_id"
CSSRS_VISIT_COL      <- NA_character_  # only if the file has several visits: the column that identifies the visit ...
CSSRS_BASELINE_VALUE <- NA_character_  # ... and the value that marks baseline. NA = use every row (history = ever yes)
CSSRS_IDEATION_COLS       <- character(0)  # empty = detect (si1l-si5l or cssrs_base_01-05). Otherwise 5 yes/no items, ordered 1 (wish to be dead) to 5 (plan and intent)
CSSRS_ATTEMPT_COLS        <- character(0)  # empty = detect (sb1l, cssrs_base_12, ...). Otherwise the yes/no item(s) for an actual attempt
CSSRS_OTHER_BEHAVIOR_COLS <- character(0)  # empty = detect (sb3l-sb5l or cssrs_base_13-15): interrupted / aborted / preparatory behaviour
# Baseline depression severity: PHQ-9 total. Give either the 9 item columns (items 1 to 9, in order) or a total column.
BASELINE_DEP_FILE       <- local({ c1 <- as.vector(outer(c("PHQ9_modified.csv", "PHQ9.csv"), c("2. Raw data", "."), function(f, d) file.path(PROJECT_DIR, d, f)))
if (any(file.exists(c1))) c1[file.exists(c1)][1] else c1[1] })   # the modified file if present
BASELINE_DEP_ID_COL     <- "src_subject_id"
BASELINE_DEP_VIA_PID    <- TRUE            # TRUE: linked to uid through Participant_ID, like the C-SSRS. FALSE: the ID column is MLife_ID
BASELINE_DEP_VISIT_COL  <- NA_character_   # optional, as for the C-SSRS
BASELINE_DEP_VISIT_VALUE <- NA_character_
BASELINE_PHQ_ITEM_COLS  <- character(0)   # empty = detect a complete 9-item set (0-3 scoring). Otherwise 9 columns, suicide item LAST
BASELINE_PHQ_TOTAL_COL  <- NA_character_
# treatmentgroup is the STUDY CONDITION (e.g. intervention vs usual care), not clinical treatment status; label it that way in the characteristics table.
CHK_DIR   <- file.path(OUT_DIR, "checkpoints")
for (d in c(OUT_DIR, CHK_DIR, OUT, OUT_TABLES)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# design
OUTCOME        <- "q9"
OUTCOME_MAX    <- 100
TRAIN_DAYS     <- 45
VAL_DAYS       <- 10
RUNIN_DAYS     <- 6       # OSF schedule: study days 1-6 are excluded, so the first training day is day 7
STUDY_DAYS     <- 90      # OSF: 90 days of data. Rows after study day 90 are dropped, so the last test day is day 90
# => train 7-51 / validate 52-61 / test 62, ..., train 35-79 / validate 80-89 / test 90 = 29 test days
MIN_TRAIN_ROWS <- 60      # a fold needs at least this many training rows
MIN_VAL_ROWS   <- 10
MIN_NONMISSING <- 0.50    # a predictor needs this share of non-missing values in the training fold (sensitivity: 0.8)
SI_CUT         <- 0       # SI-positive = item 9 > SI_CUT. Keep 0 for the main analysis. A higher cutoff is a sensitivity
# analysis of OCCURRENCE (AUC, sensitivity/specificity); MAE-type metrics assume SI_CUT = 0.
MAX_PREDICTORS <- NA_integer_   # NA = off. A number K keeps the K predictors most related to item 9 (|Spearman|) in the
# TRAINING rows of each fold (no leakage) before fitting: a feature-reduction sensitivity analysis.
MIN_SI_TRAIN   <- 10      # positives needed to fit the severity stage
MTRY_GRID      <- c(2, 5, 10)    # RF only; tuned on the validation split
NTREE          <- 500
GBM_CLS <- data.frame(n.trees = 100, interaction.depth = 3, shrinkage = 0.1, n.minobsinnode = 10)
GBM_REG <- data.frame(n.trees = 100, interaction.depth = 3, shrinkage = 0.1, n.minobsinnode = 5)
# First run = the two models named in the OSF registration: random forest, and a two-part hurdle model with gradient boosting.
# For the same-learner hurdle-vs-continuous contrasts (exploratory) add Continuous-GBM and Hurdle-RF later: checkpoints mean
# only the new combinations are fitted.
COMBOS  <- data.frame(framework = c("Continuous", "Hurdle"), learner = c("RF", "GBM"), stringsAsFactors = FALSE)
# COMBOS <- data.frame(framework = c("Continuous", "Continuous", "Hurdle", "Hurdle"), learner = c("RF", "GBM", "RF", "GBM"), stringsAsFactors = FALSE)
CONFIGS <- c("Baseline", "EMA", "EMA (dynamic)", "Passive", "Passive (dynamic)",
             "EMA + Passive", "EMA + Passive (dynamic)")

# computation / inference
SEED    <- 12361488L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # pin the random-number generator: the sample() default changed in R 3.6.0, so an unpinned seed can give different draws in another R version
set.seed(SEED)                                              # every stochastic step below also reseeds itself from SEED, so each part reproduces on its own
N_CORES <- 10        # your value from the previous version of this script
data.table::setDTthreads(1)
B_BOOT  <- 5000            # bootstrap resamples for contrast CIs
PRIMARY_MODELS <- c("Hurdle-GBM", "Continuous-RF")   # permutation tests run on these (+ references)

# OSF preregistration: H1 = each individualised model (EMA, passive, combined) beats the training-mean baseline;
# H2 = the combined model beats the EMA-only and passive-only models. Metrics: R2, Pearson r, MAE, RMSE. Main model: random forest.
# Everything else (AUC, hurdle, other learners, persistence, window and dynamic-vs-mean questions) is EXPLORATORY and is labelled so.
OSF_MAIN_MODEL  <- "Continuous-RF"
OSF_FEATURE_SET <- "dynamic"      # the OSF lists four metrics per person (mean, MAD, sample entropy, P90-P10), i.e. our "(dynamic)"
# configurations. Use "mean" if your team reads the registered models as mean-level only.


# ==============================================================================
# STAGE 1: MATCHED MODEL FITTING
# ==============================================================================
if (RUN_STAGE1) {
  # ------------------------------------------------------------------------------
  # S1.1 LOAD DATA
  # ------------------------------------------------------------------------------
  parse_date2 <- function(x) {
    x <- as.character(x)
    d <- suppressWarnings(as.Date(x))
    if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%m/%d/%Y"))
    d
  }
  parse_time2 <- function(x) {
    if (inherits(x, "POSIXct")) return(as.numeric(x))
    if (is.numeric(x)) return(as.numeric(x))
    t <- tryCatch(as.numeric(as.POSIXct(as.character(x), tz = "UTC")),
                  error = function(e) rep(NA_real_, length(x)))
    if (all(is.na(t)))
      t <- suppressWarnings(as.numeric(lubridate::ymd_hms(as.character(x), tz = "UTC")))
    t
  }
  
  dat <- data.table::fread(DATA_FILE, colClasses = list(character = "uid")) %>% as.data.frame()
  # The cleaned file also contains rows for MISSED prompts (q9 = NA and ema_day = NA). Drop them first;
  # only answered prompts have a response time, so the timestamp checks below apply to those rows only.
  n_all <- nrow(dat)
  dat <- dat[!is.na(dat[[OUTCOME]]), , drop = FALSE]
  cat(sprintf("Rows with an answered %s: %d of %d\n", OUTCOME, nrow(dat), n_all))
  dat$date_cal <- parse_date2(dat$date_rep)
  dat$ema_num  <- parse_time2(dat$ema_day)
  if (mean(is.na(dat$date_cal)) > 0.01) stop("Could not parse date_rep.")
  if (mean(is.na(dat$ema_num)) > 0.05) {
    cat("Example ema_day values:\n"); print(head(unique(as.character(dat$ema_day[is.na(dat$ema_num)])), 5))
    stop("Could not parse ema_day for ", round(100 * mean(is.na(dat$ema_num)), 1),
         "% of answered prompts (needed for the Persistence reference).")
  }
  # Guard against a silent date-only parse (every timestamp at midnight)
  if (mean(dat$ema_num %% 86400 == 0, na.rm = TRUE) > 0.5)
    stop("ema_day looks date-only after parsing (most values are exactly midnight); check its format.")
  
  dat <- dat %>%
    mutate(.ord = match(ema_time, c("morning", "afternoon", "evening"))) %>%
    arrange(uid, date_cal, .ord) %>%
    select(-.ord) %>%
    filter(!is.na(.data[[OUTCOME]])) %>%
    group_by(uid) %>%
    mutate(day_index = as.integer(date_cal - min(date_cal))) %>%
    ungroup() %>%
    filter(day_index < STUDY_DAYS) %>%                      # study day = day_index + 1; keep days 1-90 (OSF)
    as.data.frame()
  dat$row_id <- seq_len(nrow(dat))
  DAT  <- dat
  UIDS <- unique(DAT$uid)
  if (!is.na(DEBUG_N_UIDS)) { UIDS <- head(UIDS, DEBUG_N_UIDS); message("TEST MODE: ", length(UIDS), " participants") }
  
  # Fingerprint of everything that changes the fitted models. It is part of every checkpoint file name, so a checkpoint made with
  # different settings (e.g. another run-in length) or a different cleaned data file is NEVER reused silently.
  SIG <- local({
    txt <- paste(RUNIN_DAYS, STUDY_DAYS, TRAIN_DAYS, VAL_DAYS, MIN_TRAIN_ROWS, MIN_VAL_ROWS, MIN_NONMISSING, MIN_SI_TRAIN, SI_CUT,
                 MAX_PREDICTORS, NTREE, paste(MTRY_GRID, collapse = "-"), SEED, DEBUG_N_UIDS, nrow(DAT), file.size(DATA_FILE),
                 paste(unlist(GBM_CLS), collapse = ","), paste(unlist(GBM_REG), collapse = ","), sep = "|")
    h <- 0
    for (ch in utf8ToInt(txt)) h <- (h * 131 + ch) %% 2147483629          # polynomial hash; stays exact in double precision
    sprintf("%08X", as.integer(h))
  })
  cat("Settings fingerprint for checkpoints:", SIG, "\n")
  cat(sprintf("Loaded %d rows | %d participants\n", nrow(DAT), length(UIDS)))
  
  # ------------------------------------------------------------------------------
  # S1.2 FEATURE CLASSIFICATION AND CANDIDATE SETS (same logic as the original)
  # ------------------------------------------------------------------------------
  non_feature <- c("uid", "date_rep", "ema_time", "mlife_id", "compliance", "ema_timestamp",
                   "ema_day", paste0("q", 1:9), "g1", "g2",
                   "date_cal", "ema_num", "day_index", "row_id")
  all_cols <- setdiff(colnames(DAT), non_feature)
  all_cols <- all_cols[!grepl("(_coverage|_n_obs)$", all_cols)]
  EMA_WINDOWED_PATTERN   <- "^(q[1-9]|g[12])_(mean|mad)_d(1|7)$"
  SLEEP_PATTERN          <- "^sleep_(mean|mad)_(total_sleep|rem_sleep)_(day1|day7)$"
  PASSIVE_WINDOW_PATTERN <- "(hr1|hr4|day1|day7)"
  col_class <- tibble(col = all_cols) %>%
    mutate(
      is_ema_windowed = str_detect(col, regex(EMA_WINDOWED_PATTERN, ignore_case = TRUE)),
      is_sleep        = str_detect(col, regex(SLEEP_PATTERN, ignore_case = TRUE)),
      is_passive      = str_detect(col, regex(PASSIVE_WINDOW_PATTERN, ignore_case = TRUE)) &
        !is_ema_windowed & !is_sleep,
      source = case_when(is_ema_windowed ~ "ema", is_sleep ~ "passive",
                         is_passive ~ "passive", TRUE ~ "unclassified"),
      type = case_when(
        is_ema_windowed & str_detect(col, "_mean_d[17]$") ~ "mean",
        is_ema_windowed & str_detect(col, "_mad_d[17]$")  ~ "dynamic",
        is_sleep & str_detect(col, "^sleep_mean_") ~ "mean",
        is_sleep & str_detect(col, "^sleep_mad_")  ~ "dynamic",
        is_passive & str_detect(col, regex("(_mean_val|_sum_val)$", ignore_case = TRUE)) ~ "mean",
        is_passive & str_detect(col, regex("(_mad_val|_entropy|_p90p10)$", ignore_case = TRUE)) ~ "dynamic",
        TRUE ~ NA_character_),
      window = case_when(
        is_ema_windowed ~ str_replace(str_extract(col, "d[17]$"), "^d", "day"),
        is_sleep        ~ str_extract(col, "(day1|day7)$"),
        is_passive      ~ str_extract(col, regex(PASSIVE_WINDOW_PATTERN, ignore_case = TRUE)),
        TRUE ~ NA_character_)
    ) %>%
    filter(source != "unclassified", !is.na(type), !is.na(window))
  windows_standard <- col_class %>% filter(source == "passive") %>% pull(window) %>% unique()
  windows_ema      <- col_class %>% filter(source == "ema")     %>% pull(window) %>% unique()
  cat(sprintf("EMA windows: %s | passive windows: %s\n",
              paste(windows_ema, collapse = ", "), paste(windows_standard, collapse = ", ")))
  
  # ---- Feature inventory: what the data file offers and what is kept as a candidate predictor ------------------------
  # The number of observations (_n_obs) and the coverage (_coverage) of each sensor window are data-quality columns: they are removed
  # above and are NEVER predictors. This block does that exclusion explicitly, builds the final feature list and counts what is left,
  # so the Methods can quote exact numbers (printed on every run and saved to feature_inventory.csv; nothing is refit).
  stream_of  <- function(x) str_extract(x, "^(gps_home|gps|acc|hr|steps|battery|convo|sleep|sms|calls)(?=_)")
  in_file    <- setdiff(colnames(DAT), non_feature)                              # every column except identifiers and the outcome items
  qc_removed <- in_file[grepl("(_coverage|_n_obs)$", in_file)]                   # data-quality columns, excluded
  not_class  <- setdiff(all_cols, col_class$col)                                 # neither an EMA window feature nor a passive window feature
  all_feats  <- unique(col_class$col)                                            # the candidate predictors: EMA window features + passive stream features
  inv <- bind_rows(
    data.frame(group = "file", item = c("columns in the file (identifiers and outcome items excluded)", "removed: number-of-observations and coverage columns",
                                        "removed: columns that are neither EMA nor passive window features", "kept as candidate predictors"),
               n = c(length(in_file), length(qc_removed), length(not_class), length(all_feats)), stringsAsFactors = FALSE),
    col_class %>% count(source, type, name = "n") %>% transmute(group = "source and type", item = paste(source, type), n = n),
    col_class %>% filter(source == "passive") %>% mutate(stream = ifelse(is.na(stream_of(col)), "other", stream_of(col))) %>%
      count(stream, name = "n") %>% transmute(group = "passive stream", item = stream, n = n))
  # Columns the cleaning script dropped before the file was written (so they are not in the file and not counted above): this is what reconciles
  # the number of passive columns in the file with any larger number built earlier (more than 50% missing, or no variability).
  cnt_rows <- function(f) if (file.exists(f)) nrow(read.csv(f, stringsAsFactors = FALSE)) else NA_integer_
  rep_dir  <- file.path(PROJECT_DIR, "4. Results", "report")
  inv <- bind_rows(inv, data.frame(group = "dropped during cleaning (not in the file)", item = c("feature columns with more than 50% missing values", "feature columns with no variability"),
                                   n = c(cnt_rows(file.path(rep_dir, "dropped_high_missing_columns.csv")), cnt_rows(file.path(rep_dir, "dropped_no_variability_columns.csv"))), stringsAsFactors = FALSE))
  cat("\nFeature inventory:\n"); print(inv, row.names = FALSE)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE); write.csv(inv, file.path(OUT_DIR, "feature_inventory.csv"), row.names = FALSE)
  
  cols_for <- function(src, win, dyn) {
    types <- if (dyn) c("mean", "dynamic") else "mean"
    col_class %>% filter(source == src, window == win, type %in% types) %>% pull(col)
  }
  make_candidates <- function(config, force_window = NULL) {
    out <- list()
    add <- function(vars, ema_window, window) {
      vars <- unique(vars)
      if (length(vars) >= 2)
        out[[length(out) + 1L]] <<- list(vars = vars, ema_window = ema_window, window = window)
    }
    if (config == "Baseline") return(out)
    uses_ema     <- grepl("EMA", config, fixed = TRUE)
    uses_passive <- grepl("Passive", config, fixed = TRUE)
    dyn          <- grepl("dynamic", config, fixed = TRUE)
    if (uses_ema && !uses_passive) {
      for (we in windows_ema) add(cols_for("ema", we, dyn), we, NA_character_)
    } else if (!uses_ema && uses_passive) {
      for (w in windows_standard) add(cols_for("passive", w, dyn), NA_character_, w)
    } else {
      ema_block <- unlist(lapply(windows_ema, function(we) cols_for("ema", we, dyn)))
      for (w in windows_standard) add(c(ema_block, cols_for("passive", w, dyn)), "all", w)
    }
    if (!is.null(force_window)) {
      out <- Filter(function(cd) {
        w <- if (!is.na(cd$window)) cd$window else cd$ema_window
        identical(w, force_window)
      }, out)
    }
    out
  }
  
  # Exact size of each candidate pool (comment 19). Each model uses ONE candidate set per fold (the selected window);
  # n_pred in matched_predictions.csv records how many of these actually entered each fold's model.
  cand_counts <- dplyr::bind_rows(lapply(setdiff(CONFIGS, "Baseline"), function(cfg)
    dplyr::bind_rows(lapply(make_candidates(cfg), function(cd)
      data.frame(Configuration = cfg, window = if (!is.na(cd$window)) cd$window else cd$ema_window,
                 n_candidate_features = length(cd$vars), stringsAsFactors = FALSE)))))
  write.csv(cand_counts, file.path(OUT_DIR, "candidate_feature_counts.csv"), row.names = FALSE)
  cat("\nCandidate feature counts per configuration and window:\n"); print(cand_counts, row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # S1.3 CALENDAR-DAY FOLDS WITH A REAL RUN-IN
  # ------------------------------------------------------------------------------
  make_day_folds <- function(pt) {
    pt <- pt[pt$day_index >= RUNIN_DAYS, , drop = FALSE]
    if (!nrow(pt)) return(NULL)
    pt$dk <- pt$day_index - RUNIN_DAYS
    first_test <- TRAIN_DAYS + VAL_DAYS
    max_dk <- max(pt$dk)
    if (max_dk < first_test) return(NULL)
    folds <- list()
    for (t in first_test:max_dk) {
      te <- which(pt$dk == t)
      if (!length(te)) next
      tr <- which(pt$dk >= t - first_test & pt$dk < t - VAL_DAYS)
      vl <- which(pt$dk >= t - VAL_DAYS  & pt$dk < t)
      if (length(tr) < MIN_TRAIN_ROWS || length(vl) < MIN_VAL_ROWS) next
      folds[[length(folds) + 1L]] <- list(fold = t - first_test + 1L, test_day = t,
                                          tr = tr, vl = vl, te = te)
    }
    list(data = pt, folds = folds, max_possible = max_dk - first_test + 1L)
  }
  
  # ------------------------------------------------------------------------------
  # S1.4 PREPROCESSING HELPERS (same as the original)
  # ------------------------------------------------------------------------------
  usable_in <- function(cols, train_df, min_nonmissing = MIN_NONMISSING) {
    cols <- intersect(cols, names(train_df))
    if (!length(cols)) return(character(0))
    keep <- vapply(cols, function(cc) {
      x <- train_df[[cc]]
      is.numeric(x) && mean(!is.na(x)) >= min_nonmissing && sum(!is.na(x)) >= 2 &&
        stats::var(x, na.rm = TRUE) > 0
    }, logical(1))
    cols[keep]
  }
  drop_near_zero_variance <- function(cols, train_df, freqCut = 19, uniqueCut = 10) {
    cols <- intersect(cols, names(train_df))
    if (length(cols) < 1) return(cols)
    nzv <- tryCatch(caret::nearZeroVar(train_df[, cols, drop = FALSE],
                                       freqCut = freqCut, uniqueCut = uniqueCut),
                    error = function(e) integer(0))
    if (length(nzv)) cols <- cols[-nzv]
    cols
  }
  # NOTE: caret::preProcess(knnImpute) also CENTERS and SCALES (z-score using the training
  # fold). The old Methods text said min-max; describe what is actually done.
  fit_preprocess <- function(train_df, cols) {
    pp <- tryCatch(caret::preProcess(train_df[, cols, drop = FALSE], method = c("zv", "knnImpute")),
                   error = function(e) NULL)
    if (is.null(pp))
      pp <- tryCatch(caret::preProcess(train_df[, cols, drop = FALSE], method = c("zv", "medianImpute")),
                     error = function(e) NULL)
    pp
  }
  # Optional feature reduction (MAX_PREDICTORS): keep the k predictors most related to the outcome, judged on the
  # TRAINING rows of the fold only, so nothing from validation or test leaks into the choice.
  screen_predictors <- function(df, cols, y, k) {
    r <- vapply(cols, function(cc) abs(suppressWarnings(
      stats::cor(df[[cc]], y, method = "spearman", use = "pairwise.complete.obs"))), numeric(1))
    r[is.na(r)] <- 0
    cols[order(r, decreasing = TRUE)][seq_len(min(k, length(cols)))]
  }
  choose_youden_threshold <- function(true_bin, p_si, thresholds = seq(0.01, 0.99, by = 0.01)) {
    pos <- p_si[true_bin == 1]; neg <- p_si[true_bin == 0]
    if (!length(pos) || !length(neg)) return(0.5)
    tp <- colSums(outer(pos, thresholds, ">="))
    tn <- colSums(outer(neg, thresholds, "<"))
    j <- tp / length(pos) + tn / length(neg) - 1
    stats::median(thresholds[which(j == max(j, na.rm = TRUE))], na.rm = TRUE)
  }
  
  # ------------------------------------------------------------------------------
  # S1.5 LEARNERS (identical interface for both frameworks)
  # ------------------------------------------------------------------------------
  reg_fit_predict <- function(learner, x, y, x_ev, mt) {
    if (learner == "RF") {
      fit <- randomForest::randomForest(x = x, y = y, mtry = max(1, mt), ntree = NTREE)
      as.numeric(predict(fit, newdata = x_ev))
    } else {
      fit <- suppressWarnings(caret::train(
        x = as.data.frame(x), y = y, method = "gbm",
        trControl = caret::trainControl(method = "none"), tuneGrid = GBM_REG, verbose = FALSE))
      as.numeric(predict(fit, newdata = as.data.frame(x_ev)))
    }
  }
  cls_fit_predict <- function(learner, x, bin, x_ev, mt) {
    if (learner == "RF") {
      fit <- randomForest::randomForest(x = x, y = bin, mtry = max(1, mt), ntree = NTREE)
      as.numeric(predict(fit, newdata = x_ev, type = "prob")[, "SI"])
    } else {
      fit <- suppressWarnings(caret::train(
        x = as.data.frame(x), y = bin, method = "gbm",
        trControl = caret::trainControl(method = "none", classProbs = TRUE),
        tuneGrid = GBM_CLS, verbose = FALSE))
      as.numeric(predict(fit, newdata = as.data.frame(x_ev), type = "prob")[, "SI"])
    }
  }
  fit_predict <- function(framework, learner, x_tr, y_tr, x_ev, mtry = NA) {
    mt <- if (learner == "RF") min(mtry, ncol(x_tr)) else NA
    n_ev <- nrow(x_ev)
    if (framework == "Continuous") {
      p <- tryCatch(reg_fit_predict(learner, x_tr, y_tr, x_ev, mt), error = function(e) NULL)
      if (is.null(p)) return(NULL)
      return(list(p_si = rep(NA_real_, n_ev), sev = rep(NA_real_, n_ev),
                  pred = pmax(0, pmin(OUTCOME_MAX, p))))
    }
    bin <- factor(ifelse(y_tr > SI_CUT, "SI", "noSI"), levels = c("noSI", "SI"))
    if (length(unique(bin)) < 2) {
      p_si <- rep(mean(y_tr > SI_CUT), n_ev)                       # fallback: prevalence
    } else {
      p_si <- tryCatch(cls_fit_predict(learner, x_tr, bin, x_ev, mt), error = function(e) NULL)
      if (is.null(p_si)) return(NULL)
    }
    pos <- which(y_tr > SI_CUT)
    sev <- rep(0, n_ev)
    if (length(pos) >= MIN_SI_TRAIN && stats::var(y_tr[pos]) > 0) {
      s <- tryCatch(reg_fit_predict(learner, x_tr[pos, , drop = FALSE], y_tr[pos], x_ev, mt),
                    error = function(e) NULL)
      sev <- if (is.null(s)) rep(mean(y_tr[pos]), n_ev) else as.numeric(s)
    } else if (length(pos) > 0) {
      sev <- rep(mean(y_tr[pos]), n_ev)                       # fallback: mean positive value
    }
    list(p_si = p_si, sev = sev, pred = pmax(0, pmin(OUTCOME_MAX, p_si * sev)))
  }
  
  # ------------------------------------------------------------------------------
  # S1.6 ONE PARTICIPANT, ONE CONFIGURATION
  # ------------------------------------------------------------------------------
  run_participant <- function(uid_i, framework, learner, config, candidates) {
    tryCatch({
      set.seed(SEED + as.integer(uid_i))
      pt_all <- DAT[DAT$uid == uid_i, , drop = FALSE]
      if (is.na(stats::var(pt_all[[OUTCOME]])) || stats::var(pt_all[[OUTCOME]]) == 0)
        return(list(uid = uid_i, status = "zero_outcome_variance", preds = NULL, log = NULL))
      fo <- make_day_folds(pt_all)
      if (is.null(fo) || !length(fo$folds))
        return(list(uid = uid_i, status = "no_valid_folds", preds = NULL, log = NULL))
      pt <- fo$data
      usable <- names(pt)[vapply(pt, function(z) is.numeric(z) && !all(is.na(z)), logical(1))]
      preds <- list(); log <- list()
      mt_grid <- if (learner == "RF") MTRY_GRID else NA_integer_
      
      for (f in fo$folds) {
        tr <- pt[f$tr, , drop = FALSE]; vl <- pt[f$vl, , drop = FALSE]; te <- pt[f$te, , drop = FALSE]
        y_tr <- tr[[OUTCOME]]; y_vl <- vl[[OUTCOME]]; y_te <- te[[OUTCOME]]
        tv <- rbind(tr, vl); y_tv <- tv[[OUTCOME]]
        base <- data.frame(fold = f$fold, uid = uid_i, row_id = te$row_id,
                           true = y_te, true_bin = as.integer(y_te > SI_CUT),
                           n_pos_tv = sum(y_tv > SI_CUT), stringsAsFactors = FALSE)   # positives available to fit this fold
        
        if (config == "Baseline") {                 # rolling train+validation mean
          m   <- mean(y_tv)
          sev <- if (any(y_tv > SI_CUT)) mean(y_tv[y_tv > SI_CUT]) else 0
          p   <- mean(y_tv > SI_CUT)
          preds[[length(preds) + 1L]] <- cbind(base, data.frame(
            p_si = if (framework == "Hurdle") p else NA_real_,
            sev = if (framework == "Hurdle") sev else NA_real_,
            pred = m, threshold = NA_real_, sel_window = NA_character_, mtry = NA_integer_, n_pred = NA_integer_, pct_na_tv = NA_real_, pct_na_te = NA_real_,
            stringsAsFactors = FALSE))
          next
        }
        
        fold_ok <- drop_near_zero_variance(usable_in(usable, tr), tr)
        best <- list(mae = Inf)
        for (cand in candidates) {
          pv <- intersect(cand$vars, fold_ok)
          if (length(pv) < 2) next
          if (!is.na(MAX_PREDICTORS) && length(pv) > MAX_PREDICTORS)
            pv <- screen_predictors(tr, pv, y_tr, MAX_PREDICTORS)       # training rows only
          pp <- fit_preprocess(tr, pv); if (is.null(pp)) next
          x_tr <- tryCatch(predict(pp, tr[, pv, drop = FALSE]), error = function(e) NULL)
          x_vl <- tryCatch(predict(pp, vl[, pv, drop = FALSE]), error = function(e) NULL)
          if (is.null(x_tr) || is.null(x_vl)) next
          for (mt in mt_grid) {
            out <- fit_predict(framework, learner, x_tr, y_tr, x_vl, mt)
            if (is.null(out)) next
            mae <- mean(abs(out$pred - y_vl))
            if (is.finite(mae) && mae < best$mae)
              best <- list(mae = mae, vars = pv, cand = cand, mtry = mt, p_vl = out$p_si)
          }
        }
        if (!is.finite(best$mae)) {
          log[[length(log) + 1L]] <- data.frame(uid = uid_i, fold = f$fold, reason = "no_candidate_fit")
          next
        }
        thr <- if (framework == "Hurdle") choose_youden_threshold(as.integer(y_vl > SI_CUT), best$p_vl) else NA_real_
        
        bv <- drop_near_zero_variance(usable_in(best$vars, tv), tv)
        if (length(bv) < 2) {
          log[[length(log) + 1L]] <- data.frame(uid = uid_i, fold = f$fold, reason = "too_few_predictors_refit")
          next
        }
        pp2  <- fit_preprocess(tv, bv)
        x_tv <- if (!is.null(pp2)) tryCatch(predict(pp2, tv[, bv, drop = FALSE]), error = function(e) NULL)
        x_te <- if (!is.null(pp2)) tryCatch(predict(pp2, te[, bv, drop = FALSE]), error = function(e) NULL)
        if (is.null(x_tv) || is.null(x_te)) {
          log[[length(log) + 1L]] <- data.frame(uid = uid_i, fold = f$fold, reason = "preprocess_refit_failed")
          next
        }
        out <- fit_predict(framework, learner, x_tv, y_tv, x_te, best$mtry)
        if (is.null(out)) {
          log[[length(log) + 1L]] <- data.frame(uid = uid_i, fold = f$fold, reason = "final_fit_failed")
          next
        }
        sel <- if (!is.na(best$cand$window)) best$cand$window else best$cand$ema_window
        preds[[length(preds) + 1L]] <- cbind(base, data.frame(
          p_si = out$p_si, sev = out$sev, pred = out$pred, threshold = thr,
          sel_window = sel, mtry = best$mtry, n_pred = length(bv),
          pct_na_tv = mean(is.na(tv[, bv, drop = FALSE])),      # share of predictor cells that had to be imputed
          pct_na_te = mean(is.na(te[, bv, drop = FALSE])), stringsAsFactors = FALSE))
      }
      if (!length(preds))
        return(list(uid = uid_i, status = "no_fold_succeeded",
                    preds = NULL, log = if (length(log)) dplyr::bind_rows(log) else NULL))
      list(uid = uid_i, status = "ok", preds = dplyr::bind_rows(preds),
           log = if (length(log)) dplyr::bind_rows(log) else NULL)
    }, error = function(e) list(uid = uid_i, status = paste0("error: ", conditionMessage(e)),
                                preds = NULL, log = NULL))
  }
  
  # Reference predictors (no model): Persistence, previous prompt, window-mean features
  run_reference_participant <- function(uid_i) {
    tryCatch({
      pt_all <- DAT[DAT$uid == uid_i, , drop = FALSE]
      fo <- make_day_folds(pt_all)
      if (is.null(fo) || !length(fo$folds)) return(list(uid = uid_i, status = "no_valid_folds", preds = NULL))
      pt <- fo$data
      out <- list()
      for (f in fo$folds) {
        te <- pt[f$te, , drop = FALSE]
        tv <- pt[c(f$tr, f$vl), , drop = FALSE]
        m  <- mean(tv[[OUTCOME]])
        pers <- vapply(te$ema_num, function(tt) {
          j <- which(pt_all$ema_num <= tt - 4 * 3600)
          if (!length(j)) return(NA_real_)
          pt_all[[OUTCOME]][j[which.max(pt_all$ema_num[j])]]
        }, numeric(1))
        pos  <- match(te$row_id, pt_all$row_id)
        prev <- ifelse(pos > 1, pt_all[[OUTCOME]][pmax(pos - 1, 1)], NA_real_)
        d1 <- if ("q9_mean_d1" %in% names(te)) te$q9_mean_d1 else rep(NA_real_, nrow(te))
        d7 <- if ("q9_mean_d7" %in% names(te)) te$q9_mean_d7 else rep(NA_real_, nrow(te))
        fill <- function(v) ifelse(is.na(v), m, v)
        mk <- function(cfg, v) data.frame(
          fold = f$fold, uid = uid_i, row_id = te$row_id, true = te[[OUTCOME]],
          true_bin = as.integer(te[[OUTCOME]] > SI_CUT), p_si = NA_real_, sev = NA_real_,
          pred = fill(v), threshold = NA_real_, sel_window = NA_character_, mtry = NA_integer_,
          Configuration = cfg, stringsAsFactors = FALSE)
        out[[length(out) + 1L]] <- rbind(mk("Persistence", pers), mk("PrevPrompt", prev),
                                         mk("MeanD1", d1), mk("MeanD7", d7),
                                         mk("AlwaysZero", rep(0, nrow(te))))   # shows the mechanical MAE-by-outcome pattern
      }
      list(uid = uid_i, status = "ok", preds = dplyr::bind_rows(out))
    }, error = function(e) list(uid = uid_i, status = paste0("error: ", conditionMessage(e)), preds = NULL))
  }
  
  # ------------------------------------------------------------------------------
  # S1.7 RUNNERS WITH CHECKPOINTS
  # ------------------------------------------------------------------------------
  collect_run <- function(res, framework, learner, config, variant) {
    bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z$status), logical(1))
    res <- res[!bad]
    tag <- function(df) if (is.null(df) || !nrow(df)) NULL else
      cbind(df, framework = framework, learner = learner, Configuration = config, variant = variant,
            stringsAsFactors = FALSE)
    preds <- dplyr::bind_rows(lapply(res, function(z) tag(z$preds)))
    log   <- dplyr::bind_rows(lapply(res, function(z) tag(z$log)))
    status <- data.frame(uid = vapply(res, `[[`, "", "uid"), status = vapply(res, `[[`, "", "status"),
                         framework = framework, learner = learner, Configuration = config,
                         variant = variant, stringsAsFactors = FALSE)
    list(preds = preds, log = log, status = status)
  }
  run_one <- function(framework, learner, config, force_window = NULL, variant = "adaptive") {
    tag <- paste(framework, learner, gsub("[^A-Za-z0-9]+", "_", config), variant, SIG, sep = "__")
    chk <- file.path(CHK_DIR, paste0(tag, ".rds"))
    if (file.exists(chk)) { message("[checkpoint] ", tag); return(readRDS(chk)) }
    cand <- if (config == "Baseline") list() else make_candidates(config, force_window)
    if (config != "Baseline" && !length(cand)) { message("[no candidates] ", tag); return(NULL) }
    t0 <- Sys.time()
    res <- parallel::mclapply(UIDS, function(u) run_participant(u, framework, learner, config, cand),
                              mc.cores = N_CORES, mc.preschedule = FALSE)
    out <- collect_run(res, framework, learner, config, variant)
    saveRDS(out, chk)
    message(sprintf("[done] %s | %d/%d participants ok | %.1f min", tag,
                    sum(out$status$status == "ok"), length(UIDS),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    out
  }
  run_reference <- function() {
    chk <- file.path(CHK_DIR, paste0("Reference__", SIG, ".rds"))
    if (file.exists(chk)) { message("[checkpoint] Reference"); return(readRDS(chk)) }
    res <- parallel::mclapply(UIDS, run_reference_participant, mc.cores = N_CORES, mc.preschedule = FALSE)
    bad <- vapply(res, function(z) inherits(z, "try-error") || is.null(z$status), logical(1))
    res <- res[!bad]
    preds <- dplyr::bind_rows(lapply(res, `[[`, "preds"))
    preds$framework <- "Reference"; preds$learner <- "none"; preds$variant <- "adaptive"
    status <- data.frame(uid = vapply(res, `[[`, "", "uid"), status = vapply(res, `[[`, "", "status"),
                         framework = "Reference", learner = "none", Configuration = "Reference",
                         variant = "adaptive", stringsAsFactors = FALSE)
    out <- list(preds = preds, log = NULL, status = status)
    saveRDS(out, chk)
    out
  }
  
  # ------------------------------------------------------------------------------
  # S1.8 RUN EVERYTHING
  # ------------------------------------------------------------------------------
  all_res <- list()
  for (fw in c("Continuous", "Hurdle"))                              # baselines (learner-free)
    all_res[[paste(fw, "Baseline")]] <- run_one(fw, "none", "Baseline")
  all_res[["Reference"]] <- run_reference()
  for (i in seq_len(nrow(COMBOS)))
    for (cfg in setdiff(CONFIGS, "Baseline"))
      all_res[[paste(COMBOS$framework[i], COMBOS$learner[i], cfg)]] <-
    run_one(COMBOS$framework[i], COMBOS$learner[i], cfg)
  
  if (RUN_FIXED_WINDOWS) {
    for (i in seq_len(nrow(COMBOS)))
      for (cfg in setdiff(CONFIGS, "Baseline")) {
        wins <- if (grepl("Passive", cfg, fixed = TRUE)) windows_standard else windows_ema
        for (w in wins)
          all_res[[paste(COMBOS$framework[i], COMBOS$learner[i], cfg, w)]] <-
            run_one(COMBOS$framework[i], COMBOS$learner[i], cfg, force_window = w,
                    variant = paste0("fixed_", w))
      }
  }
  
  # ------------------------------------------------------------------------------
  # S1.9 WRITE OUTPUTS
  # ------------------------------------------------------------------------------
  all_res <- Filter(Negate(is.null), all_res)
  preds  <- dplyr::bind_rows(lapply(all_res, `[[`, "preds"))
  status <- dplyr::bind_rows(lapply(all_res, `[[`, "status"))
  flog   <- dplyr::bind_rows(lapply(all_res, `[[`, "log"))
  write.csv(preds,  file.path(OUT_DIR, "matched_predictions.csv"), row.names = FALSE)
  write.csv(status, file.path(OUT_DIR, "matched_status.csv"),      row.names = FALSE)
  write.csv(flog,   file.path(OUT_DIR, "matched_fold_log.csv"),    row.names = FALSE)
  
  fa <- do.call(rbind, lapply(UIDS, function(u) {
    pt <- DAT[DAT$uid == u, ]
    fo <- make_day_folds(pt)
    data.frame(uid = u, n_rows = nrow(pt), n_days_span = max(pt$day_index) + 1L,
               n_rows_after_runin = sum(pt$day_index >= RUNIN_DAYS),
               max_possible_folds = if (is.null(fo)) 0L else fo$max_possible,
               folds_with_enough_data = if (is.null(fo)) 0L else length(fo$folds))
  }))
  write.csv(fa, file.path(OUT_DIR, "fold_arithmetic.csv"), row.names = FALSE)
  
  cat("\nDone. Participants per model/config (status == ok):\n")
  print(status %>% group_by(framework, learner, Configuration, variant) %>%
          summarise(n_ok = sum(status == "ok"), n_total = n(), .groups = "drop"), n = 100)
  cat("\nSkipped-fold reasons:\n")
  if (nrow(flog)) print(table(flog$reason))
  cat("\nMax possible folds per participant (OSF schedule: test days 62-90 = 29 when study days 7-90 are available):\n")
  print(summary(fa$max_possible_folds))
  
}  # end STAGE 1

# ==============================================================================
# STAGE 2: POST-HOC STATISTICS
# ==============================================================================
if (RUN_STAGE2) {
  options(width = 200)                          # wide tables print on one line instead of being split
  save_csv <- function(df, name) { write.csv(df, file.path(OUT, name), row.names = FALSE); invisible(df) }
  fmt_iqr <- function(x, d = 3) {
    x <- x[is.finite(x)]; if (!length(x)) return(NA_character_)
    q <- stats::quantile(x, c(.25, .5, .75))
    sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), q[2], q[1], q[3])
  }
  
  # ---- metric helpers (no external packages) ----
  auc_fixed <- function(y, s) {                       # higher score = more SI; ties = 0.5
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(s, ties.method = "average")
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  avg_precision <- function(y, s) {
    if (sum(y) == 0) return(NA_real_)
    o <- order(s, decreasing = TRUE); y <- y[o]
    prec <- cumsum(y) / seq_along(y)
    mean(prec[y == 1])
  }
  
  if (isTRUE(TABLE2_ONLY)) P <- read.csv(file.path(NEW_DIR, "matched_predictions.csv"), colClasses = c(uid = "character"), stringsAsFactors = FALSE)[, "uid", drop = FALSE]
  if (!isTRUE(TABLE2_ONLY)) {                    # ---- Parts A, B and C1-C3 (C4 and C5 always run) ----
    # ==============================================================================
    # PART A. ROC-direction bias on the OLD outputs
    # ==============================================================================
    f_cont <- file.path(OLD_DIR, "continuous_predictions.csv")
    f_hurd <- file.path(OLD_DIR, "hurdle_predictions.csv")
    if (RUN_PART_A && file.exists(f_cont) && file.exists(f_hurd) && requireNamespace("pROC", quietly = TRUE)) {
      cat("\n===== PART A: auto vs fixed ROC direction (old outputs) =====\n")
      auc_auto <- function(y, s) {
        if (length(unique(y)) < 2) return(NA_real_)
        as.numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE)))
      }
      oc <- read.csv(f_cont, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
      oh <- read.csv(f_hurd, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
      pa <- bind_rows(
        oc %>% group_by(Configuration, uid) %>%
          summarise(AUC_auto = auc_auto(true_bin, pred), AUC_fixed = auc_fixed(true_bin, pred), .groups = "drop") %>%
          mutate(framework = "Continuous"),
        oh %>% group_by(Configuration, uid) %>%
          summarise(AUC_auto = auc_auto(true_bin, p_si), AUC_fixed = auc_fixed(true_bin, p_si), .groups = "drop") %>%
          mutate(framework = "Hurdle"))
      A <- pa %>% group_by(framework, Configuration) %>%
        summarise(n = sum(!is.na(AUC_fixed)),
                  median_auto = median(AUC_auto, na.rm = TRUE), median_fixed = median(AUC_fixed, na.rm = TRUE),
                  pct_gt_half_auto = mean(AUC_auto > .5, na.rm = TRUE), pct_gt_half_fixed = mean(AUC_fixed > .5, na.rm = TRUE),
                  n_flipped_up = sum(AUC_auto > .5 & AUC_fixed < .5, na.rm = TRUE), .groups = "drop")
      print(as.data.frame(A), digits = 3)
      save_csv(A, "A_roc_direction_bias.csv")
    } else cat("\nPART A skipped (old prediction files or pROC not found).\n")
    
    # ==============================================================================
    # PART B. Matched models
    # ==============================================================================
    f_new <- file.path(NEW_DIR, "matched_predictions.csv")
    if (!file.exists(f_new)) stop("matched_predictions.csv not found: run si_matched_models.R first.")
    cat("\n===== PART B: matched models =====\n")
    
    K <- function(model, config, variant = "adaptive") paste(model, config, variant, sep = " | ")
    P <- read.csv(f_new, colClasses = c(uid = "character"), stringsAsFactors = FALSE) %>%
      mutate(Model = paste0(framework, "-", learner),
             score_occ = ifelse(framework == "Hurdle", p_si, pred),    # occurrence score used for AUC/AP
             key = K(Model, Configuration, variant),
             uf = paste(uid, fold))
    
    # ---- B1. common support --------------------------------------------------------
    adaptive_keys <- unique(P$key[P$variant == "adaptive"])
    sets   <- split(P$uf, P$key)
    common <- Reduce(intersect, sets[adaptive_keys])
    ret <- data.frame(key = adaptive_keys,
                      n_folds = vapply(adaptive_keys, function(k) length(unique(sets[[k]])), 1L),
                      stringsAsFactors = FALSE)
    ret$n_common <- length(common); ret$pct_kept <- round(100 * ret$n_common / ret$n_folds, 1)
    print(ret, row.names = FALSE)
    save_csv(ret, "B1_common_support.csv")
    if (!length(common)) stop("No (uid, fold) is shared by all adaptive models; drop a model from the run or relax support.")
    if (min(ret$pct_kept) < 80) warning("Common support keeps <80% of folds for some model; see B1_common_support.csv")
    Pc <- P %>% filter(uf %in% common)
    
    # ---- B2. person-level metrics --------------------------------------------------
    pm <- function(d) {
      y <- d$true_bin; s <- d$score_occ; t <- d$true; p <- d$pred
      sst <- sum((t - mean(t))^2)
      # threshold-based metrics: hurdle models only, threshold chosen on the validation split (Youden).
      # Baseline and continuous rows have no threshold, so these are NA for them.
      has_thr <- all(!is.na(d$p_si)) && all(!is.na(d$threshold))
      yhat <- if (has_thr) as.integer(d$p_si >= d$threshold) else rep(NA_integer_, nrow(d))
      sens <- if (has_thr && any(y == 1)) mean(yhat[y == 1] == 1) else NA_real_
      spec <- if (has_thr && any(y == 0)) mean(yhat[y == 0] == 0) else NA_real_
      # severity stage alone (hurdle models): accuracy of the conditional severity prediction on SI-positive days
      sev_ok  <- !is.null(d$sev) && all(!is.na(d$sev)) && any(y == 1)
      mae_sev <- if (sev_ok) mean(abs(d$sev[y == 1] - t[y == 1])) else NA_real_
      r_sev   <- if (sev_ok && sum(y == 1) >= 3 && sd(d$sev[y == 1]) > 0 && sd(t[y == 1]) > 0)
        suppressWarnings(cor(d$sev[y == 1], t[y == 1])) else NA_real_
      med <- function(x) if (is.null(x) || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
      avg <- function(x) if (is.null(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      tibble(
        n_obs = nrow(d), n_SI = sum(y == 1), n_nonSI = sum(y == 0), n_folds = n_distinct(d$fold),
        prevalence = mean(y),
        AUC = auc_fixed(y, s), AUC_comb = auc_fixed(y, p), AP = avg_precision(y, s),
        MAE = mean(abs(p - t)),
        MAE_zero = if (any(y == 0)) mean(abs(p[y == 0] - t[y == 0])) else NA_real_,
        MAE_pos  = if (any(y == 1)) mean(abs(p[y == 1] - t[y == 1])) else NA_real_,
        RMSE = sqrt(mean((p - t)^2)),
        R2 = if (sst == 0) NA_real_ else 1 - sum((t - p)^2) / sst,
        r  = if (sd(t) == 0) NA_real_ else if (sd(p) == 0) 0 else suppressWarnings(cor(p, t)),   # constant OUTCOME: r undefined. Constant PREDICTIONS are a model failure: r = 0 and the participant stays in
        Brier = if (all(!is.na(d$p_si))) mean((d$p_si - y)^2) else NA_real_,
        Sens = sens, Spec = spec,
        BalAcc = if (is.na(sens) || is.na(spec)) NA_real_ else (sens + spec) / 2,
        pred_sd = sd(p),                       # spread of predictions: flatter predictions can raise R2 without signal
        MAE_sev = mae_sev, r_sev = r_sev,
        n_pos_tv_med = med(d$n_pos_tv), n_pred_med = med(d$n_pred), pct_na_tv = avg(d$pct_na_tv),
        const_score = sd(s) == 0)
    }
    M <- Pc %>% group_by(key, Model, Configuration, variant, uid) %>% group_modify(~ pm(.x)) %>% ungroup()
    save_csv(M, "B2_person_metrics.csv")
    
    # ---- B2c. adaptive vs fixed-window pairs scored on the SAME test days ---------------------------------
    # The common support above is built from the adaptive models. A fixed window can fail in folds where the adaptive model
    # succeeded (e.g. the 1-day EMA window), so each H4 pair is re-scored on the (uid, fold) cells that BOTH members produced.
    .pm_cache <- new.env()
    pair_metrics <- function(ka, kv) {
      id <- paste(ka, kv, sep = " || ")
      if (!is.null(.pm_cache[[id]])) return(.pm_cache[[id]])
      cv  <- Reduce(intersect, list(sets[[ka]], sets[[kv]], common))
      out <- Pc %>% filter(key %in% c(ka, kv), uf %in% cv) %>%
        group_by(key, Model, Configuration, variant, uid) %>% group_modify(~ pm(.x)) %>% ungroup()
      chk <- out %>% group_by(uid) %>% summarise(k = n_distinct(n_obs), .groups = "drop")
      stopifnot(all(chk$k == 1))                      # both members of a pair must be scored on identical observations
      assign(id, out, envir = .pm_cache)
      out
    }
    
    # ---- B2b. pooled calibration of the occurrence probability (hurdle models) ---------------
    # Bins of predicted P(SI) vs the observed rate. Pooled across people, so it shows calibration-in-the-large only.
    CAL <- Pc %>% filter(framework == "Hurdle", !is.na(p_si)) %>%
      mutate(bin = cut(p_si, breaks = c(0, .05, .1, .2, .3, .4, .5, .7, 1), include.lowest = TRUE)) %>%
      group_by(Model, Configuration, variant, bin) %>%
      summarise(n_obs = n(), n_people = n_distinct(uid), mean_predicted = mean(p_si),
                observed_rate = mean(true_bin), .groups = "drop")
    save_csv(CAL, "B2b_calibration_bins.csv")
    
    # ---- B3. summary ---------------------------------------------------------------
    S <- M %>% group_by(Model, Configuration, variant) %>%
      summarise(n_people = n(),
                pct_AUC_gt_half = round(100 * mean(AUC > .5, na.rm = TRUE)),   # before AUC is overwritten
                pct_AUC_ge_70   = round(100 * mean(AUC >= .70, na.rm = TRUE)),
                AUC = fmt_iqr(AUC),
                AP = fmt_iqr(AP), MAE = fmt_iqr(MAE, 2), MAE_zero = fmt_iqr(MAE_zero, 2),
                MAE_pos = fmt_iqr(MAE_pos, 2), R2 = fmt_iqr(R2), r = fmt_iqr(r),
                Brier = fmt_iqr(Brier), Sens = fmt_iqr(Sens), Spec = fmt_iqr(Spec), BalAcc = fmt_iqr(BalAcc),
                pred_sd = fmt_iqr(pred_sd, 2), MAE_sev = fmt_iqr(MAE_sev, 2), r_sev = fmt_iqr(r_sev),
                n_pred = fmt_iqr(n_pred_med, 0), n_pos_tv = fmt_iqr(n_pos_tv_med, 0), pct_imputed = fmt_iqr(pct_na_tv, 3),
                n_const_score = sum(const_score), .groups = "drop") %>%
      arrange(Model, match(Configuration, c(CONFIGS, "Persistence", "PrevPrompt", "MeanD1", "MeanD7", "AlwaysZero")), variant)
    save_csv(S, "B3_summary.csv")
    print(as.data.frame(S %>% filter(variant == "adaptive") %>% select(Model, Configuration, n_people, AUC, MAE, R2)),
          row.names = FALSE)
    
    # ---- B4. paired contrasts ------------------------------------------------------
    set.seed(SEED)
    contrast <- function(M, a, b, metric, family) {
      A <- M[M$key == a, c("uid", metric)]; Bb <- M[M$key == b, c("uid", metric)]
      if (!nrow(A) || !nrow(Bb)) return(NULL)
      names(A)[2] <- "a"; names(Bb)[2] <- "b"
      w <- merge(A, Bb, by = "uid"); w <- w[is.finite(w$a) & is.finite(w$b), ]
      d <- w$a - w$b
      if (length(d) < 5) return(NULL)
      boot <- replicate(B_BOOT, median(sample(d, replace = TRUE)))
      p <- tryCatch(suppressWarnings(wilcox.test(d)$p.value), error = function(e) NA_real_)
      data.frame(family = family, a = a, b = b, metric = metric, n = length(d),
                 median_diff = median(d), ci_lo = unname(quantile(boot, .025)),
                 ci_hi = unname(quantile(boot, .975)), p = p, stringsAsFactors = FALSE)
    }
    models <- sort(unique(M$Model[!M$Model %in% c("Continuous-none", "Hurdle-none", "Reference-none")]))
    base_key <- function(m) K(sub("-.*$", "-none", m), "Baseline")
    pers_key <- K("Reference-none", "Persistence")
    specs <- list()
    add <- function(a, b, family) specs[[length(specs) + 1L]] <<- data.frame(a = a, b = b, family = family, stringsAsFactors = FALSE)
    for (m in models) {
      bk <- base_key(m)
      add(K(m, "EMA + Passive"), K(m, "EMA"),      "H1_vs_EMA")
      add(K(m, "EMA + Passive"), K(m, "Passive"),  "H1_vs_Passive")
      add(K(m, "EMA + Passive"), bk,               "H1_vs_Baseline")
      add(K(m, "EMA + Passive (dynamic)"), K(m, "EMA (dynamic)"),     "H1_vs_EMA_dynamic")
      add(K(m, "EMA + Passive (dynamic)"), K(m, "Passive (dynamic)"), "H1_vs_Passive_dynamic")
      add(K(m, "EMA (dynamic)"), K(m, "EMA"),                         "H2_dynamic")
      add(K(m, "Passive (dynamic)"), K(m, "Passive"),                 "H2_dynamic")
      add(K(m, "EMA + Passive (dynamic)"), K(m, "EMA + Passive"),     "H2_dynamic")
      for (cfg in CONFIGS[-1]) { add(K(m, cfg), bk, "vs_baseline"); add(K(m, cfg), pers_key, "vs_persistence") }
    }
    for (fw in c("Continuous", "Hurdle")) add(K(paste0(fw, "-none"), "Baseline"), pers_key, "baseline_vs_persistence")
    learners <- unique(sub("^[A-Za-z]+-", "", models))
    for (L in learners) for (cfg in CONFIGS)       # H3: same learner, different framework
      add(K(paste0("Hurdle-", L), cfg), K(paste0("Continuous-", L), cfg), "H3_hurdle_vs_continuous_same_learner")
    if (all(c("RF", "GBM") %in% learners)) for (fw in c("Continuous", "Hurdle")) for (cfg in CONFIGS[-1])
      add(K(paste0(fw, "-GBM"), cfg), K(paste0(fw, "-RF"), cfg), "learner_effect_GBM_vs_RF")
    for (v in setdiff(unique(M$variant), "adaptive")) for (m in models) for (cfg in CONFIGS[-1])
      add(K(m, cfg, "adaptive"), K(m, cfg, v), "H4_adaptive_vs_fixed_window")
    specs <- bind_rows(specs) %>% filter(a %in% M$key, b %in% M$key)
    
    METRICS <- c("AUC", "AUC_comb", "AP", "MAE", "MAE_zero", "MAE_pos", "BalAcc", "R2", "MAE_sev", "r", "RMSE")
    set.seed(SEED)                                              # B4 bootstraps: the same draws on every run, whatever ran before
    res <- bind_rows(lapply(seq_len(nrow(specs)), function(i) {
      Msrc <- if (specs$family[i] == "H4_adaptive_vs_fixed_window") pair_metrics(specs$a[i], specs$b[i]) else M   # H4 pairs: same participants AND same test days
      bind_rows(lapply(METRICS, function(mt) contrast(Msrc, specs$a[i], specs$b[i], mt, specs$family[i])))
    }))
    # ---- B8b. adaptive vs fixed windows on the SAME test days (the numbers behind H4) ----------------------
    W8b <- bind_rows(lapply(which(specs$family == "H4_adaptive_vs_fixed_window"), function(i) {
      pair_metrics(specs$a[i], specs$b[i]) %>% group_by(Model, Configuration, variant) %>%
        summarise(n_people = n(), n_obs = sum(n_obs), AUC = median(AUC, na.rm = TRUE), MAE = median(MAE, na.rm = TRUE),
                  R2 = median(R2, na.rm = TRUE), r = median(r, na.rm = TRUE), .groups = "drop") %>%
        mutate(pair = paste(specs$a[i], "vs", specs$b[i]))
    }))
    if (nrow(W8b)) save_csv(W8b, "B8b_window_pairs_same_days.csv")
    if (!nrow(res)) {
      cat("\nB4: fewer than 5 people have usable metrics, so no contrasts were computed (normal in a small test run).\n")
      res <- data.frame(family = character(), a = character(), b = character(), metric = character(), n = integer(),
                        median_diff = numeric(), ci_lo = numeric(), ci_hi = numeric(), p = numeric(), stringsAsFactors = FALSE)
    }
    res <- res %>% group_by(family, metric) %>% mutate(q_BH = p.adjust(p, "BH")) %>% ungroup()
    # Confirmatory = the two preregistered hypotheses, for the preregistered main model and metrics. Everything else is exploratory.
    mod_a <- vapply(strsplit(res$a, " | ", fixed = TRUE), function(z) z[1], "")
    cfg_a <- vapply(strsplit(res$a, " | ", fixed = TRUE), function(z) z[2], "")
    osf_cfgs <- if (OSF_FEATURE_SET == "dynamic") c("EMA (dynamic)", "Passive (dynamic)", "EMA + Passive (dynamic)") else c("EMA", "Passive", "EMA + Passive")
    osf_h2_fam <- if (OSF_FEATURE_SET == "dynamic") c("H1_vs_EMA_dynamic", "H1_vs_Passive_dynamic") else c("H1_vs_EMA", "H1_vs_Passive")
    osf_h1 <- res$family == "vs_baseline" & cfg_a %in% osf_cfgs
    osf_h2 <- res$family %in% osf_h2_fam
    conf <- mod_a == OSF_MAIN_MODEL & (osf_h1 | osf_h2) & res$metric %in% c("R2", "r", "MAE", "RMSE")
    res$tier <- ifelse(conf, "confirmatory (OSF)", "exploratory")
    res$osf_hypothesis <- ifelse(conf, ifelse(osf_h1, "H1", "H2"), NA_character_)
    res$better <- ifelse(res$metric %in% c("MAE", "RMSE", "MAE_zero", "MAE_pos", "MAE_sev"), "lower", "higher")
    res$q_confirmatory <- NA_real_                          # one BH correction across all confirmatory tests of the same metric
    for (mt in unique(res$metric[conf])) { i <- conf & res$metric == mt; res$q_confirmatory[i] <- p.adjust(res$p[i], "BH") }
    save_csv(res, "B4_contrasts.csv")
    cat("\nCONFIRMATORY (OSF H1/H2) contrasts for", OSF_MAIN_MODEL, "(diff = a - b; 'better' says which direction favours a):\n")
    print(as.data.frame(res %>% filter(tier != "exploratory") %>%
                          select(osf_hypothesis, a, b, metric, better, n, median_diff, ci_lo, ci_hi, q_confirmatory)), digits = 3, row.names = FALSE)
    
    # ---- B5. per-person reliability test: circular-shift permutation ------------------------
    # Question: is this person's AUC higher than expected from two unrelated series with the same
    # persistence? Null: rotate the outcome series circularly against the score series (every
    # non-zero shift, in time order). Unlike shuffling days at random, a rotation KEEPS the
    # autocorrelation of both series, so a slowly drifting score (e.g. a rolling mean) is not
    # rewarded for matching an episodic outcome by chance. (A random day-shuffle flagged 17-25% of
    # unrelated persistent series at a nominal 5% in simulation.) One-sided: AUC > chance.
    # p = (1 + #{shifted AUC >= observed}) / n_rows, so its resolution is 1/n_rows; it is deterministic.
    perm_p <- function(d) {
      d <- d[order(d$row_id), , drop = FALSE]                  # time order (row_id increases with time)
      y <- d$true_bin; sc <- d$score_occ
      n <- length(y); n1 <- sum(y == 1); n0 <- n - n1
      if (n < 10 || n1 == 0 || n0 == 0) return(data.frame(AUC = NA_real_, p_perm = NA_real_))
      R <- rank(sc, ties.method = "average"); base <- n1 * (n1 + 1) / 2
      auc_of <- function(yv) (sum(R[yv == 1]) - base) / (n1 * n0)
      obs <- auc_of(y)
      shifted <- vapply(seq_len(n - 1), function(k) auc_of(y[((seq_len(n) - 1 + k) %% n) + 1]), numeric(1))
      data.frame(AUC = obs, p_perm = (1 + sum(shifted >= obs - 1e-12)) / n)
    }
    perm_keys <- unique(M$key[M$Model %in% c(PRIMARY_MODELS, "Reference-none", "Hurdle-none", "Continuous-none") &
                                M$variant == "adaptive"])
    set.seed(SEED)
    perm <- bind_rows(lapply(perm_keys, function(k) {
      dk <- Pc %>% filter(key == k)
      out <- dk %>% group_by(uid) %>% group_modify(~ perm_p(.x)) %>% ungroup()
      out$key <- k; out
    })) %>% group_by(key) %>% mutate(q_BH = p.adjust(p_perm, "BH")) %>% ungroup()
    save_csv(perm, "B5_person_permutation.csv")
    PS <- perm %>% group_by(key) %>%
      summarise(n = sum(!is.na(p_perm)), pct_p_lt_05 = round(100 * mean(p_perm < .05, na.rm = TRUE), 1),
                n_q_lt_05 = sum(q_BH < .05, na.rm = TRUE), .groups = "drop")
    cat("\nPeople reliably above chance (expect about 5% with p<.05 by chance alone):\n")
    print(as.data.frame(PS), row.names = FALSE)
    save_csv(PS, "B5_permutation_summary.csv")
    
    # ---- B6. window-selection stability ----------------------------------------------
    K_WIN <- function(cfg) if (grepl("Passive", cfg, fixed = TRUE)) 4L else 2L
    wf <- P %>% filter(variant == "adaptive", !is.na(sel_window)) %>%
      distinct(Model, Configuration, uid, fold, sel_window) %>% arrange(Model, Configuration, uid, fold)
    stab <- wf %>% group_by(Model, Configuration, uid) %>%
      summarise(n_folds = n(), modal_window = names(which.max(table(sel_window))),
                modal_share = max(table(sel_window)) / n(),
                switch_rate = if (n() > 1) mean(sel_window[-1] != sel_window[-n()]) else NA_real_,
                share_long = mean(sel_window %in% c("day1", "day7")),     # share of folds that chose a 1- or 7-day window
                .groups = "drop")
    null_cache <- list()
    null95 <- function(n, k) {
      nm <- paste(n, k)
      if (is.null(null_cache[[nm]])) null_cache[[nm]] <<-
          unname(quantile(replicate(2000, max(tabulate(sample.int(k, n, TRUE), k)) / n), .95))
      null_cache[[nm]]
    }
    set.seed(SEED)
    stab$null95 <- mapply(function(n, cfg) null95(n, K_WIN(cfg)), stab$n_folds, stab$Configuration)
    stab$exceeds_iid_null <- stab$modal_share > stab$null95
    save_csv(stab, "B6_window_stability_by_person.csv")
    SW <- stab %>% group_by(Model, Configuration) %>%
      summarise(n_people = n(), modal_share = fmt_iqr(modal_share, 2), switch_rate = fmt_iqr(switch_rate, 2),
                share_long_window = fmt_iqr(share_long, 2), pct_exceed_iid_null = round(100 * mean(exceeds_iid_null), 1), .groups = "drop")
    save_csv(SW, "B6_window_stability_summary.csv")
    cat("\nWindow stability (the i.i.d. null is lenient: neighbouring folds share ~97% of training data, so\n",
        "stable choices are expected even without person-specific windows. The real H4 test is adaptive vs\n",
        "fixed-window performance, family H4_adaptive_vs_fixed_window in B4_contrasts.csv).\n", sep = "")
    print(as.data.frame(SW), row.names = FALSE)
    save_csv(as.data.frame(table(wf$Model, wf$Configuration, wf$sel_window)) %>% rename(Model = Var1, Configuration = Var2, window = Var3),
             "B6_window_choice_counts.csv")
    
    # ---- B8. performance by window: adaptive vs each fixed window (needs RUN_FIXED_WINDOWS = TRUE) ----------
    # Evidence for or against 'passive sensing favours longer windows' and for H4 (personalised vs single window).
    W8 <- M %>% filter(Configuration != "Baseline", !grepl("-none$", Model)) %>%
      group_by(Model, Configuration, variant) %>%
      summarise(n_people = n(), median_AUC = round(median(AUC, na.rm = TRUE), 3),
                median_MAE = round(median(MAE, na.rm = TRUE), 2), .groups = "drop") %>%
      arrange(Model, match(Configuration, CONFIGS), match(variant, c("adaptive", "fixed_hr1", "fixed_hr4", "fixed_day1", "fixed_day7")))
    save_csv(W8, "B8_window_vs_performance.csv")
    if (n_distinct(W8$variant) > 1) { cat("\nPerformance by window (adaptive vs fixed):\n"); print(as.data.frame(W8), row.names = FALSE) } else
      cat("\nB8: only the adaptive variant exists. Set RUN_FIXED_WINDOWS <- TRUE and rerun Stage 1 to compare fixed windows.\n")
    
    # ---- B7. variance partition ---------------------------------------------------------
    if (requireNamespace("lme4", quietly = TRUE)) {
      vp <- bind_rows(lapply(PRIMARY_MODELS, function(m) {
        d <- M %>% filter(Model == m, variant == "adaptive", is.finite(AUC))
        if (n_distinct(d$Configuration) < 2) return(NULL)
        fit <- tryCatch(lme4::lmer(AUC ~ 1 + (1 | uid) + (1 | Configuration), data = d), error = function(e) NULL)
        if (is.null(fit)) { cat("  B7: variance partition could not be fitted for", m, "(too few people).\n"); return(NULL) }
        vc <- as.data.frame(lme4::VarCorr(fit))
        data.frame(Model = m, component = vc$grp, variance = vc$vcov, pct = round(100 * vc$vcov / sum(vc$vcov), 1))
      }))
      print(vp); save_csv(vp, "B7_variance_partition.csv")
    } else cat("\nB7 skipped (lme4 not installed).\n")
    
    # ==============================================================================
    # PART C. Descriptives and accounting
    # ==============================================================================
    cat("\n===== PART C =====\n")
    # ---- C1. participant flow + included vs excluded -------------------------------------
    if (file.exists(INCL_FILE)) {
      rep <- read.csv(INCL_FILE, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
      T_ <- function(x) x %in% TRUE                       # NA-safe
      s0 <- rep$uid
      s1 <- rep$uid[T_(rep$pass_ema)]
      s2 <- rep$uid[T_(rep$pass_ema) & T_(rep$pass_passive)]
      s3 <- rep$uid[T_(rep$pass_ema) & T_(rep$pass_passive) & T_(rep$pass_si)]
      s4 <- rep$uid[T_(rep$pass_ema) & T_(rep$pass_passive) & T_(rep$pass_si) & T_(rep$pass_var)]
      s5 <- rep$uid[T_(rep$pass_ema) & T_(rep$pass_passive) & T_(rep$pass_si) & T_(rep$pass_var) & T_(rep$pass_folds)]
      flow <- data.frame(step = c("In inclusion report", "EMA compliance", "+ passive coverage",
                                  "+ >=10 SI-positive and >=20 SI-negative", "+ non-zero outcome variance",
                                  "+ >=195 answered prompts (legacy rule; not used by the v2 analysis)"),
                         n = c(length(s0), length(s1), length(s2), length(s3), length(s4), length(s5)))
      flow <- rbind(flow, data.frame(step = "Analysed in every model (>=1 shared valid test day)", n = dplyr::n_distinct(M$uid)))
      print(flow, row.names = FALSE); save_csv(flow, "C1_participant_flow.csv")
      f_stat1 <- file.path(NEW_DIR, "matched_status.csv")           # configuration-specific losses
      if (file.exists(f_stat1)) {
        st1 <- read.csv(f_stat1, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
        ok1 <- st1 %>% filter(variant == "adaptive") %>% group_by(framework, learner, Configuration) %>%
          summarise(n_attempted = n(), n_with_valid_folds = sum(status == "ok"), .groups = "drop")
        save_csv(ok1, "C1_valid_folds_by_model.csv")
        cat("\nParticipants with at least one valid fold, by model and configuration:\n"); print(as.data.frame(ok1), row.names = FALSE)
      }
      rep$analysed <- rep$uid %in% s5
      cmp <- bind_rows(lapply(intersect(c("ema_compliance", "mean_cov", "n_pos", "n_neg", "q9_var", "n_obs"), names(rep)),
                              function(v) {
                                a <- rep[[v]][rep$analysed]; b <- rep[[v]][!rep$analysed & T_(rep$pass_ema)]
                                data.frame(variable = v, included = fmt_iqr(a, 2), excluded_after_EMA = fmt_iqr(b, 2),
                                           p_wilcoxon = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_))
                              }))
      print(cmp, row.names = FALSE); save_csv(cmp, "C1_included_vs_excluded.csv")   # NOTE: the comparison group here is everyone who passed the EMA criterion but was not analysed (227 - 143 = 84, incl. 2 who failed passive coverage)
      # ---- C1b. the participants excluded ONLY by the outcome-count rule (passed EMA and passive coverage, failed >= 10 SI-positive and >= 20 SI-negative)
      s_si <- rep$uid[T_(rep$pass_ema) & T_(rep$pass_passive) & !T_(rep$pass_si)]
      cat(sprintf("\nC1b: %d participants were excluded only by the outcome-count rule\n", length(s_si)))
      cmp2 <- bind_rows(lapply(intersect(c("ema_compliance", "mean_cov", "n_pos", "n_neg", "q9_var", "n_obs"), names(rep)), function(v) {
        a <- rep[[v]][rep$analysed]; b <- rep[[v]][rep$uid %in% s_si]
        data.frame(variable = v, included = fmt_iqr(a, 2), excluded_by_outcome_rule = fmt_iqr(b, 2),
                   p_wilcoxon = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_))
      }))
      print(cmp2, row.names = FALSE); save_csv(cmp2, "C1_included_vs_excluded_by_outcome_rule.csv")
      if (all(c("n_pos", "n_neg") %in% names(rep))) {                  # who are they? almost never vs almost always endorsing item 9
        sp <- rep[rep$uid %in% s_si, ]
        never <- sp$n_pos < 10 & sp$n_neg >= 20; always <- sp$n_neg < 20 & sp$n_pos >= 10; both <- sp$n_pos < 10 & sp$n_neg < 20
        split_tab <- data.frame(group = c("excluded only by the outcome-count rule", "almost never endorsed item 9 (< 10 SI-positive, >= 20 SI-negative)",
                                          "almost always endorsed it (>= 10 SI-positive, < 20 SI-negative)", "both (< 10 SI-positive and < 20 SI-negative)", "none of these (should be 0)"),
                                n = c(nrow(sp), sum(never, na.rm = TRUE), sum(always, na.rm = TRUE), sum(both, na.rm = TRUE), sum(!(never | always | both), na.rm = TRUE)))
        print(split_tab, row.names = FALSE); save_csv(split_tab, "C1_excluded_split_never_always.csv")
      } else cat("C1b: n_pos / n_neg are not in the inclusion report, so the almost-never / almost-always split was skipped.\n")
      cat("Join demographics, baseline severity and C-SSRS history onto the report by uid to extend this table.\n")
    } else cat("C1 skipped: inclusion report not found.\n")
    
    # ---- C2. outcome descriptives (for the characteristics table) ---------------------------------------------
    dat <- data.table::fread(DATA_FILE, select = c("uid", "date_rep", "q9"), colClasses = list(character = "uid")) %>%
      as.data.frame() %>% filter(!is.na(q9))
    dd <- suppressWarnings(as.Date(as.character(dat$date_rep)))
    if (all(is.na(dd))) dd <- as.Date(as.character(dat$date_rep), format = "%m/%d/%Y")
    dat$day_index <- ave(as.integer(dd), dat$uid, FUN = function(z) z - min(z))
    dat <- dat[dat$day_index < STUDY_DAYS, , drop = FALSE]        # same study days (1-90) as the models
    test_start <- RUNIN_DAYS + TRAIN_DAYS + VAL_DAYS
    desc <- dat %>% group_by(uid) %>%
      summarise(n_obs = n(), prop_nonzero = mean(q9 > 0), mean_q9 = mean(q9), sd_q9 = sd(q9),
                prop_nonzero_pretest = mean(q9[day_index < test_start] > 0),
                prop_nonzero_test_period = if (any(day_index >= test_start)) mean(q9[day_index >= test_start] > 0) else NA_real_,
                .groups = "drop")
    save_csv(desc, "C2_outcome_descriptives_by_person.csv")
    cat("\nOutcome descriptives across people (median [IQR]):\n")
    print(data.frame(variable = c("n_obs", "prop_nonzero", "mean_q9", "sd_q9", "prop_nonzero_pretest", "prop_nonzero_test_period"),
                     summary = vapply(c("n_obs", "prop_nonzero", "mean_q9", "sd_q9", "prop_nonzero_pretest", "prop_nonzero_test_period"),
                                      function(v) fmt_iqr(desc[[v]], 2), "")), row.names = FALSE)
    cat(sprintf("People with no SI in the test period: %d of %d\n",
                sum(desc$prop_nonzero_test_period == 0, na.rm = TRUE), nrow(desc)))
    
    # ---- C2b. how zero-heavy is the outcome, and what did the inclusion criteria do to it? -------------
    if (file.exists(INCL_FILE)) {
      ir <- read.csv(INCL_FILE, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
      ir <- ir[!is.na(ir$n_pos) & !is.na(ir$n_neg), , drop = FALSE]
      ir$prop_nonzero <- ir$n_pos / (ir$n_pos + ir$n_neg)
      ir$analysed <- ir$uid %in% unique(P$uid)
      grp <- function(sel, label) data.frame(
        group = label, n_people = sum(sel),
        pooled_prop_nonzero = round(sum(ir$n_pos[sel]) / sum(ir$n_pos[sel] + ir$n_neg[sel]), 3),
        median_person_prop_nonzero = fmt_iqr(ir$prop_nonzero[sel], 2),
        pct_people_never_nonzero = round(100 * mean(ir$n_pos[sel] == 0), 1),
        pct_people_mostly_nonzero = round(100 * mean(ir$prop_nonzero[sel] > .5), 1))
      zi <- rbind(grp(rep(TRUE, nrow(ir)), "All participants with outcome data"), grp(ir$analysed, "Analysed participants"))
      cat("\nHow zero-heavy is item 9, before and after the inclusion criteria?\n"); print(zi, row.names = FALSE)
      save_csv(zi, "C2b_zero_inflation_by_inclusion.csv")
    }
    pos_q9 <- dat$q9[dat$q9 > 0]
    cat(sprintf("Among nonzero responses (analysed participants): median %.0f [IQR %.0f to %.0f] on the 0-100 scale, n = %d responses.\n",
                median(pos_q9), quantile(pos_q9, .25), quantile(pos_q9, .75), length(pos_q9)))
    
    # ---- C3. skipped participants / folds ---------------------------------------------------
    f_stat <- file.path(NEW_DIR, "matched_status.csv"); f_log <- file.path(NEW_DIR, "matched_fold_log.csv")
    if (file.exists(f_stat)) {
      st <- read.csv(f_stat, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
      stt <- st %>% group_by(framework, learner, Configuration, variant, status) %>% summarise(n = n(), .groups = "drop")
      save_csv(stt, "C3_participant_status.csv"); print(as.data.frame(stt %>% filter(status != "ok")), row.names = FALSE)
    }
    if (file.exists(f_log) && file.size(f_log) > 5) {
      fl <- read.csv(f_log, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
      flt <- fl %>% group_by(framework, learner, Configuration, variant, reason) %>% summarise(n = n(), .groups = "drop")
      save_csv(flt, "C3_fold_skip_reasons.csv")
    }
    # constant-prediction people (a model failure, not an exclusion)
    cp <- M %>% filter(const_score) %>% count(Model, Configuration, variant, name = "n_people_constant_score")
    save_csv(cp, "C3_constant_score_people.csv")
    
  }                                              # ---- end of Parts A to C3 ----
  
  # ---- C4. availability of each data stream among the analysed participants (for Table 1) ---------
  streams <- c("acc", "hr", "steps", "battery", "convo", "gps", "gps_home", "sleep", "sms", "calls")
  Dall <- data.table::fread(DATA_FILE, colClasses = list(character = "uid")) %>% as.data.frame()
  Dall <- Dall[Dall$uid %in% unique(P$uid) & !is.na(Dall$q9), , drop = FALSE]
  has_by_stream <- lapply(streams, function(st) {
    cols <- if (st == "sleep") grep("^sleep_", names(Dall), value = TRUE)
    else grep(paste0("^", st, "_(hr1|hr4|day1|day7)_n_obs$"), names(Dall), value = TRUE)
    if (!length(cols)) return(NULL)
    has <- vapply(cols, function(cc) { v <- Dall[[cc]]; if (st == "sleep") !is.na(v) else (!is.na(v) & v > 0) },
                  logical(nrow(Dall)))
    tapply(rowSums(matrix(has, nrow = nrow(Dall))) > 0, Dall$uid, mean)     # share of prompts with data, per person
  })
  names(has_by_stream) <- streams
  avail <- bind_rows(lapply(streams, function(st) {
    per <- has_by_stream[[st]]
    if (is.null(per)) return(data.frame(stream = st, n_people_with_data = 0L, pct_people = 0, median_pct_prompts_with_data = NA_real_))
    data.frame(stream = st, n_people_with_data = sum(per > 0), pct_people = round(100 * mean(per > 0), 1),
               median_pct_prompts_with_data = round(100 * median(per), 1))
  }))
  ev_people <- vapply(c("sms", "calls", "convo"), function(st) {
    cols <- grep(paste0("^", st, "_(hr1|hr4|day1|day7)_sum_val$"), names(Dall), value = TRUE)
    if (!length(cols)) return(NA_integer_)
    anyv <- rowSums(do.call(cbind, lapply(cols, function(cc) !is.na(Dall[[cc]]) & Dall[[cc]] > 0))) > 0
    as.integer(sum(tapply(anyv, Dall$uid, any)))
  }, integer(1))
  avail$n_people_with_logged_events <- unname(ev_people[match(avail$stream, names(ev_people))])
  cat("(For SMS, calls and conversation, 'with data' includes zero-filled windows; 'logged_events' counts people with at least one real event.)\n")
  cat(sprintf("\nData availability among the %d analysed participants:\n", dplyr::n_distinct(Dall$uid)))
  print(avail, row.names = FALSE); save_csv(avail, "C4_stream_availability.csv")
  pp <- do.call(cbind, lapply(has_by_stream[!vapply(has_by_stream, is.null, logical(1))], function(per) as.integer(per > 0)))
  cat(sprintf("Streams with data per person (of %d): median [IQR] %s.\n", ncol(pp), fmt_iqr(rowSums(pp), 0)))
  if (all(c("sms", "calls") %in% colnames(pp)))
    cat(sprintf("People with neither SMS nor call data: %d of %d (their passive and combined models simply have no SMS/call features).\n",
                sum(pp[, "sms"] == 0 & pp[, "calls"] == 0), nrow(pp)))
  f_comp <- file.path(PROJECT_DIR, "4. Results", "report", "raw_passive_compliance_by_uid_sensor.csv")
  if (file.exists(f_comp)) {                       # the quantity behind the 75% passive-coverage criterion
    cp_ <- read.csv(f_comp, colClasses = c(uid = "character"), stringsAsFactors = FALSE)
    cp_ <- cp_[cp_$uid %in% unique(P$uid), , drop = FALSE]
    cps <- cp_ %>% group_by(sensor) %>%
      summarise(n_people = n(), median_iqr_compliance = fmt_iqr(passive_compliance, 2),
                pct_people_ge_75 = round(100 * mean(passive_compliance >= .75, na.rm = TRUE), 1), .groups = "drop")
    cat("\nPassive compliance by sensor (the 75% criterion is the mean of these six):\n"); print(as.data.frame(cps), row.names = FALSE)
    save_csv(cps, "C4_passive_compliance_by_sensor.csv")
  } else cat("(Per-sensor passive compliance file not found; only the availability table above was written.)\n")
  
  # ---- C5. Characteristics-table variables (Table 3 in the manuscript): demographics, C-SSRS history, baseline depression ---------------------------
  # Joins everything on uid and DETECTS the columns it needs from the NDA-style names (the codebooks list far more
  # columns than any one file has). It prints which columns it used; override them in the settings if they are wrong.
  as_yes <- function(x) {                          # TRUE / FALSE / NA. 1 = yes; 0 or 2 = no; negatives and >= 8 are missing codes
    z <- tolower(trimws(as.character(x)))
    num <- suppressWarnings(as.numeric(z)); isnum <- !is.na(num)
    out <- rep(NA, length(z))
    out[isnum] <- ifelse(num[isnum] < 0 | num[isnum] >= 8, NA, num[isnum] == 1)
    txt <- !isnum & z %in% c("yes", "y", "true", "t", "no", "n", "false", "f")
    out[txt] <- z[txt] %in% c("yes", "y", "true", "t")
    as.logical(out)
  }
  score03 <- function(x) { v <- suppressWarnings(as.numeric(as.character(x))); ifelse(v %in% 0:3, v, NA_real_) }   # 888/999/-9 -> NA
  norm_key <- function(x) {                        # compare IDs ignoring case, spaces and leading zeros
    z <- trimws(as.character(x)); z[z %in% c("", "NA")] <- NA
    isnum <- !is.na(z) & grepl("^[0-9]+$", z)
    z[isnum] <- sub("^0+(?=.)", "", z[isnum], perl = TRUE)
    tolower(z)
  }
  norm_uid4 <- function(x) {
    d <- gsub("\\D", "", sub("@.*$", "", as.character(x)))
    ifelse(nzchar(d), stringr::str_pad(d, 4, pad = "0"), NA_character_)
  }
  show_cols <- function(df, what) cat("  Columns in", what, ":", paste(head(names(df), 80), collapse = ", "), "\n")
  first_full_set <- function(sets, nm) { for (s in sets) if (all(s %in% nm)) return(s); NULL }
  
  # Keep the baseline rows if the file has several visits (explicit setting first, then common conventions).
  pick_baseline <- function(df, visit_col, visit_value, label) {
    n0 <- nrow(df); how <- "all rows"
    if (!is.na(visit_col) && visit_col %in% names(df) && !is.na(visit_value)) {
      df <- df[as.character(df[[visit_col]]) == as.character(visit_value) & !is.na(df[[visit_col]]), , drop = FALSE]
      how <- sprintf("%s == %s", visit_col, visit_value)
    } else {
      auto <- NULL
      num0 <- function(cc) suppressWarnings(as.numeric(as.character(df[[cc]]))) == 0
      if ("month_interview" %in% names(df) && any(num0("month_interview"), na.rm = TRUE)) {
        auto <- num0("month_interview"); how <- "month_interview == 0"
      } else if ("visnum" %in% names(df) && any(num0("visnum"), na.rm = TRUE)) {
        auto <- num0("visnum"); how <- "visnum == 0"
      } else {
        for (cc in intersect(c("visit", "timepoint_label", "version_form", "week"), names(df))) {
          m <- grepl("baseline|screen|^bl$|^t0$|month ?0|week ?0|^0$", tolower(as.character(df[[cc]])))
          if (any(m)) { auto <- m; how <- paste0(cc, " looks like baseline"); break }
        }
      }
      if (!is.null(auto)) df <- df[auto %in% TRUE, , drop = FALSE]
    }
    cat(sprintf("  %s rows: %s (%d of %d)\n", label, how, nrow(df), n0))
    df
  }
  
  CSSRS_TF_SETS <- list(                           # candidate columns per timeframe; ideation items run 1 (wish to be dead) to 5 (plan and intent)
    "3mo" = list(label = "last 3 months", suffix = "3mo",
                 ideation = list(c("si1p3", "si2p3", "si3p3", "si4p3", "si5p3")),                  # 0 = no, 1 = yes
                 attempt = c("sbaap3m", "cssrs_sa1_3mo"),                                          # 1 = yes, 2 = no
                 other = list(c("sbiap3m", "sbasiap3m", "sbpabp3m")), other_partial = TRUE,      # interrupted, aborted, preparatory
                 severity = character(0)),
    "ever" = list(label = "lifetime", suffix = "ever",
                  ideation = list(c("si1l", "si2l", "si3l", "si4l", "si5l"),
                                  c("cssrs_base_01", "cssrs_base_02", "cssrs_base_03", "cssrs_base_04", "cssrs_base_05")),
                  attempt = c("sb1l", "cssrs_base_12", "bl_c_cssrs_lifetime_sa", "cssrslifeatt_1", "cssrs_sa1"),
                  other = list(c("sb3l", "sb4l", "sb5l"), c("cssrs_base_13", "cssrs_base_14", "cssrs_base_15")), other_partial = FALSE,
                  severity = c("iilmsi", "cssrs_suicidea", "type_1_lifetime")))
  PHQ_SETS <- list(                                # 9 items scored 0-3; the 9th is always the suicide item
    list(label = "past 2 weeks (anhedonia ... si_hi)", cols = c("anhedonia", "down", "sleep_trouble", "phq_tired", "appetite",
                                                                "failure", "concentration_problems", "psychomotor_retardation", "si_hi")),
    list(label = "phq_v8_1-8 + phq_v9_9", cols = c(paste0("phq_v8_", 1:8), "phq_v9_9")),
    list(label = "phq2a-i_t1", cols = paste0("phq2", letters[1:9], "_t1")),
    list(label = "smhint_phq9_01-09", cols = sprintf("smhint_phq9_%02d", 1:9)),
    list(label = "phq_month_1-9 (past month)", cols = paste0("phq_month_", 1:9)))
  PHQ_TOTAL_CANDIDATES <- c("phq_9_total", "phq9_current")
  
  if (!is.na(DEMO_FILE) && file.exists(DEMO_FILE)) {
    cat("\n===== C5: characteristics-table variables =====\n")
    demo <- read.csv(DEMO_FILE, stringsAsFactors = FALSE, check.names = FALSE)
    fix_enc <- function(x) { if (!is.character(x)) return(x); bad <- !is.na(x) & !validUTF8(x); x[bad] <- iconv(x[bad], from = "CP1252", to = "UTF-8"); x }
    demo[] <- lapply(demo, fix_enc)                # e.g. "Bachelor\x92s" -> "Bachelor's"
    demo <- demo[, nzchar(names(demo)), drop = FALSE]   # the demographics file has a few unnamed (empty) columns; drop them
    if (!DEMO_UID_COL %in% names(demo)) { show_cols(demo, "the demographics file"); stop("C5: DEMO_UID_COL '", DEMO_UID_COL, "' not found.") }
    demo$uid <- norm_uid4(demo[[DEMO_UID_COL]])
    demo <- demo[!is.na(demo$uid) & !duplicated(demo$uid), , drop = FALSE]
    derived <- character(0)
    pid_key <- NULL                                # Participant_ID -> uid, used for every file keyed by src_subject_id
    if (DEMO_PID_COL %in% names(demo)) {
      pid_key <- stats::setNames(demo$uid, norm_key(demo[[DEMO_PID_COL]]))
      pid_key <- pid_key[!is.na(names(pid_key))]
      if (anyDuplicated(names(pid_key))) cat("  WARNING: Participant_ID is not unique in the demographics file; the first match is used.\n")
      pid_key <- pid_key[!duplicated(names(pid_key))]
    }
    cat(sprintf("Demographics: %d participants; %d of the %d analysed are present.\n", nrow(demo),
                sum(unique(P$uid) %in% demo$uid), dplyr::n_distinct(P$uid)))
    
    # ---- age in years (NDA's interview_age is in months) and automatic demographic columns ----
    if ("interview_age" %in% names(demo)) {
      a <- suppressWarnings(as.numeric(demo$interview_age)); demo$age_years <- if (median(a, na.rm = TRUE) > 150) a / 12 else a
    } else {
      agecol <- names(demo)[grepl("^age($|[^a-z])", tolower(names(demo)))][1]      # e.g. "Age at Study Start"
      if (!is.na(agecol)) demo$age_years <- suppressWarnings(as.numeric(demo[[agecol]]))
    }
    auto_demo <- unique(c(if ("age_years" %in% names(demo)) "age_years", names(demo)[grepl("(^|_)(sex|gender|race|ethnic|hisp|educ)", tolower(names(demo))) & !grepl("(_other|_spec)$", tolower(names(demo)))]))
    auto_demo <- setdiff(auto_demo, c(DEMO_UID_COL, DEMO_PID_COL, "uid", "interview_age"))
    
    # ---- C-SSRS history (Participant_ID <-> src_subject_id) ----
    if (!is.na(CSSRS_FILE) && file.exists(CSSRS_FILE)) {
      cs <- read.csv(CSSRS_FILE, stringsAsFactors = FALSE, check.names = FALSE)
      if (is.null(pid_key)) { show_cols(demo, "the demographics file"); cat("  C-SSRS skipped: DEMO_PID_COL not found.\n")
      } else if (!CSSRS_ID_COL %in% names(cs)) { show_cols(cs, "the C-SSRS file"); cat("  C-SSRS skipped: CSSRS_ID_COL not found.\n")
      } else {
        cs$uid <- unname(pid_key[norm_key(cs[[CSSRS_ID_COL]])])
        cat(sprintf("C-SSRS: %d rows, %d matched to a uid (%d unmatched), %d participants.\n", nrow(cs), sum(!is.na(cs$uid)),
                    sum(is.na(cs$uid)), dplyr::n_distinct(cs$uid[!is.na(cs$uid)])))
        cs <- cs[!is.na(cs$uid), , drop = FALSE]
        if (isTRUE(ONE_ROW_PER_PERSON)) {
          cat(sprintf("  C-SSRS: one row per person expected; %d participants appear more than once%s\n", sum(duplicated(cs$uid)),
                      if (any(duplicated(cs$uid))) " (history = any 'yes' across their rows)." else "."))
        } else cs <- pick_baseline(cs, CSSRS_VISIT_COL, CSSRS_BASELINE_VALUE, "C-SSRS")
        nm <- names(cs)
        if (identical(CSSRS_TIMEFRAME, "screen")) {
          # Screening C-SSRS. Items 1-5 (wish to be dead ... plan and intent) are coded 0 = no, 1 = yes, 2 = yes in the past month,
          # 3 = yes but not currently. cssrs_add_6 = suicidal behaviour ever; cssrs_add_7 = in the past 3 months (asked if add_6 is yes).
          items <- paste0("cssrs_add_", 1:5)
          if (!all(items %in% nm)) {
            show_cols(cs, "the C-SSRS file"); cat("  C-SSRS screener items cssrs_add_1 to cssrs_add_5 not found.\n")
          } else {
            code <- function(x) { v <- suppressWarnings(as.numeric(as.character(x))); ifelse(v %in% 0:3, v, NA_real_) }
            cat("  C-SSRS screener codes (0 = no, 1 = yes, 2 = yes in the past month, 3 = yes but not currently):\n")
            for (cc in c(items, intersect(c("cssrs_add_6", "cssrs_add_7"), nm))) {
              tb <- table(code(cs[[cc]]), useNA = "ifany")
              cat(sprintf("    %s: %s\n", cc, paste(names(tb), as.integer(tb), sep = "=", collapse = ", ")))
            }
            yn <- function(v) ifelse(is.na(v), NA_character_, ifelse(v, "Yes", "No"))
            by_uid <- function(v, f) { r <- tapply(v, cs$uid, function(z) if (all(is.na(z))) NA else f(z)); unname(r[cdf$uid]) }
            cdf <- data.frame(uid = sort(unique(cs$uid)), stringsAsFactors = FALSE)
            endorsed <- do.call(cbind, lapply(items, function(cc) code(cs[[cc]]) >= 1))
            sev_row <- ifelse(rowSums(!is.na(endorsed)) == 0, NA_real_,
                              apply(endorsed, 1, function(r) { w <- which(r %in% TRUE); if (length(w)) max(w) else 0 }))
            sev <- as.numeric(by_uid(sev_row, function(z) max(z, na.rm = TRUE)))
            # The screening C-SSRS was given to people with a non-zero PHQ-9 item 9 at screening. A person with NO valid answers is therefore
            # "not asked" if their item 9 was 0, or a gap in the record if it was above 0. By decision both are counted as No (recoded after the PHQ-9 step).
            administered <- as.logical(by_uid(rowSums(!is.na(endorsed)) > 0, function(z) any(z, na.rm = TRUE)))
            administered[is.na(administered)] <- FALSE
            cat(sprintf("  C-SSRS answered by %d of %d people; the others are counted as No below.\n",
                        sum(administered), length(administered)))
            type_lab <- c("0 None", "1 Wish to be dead", "2 Non-specific active suicidal thoughts", "3 Thoughts with method",
                          "4 Thoughts with some intent", "5 Plan with intent")
            cdf$cssrs_ideation_screen <- ifelse(!administered, "Not answered", yn(sev >= 1))
            cdf$cssrs_max_ideation_type_screen <- ifelse(!administered, "Not answered", type_lab[sev + 1])
            if ("cssrs_add_6" %in% nm) {
              b6 <- code(cs$cssrs_add_6); b7 <- if ("cssrs_add_7" %in% nm) code(cs$cssrs_add_7) else rep(NA_real_, nrow(cs))
              cdf$cssrs_behavior_lifetime <- ifelse(!administered, "Not answered", yn(as.logical(by_uid(b6 >= 1, function(z) any(z, na.rm = TRUE)))))
              b3 <- ifelse(b6 == 0, 0, b7)                                   # not asked when there was no behaviour: counts as no
              cdf$cssrs_behavior_3mo <- ifelse(!administered, "Not answered", yn(as.logical(by_uid(b3 >= 1, function(z) any(z, na.rm = TRUE)))))
            }
            demo <- merge(demo, cdf, by = "uid", all.x = TRUE)
            derived <- c(derived, setdiff(names(cdf), "uid"))
            cat("  C-SSRS variables: ideation = any item 1-5 endorsed (any timeframe); max type 0-5; behaviour ever; behaviour in the past 3 months.\n")
          }
        } else {
          tf <- CSSRS_TF_SETS[[CSSRS_TIMEFRAME]]
          if (is.null(tf)) stop("C5: CSSRS_TIMEFRAME must be one of: ", paste(names(CSSRS_TF_SETS), collapse = ", "))
          suf <- tf$suffix
          cat(sprintf("  C-SSRS timeframe: %s (CSSRS_TIMEFRAME = '%s')\n", tf$label, CSSRS_TIMEFRAME))
          ide <- if (length(CSSRS_IDEATION_COLS)) CSSRS_IDEATION_COLS else first_full_set(tf$ideation, nm)
          att <- if (length(CSSRS_ATTEMPT_COLS)) CSSRS_ATTEMPT_COLS else head(intersect(tf$attempt, nm), 1)
          oth <- if (length(CSSRS_OTHER_BEHAVIOR_COLS)) CSSRS_OTHER_BEHAVIOR_COLS else
            if (isTRUE(tf$other_partial)) unique(unlist(lapply(tf$other, function(x) intersect(x, nm)))) else first_full_set(tf$other, nm)
          sevcol <- head(intersect(tf$severity, nm), 1)
          cat("  C-SSRS columns used: ideation =", if (length(ide)) paste(ide, collapse = ", ") else "(none)",
              "| attempt =", if (length(att)) paste(att, collapse = ", ") else "(none)",
              "| other behaviour =", if (length(oth)) paste(oth, collapse = ", ") else "(none)",
              "| severity column =", if (length(sevcol)) sevcol else "(none)", "\n")
          if (!length(ide) && !length(att) && !length(sevcol)) {
            show_cols(cs, "the C-SSRS file")
            cat("  C-SSRS: no known ideation or attempt columns found. Set CSSRS_IDEATION_COLS and CSSRS_ATTEMPT_COLS.\n")
          } else {
            cdf <- data.frame(uid = sort(unique(cs$uid)), stringsAsFactors = FALSE)
            by_uid <- function(v, f) { r <- tapply(v, cs$uid, function(z) if (all(is.na(z))) NA else f(z)); unname(r[cdf$uid]) }
            row_any <- function(cols) {
              m <- do.call(cbind, lapply(cols, function(cc) as_yes(cs[[cc]])))
              ifelse(rowSums(!is.na(m)) == 0, NA, rowSums(m, na.rm = TRUE) > 0)
            }
            yn <- function(v) ifelse(is.na(v), NA_character_, ifelse(v, "Yes", "No"))
            sev <- NULL
            if (length(ide)) {
              m <- do.call(cbind, lapply(ide, function(cc) as_yes(cs[[cc]])))
              sev_row <- ifelse(rowSums(!is.na(m)) == 0, NA_real_,
                                apply(m, 1, function(r) { w <- which(r %in% TRUE); if (length(w)) max(w) else 0 }))
              sev <- by_uid(sev_row, function(z) max(z, na.rm = TRUE))
            } else if (length(sevcol)) {
              v <- suppressWarnings(as.numeric(as.character(cs[[sevcol]]))); v[v < 0 | v > 5] <- NA
              sev <- by_uid(v, function(z) max(z, na.rm = TRUE))
            }
            if (!is.null(sev)) {
              sevv <- as.numeric(sev)
              cdf[[paste0("cssrs_max_ideation_severity_", suf)]] <- sevv
              cdf[[paste0("cssrs_ideation_", suf)]] <- yn(sevv >= 1)
              cdf[[paste0("cssrs_active_ideation_", suf)]] <- yn(sevv >= 2)
            }
            if (length(att)) cdf[[paste0("cssrs_attempt_", suf)]] <- yn(as.logical(by_uid(row_any(att), function(z) any(z, na.rm = TRUE))))
            if (length(att) || length(oth))
              cdf[[paste0("cssrs_any_behavior_", suf)]] <- yn(as.logical(by_uid(row_any(c(att, oth)), function(z) any(z, na.rm = TRUE))))
            demo <- merge(demo, cdf, by = "uid", all.x = TRUE)
            derived <- c(derived, setdiff(names(cdf), "uid"))
            cat(sprintf("  C-SSRS variables are for the %s; 'any behaviour' = attempt, interrupted, aborted or preparatory. Negative and 8+ codes are missing.\n", tf$label))
          }
        }
      }
    } else cat("C-SSRS file not found; suicide history skipped.\n")
    
    # ---- baseline depression: PHQ-9 total, and the total without item 9 (which is the outcome) ----
    if (!is.na(BASELINE_DEP_FILE) && file.exists(BASELINE_DEP_FILE)) {
      dp <- read.csv(BASELINE_DEP_FILE, stringsAsFactors = FALSE, check.names = FALSE)
      if (!BASELINE_DEP_ID_COL %in% names(dp)) { show_cols(dp, "the PHQ-9 file"); cat("  PHQ-9 skipped: BASELINE_DEP_ID_COL not found.\n")
      } else {
        if (isTRUE(BASELINE_DEP_VIA_PID)) {
          if (is.null(pid_key)) stop("C5: BASELINE_DEP_VIA_PID is TRUE but DEMO_PID_COL was not found in the demographics file.")
          dp$uid <- unname(pid_key[norm_key(dp[[BASELINE_DEP_ID_COL]])])
          cat(sprintf("PHQ-9 file: %d rows, %d matched to a uid (%d unmatched).\n", nrow(dp), sum(!is.na(dp$uid)), sum(is.na(dp$uid))))
        } else dp$uid <- norm_uid4(dp[[BASELINE_DEP_ID_COL]])
        dp <- dp[!is.na(dp$uid), , drop = FALSE]
        if (isTRUE(ONE_ROW_PER_PERSON)) {
          cat(sprintf("  PHQ-9: one row per person expected; %d participants appear more than once.\n", sum(duplicated(dp$uid))))
        } else dp <- pick_baseline(dp, BASELINE_DEP_VISIT_COL, BASELINE_DEP_VISIT_VALUE, "PHQ-9")
        ord <- tryCatch({
          if ("days_baseline" %in% names(dp) && any(!is.na(suppressWarnings(as.numeric(dp$days_baseline))))) suppressWarnings(as.numeric(dp$days_baseline))
          else if ("interview_date" %in% names(dp)) as.numeric(as.Date(dp$interview_date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")))
          else seq_len(nrow(dp))
        }, error = function(e) seq_len(nrow(dp)))
        dp <- dp[order(dp$uid, ord), , drop = FALSE]
        if (any(duplicated(dp$uid))) cat(sprintf("  %d participants have several PHQ-9 rows; the earliest is used.\n", sum(duplicated(dp$uid))))
        dp <- dp[!duplicated(dp$uid), , drop = FALSE]
        items <- NULL
        if (length(BASELINE_PHQ_ITEM_COLS) == 9 && all(BASELINE_PHQ_ITEM_COLS %in% names(dp))) {
          items <- BASELINE_PHQ_ITEM_COLS; cat("  PHQ-9 items: set in BASELINE_PHQ_ITEM_COLS\n")
        } else {
          st <- NULL; for (s0 in PHQ_SETS) if (all(s0$cols %in% names(dp))) { st <- s0; break }
          if (!is.null(st)) { items <- st$cols; cat("  PHQ-9 items detected:", st$label, "\n") }
        }
        if (!is.null(items)) {
          it <- matrix(sapply(items, function(cc) score03(dp[[cc]])), nrow = nrow(dp))
          dp$phq9_total <- ifelse(rowSums(is.na(it)) == 0, rowSums(it), NA_real_)
          dp$phq9_items1to8 <- ifelse(rowSums(is.na(it[, 1:8, drop = FALSE])) == 0, rowSums(it[, 1:8, drop = FALSE]), NA_real_)
          derived <- c(derived, "phq9_total", "phq9_items1to8")
        } else {
          tot <- if (!is.na(BASELINE_PHQ_TOTAL_COL)) BASELINE_PHQ_TOTAL_COL else head(intersect(PHQ_TOTAL_CANDIDATES, names(dp)), 1)
          if (length(tot) && tot %in% names(dp)) {
            v <- suppressWarnings(as.numeric(as.character(dp[[tot]]))); v[v < 0 | v > 27] <- NA
            dp$phq9_total <- v; derived <- c(derived, "phq9_total"); cat("  PHQ-9 total column used:", tot, "(no complete item set found)\n")
          } else {
            show_cols(dp, "the PHQ-9 file")
            cat("  PHQ-9: no complete 9-item set and no total column found. Set BASELINE_PHQ_ITEM_COLS (9 columns, suicide item last) or BASELINE_PHQ_TOTAL_COL.\n")
          }
        }
        keep <- intersect(c("uid", "phq9_total", "phq9_items1to8"), names(dp))
        if (length(keep) > 1) demo <- merge(demo, dp[, keep, drop = FALSE], by = "uid", all.x = TRUE)
      }
    } else {
      cat("PHQ-9 file not found at", BASELINE_DEP_FILE, "; baseline depression skipped.\n")
      cand <- list.files(dirname(BASELINE_DEP_FILE), pattern = "phq|depress|baseline|idas|dass", ignore.case = TRUE)
      cat("  Files in that folder that might be it:", if (length(cand)) paste(cand, collapse = ", ") else "(none matched)",
          "\n  Set BASELINE_DEP_FILE to the right file name.\n")
    }
    
    # ---- C-SSRS people with no answers are counted as No (by decision) ----
    # 150 people had a PHQ-9 item 9 of 0 and were not given the C-SSRS. A further 10 early participants have a non-zero item 9 but no
    # C-SSRS record. Both groups are scored No / "0 None"; the console reports how many fall in each so the footnote can say so.
    cs_vars <- grep("^cssrs_(ideation|max_ideation|behavior)", names(demo), value = TRUE)
    if (length(cs_vars)) {
      na_ans <- demo[[cs_vars[1]]] %in% "Not answered"
      if (exists("dp") && !is.null(items) && items[9] %in% names(dp)) {
        i9 <- suppressWarnings(as.numeric(as.character(dp[[items[9]]][match(demo$uid, dp$uid)])))
        cat(sprintf("  C-SSRS: %d people had no valid answers: PHQ-9 item 9 at screening was 0 for %d (not asked) and above 0 or unknown for %d (no C-SSRS record). All are counted as No.\n",
                    sum(na_ans), sum(na_ans & i9 %in% 0), sum(na_ans & !(i9 %in% 0))))
      } else cat(sprintf("  C-SSRS: %d people had no valid answers and are counted as No.\n", sum(na_ans)))
      for (v in cs_vars) {
        na_v <- demo[[v]] %in% "Not answered"
        demo[[v]][na_v] <- if (v == "cssrs_max_ideation_type_screen") "0 None" else "No"
      }
    }
    
    # ---- study treatment condition (treatmentgroup): demographics file, else PHQ-9 file, else C-SSRS file ----
    if (!is.na(TREATMENT_COL)) {
      got <- FALSE
      if (TREATMENT_COL %in% names(demo)) {
        demo$treatment_group <- as.character(demo[[TREATMENT_COL]]); got <- TRUE; cat("  treatment_group: from the demographics file\n")
      } else {
        srcs <- list(list(df = if (exists("dp")) dp, label = "the PHQ-9 file"), list(df = if (exists("cs")) cs, label = "the C-SSRS file"))
        for (sc in srcs) {
          if (is.null(sc$df) || !TREATMENT_COL %in% names(sc$df) || !"uid" %in% names(sc$df)) next
          v <- sc$df[!is.na(sc$df[[TREATMENT_COL]]) & nzchar(trimws(as.character(sc$df[[TREATMENT_COL]]))), c("uid", TREATMENT_COL), drop = FALSE]
          v <- v[!duplicated(v$uid), , drop = FALSE]; names(v)[2] <- "treatment_group"; v$treatment_group <- as.character(v$treatment_group)
          demo <- merge(demo, v, by = "uid", all.x = TRUE); got <- TRUE
          cat(sprintf("  treatment_group: from %s (%d participants)\n", sc$label, nrow(v))); break
        }
      }
      if (got) derived <- c(derived, "treatment_group") else cat("  treatment_group: column '", TREATMENT_COL, "' not found in any file; skipped.\n", sep = "")
    }
    
    vars <- intersect(unique(c(if (length(DEMO_VARS)) DEMO_VARS else auto_demo, derived)), names(demo))
    cat("Characteristics-table variables:", paste(vars, collapse = ", "), "\n")
    if (!length(DEMO_VARS)) { cat("  (demographic columns were chosen automatically; set DEMO_VARS to override.) "); show_cols(demo, "the demographics file") }
    if (!length(vars)) {
      cat("C5: nothing to summarise yet.\n")
    } else {
      demo$analysed <- demo$uid %in% unique(P$uid)
      order_levels <- function(v, x) {
        f <- droplevels(factor(x)); lev <- levels(f); key <- NULL
        if (v == "Income")
          key <- ifelse(grepl("^Less than", lev), -1, ifelse(grepl("or more", lev), 1e9,
                                                             suppressWarnings(as.numeric(gsub("[^0-9]", "", sub(" to .*", "", lev))))))
        if (v == "Highest_Level_Education") {
          ord <- c("Less than high", "High school", "Some college", "Trade", "Associate", "Bachelor", "Master", "Doctoral")
          key <- vapply(lev, function(l) { i <- which(vapply(ord, function(o) grepl(o, l, fixed = TRUE), logical(1))); if (length(i)) i[1] else 99 }, numeric(1))
        }
        if (grepl("^cssrs_(ideation|behavior)", v)) key <- match(lev, c("Yes", "No"))
        if (is.null(key) || anyNA(key)) return(f)
        factor(f, levels = lev[order(key)])
      }
      summ_var <- function(v) {
        x <- demo[[v]]; miss <- sum(is.na(x[demo$analysed]))
        if (is.numeric(x)) {
          a <- x[demo$analysed]; b <- x[!demo$analysed]
          data.frame(variable = v, level = "median [IQR]", analysed = fmt_iqr(a, 2), not_analysed = fmt_iqr(b, 2),
                     p = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_),
                     n_missing_analysed = miss)
        } else {
          tab <- table(factor(demo$analysed, levels = c(TRUE, FALSE)), order_levels(v, x))
          pv <- tryCatch(suppressWarnings(fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value), error = function(e) NA_real_)
          data.frame(variable = v, level = colnames(tab),
                     analysed = sprintf("%d (%.0f%%)", tab[1, ], 100 * tab[1, ] / sum(tab[1, ])),
                     not_analysed = sprintf("%d (%.0f%%)", tab[2, ], 100 * tab[2, ] / sum(tab[2, ])),
                     p = c(pv, rep(NA_real_, ncol(tab) - 1)), n_missing_analysed = miss)
        }
      }
      set.seed(SEED)                              # Fisher p values use simulation: fixed seed = same p values on every rerun
      D5 <- dplyr::bind_rows(lapply(vars, summ_var))
      D5 <- dplyr::bind_rows(data.frame(variable = "Participants", level = "n", analysed = as.character(sum(demo$analysed)), not_analysed = as.character(sum(!demo$analysed)), p = NA_real_, n_missing_analysed = NA_integer_, stringsAsFactors = FALSE), D5)
      save_csv(D5, "C5_table2_and_analysed_vs_not.csv")
      cat(sprintf("\nCharacteristics-table variables: analysed n = %d, not analysed n = %d (all other participants in the demographics file)\n",
                  sum(demo$analysed), sum(!demo$analysed)))
      print(D5, row.names = FALSE)
      save_csv(demo[demo$analysed, c("uid", vars), drop = FALSE], "C5_table2_participant_level.csv")
      # ---- C5b. the same table for the 143 vs only the participants excluded by the outcome-count rule
      if (file.exists(INCL_FILE)) {
        ir <- read.csv(INCL_FILE, colClasses = c(uid = "character"), stringsAsFactors = FALSE); Tf <- function(x) x %in% TRUE
        excl_si <- ir$uid[Tf(ir$pass_ema) & Tf(ir$pass_passive) & !Tf(ir$pass_si)]; demo$excl_si <- demo$uid %in% excl_si
        summ_var2 <- function(v) {
          x <- demo[[v]]; keep <- demo$analysed | demo$excl_si
          if (is.numeric(x)) {
            a <- x[demo$analysed]; b <- x[demo$excl_si]
            data.frame(variable = v, level = "median [IQR]", analysed = fmt_iqr(a, 2), excluded_by_outcome_rule = fmt_iqr(b, 2),
                       p = tryCatch(suppressWarnings(wilcox.test(a, b)$p.value), error = function(e) NA_real_))
          } else {
            tab <- table(factor(demo$analysed[keep], levels = c(TRUE, FALSE)), order_levels(v, x[keep]))
            pv <- tryCatch(suppressWarnings(fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value), error = function(e) NA_real_)
            data.frame(variable = v, level = colnames(tab), analysed = sprintf("%d (%.0f%%)", tab[1, ], 100 * tab[1, ] / sum(tab[1, ])),
                       excluded_by_outcome_rule = sprintf("%d (%.0f%%)", tab[2, ], 100 * tab[2, ] / sum(tab[2, ])), p = c(pv, rep(NA_real_, ncol(tab) - 1)))
          }
        }
        set.seed(SEED); D5b <- dplyr::bind_rows(lapply(vars, summ_var2)); D5b <- dplyr::bind_rows(data.frame(variable = "Participants", level = "n", analysed = as.character(sum(demo$analysed)), excluded_by_outcome_rule = as.character(sum(demo$excl_si)), p = NA_real_, stringsAsFactors = FALSE), D5b); save_csv(D5b, "C5b_analysed_vs_excluded_by_outcome_rule.csv")
        cat(sprintf("\nC5b: analysed n = %d vs excluded only by the outcome-count rule n = %d (present in the demographics file: %d)\n", sum(demo$analysed), length(excl_si), sum(demo$excl_si)))
        if (all(c("n_pos", "n_neg") %in% names(ir))) {                  # C5c: baseline clinical measures for the analysed and the two kinds of excluded participants (almost never vs almost always endorsing item 9)
          set.seed(SEED)
          ty <- ir[ir$uid %in% excl_si, c("uid", "n_pos", "n_neg")]
          ty$grp <- ifelse(ty$n_pos < 10 & ty$n_neg >= 20, "almost never", ifelse(ty$n_neg < 20 & ty$n_pos >= 10, "almost always", "both"))
          d3 <- demo; d3$grp <- ifelse(d3$analysed, "analysed", NA_character_); mm <- match(ty$uid, d3$uid); d3$grp[mm[!is.na(mm)]] <- ty$grp[!is.na(mm)]
          lv <- c("analysed", "almost never", "almost always", "both"); d3 <- d3[!is.na(d3$grp), , drop = FALSE]; g3 <- factor(d3$grp, levels = lv)
          clin <- vars[grepl("^(phq9|cssrs)", vars)]
          D5c <- dplyr::bind_rows(lapply(clin, function(v) {
            x <- d3[[v]]
            if (is.numeric(x)) {
              out <- data.frame(variable = v, level = "median [IQR]", stringsAsFactors = FALSE); for (l in lv) out[[l]] <- fmt_iqr(x[g3 == l], 2)
              k3 <- g3 %in% lv[1:3] & !is.na(x)
              out$p <- if (length(unique(g3[k3])) > 1) tryCatch(kruskal.test(x[k3], droplevels(g3[k3]))$p.value, error = function(e) NA_real_) else NA_real_
              out
            } else {
              ol <- order_levels(v, x); tab <- table(g3, ol); out <- data.frame(variable = v, level = colnames(tab), stringsAsFactors = FALSE)
              for (l in lv) out[[l]] <- sprintf("%d (%.0f%%)", tab[l, ], 100 * tab[l, ] / max(1, sum(tab[l, ])))
              k3 <- g3 %in% lv[1:3]; pv <- tryCatch(fisher.test(table(droplevels(g3[k3]), ol[k3]), simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA_real_)
              out$p <- c(pv, rep(NA_real_, ncol(tab) - 1)); out
            }
          }))
          D5c <- dplyr::bind_rows(data.frame(variable = "Participants", level = "n", analysed = as.character(sum(g3 == "analysed")), `almost never` = as.character(sum(g3 == "almost never")),
                                             `almost always` = as.character(sum(g3 == "almost always")), both = as.character(sum(g3 == "both")), p = NA_real_, check.names = FALSE, stringsAsFactors = FALSE), D5c)
          save_csv(D5c, "C5c_baseline_by_excluded_type.csv")
        }
      }
    }
  } else cat("\nC5 skipped: DEMO_FILE not found (", DEMO_FILE, "). Edit the path in the settings.\n")
  
  cat("\nAll outputs written to: ", OUT, "\n")
  
}  # end STAGE 2

# ==============================================================================
# STAGE 3: REPORTING ITEMS FOR THE MANUSCRIPT (reads the Stage 1 and Stage 2 files; nothing is refit)
#   D1  participant accounting: who is in every analysis, and why not (constant outcome, constant predictions, failed fits, skipped folds)
#   D2  prompt times and the 4-hour gap
#   D3  per-person comparison with the training-mean baseline (paired circular shift)
#   D4  the 12 models against the training-mean baseline and the previous-day mean (Supplementary Table S7)
#   D5  correlation and balanced accuracy by predictor (Supplementary Table S6)
#   D6  bootstrap 95% CIs for the Table 4 medians
#   D7  higher SI-positive cutoffs, predictions re-scored (not refit)
#   D8  higher-cutoff runs that were refit (SI_CUT / RUN_TAG): comparison with the main analysis
#   D9  source for the feature table (Supplementary Table S2)
#   D10 Figure 1
# Every block is independent: an error in one is printed and the others still run. Files go to <posthoc>/working_files/reporting.
# ==============================================================================
if (RUN_STAGE3) {
  options(width = 200)
  OUT3 <- file.path(OUT, "reporting"); dir.create(OUT3, showWarnings = FALSE, recursive = TRUE)
  ANS3 <- character(0)
  say3  <- function(...) { txt <- sprintf(...); cat(txt, "\n"); ANS3 <<- c(ANS3, txt) }
  blk3  <- function(title, expr) { cat("\n==========", title, "==========\n"); ANS3 <<- c(ANS3, "", paste0("## ", title))        # each block runs in its own environment, so its working variables do not leak
  tryCatch(eval(substitute(expr), new.env(parent = globalenv())), error = function(e) say3("ERROR in this block: %s", conditionMessage(e))) }
  save3 <- function(df, name) { write.csv(df, file.path(OUT3, name), row.names = FALSE); invisible(df) }
  iqr3  <- function(x, d = 3) { x <- x[is.finite(x)]; if (!length(x)) return(NA_character_); q <- stats::quantile(x, c(.25, .5, .75)); sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), q[2], q[1], q[3]) }
  hhmm3 <- function(h) { m <- round(h * 60); sprintf("%02d:%02d", (m %/% 60) %% 24, m %% 60) }
  auc3  <- function(y, s) { n1 <- sum(y == 1); n0 <- sum(y == 0); if (n1 == 0 || n0 == 0) return(NA_real_); r <- rank(s, ties.method = "average"); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
  binom3 <- function(k, n, p0 = 0.05) if (n > 0) binom.test(k, n, p0, alternative = "greater")$p.value else NA_real_
  rd3   <- function(f) read.csv(f, colClasses = c(uid = "character"), stringsAsFactors = FALSE, check.names = FALSE)
  need3 <- function(f) { if (!file.exists(f)) stop("file not found: ", f, ". Run Stage 1 and Stage 2 first."); f }
  cnt3  <- function(x) paste(sprintf("%s %d", names(table(x)), as.integer(table(x))), collapse = ", ")
  b2_3  <- function() { b <- rd3(need3(file.path(OUT, "B2_person_metrics.csv"))); b$key <- paste(b$Model, b$Configuration, sep = " | "); b }
  own_base3 <- function(k) paste0(sub("-.*$", "-none", sub(" \\|.*$", "", k)), " | Baseline")      # the training-mean baseline of that model's framework
  MAIN3 <- list(model = paste(COMBOS$framework[1], COMBOS$learner[1], sep = "-"), config = "EMA + Passive (dynamic)")
  P3 <- NULL                                                   # Stage 1 predictions on the common support, read once when a block needs them
  get_P3 <- function() {
    if (!is.null(P3)) return(P3)
    p <- fread(need3(file.path(NEW_DIR, "matched_predictions.csv")), select = c("uid", "fold", "row_id", "true", "true_bin", "p_si", "pred", "framework", "learner", "Configuration", "variant"),
               colClasses = list(character = "uid")) %>% as.data.frame()
    p <- p[p$variant == "adaptive", , drop = FALSE]
    p$Model <- paste0(p$framework, "-", p$learner); p$key <- paste(p$Model, p$Configuration, sep = " | "); p$uf <- paste(p$uid, p$fold)
    p$score_occ <- ifelse(p$framework == "Hurdle", p$p_si, p$pred)
    common <- Reduce(intersect, split(p$uf, p$key)); p <- p[p$uf %in% common, , drop = FALSE]      # (participant, test day) produced by every adaptive model
    P3 <<- p; p
  }
  
  # ---- D1 ----------------------------------------------------------------------
  blk3("D1. Participant accounting: who is in every analysis, and why not", {
    b2 <- b2_3(); ad <- b2[b2$variant == "adaptive", ]; shared <- sort(unique(ad$uid))
    dat_u <- tryCatch({ d <- fread(DATA_FILE, select = c("uid", "q9"), colClasses = list(character = "uid")); sort(unique(d$uid[!is.na(d$q9)])) }, error = function(e) shared)
    st <- if (file.exists(f1 <- file.path(NEW_DIR, "matched_status.csv"))) rd3(f1) else NULL
    fl <- if (file.exists(f2 <- file.path(NEW_DIR, "matched_fold_log.csv"))) rd3(f2) else NULL
    miss <- setdiff(dat_u, shared)
    say3("Analysed participants (cleaned file): %d. In every model: %d. Not in every model: %d%s.", length(dat_u), length(shared), length(miss), if (length(miss)) paste0(" (", paste(miss, collapse = ", "), ")") else "")
    for (u in miss) {
      if (!is.null(st)) {
        s <- st[st$uid == u & st$variant == "adaptive", ]; bad <- s[s$status != "ok", ]
        say3("Participant %s: %d of %d adaptive fits not ok: %s", u, nrow(bad), nrow(s),
             if (nrow(bad)) paste(unique(paste0(bad$framework, "-", bad$learner, " ", bad$Configuration, " [", bad$status, "]")), collapse = "; ") else "none")
      }
      if (!is.null(fl)) {
        f <- fl[fl$uid == u & fl$variant == "adaptive", ]
        if (nrow(f)) {
          tt <- as.data.frame(table(f$Configuration, f$reason), stringsAsFactors = FALSE); tt <- tt[tt$Freq > 0, ]
          say3("  skipped folds: %s", paste(sprintf("%s %s x %d", tt$Var1, tt$Var2, tt$Freq), collapse = "; "))
        } else say3("  no skipped-fold entries: the failure happened before any fold (status no_valid_folds or zero_outcome_variance).")
      }
    }
    ref_key <- if ("Continuous-RF | EMA + Passive (dynamic)" %in% ad$key) "Continuous-RF | EMA + Passive (dynamic)" else ad$key[1]
    one <- ad[ad$key == ref_key, ]; no_si <- one$uid[one$n_SI == 0]; all_si <- one$uid[one$n_nonSI == 0]
    say3("Outcome in the test period: no ideation at all %d (%s); ideation at every test prompt %d (%s); both classes %d.", length(no_si), paste(no_si, collapse = ", "), length(all_si), paste(all_si, collapse = ", "), sum(one$n_SI > 0 & one$n_nonSI > 0))
    say3("Fewer than 10 test days with ideation (AUC defined but imprecise): %d. Participants with fewer than 29 test days: %d (minimum %d).", sum(one$n_SI > 0 & one$n_nonSI > 0 & one$n_SI < 10), sum(one$n_folds < 29), min(one$n_folds))
    if ("n_pos_tv_med" %in% names(one)) for (u in c(no_si, all_si)) { r <- one[one$uid == u, ]; say3("  %s: median %.0f SI-positive prompts in the 55-day training window; test period %d of %d prompts with ideation.", u, r$n_pos_tv_med, r$n_SI, r$n_obs) }
    has_var <- setdiff(one$uid, no_si)
    av <- do.call(rbind, lapply(split(ad, ad$key), function(g) { cp <- is.na(g$pred_sd) | g$pred_sd == 0
    data.frame(predictor = g$key[1], n = nrow(g), AUC = sum(!is.na(g$AUC)), R2 = sum(!is.na(g$R2)), r = sum(!is.na(g$r)), MAE = sum(!is.na(g$MAE)), RMSE = sum(!is.na(g$RMSE)),
               constant_predictions = sum(cp), constant_among_outcome_variance = sum(cp & g$uid %in% has_var), stringsAsFactors = FALSE) }))
    save3(av, "D1_metric_availability_by_predictor.csv"); refs <- grepl("^Reference", av$predictor); mods <- av[!refs, ]
    say3("Models and baselines: participants with AUC %s, R2 %s, r %s, MAE %s, RMSE %s; constant predictions in %d participant-model cells.",
         paste(unique(mods$AUC), collapse = "/"), paste(unique(mods$R2), collapse = "/"), paste(unique(mods$r), collapse = "/"), paste(unique(mods$MAE), collapse = "/"), paste(unique(mods$RMSE), collapse = "/"), sum(mods$constant_predictions))
    for (i in which(refs)) say3("  %-34s constant predictions among participants with outcome variance: %d", av$predictor[i], av$constant_among_outcome_variance[i])
    if (!is.null(fl)) { f <- fl[fl$variant == "adaptive", ]; say3("Skipped participant-days in the adaptive fits: %d (%s); participants affected: %d.", nrow(f), cnt3(f$reason), length(unique(f$uid))) }
    if (!is.null(st)) { s <- st[st$variant == "adaptive" & st$status != "ok", ]; if (nrow(s)) say3("Adaptive fits with a status other than ok: %s.", cnt3(s$status)) }
    nu <- tapply(ad$n_obs, ad$uid, function(z) length(unique(z)))
    say3("Same test days: adaptive models differ in the number of test prompts for %d participants (prompts per model: %s).", sum(nu > 1), paste(unique(tapply(ad$n_obs, ad$key, sum)), collapse = ", "))
    a0 <- ad[, c("key", "uid", "n_obs")]; names(a0)[3] <- "n_adaptive"
    for (v in setdiff(unique(b2$variant), "adaptive")) {
      f <- b2[b2$variant == v, c("key", "uid", "n_obs")]; m <- merge(a0, f, by = c("key", "uid"))
      say3("  %s: participants per model %d to %d; participants whose test days differ from the adaptive model's (summed over models): %d.", v, min(table(f$key)), max(table(f$key)), sum(m$n_adaptive != m$n_obs))
    }
    fx <- b2[b2$variant == "fixed_day1" & b2$Configuration %in% c("EMA", "EMA (dynamic)"), ]; nof <- setdiff(shared, fx$uid)
    say3("Participants in every adaptive model with no 1-day EMA-only fixed-window result: %d.", length(nof))
    if (!is.null(st) && length(nof)) { s <- st[st$variant == "fixed_day1" & st$Configuration %in% c("EMA", "EMA (dynamic)") & st$uid %in% nof, ]; if (nrow(s)) say3("  their status: %s.", cnt3(s$status)) }
    if (!is.null(fl) && length(nof)) {
      f <- fl[fl$variant == "fixed_day1" & fl$Configuration %in% c("EMA", "EMA (dynamic)") & fl$uid %in% nof, ]
      if (nrow(f)) { t <- sort(table(f$reason), decreasing = TRUE); say3("  skipped-fold reasons for the 1-day EMA-only fits: %s (most common: %s).", cnt3(f$reason), names(t)[1]) } else say3("  no skipped folds logged: the failures happened before any fold; see the status above.")
    }
  })
  
  # ---- D2 ----------------------------------------------------------------------
  blk3("D2. Prompt times and the 4-hour gap", {
    d <- fread(DATA_FILE, select = c("uid", "date_rep", "ema_time", "ema_day", "q9"), colClasses = list(character = c("uid", "ema_day"))) %>% as.data.frame()
    d <- d[!is.na(d$q9), , drop = FALSE]
    d$date <- suppressWarnings(as.Date(as.character(d$date_rep))); if (all(is.na(d$date))) d$date <- suppressWarnings(as.Date(as.character(d$date_rep), format = "%m/%d/%Y"))
    d$day_index <- ave(as.integer(d$date), d$uid, FUN = function(z) z - min(z))
    d <- d[d$day_index < STUDY_DAYS, , drop = FALSE]
    keep <- tryCatch({ ri <- rd3(INCL_FILE); ri$uid[ri$pass_ema %in% TRUE & ri$pass_passive %in% TRUE & ri$pass_si %in% TRUE & ri$pass_var %in% TRUE] }, error = function(e) NULL)
    if (length(keep)) d <- d[d$uid %in% keep, , drop = FALSE]
    tt <- suppressWarnings(as.POSIXct(d$ema_day, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y/%m/%d %H:%M:%OS", "%m/%d/%Y %H:%M:%OS", "%m/%d/%Y %H:%M", "%Y-%m-%dT%H:%M:%OS")))
    if (mean(is.na(tt)) > .05) stop("could not parse ema_day; examples: ", paste(head(unique(d$ema_day), 3), collapse = " | "))
    d$t <- as.numeric(tt); d$hour <- (d$t %% 86400) / 3600
    pt <- d[!is.na(d$hour), ] %>% group_by(ema_time) %>%
      summarise(n_prompts = n(), n_people = n_distinct(uid), median = hhmm3(median(hour)), q25 = hhmm3(quantile(hour, .25)), q75 = hhmm3(quantile(hour, .75)), .groups = "drop")
    pt <- pt[order(match(pt$ema_time, c("morning", "afternoon", "evening"))), ]; print(as.data.frame(pt), row.names = FALSE); save3(pt, "D2_prompt_times.csv")
    g <- d[!is.na(d$t), ]; g <- g[order(g$uid, g$date, g$t), ]; g$gap_h <- ave(g$t, g$uid, g$date, FUN = function(z) c(NA, diff(z))) / 3600; g <- g[!is.na(g$gap_h), ]
    say3("Prompt times (answered prompts, %d participants, study days 1-%d), median [IQR] clock time: %s.", length(unique(d$uid)), STUDY_DAYS,
         paste(sprintf("%s %s [%s-%s] (%d prompts)", pt$ema_time, pt$median, pt$q25, pt$q75, pt$n_prompts), collapse = "; "))
    say3("Consecutive answered prompts on the same day were a median of %s hours apart; %.1f%% were less than 4 hours apart (the previous prompt then falls inside the 4-hour gap and is left out of the EMA windows).", iqr3(g$gap_h, 1), 100 * mean(g$gap_h < 4))
    if (any(pt$ema_time == "morning") && median(d$hour[d$ema_time == "morning"], na.rm = TRUE) > 12) say3("WARNING: the median morning prompt is after noon, so ema_day is probably stored in UTC; convert to local time before quoting these.")
    save3(data.frame(n_gaps = nrow(g), median_iqr_hours = iqr3(g$gap_h, 1), pct_lt_4h = round(100 * mean(g$gap_h < 4), 1)), "D2_within_day_gaps.csv")
  })
  
  # ---- D3 ----------------------------------------------------------------------
  # For each person, model and configuration: does the model beat that person's own training-mean baseline?
  # Statistics: MAE(baseline) - MAE(model) and AUC(model) - AUC(baseline) (> 0 = model better).
  # Null: the model's DEPARTURES from the baseline (model - baseline) are rotated circularly against the outcome, every non-zero shift in time
  # order, so the baseline and the outcome stay aligned and the autocorrelation of the departures is kept. p = (1 + #{null >= observed}) / n.
  # The p value says whether the departures are better ALIGNED with the outcome than random alignment, not that the model is better than the
  # baseline, so a person counts only if the gain is positive AND p < .05. With about 80 prompts per person the smallest p is 1 / n, so BH correction
  # across people cannot reach q < .05; report the share against the 5% expected by chance (exact binomial test).
  blk3("D3. Per-person comparison with the training-mean baseline", {
    Pc <- get_P3()
    shift_test <- function(d) {
      d <- d[order(d$row_id), , drop = FALSE]; n <- nrow(d); y <- d$true; yb <- d$true_bin
      na <- data.frame(n_obs = n, MAE_base = NA_real_, MAE_model = NA_real_, gain_MAE = NA_real_, p_MAE = NA_real_, AUC_base = NA_real_, AUC_model = NA_real_, gain_AUC = NA_real_, p_AUC = NA_real_)
      if (n < 10) return(na)
      rot <- function(k) ((seq_len(n) - 1 + k) %% n) + 1
      mae_b <- mean(abs(d$pred_b - y)); mae_m <- mean(abs(d$pred - y)); obs <- mae_b - mae_m; dl <- d$pred - d$pred_b
      nul <- vapply(seq_len(n - 1), function(k) mae_b - mean(abs(pmin(OUTCOME_MAX, pmax(0, d$pred_b + dl[rot(k)])) - y)), numeric(1))
      na$MAE_base <- mae_b; na$MAE_model <- mae_m; na$gain_MAE <- obs
      na$p_MAE <- if (all(dl == 0)) 1 else (1 + sum(nul >= obs - 1e-12)) / n
      if (any(yb == 1) && any(yb == 0)) {
        ab <- auc3(yb, d$score_b); am <- auc3(yb, d$score_occ); ds <- d$score_occ - d$score_b
        nul <- vapply(seq_len(n - 1), function(k) auc3(yb, d$score_b + ds[rot(k)]) - ab, numeric(1))
        na$AUC_base <- ab; na$AUC_model <- am; na$gain_AUC <- am - ab
        na$p_AUC <- if (all(ds == 0)) 1 else (1 + sum(nul >= am - ab - 1e-12)) / n
      }
      na
    }
    models <- intersect(paste(COMBOS$framework, COMBOS$learner, sep = "-"), unique(Pc$Model))
    res <- bind_rows(lapply(models, function(m) {
      fw <- sub("-.*$", "", m)
      b <- Pc[Pc$key == paste0(fw, "-none | Baseline"), c("uid", "row_id", "pred", "score_occ")]; names(b)[3:4] <- c("pred_b", "score_b")
      bind_rows(lapply(setdiff(unique(Pc$Configuration[Pc$Model == m]), "Baseline"), function(cfg) {
        a <- Pc[Pc$Model == m & Pc$Configuration == cfg, c("uid", "row_id", "true", "true_bin", "pred", "score_occ")]
        w <- merge(a, b, by = c("uid", "row_id"))
        bind_rows(lapply(split(w, w$uid), function(d) cbind(data.frame(uid = d$uid[1], Model = m, Configuration = cfg), shift_test(d))))
      }))
    }))
    if (!nrow(res)) stop("no model and baseline pairs found in matched_predictions.csv")
    res <- res %>% group_by(Model, Configuration) %>% mutate(q_MAE = p.adjust(p_MAE, "BH"), q_AUC = p.adjust(p_AUC, "BH")) %>% ungroup()
    save3(res, "D3_person_vs_baseline.csv")
    better_p <- function(gain, p) !is.na(gain) & !is.na(p) & gain > 0 & p < .05
    S <- res %>% group_by(Model, Configuration) %>%
      summarise(n_people = n(), MAE_gain = iqr3(gain_MAE, 2), pct_better_MAE = round(100 * mean(gain_MAE > 0, na.rm = TRUE), 1),
                n_sig_better_MAE = sum(better_p(gain_MAE, p_MAE)), pct_sig_better_MAE = round(100 * mean(better_p(gain_MAE, p_MAE)), 1),
                n_q05_better_MAE = sum(!is.na(q_MAE) & q_MAE < .05 & gain_MAE > 0),
                n_AUC = sum(!is.na(gain_AUC)), AUC_gain = iqr3(gain_AUC, 3), pct_better_AUC = round(100 * mean(gain_AUC > 0, na.rm = TRUE), 1),
                n_sig_better_AUC = sum(better_p(gain_AUC, p_AUC)), pct_sig_better_AUC = round(100 * sum(better_p(gain_AUC, p_AUC)) / pmax(1, n_AUC), 1),
                n_q05_better_AUC = sum(!is.na(q_AUC) & q_AUC < .05 & gain_AUC > 0), .groups = "drop") %>%
      mutate(binom_p_MAE = mapply(binom3, n_sig_better_MAE, n_people), binom_p_AUC = mapply(binom3, n_sig_better_AUC, n_AUC))
    print(as.data.frame(S), row.names = FALSE); save3(S, "D3_person_vs_baseline_summary.csv")
    r <- S[S$Model == MAIN3$model & S$Configuration == MAIN3$config, ]
    if (nrow(r)) {
      say3("%s, %s vs each person's own baseline: MAE was lower than the baseline for %.1f%% of %d people (median [IQR] gain %s points); %d people (%.1f%%) were both better and beyond the circular-shift null at p < .05 (5%% expected by chance; binomial p = %s); %d after BH correction, which cannot reach q < .05 at this resolution.",
           MAIN3$model, MAIN3$config, r$pct_better_MAE, r$n_people, r$MAE_gain, r$n_sig_better_MAE, r$pct_sig_better_MAE, format(signif(r$binom_p_MAE, 2)), r$n_q05_better_MAE)
      say3("For AUC (n = %d with both classes): higher than the baseline for %.1f%%; %d people (%.1f%%) better and p < .05 (binomial p = %s); %d after BH.",
           r$n_AUC, r$pct_better_AUC, r$n_sig_better_AUC, r$pct_sig_better_AUC, format(signif(r$binom_p_AUC, 2)), r$n_q05_better_AUC)
    }
    say3("Methods wording: 'For each person we tested whether the model's departures from the person's training-mean baseline tracked the outcome, by rotating those departures circularly against the outcome (all non-zero shifts) and recomputing the gain over the baseline; this keeps the autocorrelation of the series. A person counted only if the gain was positive and p < .05.'")
  })
  
  # ---- D4 ----------------------------------------------------------------------
  blk3("D4. The models against the training-mean baseline and the previous-day mean (Supplementary Table S7)", {
    b2 <- b2_3(); ad <- b2[b2$variant == "adaptive", ]; set.seed(SEED)
    keys_m <- sort(unique(ad$key[!grepl("-none \\|", ad$key)])); prev_key <- "Reference-none | MeanD1"
    paired3 <- function(a, b, metric, nboot = 5000) {
      m <- merge(ad[ad$key == a, c("uid", metric)], ad[ad$key == b, c("uid", metric)], by = "uid", suffixes = c(".a", ".b"))
      d <- m[[paste0(metric, ".a")]] - m[[paste0(metric, ".b")]]; d <- d[is.finite(d)]
      if (length(d) < 5) return(data.frame(n = length(d), median_diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_))
      bt <- replicate(nboot, median(sample(d, replace = TRUE)))
      data.frame(n = length(d), median_diff = median(d), ci_lo = unname(quantile(bt, .025)), ci_hi = unname(quantile(bt, .975)), p = suppressWarnings(wilcox.test(d)$p.value))
    }
    rows <- list()
    for (m in c("AUC", "R2", "r", "MAE", "RMSE")) for (k in keys_m) for (bn in c("training-mean baseline", "previous-day mean")) {
      bk <- if (bn == "previous-day mean") prev_key else own_base3(k)
      rows[[length(rows) + 1]] <- cbind(data.frame(metric = m, key = k, benchmark = bn, stringsAsFactors = FALSE), paired3(k, bk, m))
    }
    T7 <- bind_rows(rows) %>% group_by(metric, benchmark) %>% mutate(q_BH_over_models = p.adjust(p, "BH")) %>% ungroup()
    T7$better <- ifelse(T7$metric %in% c("MAE", "RMSE"), T7$median_diff < 0, T7$median_diff > 0)
    save3(T7, "D4_all_models_vs_benchmarks.csv")
    both <- T7 %>% group_by(metric, key) %>% summarise(better_than_both = sum(better, na.rm = TRUE) == 2, sig_better_than_both = sum(better & q_BH_over_models < .05, na.rm = TRUE) == 2, .groups = "drop")
    say3("Across the %d models and 5 metrics, %d of %d model-by-metric cells are better than BOTH the training-mean baseline and the previous-day mean (point estimate); %d at q < .05.", length(keys_m), sum(both$better_than_both), nrow(both), sum(both$sig_better_than_both))
    S7 <- T7[T7$key %in% paste(c(MAIN3$model, paste(COMBOS$framework[2], COMBOS$learner[2], sep = "-")), MAIN3$config, sep = " | "), ]
    S7 <- S7 %>% group_by(metric) %>% mutate(q_BH = p.adjust(p, "BH")) %>% ungroup()      # BH across the four comparisons of each metric, as in Supplementary Table S7
    save3(S7, "S7_combined_dynamic_vs_benchmarks.csv"); print(as.data.frame(S7[, c("metric", "key", "benchmark", "n", "median_diff", "ci_lo", "ci_hi", "q_BH")]), row.names = FALSE)
  })
  
  # ---- D5 ----------------------------------------------------------------------
  blk3("D5. Correlation and balanced accuracy by predictor (Supplementary Table S6)", {
    s <- rd3(need3(file.path(OUT, "B3_summary.csv"))); s <- s[s$variant == "adaptive", intersect(c("Model", "Configuration", "n_people", "r", "Sens", "Spec", "BalAcc"), names(s))]
    save3(s, "S6_r_and_balanced_accuracy.csv"); print(s, row.names = FALSE)
    say3("Balanced accuracy of the hurdle models: %s.", paste(unique(s$BalAcc[!is.na(s$BalAcc)]), collapse = "; "))
  })
  
  # ---- D6 ----------------------------------------------------------------------
  blk3("D6. Bootstrap 95% CIs for the Table 4 medians", {
    b2 <- b2_3(); ad <- b2[b2$variant == "adaptive", ]; set.seed(SEED)
    ci <- bind_rows(lapply(split(ad, ad$key), function(g) bind_rows(lapply(c("AUC", "MAE", "R2", "r", "RMSE"), function(m) {
      x <- g[[m]]; x <- x[is.finite(x)]
      if (length(x) < 5) return(data.frame(key = g$key[1], metric = m, n = length(x), median = NA_real_, lo = NA_real_, hi = NA_real_))
      bt <- replicate(5000, median(sample(x, replace = TRUE)))
      data.frame(key = g$key[1], metric = m, n = length(x), median = median(x), lo = unname(quantile(bt, .025)), hi = unname(quantile(bt, .975)))
    }))))
    save3(ci, "T4_median_CIs.csv"); say3("Bootstrap 95%% CIs (5,000 resamples of participants) for the median AUC, MAE, R2, r and RMSE of %d predictors written to T4_median_CIs.csv.", length(unique(ci$key)))
  })
  
  # ---- D7 ----------------------------------------------------------------------
  blk3("D7. Higher SI-positive cutoffs: saved predictions re-scored (NOT refit)", {
    Pc <- get_P3(); cuts <- c(0, 2, 5, 10); main_key <- paste(MAIN3$model, MAIN3$config, sep = " | ")
    pp <- bind_rows(lapply(cuts, function(ct) { d <- Pc; d$yb <- as.integer(d$true > ct)
    d %>% group_by(key, uid) %>% summarise(AUC = auc3(yb, score_occ), .groups = "drop") %>% mutate(cutoff = ct) }))
    save3(pp, "D7_person_AUC_by_cutoff.csv")                       # person-level AUC at every cutoff (Supplementary Table S9 is built from this in Stage 4)
    res <- pp %>% group_by(key, cutoff) %>% summarise(n_people = n(), n_AUC = sum(!is.na(AUC)), median_IQR_AUC = iqr3(AUC), .groups = "drop")
    one <- Pc[Pc$key == main_key, ]
    say3("Share of test prompts with item 9 above 2, 5 and 10: %s.", paste(sprintf("%.1f%%", 100 * sapply(c(2, 5, 10), function(ct) mean(one$true > ct))), collapse = ", "))
    r <- res[res$key == main_key, ]; say3("%s, %s, median [IQR] AUC when re-scored at item 9 > %s: %s (participants with an AUC: %s).", MAIN3$model, MAIN3$config, paste(r$cutoff, collapse = ", "), paste(r$median_IQR_AUC, collapse = "; "), paste(r$n_AUC, collapse = ", "))
    say3("This is a re-scoring: continuous models are unaffected by the cutoff, but the hurdle classifier stays trained at item 9 > 0. A full redo uses SI_CUT and RUN_TAG (block D8).")
  })
  
  # ---- D8 ----------------------------------------------------------------------
  blk3("D8. Refit runs at higher cutoffs (SI_CUT, RUN_TAG): comparison with the main analysis", {
    res_dir <- dirname(dirname(OUT)); cuts <- c(0, 2, 5, 10); dirs <- c("posthoc", paste0("posthoc_cut", c(2, 5, 10)))
    ok <- file.exists(file.path(res_dir, dirs, "working_files", "B2_person_metrics.csv")) & file.exists(file.path(res_dir, dirs, "working_files", "B3_summary.csv")) & file.exists(file.path(res_dir, dirs, "working_files", "B4_contrasts.csv"))
    if (!ok[1] || sum(ok) < 2) say3("SKIPPED: needs the main run and at least one refit run. Run this script again with SI_CUT <- 2, 5 or 10 and RUN_TAG <- \"_cut2\", \"_cut5\", \"_cut10\" (RUN_FIXED_WINDOWS and RUN_PART_A off), then rerun Stage 3.") else {
      med3 <- function(x) as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", x)); tabs <- list(); conc <- character(0)
      for (i in which(ok)) {
        d <- file.path(res_dir, dirs[i], "working_files"); s <- rd3(file.path(d, "B3_summary.csv")); s <- s[s$variant == "adaptive", ]; s$key <- paste(s$Model, s$Configuration, sep = " | ")
        p <- rd3(file.path(d, "B2_person_metrics.csv")); p <- p[p$variant == "adaptive", ]; p$key <- paste(p$Model, p$Configuration, sep = " | ")
        na <- tapply(!is.na(p$AUC), p$key, sum)
        tabs[[as.character(cuts[i])]] <- data.frame(key = s$key, AUC_median = med3(s$AUC), AUC = s$AUC, n_AUC = as.integer(na[s$key]), stringsAsFactors = FALSE)
        names(tabs[[as.character(cuts[i])]])[-1] <- paste0(names(tabs[[as.character(cuts[i])]])[-1], "_cut", cuts[i])
      }
      W <- Reduce(function(a, b) merge(a, b, by = "key", all = TRUE), tabs); save3(W, "D8_refit_cutoffs_AUC_by_model.csv")
      for (ct in cuts[ok][-1]) {
        j <- complete.cases(W[[paste0("AUC_median_cut", 0)]], W[[paste0("AUC_median_cut", ct)]]); rho <- suppressWarnings(cor(W[[paste0("AUC_median_cut", 0)]][j], W[[paste0("AUC_median_cut", ct)]][j], method = "spearman"))
        say3("Cutoff %d: rank order of the %d predictors vs the main analysis, Spearman rho = %.2f; participants with an AUC per model: median %d (main: %d).", ct, sum(j), rho, as.integer(median(W[[paste0("n_AUC_cut", ct)]], na.rm = TRUE)), as.integer(median(W[["n_AUC_cut0"]], na.rm = TRUE)))
      }
      b4 <- lapply(which(ok), function(i) { x <- rd3(file.path(res_dir, dirs[i], "working_files", "B4_contrasts.csv")); x <- x[x$metric == "AUC" & x$family == "vs_baseline", c("a", "median_diff", "q_BH")]; x$cut <- cuts[i]; x })
      b4 <- bind_rows(b4); w4 <- reshape(b4, idvar = "a", timevar = "cut", direction = "wide"); save3(w4, "D8_refit_cutoffs_AUC_vs_baseline.csv")
      for (ct in cuts[ok][-1]) { s0 <- w4[[paste0("q_BH.", 0)]] < .05; sc <- w4[[paste0("q_BH.", ct)]] < .05; sg <- sign(w4[[paste0("median_diff.", 0)]]) != sign(w4[[paste0("median_diff.", ct)]])
      say3("AUC vs the training-mean baseline, cutoff %d: %d comparisons; significance (q < .05) changed in %d; direction changed in %d.", ct, nrow(w4), sum(s0 != sc, na.rm = TRUE), sum(sg, na.rm = TRUE)) }
      say3("If the rank order is stable, the differences keep their direction and few significance decisions change, the AUC findings are robust to the cutoff; otherwise say so. R2, r, MAE and RMSE of the continuous models do not depend on the cutoff; the hurdle models are refit, so their error metrics do.")
    }
  })
  
  # ---- D9 ----------------------------------------------------------------------
  blk3("D9. Source for the feature table (Supplementary Table S2)", {
    cc <- col_class; cc$stream <- ifelse(cc$source == "passive", coalesce(stream_of(cc$col), "other"), NA_character_)
    cc <- cc %>% mutate(
      stream    = ifelse(source == "ema", "EMA", stream),
      measure   = case_when(source == "ema" ~ str_extract(col, "^(q[1-9]|g[12])"),
                            stream == "sleep" ~ str_extract(col, "(total_sleep|rem_sleep)"),
                            TRUE ~ str_remove(col, paste0("^", stream, "_")) %>% str_remove("_?(mean_val|sum_val|mad_val|entropy|p90p10)$") %>%
                              str_replace("(^|_)(hr1|hr4|day1|day7)(_|$)", "_") %>% str_replace_all("_+", "_") %>% str_remove_all("^_|_$")),
      measure   = ifelse(is.na(measure) | measure == "", "(total)", measure),
      statistic = case_when(source == "ema" | stream == "sleep" ~ ifelse(grepl("_mad_", col), "MAD", "mean"), TRUE ~ str_extract(col, "(mean_val|sum_val|mad_val|entropy|p90p10)$")),
      statistic = recode(statistic, mean_val = "mean", sum_val = "sum", mad_val = "MAD", entropy = "sample entropy", p90p10 = "P90-P10"))
    long <- cc %>% transmute(source, stream, measure, statistic, type, window, column = col) %>% arrange(source, stream, measure, statistic, window)
    save3(long, "D9_S2_feature_table_long.csv")
    wide <- long %>% group_by(source, stream, measure, statistic, type) %>% summarise(windows = paste(sort(unique(window)), collapse = ", "), n_columns = n(), .groups = "drop"); save3(wide, "D9_S2_feature_table_summary.csv")
    say3("Feature table source written: %d candidate features (%d EMA, %d passive) in %d rows of stream x measure x statistic. Add a plain-language description of each passive measure to finish the table.", nrow(long), sum(long$source == "ema"), sum(long$source == "passive"), nrow(wide))
  })
  
  # ---- D10 ---------------------------------------------------------------------
  blk3("D10. Figure 3: within-person AUC of every predictor", {
    b2 <- b2_3(); ad <- b2[b2$variant == "adaptive", ]
    GREY <- "#7F7F7F"; ORANGE <- "#E69F00"; BLUE <- "#0072B2"; GREEN <- "#009E73"
    cfgs <- c("EMA", "EMA (dynamic)", "Passive", "Passive (dynamic)", "EMA + Passive", "EMA + Passive (dynamic)"); labs <- c("EMA", "EMA (dynamic)", "passive", "passive (dynamic)", "EMA + passive", "EMA + passive (dynamic)")
    cats <- rbind(
      data.frame(label = c("Training-mean baseline (continuous)", "Training-mean baseline (hurdle)", "Persistence", "Previous prompt", "Previous 1-day mean", "Previous 7-day mean", "Always zero"),
                 key = c("Continuous-none | Baseline", "Hurdle-none | Baseline", "Reference-none | Persistence", "Reference-none | PrevPrompt", "Reference-none | MeanD1", "Reference-none | MeanD7", "Reference-none | AlwaysZero"),
                 col = c(GREY, GREY, rep(ORANGE, 5)), stringsAsFactors = FALSE),
      data.frame(label = paste0("RF, ", labs), key = paste0("Continuous-RF | ", cfgs), col = BLUE, stringsAsFactors = FALSE),
      data.frame(label = paste0("GBM hurdle, ", labs), key = paste0("Hurdle-GBM | ", cfgs), col = GREEN, stringsAsFactors = FALSE))
    cats <- cats[cats$key %in% ad$key, ]; K <- nrow(cats)
    main_key <- paste(MAIN3$model, MAIN3$config, sep = " | "); base_key <- own_base3(main_key)
    rf <- ad[ad$key == main_key, c("uid", "AUC", "n_SI")]; bs <- ad[ad$key == base_key, c("uid", "AUC")]; names(bs)[2] <- "AUC_base"
    d <- merge(rf, bs, by = "uid"); d <- d[is.finite(d$AUC) & is.finite(d$AUC_base), ]; d <- d[order(d$AUC), ]; x <- seq_len(nrow(d))
    opened <- FALSE
    tryCatch({
      png(file.path(OUT_TABLES, "Figure3_within_person_auc.png"), width = 6, height = 7, units = "in", res = 300); opened <- TRUE
      layout(matrix(1:3, ncol = 1), heights = c(3.9, 1.95, 0.8)); par(cex = 0.62, mgp = c(1.8, 0.5, 0), tcl = -0.25); set.seed(SEED)
      par(mar = c(3, 13, 1.5, 4.5))                                                                             # panel A: every predictor
      plot(NA, xlim = c(0, 1), ylim = c(K + 0.5, 0.35), yaxt = "n", xlab = "Within-person AUC", ylab = "", bty = "n")
      axis(2, at = seq_len(K), labels = cats$label, las = 1, lwd = 0.6); abline(h = c(2.5, 7.5, 13.5), col = "#BBBBBB", lwd = 0.5)
      for (k in seq_len(K)) {
        v <- ad$AUC[ad$key == cats$key[k]]; v <- v[is.finite(v)]; if (!length(v)) next
        points(v, k + runif(length(v), -0.3, 0.3), pch = 16, cex = 0.5, col = adjustcolor(cats$col[k], alpha.f = 0.4))
        q <- quantile(v, c(.05, .25, .5, .75, .95)); rect(q[2], k - 0.33, q[4], k + 0.33, border = "black", lwd = 0.7)
        segments(q[3], k - 0.33, q[3], k + 0.33, lwd = 2); segments(q[1], k, q[2], k, lwd = 0.6); segments(q[4], k, q[5], k, lwd = 0.6)
        text(1.01, k, sprintf("%.3f", q[3]), adj = 0, xpd = NA)
      }
      abline(v = 0.5, lty = 2, lwd = 0.7); text(1.01, 0, "Median", adj = 0, xpd = NA, font = 2)
      par(mar = c(0.4, 13, 3.2, 4.5))                                                                           # panel B: participants ordered by the main model's AUC
      plot(NA, xlim = c(0, nrow(d) + 1), ylim = c(0, 1.22), xaxt = "n", yaxt = "n", xlab = "", ylab = "Within-person AUC", bty = "n"); axis(2, at = seq(0, 1, .2), lwd = 0.6, las = 1)
      segments(x, pmin(d$AUC, d$AUC_base), x, pmax(d$AUC, d$AUC_base), col = "#BBBBBB", lwd = 0.5); abline(h = 0.5, lty = 2, lwd = 0.7)
      points(x, d$AUC_base, pch = 1, col = GREY, cex = 0.55); points(x, d$AUC, pch = 16, col = BLUE, cex = 0.5)
      legend("topleft", legend = c("Combined dynamic random forest", "Training-mean baseline, same participant", sprintf("Model above baseline for %d of %d participants", sum(d$AUC > d$AUC_base), nrow(d))),
             pch = c(16, 1, NA), col = c(BLUE, GREY, NA), bty = "n", cex = 0.95)
      par(mar = c(3.2, 13, 0.4, 4.5))                                                                           # panel C: days with ideation
      plot(NA, xlim = c(0, nrow(d) + 1), ylim = c(0, max(d$n_SI) * 1.05), xaxt = "n", xlab = sprintf("Participants, ordered by model AUC (n = %d)", nrow(d)), ylab = "Days with ideation", bty = "n", las = 1)
      rect(x - 0.45, 0, x + 0.45, d$n_SI, col = "#9AA5B1", border = NA); axis(1, at = unique(c(1, 25, 50, 75, 100, nrow(d))), lwd = 0.6)
      dev.off(); opened <- FALSE
    }, error = function(e) { if (opened) dev.off(); stop(e) })
    say3("Figure 3 written to %s (%d participants in panels B and C). The manuscript figure was drawn with matplotlib from the same B2 file; this R version has the same content but may differ in fonts and spacing.", file.path(OUT_TABLES, "Figure3_within_person_auc.png"), nrow(d))
  })
  
  writeLines(c(sprintf("Reporting items, run %s", format(Sys.time(), "%Y-%m-%d %H:%M")), ANS3), file.path(OUT_TABLES, "reporting_answers.txt"))
  cat("\nStage 3 working files:", OUT3, "\nreporting_answers.txt and Figure 3:", OUT_TABLES, "\n")
}  # end STAGE 3


# ==============================================================================
# STAGE 4: ONE CSV PER MANUSCRIPT TABLE (reads the working files of Stages 1 to 3; nothing is refit)
#   Table 1 sensor availability        Table 3 characteristics            Table 4 performance        Table 5 H1 and H2
#   Table 6 AUC by window              Table 7 MAE by outcome status      Table 8 exploratory comparisons
#   S2 feature table                   S3 candidate predictors            S6 correlation and balanced accuracy
#   S7 models vs benchmarks            S8 per-person comparison           S9 higher cutoffs
#   S10 analysed vs the 82 excluded    S11 baseline measures by excluded type
# Table 2 (model configurations), Supplement 1 (EMA wording) and Table 9 are text or come from another script: no file here.
# The files go to 4. Results/tables/ (the only CSVs you need). Everything else is working files (see KEEP_WORKING_FILES).
# Numbers are written plainly (ASCII minus sign, no thousands separators); copy them into the manuscript tables as they are.
# ==============================================================================
if (RUN_STAGE4) {
  options(width = 200)
  if (!exists("OUT3")) OUT3 <- file.path(OUT, "reporting")
  tabs4 <- character(0)
  blk4 <- function(title, expr) { cat("\n==========", title, "==========\n"); tryCatch(eval(substitute(expr), new.env(parent = globalenv())), error = function(e) cat("ERROR in this table:", conditionMessage(e), "\n")) }
  w4   <- function(df, name) { write.csv(df, file.path(OUT_TABLES, name), row.names = FALSE, na = ""); tabs4 <<- c(tabs4, name); cat(sprintf("written: %s (%d rows)\n", name, nrow(df))) }
  rd4  <- function(f, dir = OUT) { p <- file.path(dir, f); if (!file.exists(p)) stop("working file not found: ", p, " (run Stages 1 to 3 first)"); read.csv(p, colClasses = c(uid = "character"), stringsAsFactors = FALSE, check.names = FALSE) }
  sg   <- function(x, d = 3) ifelse(is.na(x), "", sprintf(paste0("%+.", d, "f"), x))                       # +0.074 or -0.079
  nm   <- function(x, d = 3) ifelse(is.na(x), "", sprintf(paste0("%.", d, "f"), x))
  fq   <- function(q) ifelse(is.na(q), "", ifelse(q < .001, "<.001", sub("^0", "", sprintf("%.3f", q))))     # .373 or <.001
  mci  <- function(m, lo, hi, d = 3) ifelse(is.na(m), "", sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), m, lo, hi))
  first_num <- function(s) suppressWarnings(as.numeric(sub("^\\s*(-?[0-9.]+).*$", "\\1", s)))
  lowp <- function(x) sub("Passive", "passive", x, fixed = TRUE)
  CFG  <- c("EMA", "EMA (dynamic)", "Passive", "Passive (dynamic)", "EMA + Passive", "EMA + Passive (dynamic)")
  SPEC4 <- rbind(
    data.frame(section = "Baseline and reference predictors",
               label = c("Training-mean baseline (continuous)", "Training-mean baseline (hurdle)", "Persistence", "Previous prompt", "Previous 1-day mean", "Previous 7-day mean", "Always zero"),
               key = c("Continuous-none | Baseline", "Hurdle-none | Baseline", "Reference-none | Persistence", "Reference-none | PrevPrompt", "Reference-none | MeanD1", "Reference-none | MeanD7", "Reference-none | AlwaysZero"), stringsAsFactors = FALSE),
    data.frame(section = "Continuous models (random forest)", label = CFG, key = paste0("Continuous-RF | ", CFG), stringsAsFactors = FALSE),
    data.frame(section = "Hurdle models (gradient boosting)", label = CFG, key = paste0("Hurdle-GBM | ", CFG), stringsAsFactors = FALSE))
  with_sections <- function(spec, rowfun) {                      # one block of rows per section, each section with a title row, as in the manuscript tables
    out <- list(); last <- ""
    for (i in seq_len(nrow(spec))) {
      r <- rowfun(i)
      if (spec$section[i] != last) { b <- r; b[1, ] <- ""; b[1, 1] <- spec$section[i]; out[[length(out) + 1]] <- b; last <- spec$section[i] }
      out[[length(out) + 1]] <- r
    }
    do.call(rbind, out)
  }
  b3_adapt <- function() { b <- rd4("B3_summary.csv"); b$key <- paste(b$Model, b$Configuration, sep = " | "); b }
  
  blk4("Table 1. Sensor modalities and how many participants have data", {
    av <- rd4("C4_stream_availability.csv"); n_all <- length(unique(DAT$uid))
    L <- data.frame(stream = c("acc", "steps", "hr", "battery", "convo", "sms", "calls", "gps", "gps_home", "sleep", "sleep"),
                    sensor = c("Accelerometer (Garmin)", "Accelerometer (Garmin)", "Optical PPG (Garmin)", "Battery (phone)", "Microphone (phone)", "SMS (phone)", "Calls (phone)", "GPS (phone)", "GPS (phone)", "Garmin sleep algorithm", "Garmin sleep algorithm"),
                    feature = c("Vector magnitude", "Step count", "Heart rate (30-220 bpm)", "Battery percentage", "Total seconds of conversation", "SMS message body lengths per epoch", "Total call duration per epoch", "Distance traveled", "Time within 40 m of home", "Daily total sleep duration", "REM sleep duration"),
                    epoch = c("10 s", "30 s", "10 s", "1 h", "5 min", "30 min", "30 min", "30 min", "30 min", "Daily", "Daily"),
                    scales = c(rep("1 h, 4 h, 1 d, 1 wk", 9), "1 d, 7 d", "1 d, 7 d"), stringsAsFactors = FALSE)
    m <- match(L$stream, av$stream)
    out <- data.frame(`Sensor / data stream` = L$sensor, `Feature / signal` = L$feature, `Analysis epoch` = L$epoch, `Time-scales` = L$scales,
                      `N with data` = av$n_people_with_data[m], `N with logged events (SMS, calls, conversation)` = av$n_people_with_logged_events[m], check.names = FALSE, stringsAsFactors = FALSE)
    out <- rbind(out, data.frame(`Sensor / data stream` = "EMA items", `Feature / signal` = "11 items", `Analysis epoch` = "Per assessment", `Time-scales` = "1 d, 1 wk", `N with data` = n_all, `N with logged events (SMS, calls, conversation)` = NA, check.names = FALSE, stringsAsFactors = FALSE))
    names(out)[5] <- sprintf("N (of %d)", n_all); w4(out, "Table1_sensor_availability.csv")
  })
  
  blk4("Table 3. Characteristics of analysed and not analysed participants", {
    d <- rd4("C5_table2_and_analysed_vs_not.csv"); d$p <- fq(as.numeric(d$p))
    names(d)[names(d) == "variable"] <- "Characteristic"; names(d)[names(d) == "level"] <- "Level"; names(d)[names(d) == "analysed"] <- "Analysed"; names(d)[names(d) == "not_analysed"] <- "Not analysed"
    w4(d, "Table3_characteristics.csv")
  })
  
  blk4("Table 4. Within-person performance", {
    b3 <- b3_adapt(); b3 <- b3[b3$variant == "adaptive", ]; ci <- read.csv(file.path(OUT3, "T4_median_CIs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    one <- function(k, m, d) { r <- ci[ci$key == k & ci$metric == m, ]; if (nrow(r) && !is.na(r$median[1])) mci(r$median[1], r$lo[1], r$hi[1], d) else "" }
    w4(with_sections(SPEC4, function(i) { k <- SPEC4$key[i]; r <- b3[b3$key == k, ]
    data.frame(Predictor = SPEC4$label[i], n = if (nrow(r)) as.character(r$n_people[1]) else "", AUC = one(k, "AUC", 3), `% > 0.5` = if (nrow(r)) as.character(r$pct_AUC_gt_half[1]) else "",
               MAE = one(k, "MAE", 2), `R2` = one(k, "R2", 3), check.names = FALSE, stringsAsFactors = FALSE) }), "Table4_performance.csv")
  })
  
  blk4("Table 5. Preregistered hypotheses H1 and H2", {
    b4 <- rd4("B4_contrasts.csv"); b4 <- b4[!is.na(b4$osf_hypothesis) & nzchar(b4$osf_hypothesis), ]
    rf <- function(cfg) paste0("Continuous-RF | ", cfg, " | adaptive"); base <- "Continuous-none | Baseline | adaptive"
    spec <- data.frame(h = c(rep("H1", 3), rep("H2", 2)), a = rf(c("EMA (dynamic)", "Passive (dynamic)", "EMA + Passive (dynamic)", "EMA + Passive (dynamic)", "EMA + Passive (dynamic)")),
                       b = c(rep(base, 3), rf("EMA (dynamic)"), rf("Passive (dynamic)")),
                       label = c("EMA (dynamic) vs baseline", "Passive (dynamic) vs baseline", "EMA + Passive (dynamic) vs baseline", "EMA + Passive (dynamic) vs EMA (dynamic)", "EMA + Passive (dynamic) vs Passive (dynamic)"), stringsAsFactors = FALSE)
    rows <- list()
    for (i in seq_len(nrow(spec))) for (m in c("R2", "r", "MAE", "RMSE")) {
      x <- b4[b4$a == spec$a[i] & b4$b == spec$b[i] & b4$metric == m, ]
      rows[[length(rows) + 1]] <- data.frame(H = spec$h[i], Comparison = spec$label[i], Metric = ifelse(m == "R2", "R2", m),
                                             `Median difference [95% CI]` = if (nrow(x)) paste0(sg(x$median_diff[1]), " [", nm(x$ci_lo[1]), ", ", nm(x$ci_hi[1]), "]") else "", q = if (nrow(x)) fq(x$q_confirmatory[1]) else "", check.names = FALSE, stringsAsFactors = FALSE)
    }
    w4(do.call(rbind, rows), "Table5_H1_H2.csv")
  })
  
  blk4("Table 6. Median AUC by window", {
    b3 <- b3_adapt(); val <- function(k, v) { r <- b3[b3$key == k & b3$variant == v, ]; if (!nrow(r)) return(NA_character_); nm(first_num(r$AUC[1])) }
    npe <- function(k, v) { r <- b3[b3$key == k & b3$variant == v, ]; if (nrow(r)) r$n_people[1] else NA }
    sp <- SPEC4[SPEC4$section != "Baseline and reference predictors", ]
    w4(with_sections(sp, function(i) { k <- sp$key[i]; ema_only <- sp$label[i] %in% c("EMA", "EMA (dynamic)")
    cell <- function(v) { x <- val(k, v); if (is.na(x)) return("-"); if (v == "fixed_day1" && !is.na(npe(k, v)) && npe(k, v) < npe(k, "adaptive")) paste0(x, "\u2020") else x }       # dagger: fitted for fewer participants
    data.frame(Configuration = sp$label[i], Adaptive = nm(first_num(b3$AUC[b3$key == k & b3$variant == "adaptive"][1])),
               `1 hour` = if (ema_only) "-" else cell("fixed_hr1"), `4 hours` = if (ema_only) "-" else cell("fixed_hr4"), `1 day` = cell("fixed_day1"), `7 days` = cell("fixed_day7"), check.names = FALSE, stringsAsFactors = FALSE) }), "Table6_AUC_by_window.csv")
  })
  
  blk4("Table 7. Mean absolute error overall, on days without ideation, and on days with ideation", {
    b3 <- b3_adapt(); b3 <- b3[b3$variant == "adaptive", ]; sp <- SPEC4[SPEC4$key != "Hurdle-none | Baseline", ]; sp$label[sp$key == "Continuous-none | Baseline"] <- "Training-mean baseline"
    w4(with_sections(sp, function(i) { r <- b3[b3$key == sp$key[i], ]
    data.frame(Predictor = sp$label[i], Overall = r$MAE[1], `Days without ideation` = r$MAE_zero[1], `Days with ideation` = r$MAE_pos[1], check.names = FALSE, stringsAsFactors = FALSE) }), "Table7_MAE_by_outcome_status.csv")
  })
  
  blk4("Table 8. Exploratory comparisons", {
    b4 <- rd4("B4_contrasts.csv")
    cmp <- data.frame(label = c("Combined vs EMA-only (mean-only)", "Combined vs passive-only (mean-only)", "Combined vs EMA-only (dynamic)", "Combined vs passive-only (dynamic)",
                                "EMA: dynamic vs mean-only", "Passive: dynamic vs mean-only", "Combined: dynamic vs mean-only"),
                      a = c("EMA + Passive", "EMA + Passive", "EMA + Passive (dynamic)", "EMA + Passive (dynamic)", "EMA (dynamic)", "Passive (dynamic)", "EMA + Passive (dynamic)"),
                      b = c("EMA", "Passive", "EMA (dynamic)", "Passive (dynamic)", "EMA", "Passive", "EMA + Passive"), stringsAsFactors = FALSE)
    rows <- list()
    for (i in seq_len(nrow(cmp))) for (lr in list(c("RF", "Continuous-RF"), c("GBM", "Hurdle-GBM"))) {
      cellf <- function(m) { x <- b4[b4$a == paste0(lr[2], " | ", cmp$a[i], " | adaptive") & b4$b == paste0(lr[2], " | ", cmp$b[i], " | adaptive") & b4$metric == m, ]
      if (nrow(x)) sprintf("%s (%s)", sg(x$median_diff[1]), fq(x$q_BH[1])) else "" }
      rows[[length(rows) + 1]] <- data.frame(Comparison = cmp$label[i], Learner = lr[1], AUC = cellf("AUC"), MAE = cellf("MAE"), `R2` = cellf("R2"), check.names = FALSE, stringsAsFactors = FALSE)
    }
    w4(do.call(rbind, rows), "Table8_exploratory_comparisons.csv")
  })
  
  blk4("Supplementary Table S2. Feature table", {
    lg <- rd4("D9_S2_feature_table_long.csv", OUT3); lg$win <- c(hr1 = "1 h", hr4 = "4 h", day1 = "1 d", day7 = "7 d")[lg$window]
    desc <- function(d) { if (!nrow(d)) return(""); paste(vapply(unique(d$statistic), function(s) paste0(s, ": ", paste(unique(d$win[d$statistic == s])[order(match(unique(d$win[d$statistic == s]), c("1 h", "4 h", "1 d", "7 d")))], collapse = ", ")), ""), collapse = "; ") }
    g <- unique(lg[, c("source", "stream")]); out <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) { d <- lg[lg$source == g$source[i] & lg$stream == g$stream[i], ]
    data.frame(Source = ifelse(g$source[i] == "ema", "EMA", "Passive"), Stream = g$stream[i], Measure = paste(sort(unique(d$measure)), collapse = ", "), `Mean-level features (window)` = desc(d[d$type == "mean", ]),
               `Dynamic features (window)` = desc(d[d$type == "dynamic", ]), Features = nrow(d), check.names = FALSE, stringsAsFactors = FALSE) }))
    out <- rbind(out, data.frame(Source = "Total", Stream = "", Measure = "", `Mean-level features (window)` = "", `Dynamic features (window)` = "",
                                 Features = sprintf("%d (%d EMA, %d passive)", nrow(lg), sum(lg$source == "ema"), sum(lg$source == "passive")), check.names = FALSE, stringsAsFactors = FALSE))
    w4(out, "TableS2_feature_table.csv")
  })
  
  blk4("Supplementary Table S3. Candidate predictors per configuration", {
    cc <- read.csv(file.path(OUT_DIR, "candidate_feature_counts.csv"), stringsAsFactors = FALSE)
    lab <- c("EMA" = "EMA-only, mean-only", "EMA (dynamic)" = "EMA-only, dynamic", "Passive" = "Passive-only, mean-only", "Passive (dynamic)" = "Passive-only, dynamic", "EMA + Passive" = "Combined, mean-only", "EMA + Passive (dynamic)" = "Combined, dynamic")
    out <- do.call(rbind, lapply(names(lab), function(cf) { n <- cc$n_candidate_features[cc$Configuration == cf]
    data.frame(`Model configuration` = lab[[cf]], `Candidate predictors (depending on the window)` = if (!length(n)) "" else if (min(n) == max(n)) as.character(min(n)) else sprintf("%d to %d", min(n), max(n)), check.names = FALSE, stringsAsFactors = FALSE) }))
    w4(out, "TableS3_candidate_predictors.csv")
  })
  
  blk4("Supplementary Table S6. Correlation and balanced accuracy", {
    b3 <- b3_adapt(); b3 <- b3[b3$variant == "adaptive", ]; ci <- read.csv(file.path(OUT3, "T4_median_CIs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    w4(with_sections(SPEC4, function(i) { k <- SPEC4$key[i]; r <- b3[b3$key == k, ]; x <- ci[ci$key == k & ci$metric == "r", ]
    hur <- grepl("^Hurdle", k) && !grepl("none", k)
    data.frame(Predictor = SPEC4$label[i], n = if (nrow(x)) as.character(x$n[1]) else "", `r, median [95% CI]` = if (nrow(x) && !is.na(x$median[1])) mci(x$median[1], x$lo[1], x$hi[1], 3) else "not defined",
               Sensitivity = if (hur && nrow(r)) r$Sens[1] else "", Specificity = if (hur && nrow(r)) r$Spec[1] else "", `Balanced accuracy` = if (hur && nrow(r)) r$BalAcc[1] else "", check.names = FALSE, stringsAsFactors = FALSE) }), "TableS6_correlation_and_balanced_accuracy.csv")
  })
  
  blk4("Supplementary Table S7. Combined dynamic models against the training-mean baseline and the previous-day mean", {
    s7 <- read.csv(file.path(OUT3, "S7_combined_dynamic_vs_benchmarks.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    models <- c("Continuous-RF | EMA + Passive (dynamic)" = "Random forest", "Hurdle-GBM | EMA + Passive (dynamic)" = "Hurdle model")
    rows <- list()
    for (mk in names(models)) for (bn in c("training-mean baseline", "previous-day mean")) for (m in c("AUC", "r", "MAE", "RMSE", "R2")) {
      x <- s7[s7$key == mk & s7$benchmark == bn & s7$metric == m, ]; if (!nrow(x)) next
      higher <- m %in% c("AUC", "r", "R2"); good <- if (higher) x$median_diff[1] > 0 else x$median_diff[1] < 0
      fav <- if (is.na(x$q_BH[1]) || x$q_BH[1] >= .05) "-" else if (good) "Model" else "Benchmark"
      rows[[length(rows) + 1]] <- data.frame(Comparison = paste(models[[mk]], "vs", bn), Metric = ifelse(m == "R2", "R2", m),
                                             `Median paired difference [95% CI]` = paste0(sg(x$median_diff[1]), " [", sg(x$ci_lo[1]), ", ", sg(x$ci_hi[1]), "]"), q = fq(x$q_BH[1]), Favours = fav, check.names = FALSE, stringsAsFactors = FALSE)
    }
    w4(do.call(rbind, rows), "TableS7_models_vs_benchmarks.csv")
  })
  
  blk4("Supplementary Table S8. Per-person comparison with the training-mean baseline", {
    s <- read.csv(file.path(OUT3, "D3_person_vs_baseline_summary.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    s <- s[order(s$Model != "Continuous-RF", match(s$Configuration, CFG)), ]
    out <- data.frame(Model = ifelse(grepl("^Continuous", s$Model), "Random forest", "Hurdle gradient boosting"), Configuration = lowp(s$Configuration),
                      `MAE gain, median [IQR]` = s$MAE_gain, `% lower error` = s$pct_better_MAE, `% better and p < .05` = s$pct_sig_better_MAE, `n q < .05` = s$n_q05_better_MAE,
                      `AUC gain, median [IQR]` = s$AUC_gain, `% higher AUC` = s$pct_better_AUC, `% better and p < .05 (AUC)` = s$pct_sig_better_AUC, `n q < .05 (AUC)` = s$n_q05_better_AUC, check.names = FALSE, stringsAsFactors = FALSE)
    w4(out, "TableS8_person_vs_baseline.csv")
  })
  
  blk4("Supplementary Table S9. Median AUC at higher SI-positive cutoffs (predictions re-scored)", {
    pp <- read.csv(file.path(OUT3, "D7_person_AUC_by_cutoff.csv"), stringsAsFactors = FALSE, check.names = FALSE); set.seed(SEED); cuts <- c(0, 2, 5, 10)
    ord <- c("Continuous-none | Baseline", "Hurdle-none | Baseline", "Reference-none | AlwaysZero", "Reference-none | MeanD1", "Reference-none | MeanD7", "Reference-none | Persistence", "Reference-none | PrevPrompt", SPEC4$key[8:19])
    lab <- setNames(c(SPEC4$label[1:7], paste0("Random forest, ", lowp(CFG)), paste0("Hurdle gradient boosting, ", lowp(CFG))), c(SPEC4$key[1:7], SPEC4$key[8:19]))
    ncut <- sapply(cuts, function(ct) sum(!is.na(pp$AUC[pp$cutoff == ct & pp$key == "Continuous-RF | EMA + Passive (dynamic)"])))
    cell <- function(k, ct) { x <- pp$AUC[pp$key == k & pp$cutoff == ct]; x <- x[is.finite(x)]; if (length(x) < 5) return(""); bt <- replicate(5000, median(sample(x, replace = TRUE))); mci(median(x), unname(quantile(bt, .025)), unname(quantile(bt, .975)), 3) }
    out <- do.call(rbind, lapply(ord, function(k) { r <- data.frame(Predictor = lab[[k]], stringsAsFactors = FALSE); for (j in seq_along(cuts)) r[[sprintf("Above %d (n = %d)", cuts[j], ncut[j])]] <- cell(k, cuts[j]); r }))
    w4(out, "TableS9_cutoff_sensitivity.csv")
  })
  
  blk4("Supplementary Table S10. Analysed participants vs those excluded by the outcome-count rule", {
    d <- rd4("C5b_analysed_vs_excluded_by_outcome_rule.csv"); d$p <- fq(as.numeric(d$p))
    names(d)[names(d) == "variable"] <- "Characteristic"; names(d)[names(d) == "level"] <- "Level"; names(d)[names(d) == "analysed"] <- "Analysed"; names(d)[names(d) == "excluded_by_outcome_rule"] <- "Excluded by the outcome-count rule"
    w4(d, "TableS10_analysed_vs_excluded.csv")
  })
  
  blk4("Supplementary Table S11. Baseline measures by kind of excluded participant", {
    d <- rd4("C5c_baseline_by_excluded_type.csv"); d$p <- fq(as.numeric(d$p))
    names(d)[names(d) == "variable"] <- "Measure"; names(d)[names(d) == "level"] <- "Level"; names(d)[names(d) == "analysed"] <- "Analysed"
    names(d)[names(d) == "almost never"] <- "Almost never endorsed item 9"; names(d)[names(d) == "almost always"] <- "Almost always endorsed item 9"; names(d)[names(d) == "both"] <- "Both (few of each)"
    w4(d, "TableS11_baseline_by_excluded_type.csv")
  })
  
  cat(sprintf("\nStage 4: %d table files written to %s\n", length(tabs4), OUT_TABLES)); for (f in tabs4) cat("   ", f, "\n")
  # ---- Reproducibility check: the tables must be identical from run to run (same seed, same data, same R and package versions) ----
  chk_file <- file.path(OUT_TABLES, "checksums_md5.txt")
  now_md5  <- tools::md5sum(file.path(OUT_TABLES, tabs4)); names(now_md5) <- tabs4
  if (file.exists(chk_file)) {
    old <- tryCatch(read.table(chk_file, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(old) && nrow(old)) {
      both <- intersect(old$file, tabs4); changed <- both[unname(now_md5[both]) != old$md5[match(both, old$file)]]
      if (!length(both)) cat("Reproducibility check: no earlier checksums for these tables.\n")
      else if (!length(changed)) cat(sprintf("Reproducibility check: all %d tables are IDENTICAL to the previous run (seed %d).\n", length(both), SEED))
      else cat("Reproducibility check: these tables DIFFER from the previous run:", paste(changed, collapse = ", "),
               "\n   (expected only if the data, SEED, settings, or the R or package versions changed; see analysis_sessionInfo.txt)\n")
    }
  } else cat("Reproducibility check: first run; the checksums are saved, and the next run is compared with them.\n")
  write.table(data.frame(file = tabs4, md5 = unname(now_md5)), chk_file, row.names = FALSE, quote = FALSE)
  if (!KEEP_WORKING_FILES) {                                    # keep only the tables, the figure and the checkpoints
    unlink(OUT, recursive = TRUE); unlink(list.files(OUT_DIR, pattern = "\\.csv$", full.names = TRUE))
    cat("Working files deleted (posthoc/working_files and the Stage 1 CSVs). The checkpoints are kept: rerun this script to rebuild everything in minutes.\n")
  }
}  # end STAGE 4


# ------------------------------------------------------------------------------
# Record the settings and environment that produced these results
# ------------------------------------------------------------------------------
try(writeLines(c(
  sprintf("Run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("SEED = %s | TRAIN/VAL/RUNIN days = %d/%d/%d | NTREE = %d | MTRY_GRID = %s",
          SEED, TRAIN_DAYS, VAL_DAYS, RUNIN_DAYS, NTREE, paste(MTRY_GRID, collapse = ",")),
  sprintf("COMBOS = %s", paste(COMBOS$framework, COMBOS$learner, sep = "-", collapse = ", ")),
  sprintf("RNGkind = %s | SEED = %d (per-participant seeds in Stage 1: SEED + participant ID)", paste(RNGkind(), collapse = " / "), SEED),
  sprintf("DEBUG_N_UIDS = %s | RUN_FIXED_WINDOWS = %s | RUN_STAGE3 = %s | RUN_STAGE4 = %s", DEBUG_N_UIDS, RUN_FIXED_WINDOWS, RUN_STAGE3, RUN_STAGE4),
  sprintf("SI_CUT = %s | MAX_PREDICTORS = %s | MIN_NONMISSING = %s | RUN_TAG = '%s'", SI_CUT, MAX_PREDICTORS, MIN_NONMISSING, RUN_TAG),
  "", capture.output(sessionInfo())),
  file.path(OUT_DIR, "analysis_sessionInfo.txt")), silent = TRUE)