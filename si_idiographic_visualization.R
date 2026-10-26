################################################################################
# SI Idiographic Prediction Project - MANUSCRIPT FIGURES (si_manuscript_figures.R)
#
# Draws the four figures of the manuscript from the files the analysis script already saved. It does NOT refit anything.
#   Figure 1  Prediction design: (A) one test day, (B) predictor windows, (C) sliding-window cross-validation
#   Figure 2  Participant flow
#   Figure 3  Within-person AUC of the baselines, reference predictors and models (median and IQR)
#   Figure 4  Paired differences for the preregistered comparisons (H1 blue, H2 green; filled = q < .05)
#
#   After the four figures: PERSONALIZED figures, one per participant (PT), and combined versions (see the last section).
# Output (4. Results/figures_manuscript), one set of files per figure, named Figure1 ... Figure4:
#   .pdf   vector, fonts embedded (Cairo)          .tiff  600 dpi, LZW compression          .png  300 dpi preview
# These follow the Cambridge University Press artwork rules (BJPsych Open is a Cambridge journal): TIFF, EPS or PDF; at least
# 600 dpi for combined text-and-colour figures; fonts embedded and used consistently; one file per figure; plain alphanumeric names.
# Check the exact column widths and the file types accepted by the journal's own "Preparing your materials" page before submitting;
# change SIZE_IN, FORMATS or FONT_FAMILY below. No titles are drawn inside the figures (the captions are in the manuscript).
#
# Inputs (written by si_idiographic_data_analysis_1.R, in 4. Results/posthoc/working_files):
#   B2_person_metrics.csv (Figures 2 and 3)   B4_contrasts.csv (Figure 4)   C1_participant_flow.csv, C1_excluded_split_never_always.csv (Figure 2)
#   model_results_matched/matched_predictions.csv (personalized figures only)
# Run:  Rscript si_manuscript_figures.R      (set the working directory or SI_PROJECT_DIR to the project folder)
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork) })
# ==============================================================================
# CONFIG
# ==============================================================================
PROJECT_DIR <- Sys.getenv("SI_PROJECT_DIR", unset = getwd())
if (basename(PROJECT_DIR) %in% c("1. Code", "2. Raw data", "3. Cleaned data", "4. Results")) PROJECT_DIR <- dirname(PROJECT_DIR)
RUN_TAG <- ""                                    # the tag of the analysis run to draw ("" = the main analysis)
RES_DIR <- file.path(PROJECT_DIR, "4. Results")
first_file <- function(...) { p <- c(...); hit <- p[file.exists(p)]; if (length(hit)) hit[1] else p[1] }
POST  <- file.path(RES_DIR, paste0("posthoc", RUN_TAG))
WORKF <- function(f) first_file(file.path(POST, "working_files", f), file.path(POST, f))
CSV <- list(metrics = WORKF("B2_person_metrics.csv"), contrasts = WORKF("B4_contrasts.csv"), flow = WORKF("C1_participant_flow.csv"), split = WORKF("C1_excluded_split_never_always.csv"))
OUT_DIR <- file.path(RES_DIR, paste0("figures_manuscript", RUN_TAG))

RUN_FIG <- c(Figure1 = TRUE, Figure2 = TRUE, Figure3 = TRUE, Figure4 = TRUE)
FORMATS <- c("pdf", "tiff", "png")               # "pdf" and "tiff" for submission, "png" to look at
FONT_FAMILY <- "sans"                            # one sans-serif family for every figure; "Arial" or "Helvetica" if installed (extrafont/systemfonts)
BASE_PT <- 8                                     # base text size in points at final size (panel letters are 10 pt bold)
SIZE_IN <- list(Figure1 = c(7.0, 8.6), Figure2 = c(7.0, 8.2), Figure3 = c(7.0, 5.6), Figure4 = c(7.0, 5.2))   # width x height in inches
MAIN_MODEL <- "Continuous-RF"; MAIN_CFG <- "EMA + Passive (dynamic)"                 # the model whose AUC defines "n with AUC defined"
# Numbers used for Figure 2 only if the flow files are missing (the manuscript text)
FLOW_FALLBACK <- list(enrolled = 300, ema = 227, passive = 225, outcome = 143, allmodels = 142, auc = 137, r2 = 139, nsi0 = 3, nsi_all = 2, never = 35, always = 47)
# Colour-blind-safe palette (Okabe-Ito); the groups also differ in marker shape so the figures stay readable in greyscale
COL_GROUP <- c("Baseline" = "#7F7F7F", "Reference predictors" = "#E69F00", "Random forest" = "#0072B2", "Gradient boosting (hurdle)" = "#009E73")
SHP_GROUP <- c("Baseline" = 15, "Reference predictors" = 17, "Random forest" = 16, "Gradient boosting (hurdle)" = 18)
COL_H <- c(H1 = "#0072B2", H2 = "#009E73")
# Personalized figures (last section)
RUN_PERSONALIZED <- TRUE                         # one figure per participant (PT) and the combined versions
PERSONALIZED_CONFIGS <- "all"                  # "all" = six panels per participant (one per feature configuration); "main" = one panel (combined dynamic) plus where the AUC falls among participants
PERSONALIZED_PER_FILE <- TRUE                    # TRUE = also one PNG per participant; the multi-page PDF and the overview are always written
SHOW_ID <- TRUE                                  # TRUE = the study ID in the title and file name; FALSE = "Participant k of N" (use FALSE for anything that leaves the team)
PARTICIPANTS <- NULL                             # NULL = every participant, or e.g. c("0300", "0011")
RUNIN_DAYS <- 6; TRAIN_DAYS <- 45; VAL_DAYS <- 10; EMA_PER_DAY <- 3          # design of the cross-validation: fold k tests study day 61 + k
CSV$preds <- file.path(RES_DIR, paste0("model_results_matched", RUN_TAG), "matched_predictions.csv")
PERSONALIZED_DIR <- file.path(OUT_DIR, "personalized")
CFG <- c("EMA", "EMA (dynamic)", "Passive", "Passive (dynamic)", "EMA + Passive", "EMA + Passive (dynamic)")
SPEC <- data.frame(
  label = c("Training-mean baseline (continuous)", "Training-mean baseline (hurdle)", "Persistence", "Previous prompt", "Previous 1-day mean", "Previous 7-day mean", "Always zero",
            paste0("RF, ", CFG), paste0("GBM hurdle, ", CFG)),
  key = c("Continuous-none | Baseline", "Hurdle-none | Baseline", "Reference-none | Persistence", "Reference-none | PrevPrompt", "Reference-none | MeanD1", "Reference-none | MeanD7",
          "Reference-none | AlwaysZero", paste0("Continuous-RF | ", CFG), paste0("Hurdle-GBM | ", CFG)),
  group = c(rep("Baseline", 2), rep("Reference predictors", 5), rep("Random forest", 6), rep("Gradient boosting (hurdle)", 6)), stringsAsFactors = FALSE)
# ==============================================================================
# HELPERS
# ==============================================================================
need_file <- function(path) { if (!file.exists(path)) stop("File not found: ", path, "\nRun the analysis script (Stages 1 to 3) first, or check PROJECT_DIR / RUN_TAG."); path }
pt2mm <- function(pt) pt / ggplot2::.pt                                           # text size in points -> the 'size' of geom_text (mm)
theme_ms <- function() theme_classic(base_size = BASE_PT, base_family = FONT_FAMILY) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.grid = element_blank(), axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.tag = element_text(face = "bold", size = 10, family = FONT_FAMILY), legend.key.size = grid::unit(3, "mm"))
save_ms <- function(p, name) {                                                   # PDF (vector), TIFF (600 dpi, LZW) and PNG preview
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE); sz <- SIZE_IN[[name]]; written <- character(0)
  one <- function(ext, ...) tryCatch({ f <- file.path(OUT_DIR, paste0(name, ".", ext)); ggplot2::ggsave(f, p, width = sz[1], height = sz[2], units = "in", bg = "white", ...); written <<- c(written, basename(f)) },
                                     error = function(e) message("  could not write ", name, ".", ext, ": ", conditionMessage(e)))
  if ("pdf"  %in% FORMATS) one("pdf", device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf")
  if ("tiff" %in% FORMATS) one("tiff", device = "tiff", dpi = 600, compression = "lzw")
  if ("png"  %in% FORMATS) one("png", dpi = 300)
  message(sprintf("%s (%.1f x %.1f in): %s", name, sz[1], sz[2], paste(written, collapse = ", ")))
  invisible(written)
}
std_metrics <- function(path) {                                                  # person-level metrics: adaptive windows only
  m <- readr::read_csv(need_file(path), show_col_types = FALSE, col_types = readr::cols(uid = readr::col_character(), .default = readr::col_guess())) |> as.data.frame()
  m <- m[m$variant == "adaptive", , drop = FALSE]; m$key <- paste(m$Model, m$Configuration, sep = " | "); m
}
# One drawing function for the hand-built panels of Figure 1: rectangles, segments, points and text given as data frames
draw_geom <- function(G, lim) {
  p <- ggplot()
  if (!is.null(G$rect)) p <- p + geom_rect(data = G$rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, alpha = alpha), colour = NA)
  if (!is.null(G$seg))  p <- p + geom_segment(data = G$seg, aes(x = x, xend = xend, y = y, yend = yend, colour = col), linewidth = 0.3)
  if (!is.null(G$pt))   p <- p + geom_point(data = G$pt, aes(x = x, y = y, fill = fill, colour = col), shape = 21, size = 2, stroke = 0.6)
  if (!is.null(G$txt))  p <- p + geom_text(data = G$txt, aes(x = x, y = y, label = label, hjust = hjust, size = size, colour = col, fontface = face), family = FONT_FAMILY, lineheight = 0.9)
  p + scale_fill_identity() + scale_colour_identity() + scale_alpha_identity() + scale_size_identity() +
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = FALSE) + theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
}
# ==============================================================================
# FIGURE 1 - prediction design (the geometry below is fixed; it does not depend on the data)
#   panel A: one test day with three prompts   panel B: predictor windows for one prompt   panel C: sliding-window cross-validation
#   Study days: fold k tests study day 61 + k (fold 1 = day 62, fold 29 = day 90), trains on the 45 days ending 11 days earlier and
#   validates on the 10 days before the test day.
# ==============================================================================
GEOM <- list(
  A = list(
    rect = data.frame(
      xmin = c(4, 5, 4, 8.5, 4, 13, 4, 13.2),
      xmax = c(5, 9, 8.5, 12.5, 13, 17, 4.55, 13.75),
      ymin = c(3.32, 3.32, 2.22, 2.22, 1.12, 1.12, -0.59, -0.59),
      ymax = c(3.68, 3.68, 2.58, 2.58, 1.48, 1.48, -0.31, -0.31),
      fill = c("#0072B2", "#D55E00", "#0072B2", "#D55E00", "#0072B2", "#D55E00", "#D55E00", "#0072B2"),
      alpha = c(0.9, 1, 0.9, 1, 0.9, 1, 1, 1),
      stringsAsFactors = FALSE),
    seg = data.frame(
      x = c(4, 4, 4, 6, 9, 12, 15, 18, 21),
      xend = c(21.5, 21.5, 21.5, 6, 9, 12, 15, 18, 21),
      y = c(3.5, 2.4, 1.3, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62),
      yend = c(3.5, 2.4, 1.3, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72),
      col = c("#BBBBBB", "#BBBBBB", "#BBBBBB", "#555555", "#555555", "#555555", "#555555", "#555555", "#555555"),
      stringsAsFactors = FALSE),
    pt = data.frame(
      x = c(9, 9, 12.5, 4.27, 13.47),
      y = c(2.4, 1.3, 1.3, -1, -1),
      fill = c("white", "#009E73", "#009E73", "#009E73", "white"),
      col = c("#009E73", "#009E73", "#009E73", "#009E73", "#009E73"),
      stringsAsFactors = FALSE),
    txt = data.frame(
      x = c(9, 5, 3.6, 12.5, 8.5, 3.6, 17, 13, 3.6, 4.1, 6, 9, 12, 15, 18, 21, 12.75, 4.8, 14, 4.8, 14),
      y = c(3.86, 3.86, 3.5, 2.76, 2.76, 2.4, 1.66, 1.66, 1.3, 3.12, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.12, -0.45, -0.45, -1, -1),
      label = c("t", "t \u2212 4 h", "Prompt 1", "t", "t \u2212 4 h", "Prompt 2", "t", "t \u2212 4 h", "Prompt 3", "earlier days", "06", "09", "12", "15", "18", "21", "Time of day (h)", "Outcome: item 9, asks about the 4 h before t", "Predictors: data completed at or before t \u2212 4 h", "Earlier response, usable (at or before t \u2212 4 h)", "Response inside the 4-h gap, not usable"),
      hjust = c(0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0),
      size = c(2.5, 2.3, 2.5, 2.5, 2.3, 2.5, 2.5, 2.3, 2.5, 2.1, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.4, 2.2, 2.2, 2.2, 2.2),
      col = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "#444444", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
      face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain"),
      stringsAsFactors = FALSE)
  ),
  B = list(
    rect = data.frame(
      xmin = c(9, 5.4, 2.7, 8.1, 7.1, 5.4, 2.7, 5.64, 6.26, 6.88, 7.5, 8.12, 8.74, 9.36),
      xmax = c(11.6, 9, 9, 9, 9, 9, 9, 6.16, 6.78, 7.4, 8.02, 8.64, 9.26, 9.88),
      ymin = c(0.5, 6.12, 5.5, 4.55, 3.93, 3.31, 2.69, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92),
      ymax = c(7.1, 6.58, 5.96, 5.01, 4.39, 3.77, 3.15, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38),
      fill = c("#E3E3E3", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#0072B2", "#009E73"),
      alpha = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9),
      stringsAsFactors = FALSE),
    seg = data.frame(
      x = c(9, 11.6),
      xend = c(9, 11.6),
      y = c(0.5, 0.5),
      yend = c(7.1, 7.1),
      col = c("#444444", "#444444"),
      stringsAsFactors = FALSE),
    pt = NULL,
    txt = data.frame(
      x = c(9, 11.6, 10.3, 10.3, -0.35, 5.55, 2.85, -0.35, 8.25, 7.25, 5.55, 2.85, -0.35, 10.05, 5.7, 6),
      y = c(7.45, 7.45, 6.7, 5.55, 6.04, 6.35, 5.73, 3.85, 4.78, 4.16, 3.54, 2.92, 1.15, 1.15, 1.71, 0.12),
      label = c("t \u2212 4 h", "t", "4-h gap", "no data\nused", "EMA items", "1 d", "1 wk", "Passive sensors", "1 h", "4 h", "1 d", "1 wk", "Sleep", "last night", "7 nights (exempt from the gap)", "Predictor windows for one prompt (not to scale)"),
      hjust = c(0.5, 0.5, 0.5, 0.5, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0.5, 0.5),
      size = c(2.4, 2.4, 2.2, 2.1, 2.5, 2.1, 2.1, 2.5, 2.1, 2.1, 2.1, 2.1, 2.5, 2.1, 2.1, 2.1),
      col = c("black", "black", "black", "black", "black", "white", "white", "black", "white", "white", "white", "white", "black", "black", "black", "#444444"),
      face = c("plain", "plain", "plain", "plain", "bold", "plain", "plain", "bold", "plain", "plain", "plain", "plain", "bold", "plain", "plain", "plain"),
      stringsAsFactors = FALSE)
  ),
  C = list(
    rect = data.frame(
      xmin = c(0.5, 6.5, 51.5, 61.5, 7.5, 52.5, 62.5, 34.5, 79.5, 89.5, 2, 30, 68),
      xmax = c(6.5, 51.5, 61.5, 62.5, 52.5, 62.5, 63.5, 79.5, 89.5, 90.5, 3.6, 31.6, 69.6),
      ymin = c(6.9, 5.7, 5.7, 5.7, 4.85, 4.85, 4.85, 3.15, 3.15, 3.15, 8, 8, 8),
      ymax = c(7.3, 6.1, 6.1, 6.1, 5.25, 5.25, 5.25, 3.55, 3.55, 3.55, 8.4, 8.4, 8.4),
      fill = c("#8C8C8C", "#0072B2", "#E69F00", "#D55E00", "#0072B2", "#E69F00", "#D55E00", "#0072B2", "#E69F00", "#D55E00", "#0072B2", "#E69F00", "#D55E00"),
      alpha = c(0.95, 0.9, 0.95, 1, 0.9, 0.95, 1, 0.9, 0.95, 1, 1, 1, 1),
      stringsAsFactors = FALSE),
    seg = data.frame(
      x = c(61.5, 61.5, 90.5, 1, 30, 60, 90, 0.5),
      xend = c(90.5, 61.5, 90.5, 1, 30, 60, 90, 90.5),
      y = c(2.6, 2.6, 2.6, 1.55, 1.55, 1.55, 1.55, 1.7),
      yend = c(2.6, 2.75, 2.75, 1.7, 1.7, 1.7, 1.7, 1.7),
      col = c("#333333", "#333333", "#333333", "#555555", "#555555", "#555555", "#555555", "#999999"),
      stringsAsFactors = FALSE),
    pt = NULL,
    txt = data.frame(
      x = c(7.6, -1.5, 63.5, -1.5, 64.5, -1.5, 91.5, 45, 76, 1, 30, 60, 90, 45, 4.4, 32.4, 70.4),
      y = c(7.1, 5.9, 5.9, 5.05, 5.05, 3.35, 3.35, 4.2, 2.28, 1.3, 1.3, 1.3, 1.3, 0.92, 8.2, 8.2, 8.2),
      label = c("run-in (days 1 to 6)", "Fold 1", "day 62", "Fold 2", "day 63", "Fold 29", "day 90", "\u22ee", "29 test days (days 62 to 90)", "1", "30", "60", "90", "Study day", "Train: 45 days", "Validate: next 10 days (select the window)", "Test: next day (model refit on all 55 days)"),
      hjust = c(0, 1, 0, 1, 0, 1, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0),
      size = c(2.2, 2.4, 2.2, 2.4, 2.2, 2.4, 2.2, 3.4, 2.2, 2.2, 2.2, 2.2, 2.2, 2.4, 2.2, 2.2, 2.2),
      col = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
      face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain"),
      stringsAsFactors = FALSE)
  )
)
LIM <- list(A = list(x = c(0.7999999999999998, 22.1), y = c(-1.35, 4.15)), B = list(x = c(-3.6, 13.6), y = c(-0.15, 7.9)), C = list(x = c(-9, 104), y = c(0.6, 8.7)))
make_fig1 <- function() {
  pA <- draw_geom(GEOM$A, LIM$A); pB <- draw_geom(GEOM$B, LIM$B); pC <- draw_geom(GEOM$C, LIM$C)
  fig <- (pA / pB / pC) + plot_layout(heights = c(1.1, 1.15, 1.2)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 10, family = FONT_FAMILY))
  save_ms(fig, "Figure1")
}
# ==============================================================================
# FIGURE 2 - participant flow (numbers read from the analysis output)
# ==============================================================================
flow_numbers <- function(metrics) {
  fn <- FLOW_FALLBACK; src <- "the manuscript text (flow files not found)"
  if (file.exists(CSV$flow)) {
    fl <- read.csv(CSV$flow, stringsAsFactors = FALSE); pick <- function(pat) { v <- fl$n[grepl(pat, fl$step)]; if (length(v)) v[1] else NA }
    fn$enrolled <- pick("^In inclusion report"); fn$ema <- pick("^EMA compliance"); fn$passive <- pick("passive coverage"); fn$outcome <- pick("SI-positive"); fn$allmodels <- pick("^Analysed in every model"); src <- CSV$flow
  }
  main <- metrics[metrics$key == paste(MAIN_MODEL, MAIN_CFG, sep = " | "), ]
  if (nrow(main)) { fn$auc <- sum(is.finite(main$AUC)); fn$r2 <- sum(is.finite(main$R2)); fn$nsi0 <- sum(main$n_SI == 0); fn$nsi_all <- sum(main$n_nonSI == 0) }
  fn$never <- NA; fn$always <- NA
  if (file.exists(CSV$split)) { sp <- read.csv(CSV$split, stringsAsFactors = FALSE)
  fn$never <- sp$n[grepl("^almost never", sp$group)][1]; fn$always <- sp$n[grepl("^almost always", sp$group)][1] }
  message("Figure 2 numbers from: ", src); fn
}
make_fig2 <- function(metrics) {
  fn <- flow_numbers(metrics); MAINX <- 3.7; SIDEX <- 9.0; MW <- 6.0; SW <- 5.0
  bx <- function(x, y, s, w, size = 2.3) { n <- nchar(s) - nchar(gsub("\n", "", s, fixed = TRUE)) + 1; h <- 0.5 + 0.34 * (n - 1)
  data.frame(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, x = x, y = y, label = s, size = size, stringsAsFactors = FALSE) }
  excl_split <- if (is.finite(fn$never) && is.finite(fn$always)) sprintf("\n%d almost never (< 10 SI-positive)\n%d almost always (< 20 SI-negative)", fn$never, fn$always) else ""
  main <- rbind(bx(MAINX, 12.0, sprintf("Enrolled participants\nn = %d", fn$enrolled), MW),
                bx(MAINX, 10.1, sprintf("Met EMA compliance criterion\n(\u2265 75%% of EMA prompts completed)\nn = %d", fn$ema), MW),
                bx(MAINX, 8.0, sprintf("Met passive-sensing coverage criterion\n(\u2265 75%% average coverage)\nn = %d", fn$passive), MW),
                bx(MAINX, 5.9, sprintf("Met outcome-count criterion\n(\u2265 10 SI-positive and \u2265 20 SI-negative prompts)\nAnalysed: n = %d", fn$outcome), MW),
                bx(MAINX, 3.9, sprintf("Included in every model comparison\nn = %d", fn$allmodels), MW))
  side <- rbind(bx(SIDEX, 11.05, sprintf("Did not meet EMA criterion\nn = %d", fn$enrolled - fn$ema), SW), bx(SIDEX, 9.05, sprintf("Below passive-sensing coverage\nn = %d", fn$ema - fn$passive), SW),
                bx(SIDEX, 6.95, sprintf("Excluded by the outcome-count criterion\nn = %d%s", fn$passive - fn$outcome, excl_split), SW),
                bx(SIDEX, 4.9, sprintf("EMA-only models could not be fitted\nn = %d", fn$outcome - fn$allmodels), SW))
  endb <- rbind(bx(2.0, 1.2, sprintf("AUC defined\nn = %d\n(%d no ideation,\n%d ideation at every prompt)", fn$auc, fn$nsi0, fn$nsi_all), 3.9),
                bx(6.0, 1.2, sprintf("R\u00b2 and r defined\nn = %d\n(%d no ideation)", fn$r2, fn$nsi0), 3.5), bx(9.7, 1.2, sprintf("MAE and RMSE\nn = %d", fn$allmodels), 3.3))
  boxes <- rbind(main, side, endb); segs <- list(); k <- 0
  add <- function(x, xend, y, yend, arrow) { k <<- k + 1; segs[[k]] <<- data.frame(x = x, xend = xend, y = y, yend = yend, arrow = arrow) }
  for (i in 1:4) { add(MAINX, MAINX, main$ymin[i], main$ymax[i + 1], 1); ym <- (main$ymin[i] + main$ymax[i + 1]) / 2; add(MAINX, side$xmin[i], ym, ym, 1) }   # steps and branches to the exclusion boxes
  bar <- 2.55; add(MAINX, MAINX, main$ymin[5], bar, 0); add(endb$x[1], endb$x[3], bar, bar, 0); for (i in 1:3) add(endb$x[i], endb$x[i], bar, endb$ymax[i], 1)
  segs <- do.call(rbind, segs)
  p <- ggplot() +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = "black", linewidth = 0.3) +
    geom_text(data = boxes, aes(x = x, y = y, label = label, size = size), family = FONT_FAMILY, lineheight = 0.95) +
    geom_segment(data = segs[segs$arrow == 1, ], aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.3, arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")) +
    geom_segment(data = segs[segs$arrow == 0, ], aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.3) +
    scale_size_identity() + coord_cartesian(xlim = c(0, 12), ylim = c(-0.2, 13), expand = FALSE) + theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  save_ms(p, "Figure2")
}
# ==============================================================================
# FIGURE 3 - within-person AUC: median (marker) and interquartile range (line); dashed line = chance
# ==============================================================================
make_fig3 <- function(m) {
  d <- m |> dplyr::inner_join(SPEC, by = "key") |> dplyr::filter(is.finite(AUC))
  n_auc <- dplyr::n_distinct(d$uid[d$key == paste(MAIN_MODEL, MAIN_CFG, sep = " | ")])
  sm <- d |> dplyr::group_by(label, group) |> dplyr::summarise(q25 = stats::quantile(AUC, .25), med = stats::median(AUC), q75 = stats::quantile(AUC, .75), .groups = "drop") |>
    dplyr::mutate(label = factor(label, levels = rev(SPEC$label)), group = factor(group, levels = names(COL_GROUP)))
  p <- ggplot(sm, aes(y = label, colour = group, shape = group)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 0.7) + geom_point(aes(x = med), size = 2) +
    geom_text(aes(x = 0.75, label = sprintf("%.3f", med)), hjust = 1, colour = "black", size = pt2mm(BASE_PT - 1), family = FONT_FAMILY) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = COL_GROUP, guide = "none") + scale_shape_manual(values = SHP_GROUP, guide = "none") +
    coord_cartesian(xlim = c(0.35, 0.75)) + scale_x_continuous(breaks = seq(0.4, 0.7, 0.1)) +
    labs(x = "Within-person AUC", y = NULL) + theme_ms() +
    theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = BASE_PT - 1),
          panel.spacing.y = grid::unit(2, "mm"), axis.line.y = element_blank(), axis.ticks.y = element_blank())
  message("Figure 3: n = ", n_auc, " participants with AUC defined (for the caption)")
  save_ms(p, "Figure3")
}
# ==============================================================================
# FIGURE 4 - paired differences, preregistered comparisons (random forest, dynamic features)
#   markers: median paired difference; lines: 95% bootstrap interval; filled = BH q < .05, open = q >= .05
#   blue: H1 (model minus training-mean baseline); green: H2 (combined minus single-source model)
# ==============================================================================
make_fig4 <- function(path) {
  b4 <- readr::read_csv(need_file(path), show_col_types = FALSE) |> as.data.frame()
  b4 <- b4[!is.na(b4$osf_hypothesis) & nzchar(b4$osf_hypothesis), , drop = FALSE]
  if (!nrow(b4)) stop("no confirmatory rows in ", path)
  lab <- function(s) sub(" \\| adaptive$", "", sub("^Continuous-RF \\| ", "", sub("^Continuous-none \\| Baseline", "baseline", s)))
  b4$cmp <- paste(lab(b4$a), "vs", lab(b4$b))
  lev <- c("EMA (dynamic) vs baseline", "Passive (dynamic) vs baseline", "EMA + Passive (dynamic) vs baseline", "EMA + Passive (dynamic) vs EMA (dynamic)", "EMA + Passive (dynamic) vs Passive (dynamic)")
  b4 <- b4[b4$cmp %in% lev & b4$metric %in% c("R2", "r", "MAE", "RMSE"), , drop = FALSE]
  b4$cmp <- factor(b4$cmp, levels = rev(lev))
  ml <- c(R2 = "R\u00b2 (higher is better)", r = "Pearson r (higher is better)", MAE = "MAE (lower is better)", RMSE = "RMSE (lower is better)")
  b4$metric_lab <- factor(ml[b4$metric], levels = unname(ml))
  b4$fill_col <- ifelse(!is.na(b4$q_confirmatory) & b4$q_confirmatory < 0.05, unname(COL_H[b4$osf_hypothesis]), "white")
  b4$q_lab <- ifelse(is.na(b4$q_confirmatory), "", ifelse(b4$q_confirmatory < 0.001, "q < .001", sub("^0", "q = ", sprintf("%.3f", b4$q_confirmatory))))
  p <- ggplot(b4, aes(y = cmp, colour = osf_hypothesis)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.6) +
    geom_point(aes(x = median_diff, fill = fill_col), shape = 21, size = 2.2, stroke = 0.7) +
    geom_text(aes(x = ci_hi, label = q_lab), hjust = -0.15, colour = "grey25", size = pt2mm(BASE_PT - 2), family = FONT_FAMILY) +
    facet_wrap(~ metric_lab, nrow = 2, scales = "free_x") +
    scale_colour_manual(values = COL_H, name = NULL, labels = c(H1 = "H1: model minus training-mean baseline", H2 = "H2: combined minus single-source model")) +
    scale_fill_identity() + scale_x_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(x = "Median paired difference (95% bootstrap interval)", y = NULL) + theme_ms() +
    theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold", size = BASE_PT), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.spacing = grid::unit(4, "mm"))
  save_ms(p, "Figure4")
}
# ==============================================================================
# PERSONALIZED FIGURES - one per participant (PT) and combined versions
#   PERSONALIZED_CONFIGS <- "all" (default): six panels per participant, one per feature configuration, each with observed item 9, the random forest, the hurdle model and the baseline.
#   PERSONALIZED_CONFIGS <- "main": one panel. Per participant: (top) observed item 9 and the predictions of the combined dynamic random forest, the combined dynamic hurdle model and the
#   training-mean baseline over the participant's test days; (bottom) where the participant's AUC falls among all participants.
#   Combined versions: Personalized_all_participants.pdf (one page per participant) and Personalized_overview.pdf/.png (every participant on one page,
#   ordered by AUC). Study IDs appear in the titles unless SHOW_ID <- FALSE.
# ==============================================================================
SERIES <- c("Random forest (combined, dynamic)" = paste(MAIN_MODEL, MAIN_CFG, sep = " | "), "Hurdle model (combined, dynamic)" = paste("Hurdle-GBM", MAIN_CFG, sep = " | "),
            "Training-mean baseline" = "Continuous-none | Baseline")
COL_SERIES <- c("Observed" = "#D55E00", "Random forest (combined, dynamic)" = "#0072B2", "Hurdle model (combined, dynamic)" = "#009E73", "Training-mean baseline" = "#7F7F7F")
LTY_SERIES <- c("Observed" = "solid", "Random forest (combined, dynamic)" = "solid", "Hurdle model (combined, dynamic)" = "solid", "Training-mean baseline" = "dashed")
read_predictions <- function() {
  keep <- c("uid", "fold", "row_id", "true", "pred", "framework", "learner", "Configuration", "variant")
  p <- data.table::fread(need_file(CSV$preds), select = keep, colClasses = list(character = "uid")) |> as.data.frame()
  p <- p[p$variant == "adaptive", , drop = FALSE]; p$key <- paste(paste0(p$framework, "-", p$learner), p$Configuration, sep = " | ")
  keys <- if (identical(PERSONALIZED_CONFIGS, "all")) unique(c(SERIES, paste0("Continuous-RF | ", CFG), paste0("Hurdle-GBM | ", CFG))) else SERIES
  p <- p[p$key %in% keys, , drop = FALSE]; p$series <- names(SERIES)[match(p$key, SERIES)]; p
}
person_series <- function(p, uid) {                                              # observed and predicted item 9 on the study-day axis (fold 1 = day 62)
  d <- p[p$uid == uid, , drop = FALSE]; rf <- d[d$series == names(SERIES)[1], , drop = FALSE]
  if (!nrow(rf)) return(NULL)
  rf <- rf[order(rf$fold, rf$row_id), , drop = FALSE]; rf$pos <- ave(rf$fold, rf$fold, FUN = seq_along); rf$n <- ave(rf$fold, rf$fold, FUN = length)
  rf$day <- RUNIN_DAYS + TRAIN_DAYS + VAL_DAYS + rf$fold + (rf$pos - 1) / pmax(rf$n, EMA_PER_DAY)
  d <- merge(d, rf[, c("fold", "row_id", "day")], by = c("fold", "row_id"))
  rbind(data.frame(day = rf$day, series = "Observed", value = rf$true, stringsAsFactors = FALSE), data.frame(day = d$day, series = d$series, value = d$pred, stringsAsFactors = FALSE))
}
fig_person <- function(uid, p, m, label) {
  d <- person_series(p, uid); if (is.null(d) || !nrow(d)) return(NULL)
  d$series <- factor(d$series, levels = names(COL_SERIES))
  mk <- function(key) m[m$uid == uid & m$key == key, , drop = FALSE]; a <- mk(SERIES[[1]]); h <- mk(SERIES[[2]]); b <- mk(SERIES[[3]])
  g <- function(x, col, f = "%.2f") if (nrow(x) && is.finite(x[[col]][1])) sprintf(f, x[[col]][1]) else "NA"
  sub <- sprintf("%s test prompts, %s with ideation\nrandom forest: AUC %s, r %s, MAE %s  |  hurdle model: AUC %s  |  training-mean baseline: AUC %s",
                 g(a, "n_obs", "%.0f"), g(a, "n_SI", "%.0f"), g(a, "AUC"), g(a, "r"), g(a, "MAE", "%.1f"), g(h, "AUC"), g(b, "AUC"))
  ymax <- max(20, ceiling(max(d$value, na.rm = TRUE) / 10) * 10)
  p1 <- ggplot(d, aes(day, value, colour = series, linetype = series)) +
    geom_line(linewidth = 0.45, alpha = 0.9) + geom_point(data = d[d$series == "Observed", , drop = FALSE], size = 0.9) +
    scale_colour_manual(values = COL_SERIES, name = NULL) + scale_linetype_manual(values = LTY_SERIES, name = NULL) +
    scale_x_continuous(limits = c(62, 91), breaks = c(62, 70, 80, 90)) + coord_cartesian(ylim = c(0, ymax)) +
    labs(x = "Study day", y = "Item 9 (0 to 100)", subtitle = sub) + theme_ms() + theme(legend.position = "top", plot.subtitle = element_text(size = BASE_PT - 1.5))
  auc_all <- m$AUC[m$key == SERIES[[1]] & is.finite(m$AUC)]; this <- if (nrow(a)) a$AUC[1] else NA_real_
  p2 <- ggplot(data.frame(AUC = auc_all), aes(AUC)) + geom_histogram(binwidth = 0.025, fill = "grey85", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey45", linewidth = 0.3) + coord_cartesian(xlim = c(0.2, 0.9)) +
    labs(x = "Within-person AUC of the combined dynamic random forest (all participants)", y = "Participants") + theme_ms()
  if (is.finite(this)) p2 <- p2 + geom_vline(xintercept = this, colour = COL_SERIES[[2]], linewidth = 0.8) +
    annotate("text", x = 0.9, y = Inf, hjust = 1, vjust = 1.4, size = pt2mm(BASE_PT - 1), family = FONT_FAMILY, colour = COL_SERIES[[2]],
             label = sprintf("this participant: AUC %.2f, rank %d of %d", this, 1L + sum(auc_all > this), length(auc_all)))
  if (nrow(b) && is.finite(b$AUC[1])) p2 <- p2 + geom_vline(xintercept = b$AUC[1], colour = "grey35", linetype = "dashed", linewidth = 0.5)
  (p1 / p2) + plot_layout(heights = c(3, 1.1)) + plot_annotation(title = label, theme = theme(plot.title = element_text(face = "bold", size = BASE_PT + 3, family = FONT_FAMILY)))
}
CFG_LABEL <- c("EMA" = "EMA mean", "EMA (dynamic)" = "EMA mean + dynamic", "Passive" = "Passive mean", "Passive (dynamic)" = "Passive mean + dynamic",
               "EMA + Passive" = "EMA + Passive mean", "EMA + Passive (dynamic)" = "EMA + Passive mean + dynamic")
COL_MODEL <- c("Observed" = "#D55E00", "Random forest" = "#0072B2", "Hurdle model" = "#009E73", "Training-mean baseline" = "#7F7F7F")
LTY_MODEL <- c("Observed" = "solid", "Random forest" = "solid", "Hurdle model" = "solid", "Training-mean baseline" = "dashed")
person_series_all <- function(p, uid) {                                          # observed and every model's prediction on the study-day axis, one block per configuration
  d <- p[p$uid == uid, , drop = FALSE]; if (!nrow(d)) return(NULL)
  ids <- unique(d[order(d$fold, d$row_id), c("fold", "row_id", "true")]); ids$pos <- ave(ids$fold, ids$fold, FUN = seq_along); ids$n <- ave(ids$fold, ids$fold, FUN = length)
  ids$day <- RUNIN_DAYS + TRAIN_DAYS + VAL_DAYS + ids$fold + (ids$pos - 1) / pmax(ids$n, EMA_PER_DAY)
  d <- merge(d, ids[, c("fold", "row_id", "day")], by = c("fold", "row_id"))
  d$model <- ifelse(grepl("^Continuous-RF", d$key), "Random forest", ifelse(grepl("^Hurdle-GBM", d$key), "Hurdle model", "Training-mean baseline"))
  pred <- d[d$model != "Training-mean baseline", c("day", "model", "Configuration", "pred")]; names(pred) <- c("day", "series", "cfg", "value")
  base <- d[d$model == "Training-mean baseline", , drop = FALSE]
  rep_cfg <- function(day, value, series) do.call(rbind, lapply(CFG, function(cf) data.frame(day = day, series = series, cfg = cf, value = value, stringsAsFactors = FALSE)))
  rbind(rep_cfg(ids$day, ids$true, "Observed"), if (nrow(base)) rep_cfg(base$day, base$pred, "Training-mean baseline"), pred)
}
fig_person_all <- function(uid, p, m, label) {
  d <- person_series_all(p, uid); if (is.null(d) || !nrow(d)) return(NULL)
  mk <- function(key) m[m$uid == uid & m$key == key, , drop = FALSE]
  g <- function(x, col, f = "%.2f") if (nrow(x) && is.finite(x[[col]][1])) sprintf(f, x[[col]][1]) else "NA"
  a <- mk(SERIES[[1]]); b <- mk(SERIES[[3]]); auc_all <- m$AUC[m$key == SERIES[[1]] & is.finite(m$AUC)]; this <- if (nrow(a)) a$AUC[1] else NA_real_
  rk <- if (is.finite(this)) sprintf("; combined dynamic random forest AUC rank %d of %d", 1L + sum(auc_all > this), length(auc_all)) else ""
  sub <- sprintf("%s test prompts, %s with ideation  |  training-mean baseline AUC %s%s", g(a, "n_obs", "%.0f"), g(a, "n_SI", "%.0f"), g(b, "AUC"), rk)
  lab <- vapply(CFG, function(cf) sprintf("%s\nrandom forest AUC %s | hurdle AUC %s", CFG_LABEL[[cf]], g(mk(paste0("Continuous-RF | ", cf)), "AUC"), g(mk(paste0("Hurdle-GBM | ", cf)), "AUC")), "")
  d$panel <- factor(lab[match(d$cfg, CFG)], levels = unname(lab)); d$series <- factor(d$series, levels = names(COL_MODEL))
  ymax <- max(20, ceiling(max(d$value, na.rm = TRUE) / 10) * 10)
  ggplot(d, aes(day, value, colour = series, linetype = series)) +
    geom_line(linewidth = 0.4, alpha = 0.9) + geom_point(data = d[d$series == "Observed", , drop = FALSE], size = 0.6) +
    facet_wrap(~ panel, ncol = 2) +
    scale_colour_manual(values = COL_MODEL, name = NULL) + scale_linetype_manual(values = LTY_MODEL, name = NULL) +
    scale_x_continuous(limits = c(62, 91), breaks = c(62, 75, 90)) + coord_cartesian(ylim = c(0, ymax)) +
    labs(title = label, subtitle = sub, x = "Study day", y = "Item 9 (0 to 100)") + theme_ms() +
    theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(size = BASE_PT - 1, face = "bold", hjust = 0), panel.spacing = grid::unit(3, "mm"),
          plot.title = element_text(face = "bold", size = BASE_PT + 3, family = FONT_FAMILY), plot.subtitle = element_text(size = BASE_PT - 1.5))
}
make_personalized <- function(m) {
  p <- read_predictions(); dir.create(PERSONALIZED_DIR, showWarnings = FALSE, recursive = TRUE)
  uids <- sort(unique(p$uid[p$series == names(SERIES)[1]])); if (!is.null(PARTICIPANTS)) uids <- intersect(uids, PARTICIPANTS)
  if (!length(uids)) stop("no participants with predictions from ", SERIES[[1]])
  lab <- function(k) if (SHOW_ID) sprintf("Participant %s", uids[k]) else sprintf("Participant %d of %d", k, length(uids))
  tag <- function(k) if (SHOW_ID) uids[k] else sprintf("%03d", k)
  message("Personalized figures: ", length(uids), " participants")
  pdf_dev <- function(f, w, h) if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(f, width = w, height = h, onefile = TRUE) else grDevices::pdf(f, width = w, height = h, onefile = TRUE)
  six <- identical(PERSONALIZED_CONFIGS, "all"); pg <- if (six) c(7, 8.2) else c(7, 5.2)
  figs <- lapply(seq_along(uids), function(k) tryCatch((if (six) fig_person_all else fig_person)(uids[k], p, m, lab(k)), error = function(e) { message("  skip ", uids[k], ": ", conditionMessage(e)); NULL }))
  keep <- which(!vapply(figs, is.null, logical(1)))
  f_all <- file.path(OUT_DIR, "Personalized_all_participants.pdf")                                                        # combined version 1: one page per participant
  pdf_dev(f_all, pg[1], pg[2]); for (k in keep) print(figs[[k]]); grDevices::dev.off()
  if (PERSONALIZED_PER_FILE) for (k in keep) tryCatch(ggplot2::ggsave(file.path(PERSONALIZED_DIR, sprintf("Personalized_PT%s.png", tag(k))), figs[[k]], width = pg[1], height = pg[2], dpi = 200, bg = "white"),
                                                      error = function(e) message("  could not write PNG for ", uids[k], ": ", conditionMessage(e)))
  message(sprintf("Personalized: %d of %d participants; combined PDF: %s", length(keep), length(uids), basename(f_all)))
  ov <- do.call(rbind, lapply(seq_along(uids), function(k) { d <- person_series(p, uids[k]); if (is.null(d)) return(NULL)       # combined version 2: everyone on one page
  d <- d[d$series %in% c("Observed", names(SERIES)[1]), , drop = FALSE]; d$k <- k; d }))
  auc <- m$AUC[match(paste(uids, SERIES[[1]]), paste(m$uid, m$key))]; ord <- order(-auc, na.last = TRUE)
  plab <- sprintf("%s  AUC %s", vapply(seq_along(uids), tag, ""), ifelse(is.finite(auc), sprintf("%.2f", auc), "NA"))
  ov$panel <- factor(plab[ov$k], levels = plab[ord]); ov$series <- factor(ov$series, levels = c("Observed", names(SERIES)[1]))
  po <- ggplot(ov, aes(day, value, colour = series)) + geom_line(linewidth = 0.25) + facet_wrap(~ panel, ncol = 12) +
    scale_colour_manual(values = COL_SERIES[c("Observed", names(SERIES)[1])], name = NULL) + coord_cartesian(xlim = c(62, 90), ylim = c(0, 100)) +
    labs(caption = "Item 9 (0 to 100) over study days 62 to 90; panels ordered by the AUC of the combined dynamic random forest") +
    theme_void(base_family = FONT_FAMILY, base_size = 5) + theme(strip.text = element_text(size = 4.5), legend.position = "top", plot.background = element_rect(fill = "white", colour = NA),
                                                                 panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.2), plot.caption = element_text(size = 5, hjust = 0))
  h_in <- max(6, 1.25 * ceiling(length(uids) / 12) + 1)
  ggplot2::ggsave(file.path(OUT_DIR, "Personalized_overview.png"), po, width = 14, height = h_in, dpi = 200, bg = "white")
  ggplot2::ggsave(file.path(OUT_DIR, "Personalized_overview.pdf"), po, width = 14, height = h_in, device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf")
  message("Personalized overview written (", length(uids), " panels)")
}
# ==============================================================================
# RUN
# ==============================================================================
step <- function(label, expr) tryCatch(expr, error = function(e) { message("[", label, "] skipped: ", conditionMessage(e)); NULL })
if (RUN_FIG[["Figure1"]]) step("Figure 1", make_fig1())
metrics <- if (any(RUN_FIG[c("Figure2", "Figure3")]) || RUN_PERSONALIZED) step("metrics", std_metrics(CSV$metrics)) else NULL
if (RUN_FIG[["Figure2"]]) step("Figure 2", make_fig2(if (is.null(metrics)) data.frame(key = character(0)) else metrics))
if (RUN_FIG[["Figure3"]] && !is.null(metrics)) step("Figure 3", make_fig3(metrics))
if (RUN_FIG[["Figure4"]]) step("Figure 4", make_fig4(CSV$contrasts))
if (RUN_PERSONALIZED && !is.null(metrics)) step("personalized figures", make_personalized(metrics))        # last: the personalized figures and the combined versions
message("Done. Files are in: ", normalizePath(OUT_DIR, mustWork = FALSE))