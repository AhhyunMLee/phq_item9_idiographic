################################################################################
# SI Idiographic Prediction Project — STANDALONE VISUALIZATION SCRIPT
#
# Writes every idiographic (within-person) figure directly from the CSVs your
# analysis script already saves. It DOES NOT re-fit any model and needs none of
# the in-memory pipeline objects (results_main, data_model_base, configs, ...).
#
# Figures produced (PNG):
#   Continuous model
#     • Figure1_within_person_performance_continuous.png   (dashboard, best cfg)
#     • timeseries_continuous/Figure2_..._PT<uid>.png      (one per participant)
#   Hurdle model
#     • Figure1_within_person_performance_hurdle.png       (dashboard, best cfg)
#     • timeseries_hurdle/Figure2_..._PT<uid>.png          (one per participant)
#
# Inputs (point DATA_DIR at your "3. Cleaned data" folder):
#   model_comparison_participant_metrics.csv   (continuous dashboard)
#   model_comparison_predictions.csv           (continuous time series)
#   hurdle_participant_metrics.csv             (hurdle dashboard)
#   hurdle_predictions.csv                     (hurdle time series)
#
# Run:  Rscript si_idiographic_visualization.R
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(stringr); library(ggplot2); library(patchwork); library(scales)
})
# ==============================================================================
# CONFIG  — edit paths / options; everything else is automatic
# ==============================================================================
DATA_DIR   <- "."
OUTPUT_DIR <- file.path(DATA_DIR, "figures")
CSV <- list(
  cont_metrics = file.path(DATA_DIR, "model_comparison_participant_metrics.csv"),
  cont_preds   = file.path(DATA_DIR, "model_comparison_predictions.csv"),
  hurdle_metrics = file.path(DATA_DIR, "hurdle_participant_metrics.csv"),
  hurdle_preds   = file.path(DATA_DIR, "hurdle_predictions.csv")
)
RUN_CONTINUOUS <- TRUE     # set FALSE to skip continuous-model figures
RUN_HURDLE     <- TRUE     # set FALSE to skip hurdle-model figures
# Best config for each dashboard / caption. "auto" = best non-baseline config
# (hurdle: highest % AUC>0.5; continuous: highest % R²>0), ties broken by median.
# Or hard-set, e.g. BEST_HURDLE_CFG <- "EMA (dynamic)".
BEST_HURDLE_CFG <- "auto"
BEST_CONT_CFG   <- "auto"
# Sliding-window CV design (only to reconstruct "day in study" for Fig 2).
TRAIN_DAYS <- 45; VAL_DAYS <- 10; EMA_PER_DAY <- 3
PARTICIPANTS <- NULL       # NULL = every participant; or c("0300","0011")
OVERWRITE    <- TRUE
DPI          <- 220
# Cosmetics (identical to the reference figures) ------------------------------
BG_CREAM <- "#fbfaf4"; BORDER_GOLD <- "#eadb9a"
COL_OBS  <- "#c4514a"; COL_PRED <- "#2f6473"
# Configuration bookkeeping ---------------------------------------------------
CONFIG_LEVELS <- c("Baseline", "EMA", "EMA (dynamic)", "Passive",
                   "Passive (dynamic)", "EMA + Passive", "EMA + Passive (dynamic)")
CONFIG_FACET_ORDER <- CONFIG_LEVELS[-1]                       # 6 non-baseline configs
FACET_LABEL <- c(
  "EMA" = "EMA mean", "EMA (dynamic)" = "EMA mean + dyn",
  "Passive" = "Passive mean", "Passive (dynamic)" = "Passive mean + dyn",
  "EMA + Passive" = "EMA + Passive mean",
  "EMA + Passive (dynamic)" = "EMA + Passive\nmean + dyn")
TITLE_LABEL <- c(
  "EMA" = "EMA — Mean", "EMA (dynamic)" = "EMA — Mean + Dynamic",
  "Passive" = "Passive — Mean", "Passive (dynamic)" = "Passive — Mean + Dynamic",
  "EMA + Passive" = "EMA + Passive — Mean",
  "EMA + Passive (dynamic)" = "EMA + Passive — Mean + Dynamic")
# ==============================================================================
# HELPERS
# ==============================================================================
pad_uid <- function(x) stringr::str_pad(as.character(x), width = 4, pad = "0")
# rank-based (Mann–Whitney) AUC — no pROC dependency
auc_fast <- function(bin, score) {
  ok <- !is.na(bin) & !is.na(score); bin <- as.integer(bin[ok]); score <- score[ok]
  if (length(unique(bin)) < 2) return(NA_real_)
  n1 <- sum(bin == 1); n0 <- sum(bin == 0); r <- rank(score)
  (sum(r[bin == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
r2_signed <- function(true, pred) {
  ss_tot <- sum((true - mean(true, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - sum((true - pred)^2, na.rm = TRUE) / ss_tot
}
disp_scale <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(1)
  if (max(abs(v)) <= 1.5) 100 else 1
}
pick_col <- function(df, ...) { for (nm in c(...)) if (nm %in% names(df)) return(nm); NA_character_ }
std_metrics <- function(path) {
  m <- readr::read_csv(path, show_col_types = FALSE) |> as.data.frame()
  m$Configuration <- m[[pick_col(m, "Configuration", "config")]]
  m$uid <- pad_uid(m[[pick_col(m, "uid", "Participant")]])
  map <- c(AUC="AUC|auc", Sensitivity="Sensitivity|sensitivity",
           Specificity="Specificity|specificity",
           BalancedAcc="BalancedAcc|balanced_acc",
           r="^r$|Pearson|cor", R2="R2|r2|R\\^2", MAE="^MAE$|^mae$",
           RMSE="RMSE|rmse", n_SI="n_SI|n_si", n_nonSI="n_nonSI|n_nonsi",
           n_obs="n_obs")
  for (new in names(map)) if (!(new %in% names(m))) {
    hit <- grep(map[[new]], names(m), value = TRUE)
    hit <- setdiff(hit, c("Configuration","uid"))
    if (length(hit)) m[[new]] <- m[[hit[1]]]
  }
  m
}
std_preds <- function(path) {
  p <- readr::read_csv(path, show_col_types = FALSE,
                       col_types = readr::cols(uid = readr::col_character(),
                                               .default = readr::col_guess())) |> as.data.frame()
  p$Configuration <- p[[pick_col(p, "Configuration", "config")]]
  p$uid <- pad_uid(p[[pick_col(p, "uid", "Participant")]])
  if (!"true_bin" %in% names(p)) p$true_bin <- as.integer(p$true > 0)
  p
}
choose_best <- function(metrics, kind) {
  hard <- if (kind == "hurdle") BEST_HURDLE_CFG else BEST_CONT_CFG
  if (!identical(hard, "auto")) return(hard)
  mcol <- if (kind == "hurdle") "AUC" else "R2"; cut <- if (kind == "hurdle") 0.5 else 0
  metrics |>
    dplyr::filter(Configuration != "Baseline", is.finite(.data[[mcol]])) |>
    dplyr::group_by(Configuration) |>
    dplyr::summarise(pct = mean(.data[[mcol]] > cut, na.rm = TRUE),
                     med = median(.data[[mcol]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(pct), dplyr::desc(med)) |>
    dplyr::slice_head(n = 1) |> dplyr::pull(Configuration)
}
finite_rows <- function(df, col) df[df$metric == col & is.finite(df$value), , drop = FALSE]
# Right-hand per-metric DENSITY panels: smooth curve + rug + median line/label.
# Density is scaled to a peak of 1 (after_stat(scaled)) so the curve always fills
# the panel instead of being clipped when a metric is tightly concentrated.
# The density layer is only added when >= 2 finite points exist, so metrics with
# too few points no longer emit "groups with fewer than two data points" warnings.
density_panels <- function(long, order, palette, digits = 2, xlim = c(0, 1.2)) {
  lapply(order, function(mm) {
    dd  <- finite_rows(long, mm)
    med <- median(dd$value, na.rm = TRUE)
    dens_layer <- if (nrow(dd) >= 2)
      geom_density(aes(y = after_stat(scaled), fill = mm), alpha = 0.85,
                   color = "#555555", linewidth = 0.3) else NULL
    ggplot(dd, aes(x = value)) +
      dens_layer +
      geom_rug(color = palette[[mm]], alpha = 0.5, length = unit(0.05, "npc")) +
      geom_vline(xintercept = med, color = "#444444", linewidth = 0.4) +
      annotate("text", x = med, y = 0.9, label = sprintf(paste0("%.", digits, "f"), med),
               size = 2.9, hjust = -0.12) +
      annotate("text", x = -Inf, y = 0.5, label = mm, hjust = -0.05, size = 3.2) +
      scale_fill_manual(values = palette, guide = "none") +
      coord_cartesian(xlim = xlim, ylim = c(0, 1.05), clip = "off") +
      labs(x = NULL, y = NULL) + theme_void(base_size = 9) +
      theme(plot.margin = ggplot2::margin(2, 8, 2, 8))
  })
}
# ==============================================================================
# FIGURE 1 (HURDLE) — dashboard with SI / non-SI observation bars
# ==============================================================================
make_fig1_hurdle <- function(metrics, best_cfg, out_dir = OUTPUT_DIR) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  best_df <- metrics |>
    dplyr::filter(Configuration == best_cfg, is.finite(AUC), is.finite(Sensitivity),
                  is.finite(Specificity), is.finite(BalancedAcc)) |>
    dplyr::arrange(AUC) |>
    dplyr::mutate(rank = dplyr::row_number(),
                  n_SI = dplyr::coalesce(as.integer(n_SI), 0L),
                  n_nonSI = dplyr::coalesce(as.integer(n_nonSI), 0L))
  if (!nrow(best_df)) { message("Fig 1 hurdle: no rows for ", best_cfg); return(invisible(NULL)) }
  config_label <- ifelse(best_cfg %in% names(TITLE_LABEL), TITLE_LABEL[[best_cfg]],
                         gsub("_", " + ", best_cfg))
  pct_auc <- mean(best_df$AUC > 0.5, na.rm = TRUE) * 100
  metric_map <- c(AUC = "AUC", Sensitivity = "Sensitivity", Specificity = "Specificity",
                  BalancedAcc = "Balanced acc. (Youden threshold)")
  pal <- c("AUC" = "grey20", "Sensitivity" = "#7fb5a3", "Specificity" = "#e7b64b",
           "Balanced acc. (Youden threshold)" = "#e66a7a")
  dot_df <- best_df |> dplyr::select(rank, dplyr::all_of(names(metric_map))) |>
    tidyr::pivot_longer(dplyr::all_of(names(metric_map)), names_to = "metric_raw", values_to = "value") |>
    dplyr::mutate(metric = factor(metric_map[metric_raw], levels = unname(metric_map)))
  p_dot <- ggplot(dot_df, aes(value, rank, color = metric)) +
    geom_hline(yintercept = seq_len(nrow(best_df)), color = "#e7e5dc", linewidth = 0.25) +
    geom_vline(xintercept = 0.5, color = "grey55", linetype = "dashed") +
    geom_point(alpha = 0.75, size = 1.2) +
    scale_color_manual(values = pal) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    annotate("text", x = 0.04, y = max(dot_df$rank) * 0.06, label = sprintf("%.0f%%", pct_auc),
             size = 4, fontface = "bold") +
    labs(x = "Value", y = "Participants", color = NULL,
         title = sprintf("Within-person performance — %s (Hurdle model)", config_label),
         subtitle = sprintf(paste0("Binary classification  |  threshold: any SI (q9 > 0)\n",
                                   "%.0f%% of participants with an AUC > 0.5\nMean AUC (within-person) = %.2f\n",
                                   "Total observations: SI = %s | Non-SI = %s  |  N = %d participants"),
                            pct_auc, mean(best_df$AUC, na.rm = TRUE),
                            format(sum(best_df$n_SI), big.mark = ","),
                            format(sum(best_df$n_nonSI), big.mark = ","), nrow(best_df))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  # Per row: gray = total observations, red = SI portion (non-SI is the remainder).
  obs_df <- best_df |> dplyr::transmute(rank, n_nonSI, n_SI, n_total = n_nonSI + n_SI)
  p_obs <- ggplot(obs_df, aes(y = rank)) +
    geom_col(aes(x = n_total, fill = "Non-SI"), width = 0.85) +
    geom_col(aes(x = n_SI, fill = "SI"), width = 0.85) +
    scale_fill_manual(values = c("Non-SI" = "#bdbbb3", "SI" = "#e34a3a"), name = NULL) +
    labs(x = "Observations", y = NULL) + theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
  dens <- density_panels(dot_df, c("Specificity", "AUC", "Sensitivity",
                                   "Balanced acc. (Youden threshold)"), pal)
  p_density <- wrap_plots(dens, ncol = 1) +
    plot_annotation(caption = sprintf("N = %d participants.", nrow(best_df)),
                    theme = theme(plot.caption = element_text(hjust = 1, size = 8, color = "grey35")))
  fig <- p_dot + p_obs + p_density + plot_layout(widths = c(2.2, 1.15, 1.35)) +
    plot_annotation(theme = theme(
      plot.background = element_rect(fill = BG_CREAM, color = BORDER_GOLD, linewidth = 1.1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)))
  out <- file.path(out_dir, "Figure1_within_person_performance_hurdle.png")
  ggsave(out, fig, width = 10, height = 8.5, dpi = DPI, bg = BG_CREAM)
  message("Saved Fig 1 (hurdle): ", out); invisible(out)
}
# ==============================================================================
# FIGURE 1 (CONTINUOUS) — dashboard ranked by R², single test-observation bar
# ==============================================================================
make_fig1_continuous <- function(metrics, best_cfg, out_dir = OUTPUT_DIR) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  df <- metrics |>
    dplyr::filter(Configuration == best_cfg, is.finite(R2), is.finite(r), is.finite(MAE)) |>
    dplyr::arrange(R2) |> dplyr::mutate(rank = dplyr::row_number())
  if (!nrow(df)) { message("Fig 1 continuous: no rows for ", best_cfg); return(invisible(NULL)) }
  metric_map <- c(R2 = "R² (signed)", r = "Pearson r", MAE = "MAE (norm)")
  pal <- c("R² (signed)" = "#252525", "Pearson r" = "#7fb5a3", "MAE (norm)" = "#e7b64b")
  long <- df |> dplyr::select(rank, dplyr::all_of(names(metric_map))) |>
    tidyr::pivot_longer(dplyr::all_of(names(metric_map)), names_to = "metric_raw", values_to = "value") |>
    dplyr::mutate(metric = factor(metric_map[metric_raw], levels = unname(metric_map)))
  pct <- mean(df$R2 > 0, na.rm = TRUE) * 100
  x_left  <- max(-0.5, stats::quantile(long$value, 0.01, na.rm = TRUE))
  # Cap the right edge at 1.2 (not a hard 1) so the 99th-percentile term actually
  # has an effect and values slightly above 1 (e.g. MAE norm, Pearson r) are not
  # silently clipped off the panel.
  x_right <- min(1.2, max(1, stats::quantile(long$value, 0.99, na.rm = TRUE)))
  config_label <- ifelse(best_cfg %in% names(TITLE_LABEL), TITLE_LABEL[[best_cfg]],
                         gsub("_", " + ", best_cfg))
  p_dot <- ggplot(long, aes(value, rank, color = metric)) +
    geom_hline(yintercept = seq_len(nrow(df)), color = "#e7e5dc", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
    geom_point(alpha = 0.86, size = 1.15) +
    annotate("text", x = x_left + 0.05, y = max(long$rank) * 0.08, label = sprintf("%.0f%%", pct),
             size = 4, fontface = "bold") +
    scale_color_manual(values = pal, name = NULL) +
    coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
    labs(x = "Value", y = "Participants",
         title = sprintf("Within-person performance — %s (Continuous model)", config_label),
         subtitle = sprintf(paste0("Continuous regression | outcome: SI severity (q9)\n",
                                   "%.0f%% of participants with R² > 0\nMedian R² = %.3f | Median MAE = %.3f\n",
                                   "Total test observations: %s  |  N = %d participants"),
                            pct, median(df$R2, na.rm = TRUE), median(df$MAE, na.rm = TRUE),
                            format(sum(df$n_obs, na.rm = TRUE), big.mark = ","), nrow(df))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  p_obs <- ggplot(df, aes(x = n_obs, y = rank)) +
    geom_col(fill = "#a8a69e", width = 0.85) + labs(x = "Test observations", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  dens <- density_panels(long, unname(metric_map), pal, digits = 3, xlim = c(x_left, x_right))
  p_density <- wrap_plots(dens, ncol = 1) +
    plot_annotation(caption = sprintf("N = %d participants.", nrow(df)),
                    theme = theme(plot.caption = element_text(hjust = 1, size = 8, color = "grey35")))
  fig <- p_dot + p_obs + p_density + plot_layout(widths = c(2.2, 1.15, 1.35)) +
    plot_annotation(theme = theme(
      plot.background = element_rect(fill = BG_CREAM, color = BORDER_GOLD, linewidth = 1.1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)))
  out <- file.path(out_dir, "Figure1_within_person_performance_continuous.png")
  ggsave(out, fig, width = 10, height = 8.5, dpi = DPI, bg = BG_CREAM)
  message("Saved Fig 1 (continuous): ", out); invisible(out)
}
# ==============================================================================
# FIGURE 2 — per-participant observed vs predicted q9 time series
#   model_type %in% c("hurdle","continuous"); facet stat metric differs.
# ==============================================================================
add_day_in_study <- function(df) {
  df |> dplyr::group_by(Configuration, fold) |>
    dplyr::mutate(.p = dplyr::row_number(), .n = dplyr::n()) |> dplyr::ungroup() |>
    dplyr::mutate(day_in_study = TRAIN_DAYS + VAL_DAYS + fold + (.p - 1) / pmax(.n, EMA_PER_DAY)) |>
    dplyr::select(-.p, -.n)
}
panel_stats <- function(d, scl, model_type) {
  r    <- suppressWarnings(stats::cor(d$pred, d$true, use = "complete.obs"))
  rmse <- sqrt(mean((d$pred - d$true)^2, na.rm = TRUE)) * scl
  third <- if (model_type == "hurdle")
    auc_fast(d$true_bin, if ("p_si" %in% names(d)) d$p_si else d$pred) else r2_signed(d$true, d$pred)
  data.frame(cor = r, rmse = rmse, third = third)
}
make_fig2_participant <- function(uid, preds, model_type, out_dir) {
  uid <- pad_uid(uid)
  tag <- if (model_type == "hurdle") "hurdle" else "continuous"
  out <- file.path(out_dir, sprintf("Figure2_%s_timeseries_PT%s.png", tag, uid))
  if (file.exists(out) && !OVERWRITE) return(invisible(out))
  pd <- preds |> dplyr::filter(uid == !!uid, Configuration %in% CONFIG_FACET_ORDER)
  if (!nrow(pd)) return(invisible(NULL))
  pd  <- add_day_in_study(pd)
  scl <- disp_scale(c(pd$true, pd$pred))
  pd  <- pd |> dplyr::mutate(observed = true * scl, predicted = pred * scl)
  third_lab <- if (model_type == "hurdle") "AUC" else "R²"
  stats_df <- pd |> dplyr::group_by(Configuration) |>
    dplyr::group_modify(~ panel_stats(.x, scl, model_type)) |> dplyr::ungroup() |>
    dplyr::mutate(facet_label = sprintf("%s\nr = %.2f | RMSE = %.2f | %s = %.2f",
                                        FACET_LABEL[Configuration], cor, rmse, third_lab, third))
  lvl <- stats_df$facet_label[match(CONFIG_FACET_ORDER, stats_df$Configuration)]; lvl <- lvl[!is.na(lvl)]
  pd <- pd |> dplyr::left_join(stats_df |> dplyr::select(Configuration, facet_label), by = "Configuration") |>
    dplyr::mutate(facet_label = factor(facet_label, levels = lvl))
  best_row <- stats_df |> dplyr::filter(is.finite(third)) |> dplyr::arrange(dplyr::desc(third)) |>
    dplyr::slice_head(n = 1)
  best_txt <- if (nrow(best_row))
    sprintf("Best config: %s (%s = %.2f)", FACET_LABEL[best_row$Configuration], third_lab, best_row$third) else
      "Best config: not available"
  model_name <- if (model_type == "hurdle") "Hurdle model" else "Continuous model"
  day_rng <- range(pd$day_in_study, na.rm = TRUE)
  y_max <- max(100, ceiling(max(c(pd$observed, pd$predicted), na.rm = TRUE) / 20) * 20)
  plot_df <- pd |> dplyr::select(facet_label, day_in_study, observed, predicted) |>
    tidyr::pivot_longer(c(observed, predicted), names_to = "series", values_to = "q9") |>
    dplyr::mutate(series = factor(dplyr::recode(series, observed = "Observed", predicted = "Predicted"),
                                  levels = c("Observed", "Predicted"))) |>
    # Drop non-finite points explicitly so geom_line/geom_point no longer emit
    # "Removed rows containing missing values" warnings.
    dplyr::filter(is.finite(day_in_study), is.finite(q9))
  p <- ggplot(plot_df, aes(day_in_study, q9, color = series, group = series)) +
    geom_line(alpha = 0.48, linewidth = 0.45) + geom_point(size = 1.1, alpha = 0.9) +
    facet_wrap(~ facet_label, ncol = 2) +
    scale_color_manual(values = c("Observed" = COL_OBS, "Predicted" = COL_PRED), name = NULL) +
    scale_y_continuous(breaks = seq(0, y_max, 20)) +
    # coord_cartesian clips the view to [0, y_max] WITHOUT dropping out-of-range
    # points, so negative continuous predictions no longer trigger a "Removed
    # rows" warning for every participant.
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = sprintf("Within-person time series — Participant %s (%s)", uid, model_name),
         subtitle = sprintf(paste0("Predicted vs observed q9 over test days %.0f-%.0f ",
                                   "(sliding-window CV requires %d prior days) | %s"),
                            day_rng[1], day_rng[2], TRAIN_DAYS + VAL_DAYS, best_txt),
         x = "Day in study", y = "q9 (SI severity, 0-100)") +
    theme_minimal(base_size = 10) +
    theme(plot.background = element_rect(fill = BG_CREAM, color = NA),
          panel.background = element_rect(fill = BG_CREAM, color = NA),
          legend.position = "top", plot.title = element_text(face = "bold", size = 13),
          strip.text = element_text(face = "bold", size = 8), panel.grid.minor = element_blank())
  ggsave(out, p, width = 14, height = 8, dpi = DPI, bg = BG_CREAM)
  invisible(out)
}
make_fig2_all <- function(preds, model_type, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  uids <- preds |> dplyr::filter(Configuration %in% CONFIG_FACET_ORDER) |>
    dplyr::pull(uid) |> unique() |> pad_uid()
  if (!is.null(PARTICIPANTS)) uids <- intersect(uids, pad_uid(PARTICIPANTS))
  uids <- sort(uids)
  message(sprintf("Fig 2 (%s): rendering %d participants -> %s", model_type, length(uids), out_dir))
  ok <- 0L
  for (u in uids) {
    res <- tryCatch(make_fig2_participant(u, preds, model_type, out_dir),
                    error = function(e) { message("  skip PT ", u, ": ", conditionMessage(e)); NULL })
    if (!is.null(res)) ok <- ok + 1L
  }
  message(sprintf("Fig 2 (%s): saved %d / %d.", model_type, ok, length(uids)))
}
# ==============================================================================
# RUN
# ==============================================================================
run_model <- function(kind) {
  mfile <- if (kind == "hurdle") CSV$hurdle_metrics else CSV$cont_metrics
  pfile <- if (kind == "hurdle") CSV$hurdle_preds   else CSV$cont_preds
  if (!file.exists(mfile) || !file.exists(pfile)) {
    message(sprintf("[%s] missing CSV(s): %s / %s — skipping.", kind, mfile, pfile)); return(invisible())
  }
  metrics <- std_metrics(mfile); preds <- std_preds(pfile)
  best <- choose_best(metrics, kind); message(sprintf("[%s] best config: %s", kind, best))
  if (kind == "hurdle") make_fig1_hurdle(metrics, best) else make_fig1_continuous(metrics, best)
  make_fig2_all(preds, kind, file.path(OUTPUT_DIR, sprintf("timeseries_%s", kind)))
}
if (RUN_CONTINUOUS) run_model("continuous")
if (RUN_HURDLE)     run_model("hurdle")
message("Done. PNGs are in: ", normalizePath(OUTPUT_DIR))