# PHQ-9 Item 9 Prediction Pipeline

**Project:** Within-person prediction of PHQ-9 item 9 / suicidal ideation  
**Author:** Ahhyun Lee  
**Last updated:** June 2026

This repository contains the final R code for the PHQ-9 item 9 prediction project. The code is organized as a single script that can rebuild cleaned features, run cached or full prediction models, generate final manuscript tables/figures, and optionally run exploratory SHAP feature-importance summaries.

## Main Script

Use the single analysis script:

```r
source("PHQ9_Item9_Predidction.R")
```

The filename intentionally matches the current project file. If you rename the file, update the `source()` command accordingly.

## Expected Folder Layout

The script expects the study folder to use this layout:

```text
SI_Prediction/
├── 1. Code/
├── 2. Raw data/
├── 3. Cleaned data/
└── 4. Results/
```

If the script is run from `1. Code`, it automatically treats the parent folder as the project root. Otherwise, set the project root explicitly before sourcing:

```r
Sys.setenv(SI_PROJECT_BASE = "/path/to/SI_Prediction")
source("PHQ9_Item9_Predidction.R")
```

## Run Flags

The script is submitted with all optional analyses turned off. To run a section, edit the `RUN_*` flags at the top of `PHQ9_Item9_Predidction.R` before sourcing the file.

For manuscript output generation using existing cleaned data and cached model results, use:

```r
RUN_DATA_CLEANING <- FALSE
RUN_MAIN_PIPELINE <- TRUE
RUN_FINAL_MANUSCRIPT_OUTPUTS <- TRUE
RUN_FEATURE_IMPORTANCE <- FALSE
RUN_ALL_PT_CONTINUOUS_TIMESERIES <- FALSE
RUN_ALL_PT_HURDLE_TIMESERIES <- FALSE
RUN_EXTRA_DIAGNOSTIC_FIGURES <- FALSE
RUN_BASE_RATE_ANALYSIS <- FALSE
```

With these defaults, the script loads model-ready files and cached model results when available, then writes the final report outputs.

## Rebuilding Data

To rebuild cleaned EMA/passive-sensing features from raw files, set:

```r
RUN_DATA_CLEANING <- TRUE
```

The cleaning pipeline:

- builds EMA features;
- computes passive-sensing dynamic indices;
- drops stress features;
- applies EMA compliance, passive sensing compliance, and SI-count inclusion criteria;
- drops all-NA feature columns after participant inclusion.

Dynamic-index columns are **not** filtered using a 75% observed or 75% zero rule. Near-zero-variance predictors are removed within training folds during modeling using `caret::nearZeroVar()`.

## Prediction Models

The script evaluates:

- continuous random forest models for SI severity;
- hurdle models for binary SI occurrence and severity conditional on SI presence;
- six feature configurations plus a baseline benchmark:
  - EMA
  - EMA (dynamic)
  - Passive
  - Passive (dynamic)
  - EMA + Passive
  - EMA + Passive (dynamic)

Models use sliding-window within-person cross-validation. Time-window selection is optimized within training/validation folds, and near-zero-variance filtering is applied within training folds to avoid leakage.

## Final Report Outputs

Main manuscript outputs are written to:

```text
4. Results/report/
```

Expected final report files include:

```text
Table1_demographics.csv
Table2_continuous_model_results_main.csv
Table2_continuous_model_results_main.docx
Table3_hurdle_model_results_main.csv
Table3_hurdle_model_results_main.docx
Table4_feature_importance_shap.csv
feature_importance_all_in_one_aggregate.csv
Figure_hurdle_model_performance_main.png
```

Optional participant-level time-series figures are written to:

```text
4. Results/report/all_PT/continuous/
4. Results/report/all_PT/hurdle/
```

Enable these with:

```r
RUN_ALL_PT_CONTINUOUS_TIMESERIES <- TRUE
RUN_ALL_PT_HURDLE_TIMESERIES <- TRUE
```

## Feature Importance

Exploratory feature-importance summaries use mean absolute approximate SHAP values among successful models:

- continuous models with `R² > 0`;
- hurdle models with `AUC > 0.5`.

To run feature importance:

```r
RUN_FEATURE_IMPORTANCE <- TRUE
```

By default, detailed SHAP working files are written to a temporary folder and deleted after the clean aggregate files are created, so the Google Drive report folder does not fill with intermediate `importance_*`, `summary_*`, `participants_*`, and `top_features_*` files. Final feature-importance outputs are:

```text
4. Results/report/Table4_feature_importance_shap.csv
4. Results/report/feature_importance_all_in_one_aggregate.csv
```

To keep the detailed diagnostic files (`importance_*`, `summary_*`, `participants_*`, `top_features_*`), set:

```r
KEEP_FEATURE_IMPORTANCE_DETAILS <- TRUE
```

Older diagnostic plots such as `F1_model_performance_main.png`, `F2_withinperson_continuous_main.png`, `F3_si_obs_counts_main.png`, `F8_demographics.png`, and `Figure_demographic_description.png` are not final report outputs and are removed during the final report build. To intentionally generate extra diagnostic figures, set:

```r
RUN_EXTRA_DIAGNOSTIC_FIGURES <- TRUE
```

## Optional Base-Rate Analysis

To quantify zero-inflation/restricted range in item 9 and test whether within-person predictability scales with the participant-level item 9 base rate, provide:

```text
4. Results/report/raw_outcome.csv
4. Results/report/perperson_performance.csv
```

Then run:

```r
RUN_BASE_RATE_ANALYSIS <- TRUE
source("PHQ9_Item9_Predidction.R")
```

This optional section writes:

```text
4. Results/report/perperson_baserate.csv
4. Results/report/baserate_predictability_results.csv
4. Results/report/Figure_base_rate_histogram.png
4. Results/report/Figure_base_rate_vs_predictability.png
```

## Privacy And GitHub Upload

This repository should contain code and documentation only. Do not upload participant-level data, raw data, cleaned data, model caches, generated result tables, or participant-level figures unless your institution has explicitly approved public sharing.

Do not commit:

- `2. Raw data/`
- `3. Cleaned data/`
- `4. Results/`
- `.RData`
- `.RDataTmp`
- `.Rhistory`
- `.Rproj.user/`
- model cache files
- participant-level figures or CSVs

Recommended public files:

```text
PHQ9_Item9_Predidction.R
PHQ9_Item9_Prediction.Rproj
README.md
```

## Reproducibility Notes

The script installs missing R packages when needed. For a clean rerun on another machine, make sure the expected folder layout exists and set `SI_PROJECT_BASE` if the script is not run from `1. Code`.
