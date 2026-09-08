# PACS Analysis Code

Analysis code evaluating a Payments for Agrobiodiversity Conservation
Services (PACS) scheme supporting in-situ conservation of enset
(*Ensete ventricosum*) landrace diversity in southwest Ethiopia.

## Results sections

Which script(s) produced the analysis behind each results section.

| Results section | Script(s) |
|---|---|
| Conservation threat assessment (landrace counts, rarity/loss, trait representation among threatened landraces) | `01`, `02` |
| PACS increased community landrace richness relative to counterfactual expectations | `03`, `04` |
| PACS contracts improved conservation outcomes for threatened landraces | `06` |
| Programme effects extended beyond contracted landraces through crowding in | `07` |
| Bid heterogeneity explains farmer opportunity costs | `08` |
| Estimating the cost of conserving range-wide enset genetic resources | `09` |

`05_clean_model_offer_selection.R` isn't tied to a single results section — it
prepares `select_offers.csv`, used as an input by `06`–`09`.

## Structure

```
PACS/
├── Code/     scripts, numbered in run order (01-09)
└── Data/     input data files
```

Scripts must be run in numeric order — later scripts depend on outputs
produced by earlier ones (see **Pipeline** below). Each script's `setwd()`
points to a local path and needs updating to your own copy of `PACS/Code`.

## Requirements

R 4.5.1

```r
install.packages(c(
  "dplyr", "readr", "tidyr", "stringr", "tibble", "janitor", "purrr",
  "ggplot2", "patchwork", "ggspatial", "ggpubr", "scales",
  "nlme", "lme4", "MuMIn", "emmeans", "ggeffects", "fixest",
  "modelsummary", "broom", "ompr", "ompr.roi", "ROI", "ROI.plugin.glpk",
  "corrplot", "Hmisc", "car",
  "sf", "terra", "rnaturalearth", "rnaturalearthdata"
))
```

## Data

Original input files (in `PACS/Data`) — not files produced by another script
in this pipeline:

| File | Description |
|---|---|
| `T1_T2_rarity_assessment_data.csv` | Community x landrace survey data, T1 baseline + T2 follow-up |
| `Use_data_final.csv` | Trait/use scores per landrace |
| `clean_control_hh_data.csv` | Control household data: T1/T2 landrace lists, coordinates, years |
| `clean_offers.csv` | Cleaned tender offers |
| `gender_balance.csv` | Community, farmers, women |
| `landrace_name_id_crosswalk.csv` | Landrace_name <-> Landrace_id_tender |
| `community_name_crosswalk.csv` | Community <-> Community_name |
| `Landrace_size_data_clean.csv` | Landrace-level plant measurements: age, maturity age, pseudostem height and circumference |
| `enset_per_farmer_by_community.csv` | Mean enset plants per farmer by community. Produced by `community_enset_per_farmer.R` — not part of this numbered set |
| `spatial_data.zip` | See below |

**`spatial_data.zip`** must be unzipped into `PACS/Data` before running
script `03`. Contains `surveyed_kebeles_PACS.shp` (+ its `.shx`/`.dbf`/`.prj`
sidecar files), `Rift_valley_boundary.kml`, and `agroecological_zones.tif`.

## Pipeline

| Script | What it does | Outputs reused by another script | Reused by |
|---|---|---|---|
| `01_clean_rarity_status.R` | Classifies landrace conservation threat/rarity status from survey data | `landrace_avg_status_landscape.csv` | 02, 05, 06 |
| `02_clean_rarity_vs_traits.R` | Summarises which use traits are represented among threatened/lost landraces | — | — |
| `03_clean_control_household_sample_prep.R` | Builds the 35-37°E control household bootstrap pool; agroecological zone map | `hh_boot_35_37E.rds` | 04 |
| `04_clean_bootstrap_analysis_code.R` | Bootstrap comparison of community-level landrace richness (per-household and total unique) against pseudo-control communities | — | — |
| `05_clean_model_offer_selection.R` | Cleans tender bid/offer data; flags which offers were selected | `select_offers.csv` | 06, 07, 08, 09 |
| `06_clean_landrace_rarity_effect_analysis.R` | DiD analysis of PACS contract effects on threatened landraces (rarity status, household adoption, plant abundance) within treated communities | `landrace_effect_analysis_df_long.csv` | 07 |
| `07_clean_crowding_in_out_analysis.R` | Tests whether programme coverage produced spillover (crowding-in) effects on non-contracted landraces | — | — |
| `08_clean_bid_heterogeneity.R` | Analyses bid heterogeneity to identify drivers of farmers' compensation demands | — | — |
| `09_clean_scaling_up_costs.R` | Projects the national cost of scaling PACS-style incentives to conserve enset genetic resources range-wide | — | — |
