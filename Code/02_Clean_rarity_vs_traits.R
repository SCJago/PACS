# ============================================================
# TRAIT REPRESENTATION AMONG THREATENED ENSET LANDRACES
# ============================================================
#
# FILES NEEDED
# ----------------------
#   - landrace_avg_status_landscape.csv   (landscape rarity status per
#                                           landrace; output of the
#                                           clean_rarity_status.R script)
#   - Use_data_final.csv                  (trait/use scores per landrace)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Joins each landrace's landscape rarity status onto its recorded trait
# scores, then quantifies how common each useful trait (e.g. drought
# tolerance, kocho/bulla production, fodder, fibre, medicine) is among
# threatened landraces (Restricted/Endangered/Critical/Lost) - both as
# "present" (score >= 1) and "preferred" (score == 2).
#
# OUTPUTS GENERATED
# ----------------------
#   threatened_landrace_trait_summary.csv   one row per trait: how many
#       threatened landraces have trait data, how many have the trait
#       present/preferred, and the corresponding percentages.
#
# ============================================================


# ----------------------------
# 0. Set working directory
# ----------------------------
setwd("~/Papers/PACS/Landscape_rarity_vs_traits/clean_code")


# ----------------------------
# 1. Load packages
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
  library(stringr)
  library(tidyr)
})


# ----------------------------
# 2. Load data
# ----------------------------
status <- read_csv("landrace_avg_status_landscape.csv") %>%
  clean_names()

traits <- read_csv("Use_data_final.csv") %>%
  clean_names()

# Inspect names if needed
names(status)
names(traits)


# ----------------------------
# 3. Identify key columns
# ----------------------------

# Landrace name column in each data set - edit manually if your column
# names differ after clean_names()
status_landrace_col <- "landrace_name"
traits_landrace_col <- "landrace"

# Rarity category column, and which categories count as "threatened".
# Expected categories: Secure, Common, Restricted, Endangered, Critical, Lost
status_category_col <- "landscape_category"
threatened_categories <- c("Restricted", "Endangered", "Critical", "Lost")


# ----------------------------
# 4. Flag threatened landraces
# ----------------------------
status_threat <- status %>%
  mutate(
    landrace = .data[[status_landrace_col]],
    landscape_status = str_to_title(.data[[status_category_col]]),
    threatened = landscape_status %in% threatened_categories
  ) %>%
  select(landrace, landscape_status, threatened)


# ----------------------------
# 5. Prepare trait variables
# ----------------------------

# Edit this vector to match trait columns in use_data_final after
# clean_names(). These are likely candidates based on the manuscript.
trait_cols <- c(
  "drought",
  "kocho",
  "bulla",
  "fodder",
  "fibre",
  "medicine"
)

# Keep only trait columns that actually exist
trait_cols <- trait_cols[trait_cols %in% names(traits)]

if (length(trait_cols) == 0) {
  stop("No trait columns found. Check names(traits) and update trait_cols.")
}

traits_clean <- traits %>%
  mutate(landrace = .data[[traits_landrace_col]]) %>%
  select(landrace, all_of(trait_cols)) %>%
  distinct()


# ----------------------------
# 6. Join status and trait data
# ----------------------------
landrace_traits <- status_threat %>%
  left_join(traits_clean, by = "landrace")


# ----------------------------
# 7. Convert trait scores to presence/preference indicators
# ----------------------------

# Assumption from manuscript:
#   0 = not used / absent
#   1 = used
#   2 = preferred
#
# Trait present: score >= 1
# Trait preferred: score == 2
trait_summary_long <- landrace_traits %>%
  pivot_longer(
    cols = all_of(trait_cols),
    names_to = "trait",
    values_to = "score"
  ) %>%
  mutate(
    trait_present = !is.na(score) & score >= 1,
    trait_preferred = !is.na(score) & score == 2
  )


# ----------------------------
# 8. Trait representation among threatened landraces
# ----------------------------
threatened_trait_summary <- trait_summary_long %>%
  filter(threatened) %>%
  group_by(trait) %>%
  summarise(
    threatened_landraces_with_data = sum(!is.na(score)),
    threatened_landraces_trait_present = sum(trait_present, na.rm = TRUE),
    threatened_landraces_trait_preferred = sum(trait_preferred, na.rm = TRUE),
    pct_threatened_trait_present = 100 * threatened_landraces_trait_present /
      threatened_landraces_with_data,
    pct_threatened_trait_preferred = 100 * threatened_landraces_trait_preferred /
      threatened_landraces_with_data,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_threatened_trait_present))

print(threatened_trait_summary)


# ----------------------------
# 9. Save outputs
# ----------------------------
write_csv(threatened_trait_summary, "threatened_landrace_trait_summary.csv")
