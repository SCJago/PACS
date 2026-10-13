# ============================================================
# PACS CROWDING-IN / CROWDING-OUT SPILLOVER ANALYSIS
# ============================================================
#
# FILES NEEDED
# ----------------------
#   - landrace_effect_analysis_df_long.csv   
#     produced by script 06
#   - select_offers.csv
#     produced by script 05
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Tests whether PACS programme coverage generated indirect (spillover)
# effects on non-contracted landraces:
#   1. Builds community-level programme intensity from accepted bids:
#      number of landraces contracted, total payment, and payment per
#      farmer.
#   2. Builds community-level programme coverage measures at
#      baseline: share of landraces contracted, number contracted,
#      and cross-checks that count against the bids file directly.
#   3. Checks whether programme coverage was systematically related
#      to baseline conservation conditions (Spearman correlations),
#      to rule out selection effects driving spillover estimates.
#   4. Restricts to non-contracted landraces, and estimates linear
#      fixed-effects models (community + landrace FE, SEs clustered
#      by community) testing whether post-implementation changes in
#      rarity status, households cultivating, plant abundance, for 
#      non-contracted landraces differed with programme coverage.
#   5. Repeats the spillover models using two alternative programme
#      intensity measures (count of contracts, GBP per farmer) as a
#      robustness check.
#
# Positive effects = crowding-in (greater PACS exposure improves
# outcomes for non-targeted landraces). Negative effects =
# crowding-out (contracts displace non-target landrace cultivation).
#
# OUTPUTS GENERATED
# ----------------------
#   COMMUNITY PAYMENTS SUMMARY: landraces contracted, total payment,
#       and payment per farmer, by community
#       -> community_payments_summary.csv
#
#   BASELINE SELECTION CHECK: Spearman correlations between programme
#       coverage and baseline conservation conditions
#       (Supplementary Table 3) -> baseline_selection_correlations.csv
#
#   PRIMARY SPILLOVER MODELS (share-based programme coverage):
#       rarity status, households, plant abundance, disappearance for
#       non-contracted landraces -> models_spillover_share_paid.html
#
#   ROBUSTNESS SPILLOVER MODELS (count-based and GBP-based programme
#       intensity) -> models_spillover_no_paid.html,
#       models_spillover_GBP_paid.html (Supplementary Table 4)
#
#   ALL SPILLOVER MODEL COEFFICIENTS: tidied coefficients for every
#       spillover model (share, count, GBP; all four outcomes)
#       -> model_coefficients_all.csv
#
# ============================================================

setwd("~/Papers/PACS/Crowding_in_out_analysis/Clean_code")

# ----------------------------
# 0. Load packages
# ----------------------------

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(broom)
library(fixest)
library(modelsummary)


# ----------------------------
# 1. Load data
# ----------------------------
landrace_panel <- read_csv("landrace_effect_analysis_df_long.csv", show_col_types = FALSE)
selected_bids <- read_csv("select_offers.csv", show_col_types = FALSE)


# ----------------------------
# 2. Check required columns are present
# ----------------------------

required_cols <- c("Community", "Landrace_ID", "Timepoint", "paid",
                   "rarity_status", "no_households", "avr_plants", "present")
missing_cols <- setdiff(required_cols, names(landrace_panel))
if (length(missing_cols) > 0) {
  stop("landrace_effect_analysis_df_long.csv is missing expected column(s): ",
       paste(missing_cols, collapse = ", "))
}

required_bid_cols <- c("Community", "Landrace_ID", "GBP", "No_farmers")
missing_bid_cols <- setdiff(required_bid_cols, names(selected_bids))
if (length(missing_bid_cols) > 0) {
  stop("select_offers.csv is missing expected column(s): ",
       paste(missing_bid_cols, collapse = ", "))
}


# ----------------------------
# 3. Clean
# ----------------------------
# - trim stray whitespace in identifier columns 
# - coerce Timepoint and paid to plain 0/1 integers, in case they were
#   read in as text, logical, or factor
# - drop rows with no landrace or community identifier, since these
#   can't be used in any grouped calculation below
landrace_panel <- landrace_panel %>%
  mutate(
    Community = str_trim(as.character(Community)),
    Landrace_ID  = str_trim(as.character(Landrace_ID)),
    Timepoint = as.integer(Timepoint),
    paid      = as.integer(paid)
  ) %>%
  filter(!is.na(Community), !is.na(Landrace_ID), Community != "", Landrace_ID != "")

cat("Rows after cleaning:", nrow(landrace_panel), "\n")
cat("Communities:", n_distinct(landrace_panel$Community), "\n")
cat("Landraces:", n_distinct(landrace_panel$Landrace_ID), "\n")
cat("Timepoints present:", paste(sort(unique(landrace_panel$Timepoint)), collapse = ", "), "\n")


# ============================================================
# COMMUNITY PROJECT INTENSITY (from accepted bids)
# ============================================================

# ----------------------------
# 4. Community project intensity
# ----------------------------
# n_landraces_selected: number of distinct landraces under contract in
#     the community, from the accepted bids.
# community_total_GBP: total PACS payment accepted in the community.
# community_GBP_per_farmer: total payment divided by the number of
#     farmers in the community - the programme intensity measure used
#     as a robustness check on programme coverage below.
community_payments <- selected_bids %>%
  mutate(
    Community = as.character(Community) %>% str_trim() %>% str_squish(),
    Landrace_ID = as.character(Landrace_ID) %>% str_trim() %>% str_squish(),
    GBP = as.numeric(GBP),
    No_farmers = as.numeric(No_farmers)
  ) %>%
  group_by(Community) %>%
  summarise(
    n_landraces_selected = n_distinct(Landrace_ID),
    community_total_GBP = sum(GBP, na.rm = TRUE),
    No_farmers = No_farmers[!is.na(No_farmers)][1],
    community_GBP_per_farmer = if_else(!is.na(No_farmers) & No_farmers > 0, community_total_GBP / No_farmers, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(desc(community_total_GBP))

print(community_payments)

write.csv(community_payments, "community_payments_summary.csv", row.names = FALSE)

# Attach community project intensity to the panel
landrace_panel <- landrace_panel %>%
  left_join(
    community_payments %>% select(Community, n_landraces_selected, community_total_GBP, community_GBP_per_farmer),
    by = "Community"
  )


# ============================================================
# COMMUNITY PROGRAMME COVERAGE MEASURES
# ============================================================

# ----------------------------
# 5. Community-level programme coverage, measured at baseline
# ----------------------------
# share_paid: proportion of landraces in the community under PACS
#     contract - the primary programme coverage measure used in text.
# no_paid: absolute number of contracted landraces in the community -
#     used as a robustness check on the coverage measure.
community_coverage <- landrace_panel %>%
  filter(Timepoint == 0) %>%
  group_by(Community) %>%
  summarise(
    share_paid = mean(paid, na.rm = TRUE),
    no_paid = sum(paid, na.rm = TRUE),
    .groups = "drop"
  )

landrace_panel_coverage <- landrace_panel %>%
  left_join(community_coverage, by = "Community")


# ----------------------------
# 5b. Consistency check: contract counts from panel vs bids file
# ----------------------------
# no_paid (from the paid flag already attached to the panel) and
# n_landraces_selected (from select_offers.csv directly) should agree
# for every community, since paid was originally derived from the
# same bids file upstream. Flags any community where they don't, so a
# mismatch (e.g. a name mismatch dropped through a join, or the bids
# file updated after the panel was last built) doesn't silently feed
# a wrong count into the COUNT-based robustness models below.
contract_count_check <- community_coverage %>%
  select(Community, no_paid) %>%
  full_join(
    community_payments %>% select(Community, n_landraces_selected),
    by = "Community"
  ) %>%
  mutate(matches = no_paid == n_landraces_selected)

contract_count_mismatches <- contract_count_check %>% filter(!matches | is.na(matches))

if (nrow(contract_count_mismatches) > 0) {
  warning(
    "Contract counts disagree between the panel (no_paid) and select_offers.csv ",
    "(n_landraces_selected) for ", nrow(contract_count_mismatches), " community/communities - ",
    "see contract_count_mismatches."
  )
  print(contract_count_mismatches)
} else {
  cat("Contract count check: no_paid matches n_landraces_selected for all communities.\n")
}


# ----------------------------
# 6. Baseline landrace diversity per community
# ----------------------------
baseline_diversity <- landrace_panel %>%
  filter(Timepoint == 0) %>%
  group_by(Community) %>%
  summarise(
    baseline_landrace_count = n_distinct(Landrace_ID[present == 1]),
    .groups = "drop"
  )

landrace_panel_coverage <- landrace_panel_coverage %>%
  left_join(baseline_diversity, by = "Community")


# ============================================================
# BASELINE SELECTION CHECK
# ============================================================
# Checks whether programme coverage was systematically related to
# baseline community characteristics, as an additional robustness
# check.

# ----------------------------
# 7. Community-level baseline conservation conditions
# ----------------------------
# Built from NON-CONTRACTED landraces only (paid == 0) at baseline,
# since these are the pre-existing conditions programme coverage is
# being checked against.
baseline_community_outcomes <- landrace_panel_coverage %>%
  filter(Timepoint == 0, paid == 0) %>%
  group_by(Community) %>%
  summarise(
    baseline_mean_rarity = mean(rarity_status, na.rm = TRUE),
    baseline_mean_hh = mean(no_households, na.rm = TRUE),
    baseline_presence_rate = mean(present, na.rm = TRUE),
    baseline_mean_plants = mean(avr_plants, na.rm = TRUE),
    baseline_landrace_count = first(baseline_landrace_count),
    community_GBP_per_farmer = first(community_GBP_per_farmer),
    share_paid = first(share_paid),
    no_paid = first(no_paid),
    .groups = "drop"
  )

print(baseline_community_outcomes)


# ----------------------------
# 8. Spearman correlations: programme coverage vs baseline conditions
# ----------------------------
# Spearman's rank correlations between the share of landraces
# contracted in each community and baseline conservation conditions -
# the test reported in the manuscript (Supplementary Table 3).
run_spearman <- function(x, y, label) {
  test <- cor.test(x, y, method = "spearman", exact = FALSE)
  tibble(
    baseline_variable = label,
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

baseline_selection_correlations <- bind_rows(
  run_spearman(baseline_community_outcomes$share_paid, baseline_community_outcomes$baseline_mean_rarity, "Mean rarity score"),
  run_spearman(baseline_community_outcomes$share_paid, baseline_community_outcomes$baseline_mean_hh, "Mean households per landrace"),
  run_spearman(baseline_community_outcomes$share_paid, baseline_community_outcomes$baseline_mean_plants, "Mean plant abundance"),
  run_spearman(baseline_community_outcomes$share_paid, baseline_community_outcomes$baseline_landrace_count, "Baseline landrace diversity")
)

print(baseline_selection_correlations)

write_csv(
  baseline_selection_correlations,
  "baseline_selection_correlations.csv"
)


# ============================================================
# SPILLOVER MODELS: NON-CONTRACTED LANDRACES ONLY
# ============================================================
# Spillover analysis, using only landraces that did not receive
# contracts within each community.

# ----------------------------
# 9. Restrict to non-contracted landraces
# ----------------------------
non_contracted_panel <- landrace_panel_coverage %>% filter(paid == 0)


# ----------------------------
# 10. Primary spillover models (share-based programme coverage)
# ----------------------------
# Tests whether post-implementation changes in rarity status,
# households cultivating, and plant abundance for non-contracted
# landraces were related to the community-level share of landraces
# under contract. Community and landrace fixed effects control for
# time-invariant differences; standard errors clustered by community.
m_crowd_rarity_SHARE <- feols(
  rarity_status ~ Timepoint:share_paid | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)

m_crowd_hh_SHARE <- feols(
  no_households ~ Timepoint:share_paid | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)

m_crowd_plants_SHARE <- feols(
  avr_plants ~ Timepoint:share_paid | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)


# ----------------------------
# 12. Export primary (share-based) spillover model table
# ----------------------------
models_spill_share <- list(
  "Rarity status" = m_crowd_rarity_SHARE,
  "Households"    = m_crowd_hh_SHARE,
  "Plants"        = m_crowd_plants_SHARE
)

modelsummary(
  models_spill_share,
  output = "models_spillover_share_paid.html",
  stars = TRUE,
  statistic = "({std.error})",
  gof_omit = "IC|Log|Adj|Within|Std.Errors|FE"
)


# ============================================================
# ROBUSTNESS: ALTERNATIVE PROGRAMME INTENSITY MEASURES
# ============================================================
# Re-estimates the spillover models using two alternative programme
# intensity measures - absolute number of contracted landraces, and
# average PACS payment per farmer - to confirm findings are not
# driven by the choice of coverage measure (Supplementary Table 4).

# ----------------------------
# 13. Robustness models: absolute count of contracts
# ----------------------------
m_crowd_rarity_COUNT <- feols(
  rarity_status ~ Timepoint:no_paid | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)

m_crowd_hh_COUNT <- feols(
  no_households ~ Timepoint:no_paid | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)

m_crowd_plants_COUNT <- feols(
  avr_plants ~ Timepoint:no_paid | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)



models_spill_count <- list(
  "Rarity status" = m_crowd_rarity_COUNT,
  "Households"    = m_crowd_hh_COUNT,
  "Plants"        = m_crowd_plants_COUNT
)

modelsummary(
  models_spill_count,
  output = "models_spillover_no_paid.html",
  stars = TRUE,
  statistic = "({std.error})",
  gof_omit = "IC|Log|Adj|Within|Std.Errors|FE"
)


# ----------------------------
# 14. Robustness models: GBP per farmer
# ----------------------------
m_crowd_rarity_GBP <- feols(
  rarity_status ~ Timepoint:community_GBP_per_farmer | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)

m_crowd_hh_GBP <- feols(
  no_households ~ Timepoint:community_GBP_per_farmer | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)

m_crowd_plants_GBP <- feols(
  avr_plants ~ Timepoint:community_GBP_per_farmer | Community + Landrace_ID,
  data = non_contracted_panel,
  cluster = ~Community
)


models_spill_gbp <- list(
  "Rarity status" = m_crowd_rarity_GBP,
  "Households"    = m_crowd_hh_GBP,
  "Plants"        = m_crowd_plants_GBP
)

modelsummary(
  models_spill_gbp,
  output = "models_spillover_GBP_paid.html",
  stars = TRUE,
  statistic = "({std.error})",
  gof_omit = "IC|Log|Adj|Within|Std.Errors|FE"
)


# ----------------------------
# 15. Export all spillover model coefficients
# ----------------------------
# Tidied coefficients for every spillover model (share, count, GBP;
# all four outcomes), for cross-checking against the manuscript text.
coefs_share <- map_dfr(models_spill_share, broom::tidy, .id = "model") %>%
  mutate(measure = "share")

coefs_count <- map_dfr(models_spill_count, broom::tidy, .id = "model") %>%
  mutate(measure = "count")

coefs_gbp <- map_dfr(models_spill_gbp, broom::tidy, .id = "model") %>%
  mutate(measure = "GBP")

coefs_all <- bind_rows(coefs_share, coefs_count, coefs_gbp)

write_csv(coefs_all, "model_coefficients_all.csv")
