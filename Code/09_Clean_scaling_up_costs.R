# ===============================================================
# SCALING UP CONSERVATION COSTS
# ===============================================================
# Author: Sophie Jago
#
# FILES NEEDED
# ----------------------
#   - T1_T2_rarity_assessment_data_share.csv 
#   - select_offers.csv                       (produced in script 05)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Estimates the household-level compensation required to scale PACS-
# style incentives to nationwide in-situ conservation of enset
# diversity, combining observed intervention data with national
# landrace richness estimates (Borrell et al. 2019).
#
#   1. Household gains: change in cultivating households per
#      community between baseline (T1) and follow-up (T2).
#   2. Community-level total paid: total compensation actually paid
#      per community, summed directly from select_offers.csv's 
#   3. Cost per additional household: total paid per community / hh
#      gain, summarised as the interquartile range (25th/50th/75th
#      percentile) across communities to capture cost heterogeneity.
#   4. National scenarios: crosses landrace-richness uncertainty
#      (Borrell et al. 2019: 475 conservative / 1270 upper) with the
#      three cost scenarios to produce a scenario matrix of total
#      national compensation required.
#
# OUTPUTS GENERATED
# ----------------------
# All saved to the working directory.
#
#   community_total_paid.csv        Total compensation paid (ETB/USD)
#       per community, derived from selected offers - step 3.
#
#   community_cost_per_household.csv  Cost per additional cultivating
#       household, per community - step 4.
#
#   cost_scenarios.csv              Interquartile cost-per-household
#       scenarios (lower/median/upper) - step 5.
#
#   scaling_estimates_rounded.csv   Full scenario matrix (landrace
#       richness x cost scenario), manuscript-ready rounding - step 9.
#
#   headline_estimates.csv          Overall lower/central/upper
#       national cost bounds - step 10.
#
# ===============================================================


# ----------------------------
# 0. Setup
# ----------------------------
setwd("~/Papers/PACS/Scaling_up_costs")

library(dplyr)
library(readr)
library(tidyr)
library(tibble)

etb_to_usd <- 0.0183


# ----------------------------
# 1. Load landrace change data
# ----------------------------
T1_T2_change_dat <- read_csv("T1_T2_rarity_assessment_data_share.csv")

T1_T2_simple <- T1_T2_change_dat %>%
  select(
    Community,
    T1_no_households,
    T2_no_households
  )


# ----------------------------
# 2. Calculate household gains by community
# ----------------------------
community_gain <- T1_T2_simple %>%
  mutate(
    change_households = T2_no_households - T1_no_households
  ) %>%
  group_by(Community) %>%
  summarise(
    hh_gain = sum(change_households, na.rm = TRUE),
    .groups = "drop"
  )


# ----------------------------
# 3. Community-level total paid, derived from actual selected offers
# ----------------------------
# select_offers.csv includes the contracted price (Birr) directly for
# each selected landrace - total paid per community is the sum of
# Birr across every landrace selected in that community. 
select_offers <- read_csv("select_offers.csv")

# Consistency check: any selected offer with a missing Birr value
# would otherwise silently understate that community's total.
missing_birr <- select_offers %>% filter(is.na(Birr))

if (nrow(missing_birr) > 0) {
  warning(
    nrow(missing_birr), " row(s) in select_offers.csv have a missing ",
    "Birr value - see missing_birr. These will be excluded from ",
    "community_total_paid."
  )
  print(missing_birr)
} else {
  cat("No missing Birr values in select_offers.csv.\n")
}

community_total_paid <- select_offers %>%
  group_by(Community) %>%
  summarise(
    total_paid_ETB = sum(Birr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(total_paid_USD = total_paid_ETB * etb_to_usd)

print(community_total_paid)

write_csv(community_total_paid, "community_total_paid.csv")


# ----------------------------
# 4. Calculate community-level cost per additional household
# ----------------------------
community_cost_per_household <- community_gain %>%
  left_join(community_total_paid, by = "Community") %>%
  filter(!is.na(total_paid_USD)) %>%
  filter(hh_gain > 0) %>%
  mutate(
    cost_per_household = total_paid_USD / hh_gain
  )

print(community_cost_per_household)

write_csv(community_cost_per_household, "community_cost_per_household.csv")


# ----------------------------
# 5. Cost uncertainty: quartiles of observed community costs
# ----------------------------
cost_scenarios <- quantile(
  community_cost_per_household$cost_per_household,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
) %>%
  as.numeric() %>%
  tibble(
    cost_scenario = c("Lower cost (Q1)", "Median cost", "Upper cost (Q3)"),
    cost_per_household_USD = .
  )

print(cost_scenarios)

write_csv(cost_scenarios, "cost_scenarios.csv")


# ----------------------------
# 6. National landrace uncertainty
# ----------------------------
# Borrell et al. (2019):
# 475 phonetic groups = conservative estimate
# 1270 vernacular names = upper estimate
landrace_scenarios <- tribble(
  ~landrace_scenario, ~total_landraces_ethiopia,
  "Conservative estimate (phonetic groups)", 475,
  "Upper estimate (vernacular names)",       1270
)


# ----------------------------
# 7. Scaling assumptions
# ----------------------------
# prop_safe: proportion of assessed landraces classified as "secure"
# or "common" at baseline - the remaining (1 - prop_safe) are assumed
# to require further conservation effort.
# safe_threshold: target number of cultivating households per
# landrace, a precautionary benchmark from conservation genetic theory
# (Drucker & Ramirez 2020; Franklin 1980; Soule 1980).
prop_safe <- 6/141
safe_threshold <- 50


# ----------------------------
# 8. Scenario matrix: landrace uncertainty x cost uncertainty
# ----------------------------
scaling_estimates <- crossing(
  landrace_scenarios,
  cost_scenarios
) %>%
  mutate(
    unsafe_landraces = total_landraces_ethiopia * (1 - prop_safe),
    households_needed = unsafe_landraces * safe_threshold,
    total_cost_USD = households_needed * cost_per_household_USD
  )

print(scaling_estimates)


# ----------------------------
# 9. Rounded manuscript-ready output
# ----------------------------
scaling_estimates_rounded <- scaling_estimates %>%
  mutate(
    unsafe_landraces = round(unsafe_landraces),
    households_needed = round(households_needed),
    cost_per_household_USD = round(cost_per_household_USD, 2),
    total_cost_USD = round(total_cost_USD, 0)
  ) %>%
  select(
    landrace_scenario,
    total_landraces_ethiopia,
    cost_scenario,
    cost_per_household_USD,
    unsafe_landraces,
    households_needed,
    total_cost_USD
  )

print(scaling_estimates_rounded)

write_csv(scaling_estimates_rounded, "scaling_estimates_rounded.csv")


# ----------------------------
# 10. Overall lower / central / upper headline values
# ----------------------------
headline_estimates <- scaling_estimates %>%
  summarise(
    lower_bound_USD = min(total_cost_USD, na.rm = TRUE),
    median_like_USD = total_cost_USD[
      landrace_scenario == "Upper estimate (vernacular names)" &
        cost_scenario == "Median cost"
    ],
    upper_bound_USD = max(total_cost_USD, na.rm = TRUE)
  ) %>%
  mutate(
    across(everything(), ~ round(.x, 0))
  )

print(headline_estimates)

write_csv(headline_estimates, "headline_estimates.csv")

