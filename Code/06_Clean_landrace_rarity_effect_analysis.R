# ============================================================
# PACS ENSET LANDRACE PES: DATA CLEANING +
# LANDRACE-LEVEL EFFECTS OF PACS CONTRACTS (DiD ANALYSIS)
# ============================================================
#
# FILES NEEDED
# ----------------------
#   - T1_T2_rarity_assessment_data_share.csv   
#   - select_offers.csv                        (produced in script 05)
#   - landrace_avg_status_landscape.csv        (produced in script 01)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Flags each community-landrace pair as contracted ("paid") if it
# appears in the selected tender offers, then reshapes the panel from
# wide to long format - one row per Community x Landrace x
# Timepoint (before/after).
#
# Restricts that panel to the primary sample - landraces classified as
# Restricted, Endangered or Critical at baseline (so treated and
# control landraces are drawn from the same conservation-risk pool) -
# and estimates difference-in-differences models (community + landrace
# fixed effects, standard errors clustered by community) for the
# primary outcome (rarity status) and two mechanism outcomes
# (households cultivating, plants per household).
#
# Exports regression tables and produces baseline-referenced DiD
# prediction plots with full confidence intervals (Fig. 4A: rarity
# status; Fig. 4B: households + plants).
#
# OUTPUTS GENERATED
# ----------------------
#   landrace_effect_analysis_df_long.csv   long-format community-landrace
#       panel (all landraces, all conservation-risk categories), with
#       clean column names and types, and landscape_category attached.
#       Input file for the crowding-in/crowding-out spillover analysis
#       script.
#
#   descriptives_rare.csv       descriptive stats for the primary
#       (rare-landrace) sample, by timepoint x contract status.
#
#   models_pacs_main.html       regression table (rarity, households,
#       plants models).
#
#   model_coefficients_all.csv  tidied coefficients for all three
#       models.
#
#   Figure2_line_plots.png      baseline-referenced DiD prediction plot
#       (Fig. 4A rarity; Fig. 4B households/plants).
#
# ============================================================


# ----------------------------
# 0. Set working directory
# ----------------------------
setwd("~/Papers/PACS/Landrace_level_effectiveness_analysis/Clean_code/Cleaning_data")


# ----------------------------
# 1. Load packages
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(fixest)
  library(modelsummary)
  library(broom)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})


# ----------------------------
# 2. Load inputs
# ----------------------------
before_after_data <- read.csv("T1_T2_rarity_assessment_data_share.csv")
selected_bids     <- read.csv("select_offers.csv")


# ----------------------------
# 3. Clean T2_participants (drop rows never followed up at T2)
# ----------------------------
# Removes rows with the literal text "NA" in T2_participants (as well
# as true NA), i.e. communities that were not revisited at T2.
before_after_data <- before_after_data %>%
  mutate(
    T2_participants = as.character(T2_participants),
    T2_participants = str_trim(T2_participants)
  ) %>%
  filter(!is.na(T2_participants), T2_participants != "NA")

print(before_after_data %>% distinct(Community) %>% nrow())


# ----------------------------
# 4. Prepare selected_bids matching keys
# ----------------------------
selected_bids_keys <- selected_bids %>%
  mutate(
    Community = str_trim(Community),
    Landrace = str_trim(Landrace_ID)
  ) %>%
  distinct(Community, Landrace_ID)

selected_key_vec <- selected_bids_keys %>%
  mutate(key = str_c(Community, "||", Landrace_ID)) %>%
  pull(key)


# ----------------------------
# 5. Flag contracted ("paid") community-landrace pairs
# ----------------------------
before_after_data <- before_after_data %>%
  mutate(
    Community = as.character(Community) %>% str_trim() %>% str_squish(),
    Landrace_ID  = as.character(Landrace_ID)  %>% str_trim() %>% str_squish(),
    key = str_c(Community, "||", Landrace_ID),
    paid = if_else(key %in% selected_key_vec, "Yes", "No")
  ) %>%
  select(-key)

print(colnames(before_after_data))


# ----------------------------
# 6. Reshape to long: Before (T1) vs After (T2)
#    One row per Community x Landrace x Timepoint
# ----------------------------
before_after_data$T2_participants <- as.integer(before_after_data$T2_participants)

before_after_long <- before_after_data %>%
  pivot_longer(
    cols = matches("^T[12]_"),
    names_to = c("Timepoint", ".value"),
    names_pattern = "^(T[12])_(.*)$"
  ) %>%
  mutate(
    Timepoint = dplyr::recode(Timepoint, "T1" = "before", "T2" = "after"),
    Timepoint = factor(Timepoint, levels = c("before", "after"))
  )


# ----------------------------
# 7. Compute derived variables (long dataframe)
# ----------------------------
before_after_long <- before_after_long %>%
  mutate(
    present = if_else(no_households > 0, 1L, 0L, missing = 0L)
  )


# ----------------------------
# 8. Build long-format panel (raw column names)
# ----------------------------
landrace_panel_raw <- before_after_long %>%
  mutate(
    paid = if_else(paid == "Yes", 1L, 0L),
    Timepoint = if_else(Timepoint == "after", 1L, 0L)
  ) %>%
  select(
    Community,
    Landrace_ID,
    paid,
    Timepoint,
    present,
    no_households,
    Landrace_avr_ind,
    rarity_status
  )


# ----------------------------
# 9. Rename and type columns for analysis
# ----------------------------
landrace_panel <- landrace_panel_raw %>%
  mutate(
    Community  = as.factor(Community),
    Landrace_ID   = as.factor(Landrace_ID),
    paid            = as.integer(paid),
    Timepoint       = as.integer(Timepoint),
    present         = as.numeric(present),
    no_households   = as.numeric(no_households),
    avr_plants      = as.numeric(Landrace_avr_ind),
    rarity_status   = as.numeric(rarity_status)
  )

# ----------------------------
# 9A. extract raw numbers for fig 2A and reporting raw averages
# ----------------------------

# Median and interquartile range, formatted as "median (Q1–Q3)"
med_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.0f (%.0f–%.0f)", q[2], q[1], q[3])
}

# Typical households and plants per household for each rarity score
by_score <- split(landrace_panel, landrace_panel$rarity_status)
score_table <- data.frame(
  rarity_status   = names(by_score),
  n_obs           = sapply(by_score, nrow),
  households      = sapply(by_score, function(d) med_iqr(d$no_households)),
  plants_per_hh   = sapply(by_score, function(d) med_iqr(d$avr_plants))
)
print(score_table, row.names = FALSE)
write.csv(score_table, "raw_score_median_table.csv")

# Mean rarity score by contract status (paid) and timepoint (0 = baseline, 1 = endline)
rarity_means <- aggregate(rarity_status ~ paid + Timepoint, data = landrace_panel,
                          FUN = function(x) round(mean(x), 2), na.action = na.omit)
print(rarity_means)
write.csv(rarity_means, "raw_rarity_means_t1t2.csv")


# ----------------------------
# 10. Join landscape rarity classification
# ----------------------------
rarity_lookup <- read_csv("landrace_avg_status_landscape.csv", show_col_types = FALSE) %>%
  select(Landrace_ID, landscape_category) %>%   # keep only needed column
  distinct() %>%                                  # ensure one row per landrace
  mutate(
    Landrace_ID = as.factor(Landrace_ID),
    landscape_category = as.factor(landscape_category)
  )

landrace_panel <- landrace_panel %>%
  left_join(rarity_lookup, by = "Landrace_ID")


# ----------------------------
# 11. Save long-format panel (clean column names, all landraces)
# ----------------------------
write_csv(landrace_panel, "landrace_effect_analysis_df_long.csv")


# ----------------------------
# 12. Restrict to primary sample: rare landraces at baseline
# ----------------------------
# Restricted, Endangered or Critical at baseline - treated and control
# landraces are drawn from the same conservation-risk pool. Epo and
# Waango are kept in regardless of their baseline category, as these are Lost 
# and Common respectively and were additional options given to farmers to select. 
eligible_categories <- c("Critical", "Endangered", "Restricted")
selected_outside_eligible <- c("LR0111", "LR0101")


landrace_panel_rare <- landrace_panel %>%
  filter(
    landscape_category %in% eligible_categories |
      Landrace_ID %in% selected_outside_eligible
  ) %>%
  droplevels()

# Sample size check
sample_check <- landrace_panel_rare %>%
  summarise(
    n_obs = n(),
    n_landraces = n_distinct(Landrace_ID),
    n_communities = n_distinct(Community)
  )
print(sample_check)



# ----------------------------
# 13. Descriptives (rare sample)
# ----------------------------
desc_rare <- landrace_panel_rare %>%
  group_by(Timepoint, paid) %>%
  summarise(
    n = n(),
    mean_hh = mean(no_households, na.rm = TRUE),
    mean_plants = mean(avr_plants, na.rm = TRUE),
    share_present = mean(present, na.rm = TRUE),
    mean_rarity = mean(rarity_status, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(desc_rare, "descriptives_rare.csv")


# ----------------------------
# 14. Primary analysis: rarity status
# ----------------------------
m_rarity_primary <- feols(
  rarity_status ~ Timepoint * paid | Community + Landrace_ID,
  data = landrace_panel_rare,
  cluster = ~Community
)


# ----------------------------
# 15. Secondary outcomes (mechanisms): households, plants
# ----------------------------
m_hh_primary <- feols(
  no_households ~ Timepoint * paid | Community + Landrace_ID,
  data = landrace_panel_rare,
  cluster = ~Community
)

m_total_primary <- feols(
  avr_plants ~ Timepoint * paid | Community + Landrace_ID,
  data = landrace_panel_rare,
  cluster = ~Community
)


# ----------------------------
# 16. Export regression tables
# ----------------------------
models_main <- list(
  "Primary: Rarity status (rare)" = m_rarity_primary,
  "Mechanism: Households"         = m_hh_primary,
  "Mechanism: Plants"             = m_total_primary
)

modelsummary(
  models_main,
  output = "models_pacs_main.html",
  stars = TRUE,
  statistic = "({std.error})",
  gof_omit = "IC|Log|Adj|Within|Std.Errors|FE"
)

coefs_all <- map_dfr(
  models_main,
  broom::tidy,
  .id = "model"
)

write_csv(coefs_all, "model_coefficients_all.csv")


# ----------------------------
# 17. DiD prediction plots (Fig. 4A rarity; Fig. 4B households/plants)
# ----------------------------
# Baseline-referenced DiD prediction plot: absolute predicted values
# start from the observed baseline (T1, non-contracted) mean, with the
# model's coefficients added on top, and standard errors 
# propagated from the model's clustered variance-covariance matrix.
#
# The baseline SE is computed by collapsing to one mean per community
# first, then taking the SE across those community-level means - this
# matches cluster = ~Community_name used in the feols() models, so all
# four points on each plot use a consistent clustering assumption

baseline_referenced_did_plot <- function(data, outcome, model, y_lab, y_limits = NULL,
                                         point_size = 3, dodge = 0) {
  
  # Community-clustered baseline mean + SE
  baseline_by_comm <- data %>%
    filter(Timepoint == 0, paid == 0) %>%
    group_by(Community) %>%
    summarise(comm_mean = mean(.data[[outcome]], na.rm = TRUE), .groups = "drop")
  
  mu0 <- mean(baseline_by_comm$comm_mean, na.rm = TRUE)
  se0 <- sd(baseline_by_comm$comm_mean, na.rm = TRUE) /
    sqrt(sum(!is.na(baseline_by_comm$comm_mean)))
  
  b <- coef(model)
  V <- vcov(model)  # already clustered by community (cluster = ~Community_name)
  
  b1 <- b["Timepoint"]
  b2 <- b["paid"]
  b3 <- b["Timepoint:paid"]
  
  pred_vals <- c(
    T1_control = mu0,
    T2_control = mu0 + b1,
    T1_treated = mu0 + b2,
    T2_treated = mu0 + b1 + b2 + b3
  )
  
  se_vals <- c(
    T1_control = se0,
    T2_control = sqrt(se0^2 + V["Timepoint", "Timepoint"]),
    T1_treated = sqrt(se0^2 + V["paid", "paid"]),
    T2_treated = sqrt(
      se0^2 +
        V["Timepoint", "Timepoint"] +
        V["paid", "paid"] +
        V["Timepoint:paid", "Timepoint:paid"] +
        2 * V["Timepoint", "paid"] +
        2 * V["Timepoint", "Timepoint:paid"] +
        2 * V["paid", "Timepoint:paid"]
    )
  )
  
  plot_df <- data.frame(
    Time = rep(c("T1", "T2"), 2),
    Group = rep(c("Non-contracted", "Contracted"), each = 2),
    estimate = pred_vals,
    se = se_vals
  )
  plot_df$lower <- plot_df$estimate - 1.96 * plot_df$se
  plot_df$upper <- plot_df$estimate + 1.96 * plot_df$se
  
  p <- ggplot(plot_df, aes(x = Time, y = estimate, color = Group, shape = Group, group = Group)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = point_size, position = position_dodge(width = dodge)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, position = position_dodge(width = dodge)) +
    scale_color_manual(values = c("Non-contracted" = "#3F8ED5", "Contracted" = "#194975")) +
    scale_shape_manual(values = c("Non-contracted" = 16, "Contracted" = 17)) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    labs(y = y_lab, x = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
  
  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits, expand = c(0, 0))
  
  p
}

# ---- Fig. 4A: rarity status ----
rarity_line <- baseline_referenced_did_plot(
  landrace_panel_rare, "rarity_status", m_rarity_primary,
  y_lab = "Within-community rarity score",
  y_limits = c(1, 6)
)

# ---- Fig. 4B (i): households cultivating ----
household_line <- baseline_referenced_did_plot(
  landrace_panel_rare, "no_households", m_hh_primary,
  y_lab = "Households cultivating landrace",
  dodge = 0.08,
  y_limits = c(0,25)
)

# ---- Fig. 4B (ii): plants per household ----
plant_line <- baseline_referenced_did_plot(
  landrace_panel_rare, "avr_plants", m_total_primary,
  y_lab = "Plants per household",
  y_limits = c(0,31)
)

# ---- Combine: Fig. 4A (top) + Fig. 4B panels (bottom row) ----
layout <- "
.A
BC
"

combined_plot <- rarity_line + plot_spacer() +
  household_line + plant_line +
  plot_layout(design = layout)

combined_plot

ggsave(
  filename = "Figure2_line_plots.png",
  plot = combined_plot,
  width = 6,
  height = 6,
  dpi = 600
)
