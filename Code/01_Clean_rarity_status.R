# ============================================================
# LANDRACE RARITY SCORING (T1 BASELINE + T1 vs T2 CHANGE)
# ============================================================
#
# FILES NEEDED
# ----------------------
#   - T1_T2_rarity_assessment_data.csv   (community x landrace survey data,
#                                          T1 baseline + T2 follow-up)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Community-level rarity scores (assigned per landrace per
# community, based on proportion of households cultivating it and its
# abundance within households - Fig 2A) are aggregated across the study
# landscape into an overall conservation status per landrace, from
# Secure to Critical to Lost (Fig 2B). This script then:
# - builds the Fig 2C rarity distribution + redundancy plots
# - builds potentially useful supplementary plots showing spread of 
#   landscape status
# - pulls headline stats for the manuscript results text
# - compares T1 vs T2 rarity composition per community (Fig 3A)
#
# OUTPUTS GENERATED
# ----------------------
#   T1_fullgrid_with_absences.csv 
#
#   landrace_avg_status_landscape.csv    one row per landrace: mean
#       landscape rarity score and landscape_category.
#
#   FIG 2C MAIN FIGURE PLOTS: landrace_rarity_bar.png (2Ci) and
#       p_comm_vs_plantsperhouse_scatt.png (2Cii), plus a combined
#       side-by-side panel (landscape_rarity_T1_panel.png).
#
#   SUPPLEMENT PLOTS: households/plants per household by rarity
#       category, as boxplots.
#
#   HEADLINE STATS CSVs: landscape_profile_T1.csv (n landraces per
#       category), threatened_share_T1.csv (share threatened/lost),
#       headline_numbers_T1.csv (key manuscript figures).
#
#   T1 vs T2 COMMUNITY CHANGE (Fig. 3A): T1_T2_stacked_bar_community_specific.png,
#       a per-community faceted comparison of rarity composition.
#
# ============================================================


# ----------------------------
# 0. Set working directory
# ----------------------------
setwd("~/Papers/PACS/01_Rarity_status/Clean_code")


# ----------------------------
# 1. Load packages
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(ggh4x)   # facet_wrap2() used in the T1 vs T2 community panel
})


# ----------------------------
# 2. Constants: labels + colours (keep consistent everywhere)
# ----------------------------

# Rarity code -> readable label
xlabel <- c(
  "1" = "Secure",
  "2" = "Common",
  "3" = "Restricted",
  "4" = "Endangered",
  "5" = "Critical",
  "6" = "Lost"
)

# Rarity code -> plot colour
cols <- c(
  "1" = "#2C7BB6",
  "2" = "#7FCDBB",
  "3" = "#F6E8A6",
  "4" = "#F4A259",
  "5" = "#D73027",
  "6" = "#A6A6A6"
)

# Same colours, keyed by category name instead of code (for fill/colour aes)
cols_by_cat <- c(
  "Secure"     = cols[["1"]],
  "Common"     = cols[["2"]],
  "Restricted" = cols[["3"]],
  "Endangered" = cols[["4"]],
  "Critical"   = cols[["5"]],
  "Lost"       = cols[["6"]]
)

# The three "threatened" categories, per Fig 2B
rare_levels <- c("Restricted", "Endangered", "Critical")

# Category name -> rarity code (inverse of xlabel), used to convert
# landscape_category back into a numeric code for plotting/counting
category_to_code <- c(
  "Secure"     = 1L,
  "Common"     = 2L,
  "Restricted" = 3L,
  "Endangered" = 4L,
  "Critical"   = 5L,
  "Lost"       = 6L
)


# ----------------------------
# 3. Load data + build T1 landscape rarity grid
# ----------------------------
T1_T2_rarity_assessment_data <- read.csv("T1_T2_rarity_assessment_data.csv")

# Master lists of every community and every landrace in the study,
# used below to build the full Community x Landrace grid
all_communities <- T1_T2_rarity_assessment_data %>% distinct(Community_name)
all_landraces   <- T1_T2_rarity_assessment_data %>% distinct(Landrace_name)

t1_fullgrid <- T1_T2_rarity_assessment_data %>%
  
  # Keep only the T1 columns needed for the rarity assessment
  select(Community_name, Landrace_name, T1_rarity_status,
         T1_no_households, T1_Landrace_avr_ind) %>%
  
  # Expand to every Community x Landrace combination. Pairs missing
  # from the raw data are added as new rows, with NA in every other column
  complete(
    Community_name = all_communities$Community_name,
    Landrace_name  = all_landraces$Landrace_name
  ) %>%
  
  mutate(
    # No record of households growing it -> 0 households
    T1_no_households = replace_na(T1_no_households, 0),
    
    # No rarity status assignable -> code 6 ("Lost"), covering both
    # confirmed absences and newly-added (never surveyed) rows
    T1_rarity_status = replace_na(as.integer(T1_rarity_status), 6L)
    
    # T1_Landrace_avr_ind is left untouched:
    #  - confirmed absence (was in raw data)          -> real 0 value
    #  - newly-added row (never mentioned as "lost")  -> stays NA
    # This stops newly-added rows being counted as real zero-plant
    # observations, which would deflate the "when present" means below
  )

landrace_avg_status_landscape <- t1_fullgrid %>%
  group_by(Landrace_name) %>%
  summarise(
    
    # Mean rarity status across ALL communities (absences = 6),
    # i.e. the "populate absences throughout" step from Fig 2B
    avg_status_landscape = mean(T1_rarity_status, na.rm = TRUE),
    
    # Number / share of communities actually growing this landrace
    n_communities_present     = sum(T1_no_households > 0, na.rm = TRUE),
    share_communities_present = mean(T1_no_households > 0, na.rm = TRUE),
    
    # Totals across all communities where present
    total_n_hh_present   = sum(T1_no_households),
    total_plants_present = sum(T1_Landrace_avr_ind * T1_no_households, na.rm = TRUE),
    
    # Averages calculated ONLY where the landrace is present (households > 0
    # and plant count not NA), so never-surveyed rows don't drag these down
    mean_hh_when_present     = mean(T1_no_households[T1_no_households > 0 & !is.na(T1_Landrace_avr_ind)], na.rm = TRUE),
    mean_plants_when_present = mean(T1_Landrace_avr_ind[T1_no_households > 0 & !is.na(T1_Landrace_avr_ind)], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  # Rank landraces from rarest to most secure
  arrange(desc(avg_status_landscape), n_communities_present) %>%
  
  mutate(
    # Bin the landscape mean into rarity categories (Fig 2B thresholds)
    landscape_category = case_when(
      avg_status_landscape == 6    ~ "Lost",
      avg_status_landscape >= 5.9  ~ "Critical",
      avg_status_landscape >= 5.6  ~ "Endangered",
      avg_status_landscape >= 5.2  ~ "Restricted",
      avg_status_landscape >= 4.0  ~ "Common",
      avg_status_landscape >= 3.0  ~ "Secure",
      TRUE                         ~ "Outside expected range"  # sanity check flag
    ),
    
    # Fix category order for consistent plotting/reporting
    landscape_category = factor(
      landscape_category,
      levels = c("Secure", "Common", "Restricted", "Endangered",
                 "Critical", "Lost", "Outside expected range")
    )
  )

print(landrace_avg_status_landscape)

write.csv(t1_fullgrid, "T1_fullgrid_with_absences.csv", row.names = FALSE)
write.csv(landrace_avg_status_landscape, "landrace_avg_status_landscape.csv", row.names = FALSE)


# ----------------------------
# 4. Fig 2C main figure plots
# ----------------------------

# Landscape-level counts: how many landraces fall in each rarity
# category across the whole study area
landrace_summary_combined <- landrace_avg_status_landscape %>%
  filter(as.character(landscape_category) %in% names(category_to_code)) %>%
  mutate(T1_rarity_status = unname(category_to_code[as.character(landscape_category)])) %>%
  count(T1_rarity_status, name = "landrace_count") %>%
  complete(T1_rarity_status = 1:6, fill = list(landrace_count = 0)) %>%
  mutate(T1_rarity_status = factor(T1_rarity_status, levels = 1:6))

# Fig 2C(i): bar chart of landscape-level rarity distribution -
# how many landraces fall into each rarity category overall
p_combined_bar <- ggplot(
  landrace_summary_combined,
  aes(x = T1_rarity_status, y = landrace_count, fill = T1_rarity_status)
) +
  geom_col(width = 0.8, alpha = 0.9, colour = "black") +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = xlabel) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Number of landraces") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("landrace_rarity_bar.png",
       p_combined_bar, width = 4, height = 4, dpi = 300)

# Fig 2C(ii): redundancy plot - how widespread (x) vs how abundant
# when present (y) each landrace is, coloured by rarity category
p_communites_vs_plantsperhouse <- ggplot(
  landrace_avg_status_landscape %>% filter(!is.na(landscape_category)),
  aes(n_communities_present, mean_plants_when_present, color = landscape_category)
) +
  geom_point(size = 2, alpha = 0.9) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = cols_by_cat, guide = "none") +
  labs(
    x = "Communities cultivating landrace",
    y = "Mean plants per household\nwhen present (log10)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(margin = margin(t = 4)),
    plot.margin = margin(t = 5.5, r = 5.5, b = 22, l = 5.5)
  )

ggsave("p_comm_vs_plantsperhouse_scatt.png",
       p_communites_vs_plantsperhouse, width = 4.2, height = 4, dpi = 300)

# Combined panel: the two main Fig 2C plots side by side
landscape_rarity_plots_T1 <- ggarrange(
  p_combined_bar,
  p_communites_vs_plantsperhouse,
  ncol = 2,
  nrow = 1
)

ggsave("landscape_rarity_T1_panel.png",
       landscape_rarity_plots_T1, width = 8, height = 4, dpi = 300)


# ----------------------------
# 5. Supplement plots (households/plants by rarity category)
# ----------------------------

# Attach landscape category to each community x landrace row (T1)
t1_obs <- T1_T2_rarity_assessment_data %>%
  left_join(landrace_avg_status_landscape %>% select(Landrace_name, landscape_category),
            by = "Landrace_name") %>%
  mutate(
    present_T1      = !is.na(T1_no_households) & T1_no_households > 0,
    T1_total_plants = as.numeric(T1_no_households) * as.numeric(T1_Landrace_avr_ind),
    landscape_category = factor(
      landscape_category,
      levels = c("Secure", "Common", "Restricted", "Endangered", "Critical", "Lost")
    )
  )

# Restrict to rows where the landrace is actually present and has a
# valid landscape category assigned
t1_present <- t1_obs %>%
  filter(present_T1, !is.na(landscape_category))

# Summary stats per landscape category (sanity check before plotting)
t1_summary <- t1_present %>%
  group_by(landscape_category) %>%
  summarise(
    mean_households    = mean(T1_no_households, na.rm = TRUE),
    mean_plants_per_hh  = mean(T1_Landrace_avr_ind, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

print(t1_summary)

# Households cultivating each landrace, by landscape rarity category
p_hh_by_landscape_cat <- ggplot(t1_present, aes(x = landscape_category, y = T1_no_households, fill = landscape_category)) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_fill_manual(values = cols_by_cat, guide = "none") +
  labs(x = NULL, y = "Households cultivating landrace\n per community") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("box_households_by_landscape_rarity_T1.png",
       p_hh_by_landscape_cat, width = 4, height = 4, dpi = 300)

# Total plants (households x mean plants), by landscape rarity category
p_plants_by_landscape_cat <- ggplot(t1_present, aes(x = landscape_category, y = T1_total_plants, fill = landscape_category)) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_fill_manual(values = cols_by_cat, guide = "none") +
  scale_y_log10(labels = comma) +
  labs(x = NULL, y = "Total plants\n (households x mean plants; log10)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("box_total_plants_by_landscape_rarity_T1.png",
       p_plants_by_landscape_cat, width = 4, height = 4, dpi = 300)

# Mean plants per household, by landscape rarity category (only where recorded)
p_meanplants_by_landscape_cat <- ggplot(
  t1_present %>% filter(!is.na(T1_Landrace_avr_ind)),
  aes(x = landscape_category, y = T1_Landrace_avr_ind, fill = landscape_category)
) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_fill_manual(values = cols_by_cat, guide = "none") +
  scale_y_log10(labels = comma) +
  labs(x = NULL, y = "Mean plants per household\n (log10)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("box_mean_plants_per_household_by_landscape_rarity_T1.png",
       p_meanplants_by_landscape_cat, width = 4, height = 4, dpi = 300)


# ----------------------------
# 6. Headline stats CSVs (T1)
# ----------------------------

# Number and share of landraces in each landscape rarity category
landscape_profile <- landrace_avg_status_landscape %>%
  filter(!is.na(landscape_category), landscape_category != "Outside expected range") %>%
  mutate(
    landscape_category = factor(
      landscape_category,
      levels = c("Secure", "Common", "Restricted", "Endangered", "Critical", "Lost")
    )
  ) %>%
  count(landscape_category, name = "n_landraces") %>%
  mutate(share = n_landraces / sum(n_landraces))

# Overall share of landraces that are "threatened" (Restricted/Endangered/
# Critical) or fully "Lost" across the landscape
threatened_share <- landrace_avg_status_landscape %>%
  filter(!is.na(landscape_category), landscape_category != "Outside expected range") %>%
  summarise(
    n_total          = n(),
    n_threatened     = sum(landscape_category %in% rare_levels),
    share_threatened = n_threatened / n_total,
    n_lost           = sum(landscape_category == "Lost"),
    share_lost       = n_lost / n_total
  )

# Headline numbers reported in the manuscript text, e.g.:
# "141 unique enset landraces, of which 64 occurred in two or less
#  communities and 36 were known but had become locally extinct"
headline_numbers <- landrace_avg_status_landscape %>%
  filter(!is.na(landscape_category), landscape_category != "Outside expected range") %>%
  summarise(
    n_unique_landraces  = n(),
    n_low_redundancy    = sum(n_communities_present %in% c(1, 2)),   # present in <=2 communities
    n_locally_extinct   = sum(landscape_category == "Lost")
  )

print(landscape_profile)
print(threatened_share)
print(headline_numbers)

write.csv(landscape_profile, "landscape_profile_T1.csv", row.names = FALSE)
write.csv(threatened_share, "threatened_share_T1.csv", row.names = FALSE)
write.csv(headline_numbers, "headline_numbers_T1.csv", row.names = FALSE)


# ----------------------------
# 7. T1 vs T2 change plot by community (Fig. 3A; only communities
#    with T2 data)
# ----------------------------

# Reshape T1 + T2 rarity status into long format (one row per
# Community x Landrace x Timepoint), then re-derive a rounded
# rarity code per Timepoint for counting/plotting
community_change_plot_data <- T1_T2_rarity_assessment_data %>%
  
  # Only communities that were actually revisited at T2
  filter(!is.na(T2_participants)) %>%
  
  pivot_longer(
    cols = c(T1_rarity_status, T2_rarity_status),
    names_to = "Timepoint_raw",
    values_to = "rarity_status"
  ) %>%
  mutate(
    Timepoint = recode(Timepoint_raw, "T1_rarity_status" = "T1", "T2_rarity_status" = "T2"),
    Timepoint = factor(Timepoint, levels = c("T2", "T1"))
  ) %>%
  
  # Average rarity status per Community x Timepoint x Landrace, then
  # bin into whole rarity codes (1-6) for counting landraces per category
  group_by(Community_name, Timepoint, Landrace_name) %>%
  summarise(avg_status = mean(rarity_status, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    T_rarity_status = case_when(
      avg_status == 6 ~ 6L,
      avg_status >= 5 ~ 5L,
      avg_status >= 4 ~ 4L,
      avg_status >= 3 ~ 3L,
      avg_status >= 2 ~ 2L,
      avg_status >= 1 ~ 1L,
      TRUE ~ NA_integer_
    ),
    T_rarity_status = factor(T_rarity_status, levels = 1:6)
  ) %>%
  
  # Count landraces per Community x Timepoint x rarity code, filling
  # any missing combination with 0 so every facet/category is represented
  count(Community_name, Timepoint, T_rarity_status, name = "landrace_count") %>%
  complete(
    Community_name, Timepoint, T_rarity_status,
    fill = list(landrace_count = 0)
  )

# Stacked bar version (Fig. 3A), excluding "Lost" (code 6) so the bars
# focus on the graded rarity categories rather than being dominated by losses
community_change_plot_data_NoLost <- community_change_plot_data %>%
  filter(T_rarity_status != 6)

p_T1T2_comm_bar <- ggplot(
  community_change_plot_data_NoLost,
  aes(x = Timepoint, y = landrace_count, fill = T_rarity_status)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(limits = c("T1", "T2")) +
  scale_y_continuous(expand = c(0), limits = c(0, 27)) +
  facet_wrap2(~Community_name, ncol = 4, axes = "all") +
  theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

p_T1T2_comm_bar

ggsave("T1_T2_stacked_bar_community_specific.png",
       p_T1T2_comm_bar, width = 6.5, height = 6, dpi = 300)
