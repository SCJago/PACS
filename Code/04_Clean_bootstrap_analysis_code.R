# ============================================================
# POOLED BOOTSTRAP PSEUDO-CONTROL ANALYSIS
# ============================================================
#
# FILES NEEDED
# ----------------------
#   - hh_boot_35_37E.rds          (control household data - produced
#                                  by script 03 but due to confidential
#                                  coordinates provided in Data folder)
#   - T1_T2_rarity_assessment_data_share.csv  (treated community data)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# For each bootstrap replicate:
#   1. For each treated community, sample a pseudo-control community
#      of the same size from the control household pool.
#   2. Calculate change in outcomes for the pseudo-control community.
#   3. Calculate treated change minus pseudo-control change - this
#      community's effect for this replicate.
#   4. Average these effects across treated communities to get the
#      pooled effect for this replicate.
#
# This produces a bootstrap distribution of the mean treatment effect
# (pooled across communities), and a bootstrap distribution of each
# individual community's own effect.
#
# OUTPUTS GENERATED
# ----------------------
#   OBSERVED TREATED CHANGES: raw mean change in treated communities
#       (T2 - T1), before any comparison to controls, with standard
#       error.
#
#   CONTROL BACKGROUND: how much the pseudo-control communities
#       themselves changed on average in the absence of the PACS
#       intervention.
#
#   POOLED TREATMENT EFFECT: the mean effect (treated change minus
#       pseudo-control change) averaged across all treated
#       communities, with a one-sided empirical p-value.
#
#   COMMUNITY-LEVEL TREATMENT EFFECTS (Fig. 3B): mean + 95% CI of the
#       effect for EACH individual treated community.
#
#   RESULTS SUMMARY: pulls the key numbers from the above
#       into one table, in the same order as the results paragraph.
#
# ============================================================


# ----------------------------
# 0. Set working directory
# ----------------------------
setwd("~/Papers/PACS/Bootstrapping_controls/clean_code")


# ----------------------------
# 1. Load packages
# ----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(patchwork)


# ----------------------------
# 2. Load data
# ----------------------------
hh_boot <- readRDS("hh_boot_35_37E.rds")
treated <- read_csv("T1_T2_rarity_assessment_data_share.csv", show_col_types = FALSE)


# ----------------------------
# 3. Summarise treated communities
# ----------------------------
treated_summary <- treated %>%
  filter(T1_participants > 0, T2_participants > 0) %>%
  group_by(Community) %>%
  summarise(
    T1_participants = first(T1_participants),
    T2_participants = first(T2_participants),
    
    avg_T1 = sum(T1_no_households, na.rm = TRUE) / T1_participants,
    avg_T2 = sum(T2_no_households, na.rm = TRUE) / T2_participants,
    change_avg = avg_T2 - avg_T1,
    
    unique_T1 = sum(T1_no_households > 0, na.rm = TRUE),
    unique_T2 = sum(T2_no_households > 0, na.rm = TRUE),
    change_unique = unique_T2 - unique_T1,
    
    N_households = T1_participants,
    .groups = "drop"
  )

print(treated_summary)


# ----------------------------
# 3b. Observed treated changes
# ----------------------------
observed_treated_summary <- tibble(
  outcome = c(
    "Average landraces per household",
    "Unique landraces per community"
  ),
  mean_observed_change = c(
    mean(treated_summary$change_avg, na.rm = TRUE),
    mean(treated_summary$change_unique, na.rm = TRUE)
  ),
  se_observed_change = c(
    sd(treated_summary$change_avg, na.rm = TRUE) / sqrt(sum(!is.na(treated_summary$change_avg))),
    sd(treated_summary$change_unique, na.rm = TRUE) / sqrt(sum(!is.na(treated_summary$change_unique)))
  ),
  n_communities_positive = c(
    sum(treated_summary$change_avg > 0, na.rm = TRUE),
    sum(treated_summary$change_unique > 0, na.rm = TRUE)
  ),
  n_communities_total = nrow(treated_summary)
)

print(observed_treated_summary)

# ----------------------------
# 2c. Average landraces per household (descriptive)
# ----------------------------

# Control households (one row per household)
control_hh_summary <- hh_boot %>%
  summarise(
    n_households = n(),
    mean_T1 = mean(n_T1, na.rm = TRUE), sd_T1 = sd(n_T1, na.rm = TRUE),
    median_T1 = median(n_T1, na.rm = TRUE),
    mean_T2 = mean(n_T2, na.rm = TRUE), sd_T2 = sd(n_T2, na.rm = TRUE),
    median_T2 = median(n_T2, na.rm = TRUE)
  )
print(control_hh_summary)

# Treated communities (community-level averages, then summarised across communities)
treated_hh_summary <- treated_summary %>%
  summarise(
    n_communities = n(),
    mean_T1 = mean(avg_T1), sd_T1 = sd(avg_T1),
    mean_T2 = mean(avg_T2), sd_T2 = sd(avg_T2),
    # household-weighted version: total household-landrace pairs / total households
    weighted_T1 = sum(avg_T1 * T1_participants) / sum(T1_participants),
    weighted_T2 = sum(avg_T2 * T2_participants) / sum(T2_participants)
  )
print(treated_hh_summary)

# ----------------------------
# 4. Bootstrap settings
# ----------------------------
set.seed(123)
B <- 5000


# ----------------------------
# 5. Function to calculate pseudo-community outcomes
# ----------------------------
compute_stats <- function(hh_sample) {
  
  avg_T1 <- mean(hh_sample$n_T1, na.rm = TRUE)
  avg_T2 <- mean(hh_sample$n_T2, na.rm = TRUE)
  
  unique_T1 <- length(unique(unlist(hh_sample$T1_clean)))
  unique_T2 <- length(unique(unlist(hh_sample$T2_clean)))
  
  c(
    change_avg = avg_T2 - avg_T1,
    change_unique = unique_T2 - unique_T1
  )
}


# ----------------------------
# 6. One pooled bootstrap replicate
# ----------------------------
one_pooled_bootstrap <- function(treated_summary, hh_boot) {
  
  effects <- treated_summary %>%
    rowwise() %>%
    mutate(
      boot_stats = list({
        samp_idx <- sample.int(
          n = nrow(hh_boot),
          size = N_households,
          replace = TRUE
        )
        compute_stats(hh_boot[samp_idx, ])
      }),
      
      boot_change_avg = boot_stats["change_avg"],
      boot_change_unique = boot_stats["change_unique"],
      
      effect_avg = change_avg - boot_change_avg,
      effect_unique = change_unique - boot_change_unique
    ) %>%
    ungroup()
  
  # Per-community effects, named by community
  community_effect_avg <- setNames(
    effects$effect_avg,
    paste0("effect_avg_", effects$Community)
  )
  community_effect_unique <- setNames(
    effects$effect_unique,
    paste0("effect_unique_", effects$Community)
  )
  
  c(
    mean_effect_avg = mean(effects$effect_avg, na.rm = TRUE),
    mean_effect_unique = mean(effects$effect_unique, na.rm = TRUE),
    
    mean_control_avg = mean(effects$boot_change_avg, na.rm = TRUE),
    mean_control_unique = mean(effects$boot_change_unique, na.rm = TRUE),
    
    community_effect_avg,
    community_effect_unique
  )
}

# ----------------------------
# 7. Run pooled bootstrap
# ----------------------------
pooled_boot <- replicate(
  B,
  one_pooled_bootstrap(treated_summary, hh_boot)
)

boot_mean_effect_avg <- pooled_boot["mean_effect_avg", ]
boot_mean_effect_unique <- pooled_boot["mean_effect_unique", ]

boot_mean_control_avg <- pooled_boot["mean_control_avg", ]
boot_mean_control_unique <- pooled_boot["mean_control_unique", ]


# ----------------------------
# 8. Summarise pooled bootstrap results
# ----------------------------
control_background <- tibble(
  outcome = c(
    "Average landraces per household",
    "Unique landraces per community"
  ),
  
  mean_control = c(
    mean(boot_mean_control_avg),
    mean(boot_mean_control_unique)
  ),
  
  lo_95 = c(
    quantile(boot_mean_control_avg, 0.025),
    quantile(boot_mean_control_unique, 0.025)
  ),
  
  hi_95 = c(
    quantile(boot_mean_control_avg, 0.975),
    quantile(boot_mean_control_unique, 0.975)
  )
)

print(control_background)


pooled_results <- tibble(
  outcome = c(
    "Average landraces per household",
    "Unique landraces per community"
  ),
  
  mean_effect = c(
    mean(boot_mean_effect_avg, na.rm = TRUE),
    mean(boot_mean_effect_unique, na.rm = TRUE)
  ),
  
  lo_95 = c(
    quantile(boot_mean_effect_avg, 0.025, na.rm = TRUE),
    quantile(boot_mean_effect_unique, 0.025, na.rm = TRUE)
  ),
  
  hi_95 = c(
    quantile(boot_mean_effect_avg, 0.975, na.rm = TRUE),
    quantile(boot_mean_effect_unique, 0.975, na.rm = TRUE)
  ),
  
  p_one_sided = c(
    (1 + sum(boot_mean_effect_avg <= 0, na.rm = TRUE)) / (B + 1),
    (1 + sum(boot_mean_effect_unique <= 0, na.rm = TRUE)) / (B + 1)
  )
)

print(pooled_results)

write_csv(
  pooled_results,
  "pooled_bootstrap_treatment_effect_results.csv"
)


# ----------------------------
# 9. Save bootstrap draws
# ----------------------------
pooled_boot_draws <- tibble(
  bootstrap_id = 1:B,
  mean_effect_avg = boot_mean_effect_avg,
  mean_effect_unique = boot_mean_effect_unique
)

write_csv(
  pooled_boot_draws,
  "pooled_bootstrap_draws.csv"
)


# ----------------------------
# 10. Community-level effects (Fig. 3B)
# ----------------------------
community_names <- treated_summary$Community

community_effects <- tibble(
  Community = rep(community_names, times = 2),
  metric = rep(
    c("Average landraces per household", "Unique landraces per community"),
    each = length(community_names)
  ),
  effect = c(
    sapply(community_names, function(cn) mean(pooled_boot[paste0("effect_avg_", cn), ], na.rm = TRUE)),
    sapply(community_names, function(cn) mean(pooled_boot[paste0("effect_unique_", cn), ], na.rm = TRUE))
  ),
  effect_lo = c(
    sapply(community_names, function(cn) quantile(pooled_boot[paste0("effect_avg_", cn), ], 0.025, na.rm = TRUE)),
    sapply(community_names, function(cn) quantile(pooled_boot[paste0("effect_unique_", cn), ], 0.025, na.rm = TRUE))
  ),
  effect_hi = c(
    sapply(community_names, function(cn) quantile(pooled_boot[paste0("effect_avg_", cn), ], 0.975, na.rm = TRUE)),
    sapply(community_names, function(cn) quantile(pooled_boot[paste0("effect_unique_", cn), ], 0.975, na.rm = TRUE))
  )
)

print(community_effects)

write_csv(
  community_effects,
  "community_level_bootstrap_effects.csv"
)


# ----------------------------
# 11. Figure: individual community-level effects
# ----------------------------
plot_df <- community_effects %>%
  mutate(
    Community = as.character(Community),
    metric = as.character(metric)
  )

community_levels <- sort(unique(plot_df$Community))

plot_df <- plot_df %>%
  mutate(
    Community = factor(Community, levels = rev(community_levels))
  )

plot_df_avg <- plot_df %>%
  filter(metric == "Average landraces per household")

plot_df_unique <- plot_df %>%
  filter(metric == "Unique landraces per community")

make_effect_plot <- function(df, xlab, x_limits = NULL, show_y_labels = TRUE) {
  
  p <- ggplot(df, aes(y = Community)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
    geom_errorbarh(
      aes(xmin = effect_lo, xmax = effect_hi),
      height = 0.2,
      linewidth = 0.35,
      colour = "#A6CAEC"
    ) +
    geom_point(aes(x = effect), size = 1.7, colour = "#215F9A") +
    theme_classic(base_size = 8) +
    labs(x = xlab, y = NULL) +
    theme(
      axis.text.y = if (show_y_labels) element_text(size = 8) else element_blank(),
      axis.ticks.y = if (show_y_labels) element_line() else element_blank(),
      axis.title.x = element_text(size = 8),
      plot.margin = margin(5.5, 20, 5.5, 5.5)
    )
  
  if (!is.null(x_limits)) {
    p <- p + scale_x_continuous(limits = x_limits)
  }
  
  p
}

p_avg <- make_effect_plot(
  plot_df_avg,
  xlab = "Community-level treatment effect for\n average landraces per household",
  x_limits = c(-1, 8),
  show_y_labels = TRUE
)

p_unique <- make_effect_plot(
  plot_df_unique,
  xlab = "Community-level treatment effect for\n total unique landraces",
  x_limits = c(-2, 15),
  show_y_labels = FALSE
)

final_plot <- p_avg + p_unique + plot_layout(ncol = 2)

final_plot

ggsave(
  "community_level_effects_individual_only.png",
  plot = final_plot,
  width = 7,
  height = 3.5,
  dpi = 300
)


# ----------------------------
# 12. Result summary table
# ----------------------------
result_summary <- tibble(
  outcome = c(
    "Average landraces per household",
    "Unique landraces per community"
  ),
  observed_treated_mean_change = observed_treated_summary$mean_observed_change,
  observed_treated_se          = observed_treated_summary$se_observed_change,
  control_mean_change   = control_background$mean_control,
  control_lo_95         = control_background$lo_95,
  control_hi_95          = control_background$hi_95,
  pooled_mean_effect     = pooled_results$mean_effect,
  pooled_effect_lo_95    = pooled_results$lo_95,
  pooled_effect_hi_95    = pooled_results$hi_95,
  p_one_sided            = pooled_results$p_one_sided
)

print(result_summary)

write_csv(
  result_summary,
  "manuscript_results_summary.csv"
)

