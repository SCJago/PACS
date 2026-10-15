# ===============================================================
# BID HETEROGENEITY ANALYSIS
# ===============================================================
#
# FILES NEEDED
# ----------------------
#   - clean_offers.csv                  
#   - select_offers.csv                   (produced in script 05)
#   - landrace_tender_id_crosswalk.csv    
#   - Use_data_final_share.csv
#   - Landrace_size_data_clean_share.csv
#   - enset_per_farmer_by_community.csv  
#   - gender_balance.csv
#
# WHAT THIS SCRIPT DOES
# ----------------------
# PART A - COMMUNITY BID PRICE PLOT
#   For each community, plots the average bid price per plant (bars + SE),
#   with individual offer-level bid prices overlaid as jittered
#   points.
#
# PART B - LANDRACE ADOPTION COMPENSATION MODEL
#   Builds an offer-level dataset with landrace use traits, landrace
#   size/growth metrics, and community characteristics attached. Fits
#   a weighted linear mixed-effects model (nlme::lme) of
#   log(compensation per plant), community random intercept, testing
#   whether compensation varied with individual use traits, growth
#   rate, number of plants offered, and community characteristics.
#   Model selection via dredge (AICc) picks a single best model (no
#   model averaging); diagnostics for that model, and effect plots
#   for its significant predictors only.
#
# OUTPUTS GENERATED
# ----------------------
# All saved to an "outputs" subfolder of the working directory
# (created automatically if it doesn't exist).
#
#   community_cost_fig_usd.png          Part A: average bid price by
#       community, with offer-level points (selected vs not).
#
#   trait_correlation_plot.png          Part B: Spearman correlations
#       among landrace use traits and growth metrics.
#
#   bid_price_range_summary.csv         Part B: n bids, n communities,
#       n landraces, and min/max/mean bid price per plant.
#
#   bid_anova_results.csv               Part B: F, df, and p-value for
#       the one-way community and landrace ANOVAs.
#
#   null_model_variance_components.csv  Part B: community/landrace/
#       residual variance breakdown from the unconditional model.
#
#   predictor_vif.csv                   Part B: variance inflation
#       factors for every candidate predictor in global_mod.
#
#   dredge_model_selection_results.csv  Part B: full model comparison
#       table (AICc ranking) from dredge().
#
#   main_model_coefficients.csv         Part B: coefficients, SEs,
#       t/p-values, and 95% CIs for best_mod (the single top-ranked
#       model from dredge) - the sole reported model.
#
#   main_model_variance_icc.csv         Part B: community/residual
#       variance, ICC, and marginal/conditional R² for best_mod.
#
#   model_diagnostics_residuals_qq.png  Part B: residuals vs fitted,
#       QQ plot of residuals, and QQ plot of community random effects.
#
#   main_model_effect_plots_updated.png Part B: predicted compensation
#       per plant for each SIGNIFICANT predictor in best_mod only.
#
# ===============================================================


# ----------------------------
# 0. Setup
# ----------------------------
setwd("~/Papers/PACS/Bid_heterogeneity_analysis/Clean_code")

# All generated outputs (figures, CSVs) are saved to this subfolder,
# separate from the many input files living in the working directory.
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(grid)
library(nlme)
library(lme4)
library(MuMIn)
library(emmeans)
library(ggeffects)
library(corrplot)
library(Hmisc)
library(car)

options(stringsAsFactors = FALSE)
options(na.action = "na.fail")

scale2 <- function(x) as.numeric(scale(x))

fmt_p <- function(x) {
  ifelse(x < 0.001, "p < 0.001", paste0("p = ", round(x, 3)))
}


# ===============================================================
# PART A: COMMUNITY BID PRICE PLOT
# ===============================================================

# ----------------------------
# 1. Load data
# ----------------------------
clean_offers <- read.csv("clean_offers.csv")
select_offers <- read.csv("select_offers.csv")
landrace_tender_crosswalk <- read.csv("landrace_tender_crosswalk_share.csv")

# ----------------------------
# 2. Link selected offers to community letter codes and landrace index
# ----------------------------
selection_join <- select_offers %>%
  left_join(landrace_tender_crosswalk, by = "Landrace_ID")

# Flag any selected offer that failed to match a community and/or
# landrace via the crosswalks (e.g. a spelling mismatch) - these
# would otherwise silently show up as "Not selected" below.
unmatched_selection <- selection_join %>%
  filter(is.na(Community) | is.na(Landrace_id_tender))

if (nrow(unmatched_selection) > 0) {
  warning(
    nrow(unmatched_selection), " row(s) in select_offers.csv could not be ",
    "matched to a  landrace via the crosswalk - see ",
    "unmatched_selection. These will NOT be flagged as selected offers below."
  )
  print(unmatched_selection)
} else {
  cat("All rows in select_offers.csv matched a community and landrace.\n")
}

selected_pairs <- selection_join %>%
  filter(!is.na(Community), !is.na(Landrace_id_tender)) %>%
  mutate(
    community_letter = str_remove(Community, "Community\\s+"),
    landrace_index = as.numeric(Landrace_id_tender)
  ) %>%
  select(community_letter, landrace_index) %>%
  distinct()


# ----------------------------
# 3. Prepare offer-level data
# ----------------------------
offers_plot <- clean_offers %>%
  mutate(
    price_per_plant_usd = (cost_etb / individuals) * 0.0183,
    community_letter = str_remove(community, "Community\\s+"),
    landrace_index = as.numeric(landrace_index)
  ) %>%
  filter(
    !is.na(community),
    !is.na(price_per_plant_usd),
    individuals > 0
  ) %>%
  left_join(
    selected_pairs %>% mutate(selected_offer_point = "Selected offer"),
    by = c("community_letter", "landrace_index")
  ) %>%
  mutate(
    selected_offer_point = if_else(
      is.na(selected_offer_point),
      "Not selected",
      selected_offer_point
    ),
    community_selected = if_else(
      community_letter %in% selected_pairs$community_letter,
      "Community with selected offer",
      "Community without selected offer"
    )
  )


# ----------------------------
# 4. Community summary (mean + SE)
# ----------------------------
community_summary <- offers_plot %>%
  group_by(community, community_selected) %>%
  dplyr::summarise(
    mean_price_usd = mean(price_per_plant_usd, na.rm = TRUE),
    n_bids = dplyr::n(),
    sd_price_usd = sd(price_per_plant_usd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se = ifelse(n_bids > 1, sd_price_usd / sqrt(n_bids), 0)
  ) %>%
  arrange(mean_price_usd)

# Order communities by mean price, low to high, and apply the same
# ordering to the offer-level data so points line up with their bars
community_levels <- community_summary$community

community_summary <- community_summary %>%
  mutate(community = factor(community, levels = community_levels))

offers_plot <- offers_plot %>%
  mutate(community = factor(community, levels = community_levels))


# ----------------------------
# 5. Plot: average bid price by community, with offer-level points
# ----------------------------
p_comm <- ggplot(
  community_summary,
  aes(x = community, y = mean_price_usd, fill = community_selected)
) +
  geom_col(width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_price_usd - se,
        ymax = mean_price_usd + se),
    width = 0.2,
    linewidth = 0.3,
    colour = "black"
  ) +
  geom_point(
    data = offers_plot,
    aes(
      x = community,
      y = price_per_plant_usd,
      colour = selected_offer_point,
      shape = selected_offer_point
    ),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.15, height = 0),
    size = 0.8,
    alpha = 1
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Community with selected offer" = "#1F4E79",
      "Community without selected offer" = "#6FA8DC"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Selected offer" = "black",
      "Not selected" = "grey40"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Selected offer" = 17,
      "Not selected" = 1
    )
  ) +
  scale_y_continuous(
    limits = c(0, 29),
    breaks = seq(0, 25, by = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Individual communities",
    y = "Average bid price per plant\n(USD, 2023 average rate: 1 ETB = $0.0183)",
    fill = "Community status",
    colour = "Landrace offer",
    shape = "Landrace offer"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = c(0.98, 0.05),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 9),
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.3
    ),
    legend.key.size = unit(0.35, "cm")
  )

p_comm


# ----------------------------
# 6. Save figure
# ----------------------------
ggsave(
  file.path(out_dir, "community_cost_fig_usd.png"),
  p_comm,
  width = 4,
  height = 3.5,
  dpi = 300,
  bg = "white"
)


# ===============================================================
# PART B: LANDRACE ADOPTION COMPENSATION MODEL
# ===============================================================

# ----------------------------
# 7. Load remaining data
# ----------------------------
# clean_offers is already loaded (step 1) - reused here rather than
# re-reading the same file.
use_data      <- read.csv("Use_data_final_share.csv")
size_data     <- read.csv("Landrace_size_data_clean_share.csv")

community_enset  <- read.csv("enset_per_farmer_by_community.csv")
community_gender <- read.csv("gender_balance.csv")

community_characteristics <- full_join(community_enset, community_gender, by = "community")
community_characteristics$prop_women <- community_characteristics$women / community_characteristics$farmers


# ----------------------------
# 7b. Consistency check: community codes match between
# community_characteristics and clean_offers
# ----------------------------
offer_communities <- clean_offers %>% distinct(community)

unmatched_offer_communities <- offer_communities %>%
  anti_join(community_characteristics, by = "community")

if (nrow(unmatched_offer_communities) > 0) {
  warning(
    nrow(unmatched_offer_communities), " community/communities in ",
    "clean_offers.csv have no matching row in community_characteristics ",
    "(enset/gender data) - see unmatched_offer_communities. These will get ",
    "NA for mean_enset_per_farmer, farmers, and prop_women after the join ",
    "in step 10."
  )
  print(unmatched_offer_communities)
} else {
  cat("All communities in clean_offers.csv matched community_characteristics.\n")
}


# ----------------------------
# 7c. Build landrace index lookup
# ----------------------------
# clean_offers only has landrace_index (numeric), not a landrace ID -
# landrace_crosswalk (Landrace_ID <-> Landrace_id_tender, already
# loaded in Part A) provides that mapping.
landrace_index_lookup <- landrace_tender_crosswalk %>%
  rename(Landrace_ID = Landrace_ID, landrace_index = Landrace_id_tender) %>%
  select(Landrace_ID, landrace_index)


# ----------------------------
# 8. Calculate landrace-level size metrics
# ----------------------------
size_summary <- size_data %>%
  rename(
    age_years         = Age.in.years,
    maturity_age      = Maturity.harvest.age,
    pseudostem_height = Pseudostem.height,
    circ_base         = Pseudostem.circumference.at.the.base,
    circ_middle       = Pseudostem.circumference.at.the.middle,
    circ_top          = Pseudostem.circumference.at.the.top
  ) %>%
  mutate(
    r_base   = circ_base / (2 * pi),
    r_middle = circ_middle / (2 * pi),
    r_top    = circ_top / (2 * pi),
    
    half_height = pseudostem_height / 2,
    
    volume_lower = (1 / 3) * pi * half_height *
      (r_base^2 + r_base * r_middle + r_middle^2),
    
    volume_upper = (1 / 3) * pi * half_height *
      (r_middle^2 + r_middle * r_top + r_top^2),
    
    pseudostem_volume = volume_lower + volume_upper,
    
    pseudostem_volume_growth = pseudostem_volume / age_years
  ) %>%
  group_by(Landrace_ID) %>%
  summarise(
    pseudostem_volume_growth_mean = mean(pseudostem_volume_growth, na.rm = TRUE),
    maturity_age_mean             = mean(maturity_age, na.rm = TRUE),
    .groups = "drop"
  )


# ----------------------------
# 9. Build landrace predictor table + trait correlation plot
# ----------------------------
landrace_predictors <- use_data %>%
  left_join(size_summary, by = "Landrace_ID")

main_traits <- landrace_predictors %>%
  select(
    Kocho,
    Bulla,
    Amicho,
    Medicine,
    Fodder,
    Fibre,
    Drought,
    Disease,
    pseudostem_volume_growth_mean,
    maturity_age_mean
  ) %>%
  mutate(across(everything(), as.numeric))

cor_results <- rcorr(as.matrix(main_traits), type = "spearman")

r_mat <- cor_results$r
p_mat <- cor_results$P

label_map <- c(
  Kocho = "Kocho",
  Bulla = "Bulla",
  Amicho = "Amicho",
  Medicine = "Medicine",
  Fodder = "Fodder",
  Fibre = "Fibre",
  Drought = "Drought",
  Disease = "Disease",
  pseudostem_volume_growth_mean = "Pseudostem volume",
  maturity_age_mean = "Maturity age"
)

rownames(r_mat) <- dplyr::recode(rownames(r_mat), !!!label_map)
colnames(r_mat) <- dplyr::recode(colnames(r_mat), !!!label_map)

rownames(p_mat) <- dplyr::recode(rownames(p_mat), !!!label_map)
colnames(p_mat) <- dplyr::recode(colnames(p_mat), !!!label_map)

n <- ncol(r_mat)

png(
  filename = file.path(out_dir, "trait_correlation_plot.png"),
  width = 9,
  height = 8,
  units = "in",
  res = 300
)

corrplot(
  r_mat,
  method = "color",
  type = "upper",
  order = "original",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.7,
  tl.offset = 1,
  diag = FALSE,
  col = colorRampPalette(c("#D55E00", "white", "#0072B2"))(100),
  cl.length = 5,
  mar = c(0, 0, 0, 0)
)

cells <- which(upper.tri(r_mat), arr.ind = TRUE)
x <- cells[, 2]
y <- n - cells[, 1] + 1

labs <- sprintf("%.2f", r_mat[cells])
fnts <- ifelse(!is.na(p_mat[cells]) & p_mat[cells] < 0.05, 2, 1)

text(
  x = x,
  y = y,
  labels = labs,
  font = fnts,
  cex = 0.8,
  col = "black",
  xpd = NA
)

dev.off()


# ----------------------------
# 10. Build offer-level dataset
# ----------------------------

# Consistency check: does every landrace_index in clean_offers have a
# match in landrace_index_lookup? Unmatched rows get NA for Landrace
# and are silently dropped by the !is.na(Landrace) filter below, so
# this makes that visible before it happens.
unmatched_landrace_index <- clean_offers %>%
  distinct(landrace_index) %>%
  anti_join(landrace_index_lookup, by = "landrace_index")

if (nrow(unmatched_landrace_index) > 0) {
  warning(
    nrow(unmatched_landrace_index), " landrace_index value(s) in ",
    "clean_offers.csv have no match in landrace_index_lookup - see ",
    "unmatched_landrace_index. These rows will be dropped below."
  )
  print(unmatched_landrace_index)
} else {
  cat("All landrace_index values in clean_offers.csv matched landrace_index_lookup.\n")
}

offers_model <- clean_offers %>%
  select(community, landrace_index, individuals, cost_etb) %>%
  left_join(
    landrace_index_lookup,
    by = "landrace_index"
  ) %>%
  mutate(
    cost_per_plant_offer_etb = cost_etb / individuals,
    cost_per_plant_offer_usd = cost_per_plant_offer_etb * 0.0183
  ) %>%
  filter(
    is.finite(individuals),
    individuals > 0,
    is.finite(cost_etb),
    is.finite(cost_per_plant_offer_etb),
    !is.na(Landrace_ID)
  )

# Consistency check: does every Landrace name here have a match in
# landrace_predictors? Unlike the check above, an unmatched Landrace
# is NOT filtered out below - it would silently carry NA through every
# use-trait column and only surface later as a na.fail error.
unmatched_landrace_predictors <- offers_model %>%
  distinct(Landrace_ID) %>%
  anti_join(landrace_predictors, by = "Landrace_ID")

if (nrow(unmatched_landrace_predictors) > 0) {
  warning(
    nrow(unmatched_landrace_predictors), " landrace(s) in the offer data ",
    "have no match in landrace_predictors (use_data/size_data) - see ",
    "unmatched_landrace_predictors. These will carry NA use-trait values ",
    "and will cause lme()/dredge() to fail (na.action = \"na.fail\") unless ",
    "resolved first."
  )
  print(unmatched_landrace_predictors)
} else {
  cat("All landraces in the offer data matched landrace_predictors.\n")
}

offers_model <- offers_model %>%
  left_join(landrace_predictors, by = "Landrace_ID") %>%
  mutate(
    community = factor(community)
  ) %>%
  left_join(community_characteristics %>%
              select(community, mean_enset_per_farmer, farmers, prop_women),
            by = "community")

# Overall bid-level range and sample sizes (not group averages) - the
# headline numbers reported in the results text (n bids, n landraces,
# n communities, min/max/mean price per plant across all bids).
bid_summary <- offers_model %>%
  summarise(
    n_bids = n(),
    n_communities = n_distinct(community),
    n_landraces = n_distinct(Landrace_ID),
    min_price_usd = min(cost_per_plant_offer_usd, na.rm = TRUE),
    max_price_usd = max(cost_per_plant_offer_usd, na.rm = TRUE),
    mean_price_usd = mean(cost_per_plant_offer_usd, na.rm = TRUE)
  )
print(bid_summary)

write_csv(bid_summary, file.path(out_dir, "bid_price_range_summary.csv"))


# ----------------------------
# 11. Prepare model variables
# ----------------------------
# Report how many offers have missing size data (from
# Landrace_size_data_clean.csv); these are mean-imputed below.
missing_size_offers <- offers_model %>%
  filter(is.na(pseudostem_volume_growth_mean) | is.na(maturity_age_mean))

missing_size_landraces <- missing_size_offers %>% distinct(Landrace_ID)

message(
  "Size data mean-imputed for ", nrow(missing_size_offers), " offer(s) ",
  "covering ", nrow(missing_size_landraces), " landrace(s) ",
  "(IDs in missing_size_landraces)."
)
# Keep unscaled versions of log_individuals and pseudostem growth so
# effect-plot axes can be shown in original units.
offers_model <- offers_model %>%
  mutate(
    across(
      c(
        pseudostem_volume_growth_mean,
        maturity_age_mean
      ),
      ~ replace(., is.na(.), mean(., na.rm = TRUE))
    )
  ) %>%
  mutate(
    log_individuals_unscaled = log(individuals),
    pseudostem_volume_growth_mean_unscaled = pseudostem_volume_growth_mean,
    prop_women_unscaled = prop_women
  ) %>%
  mutate(
    log_individuals = scale2(log_individuals_unscaled),
    pseudostem_volume_growth_mean = scale2(pseudostem_volume_growth_mean_unscaled),
    maturity_age_mean = scale2(maturity_age_mean),
    mean_enset_per_farmer = scale2(mean_enset_per_farmer),
    farmers = scale2(farmers),
    prop_women = scale2(prop_women_unscaled)
  ) %>%
  mutate(
    across(
      c(
        Kocho, Bulla, Amicho, Medicine,
        Fodder, Fibre,
        Drought, Disease
      ),
      as.numeric
    ),
    Kocho_c = as.numeric(scale(Kocho, scale = FALSE))
  )

# Means and SDs for back-transforming effect-plot x-axes
mean_log_ind <- mean(offers_model$log_individuals_unscaled)
sd_log_ind   <- sd(offers_model$log_individuals_unscaled)

mean_ps <- mean(offers_model$pseudostem_volume_growth_mean_unscaled)
sd_ps   <- sd(offers_model$pseudostem_volume_growth_mean_unscaled)

mean_wm <- mean(offers_model$prop_women_unscaled)
sd_wm   <- sd(offers_model$prop_women_unscaled)


# ----------------------------
# 11b. Variance components check
# ----------------------------
mod_null <- lmer(
  log(cost_per_plant_offer_usd + 1) ~ 1 +
    (1 | community) + (1 | Landrace_ID),
  data = offers_model
)

vc <- as.data.frame(VarCorr(mod_null))

# Named with a _null suffix so these aren't overwritten by the
# best_mod variance components computed later (step 14) - both sets
# get saved to their own CSV.
community_var_null <- vc$vcov[vc$grp == "community"]
landrace_var_null   <- vc$vcov[vc$grp == "Landrace_ID"]
resid_var_null      <- vc$vcov[vc$grp == "Residual"]

total_var_null <- community_var_null + landrace_var_null + resid_var_null

pct_community <- community_var_null / total_var_null * 100
pct_landrace  <- landrace_var_null  / total_var_null * 100
pct_residual  <- resid_var_null     / total_var_null * 100

mod_comm_aov <- aov(
  log(cost_per_plant_offer_usd + 1) ~ community,
  data = offers_model
)
summary(mod_comm_aov)

mod_landrace_aov <- aov(
  log(cost_per_plant_offer_usd + 1) ~ Landrace_ID,
  data = offers_model
)
summary(mod_landrace_aov)

# Tidied F-test results (F statistic, df, p-value) for both one-way
# ANOVAs, matching the numbers reported in the results text.
comm_aov_table <- summary(mod_comm_aov)[[1]]
landrace_aov_table <- summary(mod_landrace_aov)[[1]]

bid_anova_results <- tibble(
  test = c("community", "landrace"),
  df1 = c(comm_aov_table$Df[1], landrace_aov_table$Df[1]),
  df2 = c(comm_aov_table$Df[2], landrace_aov_table$Df[2]),
  f_value = c(comm_aov_table$`F value`[1], landrace_aov_table$`F value`[1]),
  p_value = c(comm_aov_table$`Pr(>F)`[1], landrace_aov_table$`Pr(>F)`[1])
)

write_csv(bid_anova_results, file.path(out_dir, "bid_anova_results.csv"))

null_model_variance_components <- tibble(
  source = c("community", "landrace", "residual"),
  variance = c(community_var_null, landrace_var_null, resid_var_null),
  pct_of_total = c(pct_community, pct_landrace, pct_residual)
)

write_csv(null_model_variance_components, file.path(out_dir, "null_model_variance_components.csv"))


# ----------------------------
# 11c. Variance inflation factors (VIF)
# ----------------------------
# Flags predictors whose SE is inflated by correlation with the other
# predictors (VIF > 5). 
vif_proxy_mod <- lm(
  log(cost_per_plant_offer_usd + 1) ~
    log_individuals +
    Kocho_c +
    Bulla +
    Amicho +
    Fodder +
    Fibre +
    Drought +
    Disease +
    Medicine +
    pseudostem_volume_growth_mean +
    maturity_age_mean +
    mean_enset_per_farmer +
    farmers +
    prop_women,
  data = offers_model
)

vif_raw <- car::vif(vif_proxy_mod)

if (is.matrix(vif_raw)) {
  vif_table <- as.data.frame(vif_raw) %>%
    rownames_to_column("predictor") %>%
    mutate(vif = `GVIF^(1/(2*Df))`^2) %>%
    select(predictor, vif)
} else {
  vif_table <- enframe(vif_raw, name = "predictor", value = "vif")
}

print(vif_table %>% arrange(desc(vif)))

write_csv(vif_table, file.path(out_dir, "predictor_vif.csv"))


# ----------------------------
# 12. Global model for dredge
# ----------------------------
global_mod <- lme(
  fixed = log(cost_per_plant_offer_usd + 1) ~
    log_individuals +
    Kocho_c +
    Bulla +
    Amicho +
    Fodder +
    Fibre +
    Drought +
    Disease +
    pseudostem_volume_growth_mean +
    maturity_age_mean +
    Kocho_c:pseudostem_volume_growth_mean +
    mean_enset_per_farmer +
    farmers +
    prop_women,
  random = ~ 1 | community,
  weights = varPower(form = ~ fitted(.)),
  data = offers_model,
  method = "ML"
)


# ----------------------------
# 13. Model selection - single best model, no averaging
# ----------------------------
# subset =: interaction terms are only tested in models that also
# include both of their main effects 
dredge_results <- dredge(
  global_mod,
  rank = "AICc",
  subset =
    dc("Kocho_c:pseudostem_volume_growth_mean", "Kocho_c") &&
    dc("Kocho_c:pseudostem_volume_growth_mean", "pseudostem_volume_growth_mean")
)

head(dredge_results, 10)

write_csv(as.data.frame(dredge_results), file.path(out_dir, "dredge_model_selection_results.csv"))

best_mod <- get.models(dredge_results, subset = 1)[[1]]
summary(best_mod)


# ----------------------------
# 14. Final reported model
# ----------------------------
# best_mod (the single top-ranked model from dredge, step 13) is used
# directly below 
r.squaredGLMM(best_mod)
intervals(best_mod)

# Tidied coefficient table, with confidence intervals from intervals()
best_mod_ci <- as.data.frame(intervals(best_mod)$fixed) %>%
  rownames_to_column("term") %>%
  rename(ci_lower = lower, ci_upper = upper) %>%
  select(term, ci_lower, ci_upper)

main_model_coefficients <- as.data.frame(summary(best_mod)$tTable) %>%
  rownames_to_column("term") %>%
  rename(
    estimate = Value,
    std_error = `Std.Error`,
    df = DF,
    t_value = `t-value`,
    p_value = `p-value`
  ) %>%
  left_join(best_mod_ci, by = "term")

write_csv(main_model_coefficients, file.path(out_dir, "main_model_coefficients.csv"))

vc <- VarCorr(best_mod)

community_var <- as.numeric(vc["(Intercept)", "Variance"])
resid_var <- as.numeric(vc["Residual", "Variance"])

icc <- community_var / (community_var + resid_var)
icc

r2_glmm <- r.squaredGLMM(best_mod)

main_model_variance_icc <- tibble(
  community_variance = community_var,
  residual_variance = resid_var,
  icc = icc,
  r2_marginal = r2_glmm[1, "R2m"],
  r2_conditional = r2_glmm[1, "R2c"]
)

write_csv(main_model_variance_icc, file.path(out_dir, "main_model_variance_icc.csv"))


# ----------------------------
# 15. Diagnostics
# ----------------------------
offers_model$fitted   <- fitted(best_mod)
offers_model$resid    <- resid(best_mod, type = "normalized")
offers_model$pred_etb <- exp(fitted(best_mod))

community_re <- ranef(best_mod)[, 1]

# Largest residuals - quick check for problem observations
offers_model[order(abs(offers_model$resid), decreasing = TRUE), ][1:10, ]

p_resid_fitted <- ggplot(offers_model, aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Fitted values", y = "Residuals", tag = "A")

qq_data <- data.frame(
  sample = sort(offers_model$resid),
  theoretical = qnorm(ppoints(length(offers_model$resid)))
)

p_qq <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_abline(
    intercept = mean(offers_model$resid, na.rm = TRUE),
    slope = sd(offers_model$resid, na.rm = TRUE),
    linetype = "dashed"
  ) +
  theme_classic() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", tag = "B")

# QQ plot of random effects (community intercepts) - assumption check
# alongside the residual plots above, per the methods text.
qq_re_data <- data.frame(
  sample = sort(community_re),
  theoretical = qnorm(ppoints(length(community_re)))
)

p_qq_re <- ggplot(qq_re_data, aes(x = theoretical, y = sample)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_abline(
    intercept = mean(community_re, na.rm = TRUE),
    slope = sd(community_re, na.rm = TRUE),
    linetype = "dashed"
  ) +
  theme_classic() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", tag = "C")

p_diagnostics <- p_resid_fitted + p_qq + p_qq_re + plot_layout(ncol = 3)

p_diagnostics

ggsave(
  filename = file.path(out_dir, "model_diagnostics_residuals_qq.png"),
  plot = p_diagnostics,
  width = 10,
  height = 3.5,
  dpi = 300
)


# ----------------------------
# 16. Effect plots - significant predictors only
# ----------------------------

plot_drought <- ggpredict(best_mod, terms = "Drought")
plot_fodder <- ggpredict(best_mod, terms = "Fodder")
plot_individuals <- ggpredict(best_mod, terms = "log_individuals [all]")
plot_women <- ggpredict(best_mod, terms = "prop_women [all]")

kocho_vals <- sort(unique(offers_model$Kocho_c))
plot_interaction <- ggpredict(
  best_mod,
  terms = c(
    "pseudostem_volume_growth_mean [all]",
    paste0("Kocho_c [", paste(kocho_vals, collapse = ", "), "]")
  )
)

# Relabel discrete axes / groups
plot_drought$x <- factor(
  plot_drought$x,
  levels = c(0, 1, 2),
  labels = c("Not used", "Used", "Preferred")
)

plot_fodder$x <- factor(
  plot_fodder$x,
  levels = c(0, 1, 2),
  labels = c("Not used", "Used", "Preferred")
)

plot_interaction$group <- factor(
  plot_interaction$group,
  labels = c("Not used", "Used", "Preferred")
)

use_cols <- c(
  "Not used" = "grey80",
  "Used" = "orange2",
  "Preferred" = "pink4"
)

p_drought <- ggplot(plot_drought, aes(x = x, y = predicted, colour = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, colour = x), width = 0.08) +
  scale_colour_manual(values = use_cols) +
  scale_y_continuous(limits = c(1, 6.5), breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "β = ", round(summary(best_mod)$tTable["Drought", "Value"], 3),
      "\n", fmt_p(summary(best_mod)$tTable["Drought", "p-value"])
    ),
    hjust = 1.05, vjust = 1.2, size = 3.2
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Drought tolerance", y = NULL)

p_fodder <- ggplot(plot_fodder, aes(x = x, y = predicted, colour = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, colour = x), width = 0.08) +
  scale_colour_manual(values = use_cols) +
  scale_y_continuous(limits = c(1, 6.5), breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "β = ", round(summary(best_mod)$tTable["Fodder", "Value"], 3),
      "\n", fmt_p(summary(best_mod)$tTable["Fodder", "p-value"])
    ),
    hjust = 1.05, vjust = 1.2, size = 3.2
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Fodder suitability", y = NULL)

ind_breaks_orig <- c(10, 50, 100, 250, 500, 1000)
ind_breaks_scaled <- (log(ind_breaks_orig) - mean_log_ind) / sd_log_ind

p_plants <- ggplot(plot_individuals, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  scale_x_continuous(breaks = ind_breaks_scaled, labels = ind_breaks_orig) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), expand = c(0, 0)) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "β = ", round(summary(best_mod)$tTable["log_individuals", "Value"], 3),
      "\n", fmt_p(summary(best_mod)$tTable["log_individuals", "p-value"])
    ),
    hjust = 1.05, vjust = 1.2, size = 3.2
  ) +
  theme_classic() +
  labs(x = "Number of plants offered", y = NULL)

breaks_orig_women <- pretty(offers_model$prop_women_unscaled, n = 4)
breaks_scaled_women <- (breaks_orig_women - mean_wm) / sd_wm

p_women <- ggplot(plot_women, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1, colour = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  scale_x_continuous(breaks = breaks_scaled_women, labels = breaks_orig_women) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), expand = c(0, 0)) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "β = ", round(summary(best_mod)$tTable["prop_women", "Value"], 3),
      "\n", fmt_p(summary(best_mod)$tTable["prop_women", "p-value"])
    ),
    hjust = 1.05, vjust = 1.2, size = 3.2
  ) +
  theme_classic() +
  labs(x = "Proportion of female farmers", y = NULL)

ps_breaks_orig <- pretty(offers_model$pseudostem_volume_growth_mean_unscaled, n = 4)
ps_breaks_scaled <- (ps_breaks_orig - mean_ps) / sd_ps

p_interaction <- ggplot(
  plot_interaction,
  aes(x = x, y = predicted, colour = group, fill = group)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = use_cols, name = "Kocho suitability") +
  scale_fill_manual(values = use_cols, name = "Kocho suitability") +
  scale_x_continuous(
    breaks = ps_breaks_scaled,
    labels = format(ps_breaks_orig, big.mark = ",", trim = TRUE)
  ) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5), expand = c(0, 0)) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "β = ", round(
        summary(best_mod)$tTable["Kocho_c:pseudostem_volume_growth_mean", "Value"], 3
      ),
      "\n", fmt_p(
        summary(best_mod)$tTable["Kocho_c:pseudostem_volume_growth_mean", "p-value"]
      )
    ),
    hjust = 1.05, vjust = 1.2, size = 3.2
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.01, 1.03),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(-0.15, "cm"),
    legend.key.height = unit(0.3, "cm")
  ) +
  labs(
    x = expression(atop("Pseudostem volume growth", "(" * cm^3 * " year"^{-1} * ")")),
    y = NULL
  )

p_interaction

common_margin <- theme(plot.margin = margin(t = 12, r = 8, b = 12, l = 8))

p_drought <- p_drought + common_margin
p_fodder <- p_fodder + common_margin
p_women <- p_women + common_margin
p_plants <- p_plants + common_margin
p_interaction <- p_interaction + common_margin

main_panel_updated <-
  (p_drought | p_fodder) /
  (p_women | p_plants) /
  p_interaction

final_effect_plots_updated <-
  wrap_elements(
    grid::textGrob(
      "Predicted compensation per plant (USD)",
      rot = 90,
      gp = gpar(fontsize = 12)
    )
  ) +
  main_panel_updated +
  plot_layout(widths = c(0.06, 1))

final_effect_plots_updated

ggsave(
  file.path(out_dir, "main_model_effect_plots_updated.png"),
  final_effect_plots_updated,
  width = 4.7,
  height = 7,
  dpi = 300
)

