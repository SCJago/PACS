# ============================================================
# PACS ENSET TENDER OPTIMISATION
# ============================================================
#
# FILES NEEDED
# ----------------------
#   - clean_offers.csv                      (cleaned tender
#                                             offers)
#   - gender_balance.csv                    (community, farmers, women)
#   - landrace_name_id_crosswalk.csv        (Landrace_name <-> Landrace_id_tender)
#   - landrace_avg_status_landscape.csv     (Landrace_name, landscape_category)
#   - community_name_crosswalk.csv          (Community <-> Community_name)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Takes the cleaned tender offers (individuals, cost per landrace per
# community) and joins gender balance and landscape rarity data.
#
# Solves a lexicographic MILP over the offers:
#   Stage 1: maximise number of landraces selected (within budget)
#   Stage 2: among those, minimise total cost
#   Stage 3: among those, maximise individuals
#   Stage 4: among those, maximise rarity (tie-break)
#
# OUTPUTS GENERATED
# ----------------------
#   select_offers.csv   final selected offers: community, landrace,
#       plants, cost (Birr/GBP), and number of farmers.
#
# ============================================================


# ----------------------------
# 0. Set working directory
# ----------------------------
setwd("~/Papers/PACS/Model_offer_selection/Clean_code")


# ----------------------------
# 1. Load packages
# ----------------------------
library(readr)
library(dplyr)

library(ompr)
library(ompr.roi)
library(ROI)
library(ROI.plugin.glpk)


# ----------------------------
# 2. Load cleaned offers
# ----------------------------
clean_offers <- read.csv("clean_offers.csv")


# ----------------------------
# 3. Offer summaries
# ----------------------------
offers_corrected_summary <- clean_offers %>%
  group_by(community) %>%
  summarise(n_offers = n(), .groups = "drop") %>%
  arrange(desc(n_offers))

offers_corrected_totals <- clean_offers %>% summarise(
  n_offers = n(),
  n_communities = n_distinct(community),
  n_landraces = n_distinct(landrace_index),
  total_possible_individuals = sum(individuals),
  total_possible_cost_etb = sum(cost_etb)
)

print(offers_corrected_summary, n = 22)
print(offers_corrected_totals)


# ----------------------------
# 4. Load gender balance data
# ----------------------------
gender_df <- read_csv("gender_balance.csv", show_col_types = FALSE) %>%
  mutate(
    community = as.character(community),
    farmers = as.numeric(farmers),
    women = as.numeric(women),
    pct_women = women / farmers
  ) %>%
  # Ensure it matches the cleaned offers set (Community H already removed)
  semi_join(clean_offers %>% distinct(community), by = "community")

farmers_by_comm <- setNames(gender_df$farmers, gender_df$community)
women_by_comm   <- setNames(gender_df$women,   gender_df$community)

print(gender_df, n = Inf)


# ----------------------------
# 5. Budget / sets
# ----------------------------
budget_gbp <- 9050
fx <- 68.03
budget_etb <- budget_gbp * fx

max_lr_per_comm <- 4
min_female_share <- 0.35

comms <- sort(unique(clean_offers$community))
lrs   <- sort(unique(clean_offers$landrace_index))


# ----------------------------
# 6. Base MILP builder (constraints only)
# ----------------------------
build_base_model <- function(
    offers,
    comms,
    lrs,
    budget_etb,
    max_lr_per_comm = 4,
    farmers_by_comm = NULL,
    women_by_comm   = NULL,
    min_female_share = NULL
) {
  
  stopifnot(all(c("id","community","landrace_index","cost_etb") %in% names(offers)))
  if (!is.null(min_female_share)) stopifnot(!is.null(farmers_by_comm), !is.null(women_by_comm))
  
  bigC <- length(lrs)    # big-M for community linking
  bigL <- length(comms)  # big-M for landrace linking
  
  m <- MIPModel() %>%
    add_variable(x[i], i = offers$id, type = "binary") %>%
    add_variable(y[c], c = comms,     type = "binary") %>%
    add_variable(z[l], l = lrs,       type = "binary") %>%
    
    # Budget
    add_constraint(sum_expr(offers$cost_etb[i] * x[i], i = offers$id) <= budget_etb) %>%
    
    # Link community selection y[c] <-> any chosen offers in community c
    add_constraint(sum_expr(x[i], i = offers$id[offers$community == c]) <= bigC * y[c], c = comms) %>%
    add_constraint(sum_expr(x[i], i = offers$id[offers$community == c]) >= y[c],       c = comms) %>%
    
    # Link landrace selection z[l] <-> any chosen offers for landrace l
    add_constraint(sum_expr(x[i], i = offers$id[offers$landrace_index == l]) <= bigL * z[l], l = lrs) %>%
    add_constraint(sum_expr(x[i], i = offers$id[offers$landrace_index == l]) >= z[l],       l = lrs) %>%
    
    # Max 4 selected landraces per community
    add_constraint(sum_expr(x[i], i = offers$id[offers$community == c]) <= max_lr_per_comm, c = comms) %>%
    
    # At most one selected offer per landrace (one community per landrace)
    add_constraint(sum_expr(x[i], i = offers$id[offers$landrace_index == l]) <= 1, l = lrs)
  
  # Gender constraint: overall female share across selected communities
  if (!is.null(min_female_share)) {
    m <- m %>%
      add_constraint(
        sum_expr(women_by_comm[c] * y[c],   c = comms) >=
          min_female_share * sum_expr(farmers_by_comm[c] * y[c], c = comms)
      )
  }
  
  m
}

base_model <- build_base_model(
  offers = clean_offers,
  comms = comms,
  lrs = lrs,
  budget_etb = budget_etb,
  max_lr_per_comm = max_lr_per_comm,
  farmers_by_comm = farmers_by_comm,
  women_by_comm = women_by_comm,
  min_female_share = min_female_share
)

summary(base_model)


# ----------------------------
# 7. Lexicographic optimisation (Stages 1-3)
# ----------------------------

# ---- Stage 1: maximise number of landraces selected ----
m1 <- base_model %>%
  set_objective(sum_expr(z[l], l = lrs), "max")

r1 <- solve_model(m1, with_ROI(solver = "glpk"))

max_landraces <- get_solution(r1, z[l]) %>%
  filter(value > 0.5) %>%
  nrow()

cat("Max landraces:", max_landraces, "\n")

# ---- Stage 2: minimise cost among max-landrace solutions ----
m2 <- base_model %>%
  add_constraint(sum_expr(z[l], l = lrs) == max_landraces) %>%
  set_objective(sum_expr(clean_offers$cost_etb[i] * x[i], i = clean_offers$id), "min")

r2 <- solve_model(m2, with_ROI(solver = "glpk"))

min_cost <- get_solution(r2, x[i]) %>%
  filter(value > 0.5) %>%
  left_join(clean_offers, by = c("i" = "id")) %>%
  summarise(total = sum(cost_etb)) %>%
  pull(total)

cat("Min cost at max landraces (ETB):", min_cost, "\n")
cat("Min cost at max landraces (GBP):", min_cost / fx, "\n")

# ---- Stage 3: maximise individuals among (max landraces, min cost) ----
m3 <- base_model %>%
  add_constraint(sum_expr(z[l], l = lrs) == max_landraces) %>%
  add_constraint(sum_expr(clean_offers$cost_etb[i] * x[i], i = clean_offers$id) == min_cost) %>%
  set_objective(sum_expr(clean_offers$individuals[i] * x[i], i = clean_offers$id), "max")

r3 <- solve_model(m3, with_ROI(solver = "glpk"))

selected_stage3 <- get_solution(r3, x[i]) %>%
  filter(value > 0.5) %>%
  left_join(clean_offers, by = c("i" = "id"))

max_individuals <- selected_stage3 %>%
  summarise(total = sum(individuals)) %>%
  pull(total)

cat("Max individuals at (max landraces, min cost):", max_individuals, "\n")


# ----------------------------
# 8. Stage 4 tie-break: maximise rarity
# ----------------------------
landrace_name_id_crosswalk <- read.csv("landrace_name_id_crosswalk.csv")

landscape_rarity <- read.csv("landrace_avg_status_landscape.csv")
landscape_rarity <- subset(landscape_rarity, select = c(Landrace_name, landscape_category))

landscape_rarity_w_id <- left_join(landrace_name_id_crosswalk, landscape_rarity,
                                   by = c("Landrace_name"))

rarity_tbl <- landscape_rarity_w_id %>%
  transmute(
    landrace_index = as.integer(Landrace_id_tender),
    landscape_category = as.character(landscape_category)
  ) %>%
  mutate(
    rarity_score = case_when(
      landscape_category == "Common"     ~ 1L,
      landscape_category == "Restricted" ~ 2L,
      landscape_category == "Endangered" ~ 3L,
      landscape_category == "Critical"   ~ 4L,
      landscape_category == "Lost"       ~ 5L,
      TRUE ~ 0L
    )
  ) %>%
  distinct(landrace_index, .keep_all = TRUE)

offers_rarity <- clean_offers %>%
  left_join(rarity_tbl %>% select(landrace_index, rarity_score), by = "landrace_index") %>%
  mutate(rarity_score = ifelse(is.na(rarity_score), 0L, rarity_score))

m4 <- base_model %>%
  add_constraint(sum_expr(z[l], l = lrs) == max_landraces) %>%
  add_constraint(sum_expr(offers_rarity$cost_etb[i] * x[i], i = offers_rarity$id) == min_cost) %>%
  add_constraint(sum_expr(offers_rarity$individuals[i] * x[i], i = offers_rarity$id) == max_individuals) %>%
  set_objective(sum_expr(offers_rarity$rarity_score[i] * x[i], i = offers_rarity$id), "max")

r4 <- solve_model(m4, with_ROI(solver = "glpk"))

selected_rarity <- get_solution(r4, x[i]) %>%
  filter(value > 0.5) %>%
  left_join(offers_rarity, by = c("i" = "id")) %>%
  arrange(community, landrace_index)


# ----------------------------
# 9. Output / checks
# ----------------------------
final_summary <- selected_rarity %>%
  summarise(
    selected_offers   = n(),
    n_landraces       = n_distinct(landrace_index),
    n_communities     = n_distinct(community),
    total_individuals = sum(individuals),
    total_cost_etb    = sum(cost_etb),
    total_cost_gbp    = sum(cost_etb) / fx,
    cost_per_individual = total_individuals / total_cost_gbp
  )
print(final_summary)

# Gender check on selected communities
gender_check <- selected_rarity %>%
  distinct(community) %>%
  left_join(gender_df, by = "community") %>%
  summarise(
    farmers = sum(farmers),
    women   = sum(women),
    pct_women = women / farmers
  )
print(gender_check)

# Landraces per community (should be <= 4)
print(
  selected_rarity %>%
    count(community, name = "landraces_in_comm") %>%
    arrange(desc(landraces_in_comm))
)

# Final selected pairs
sel_pairs_rarity <- selected_rarity %>%
  distinct(community, landrace_index) %>%
  arrange(community, landrace_index)

print(sel_pairs_rarity)


# ----------------------------
# 10. Build select_offers.csv (final selected offers)
# ----------------------------
community_name_crosswalk <- read.csv("community_name_crosswalk.csv")

select_offers <- selected_rarity %>%
  left_join(community_name_crosswalk, by = c("community" = "Community")) %>%
  left_join(landrace_name_id_crosswalk, by = c("landrace_index" = "Landrace_id_tender")) %>%
  left_join(gender_df %>% select(community, farmers), by = "community") %>%
  transmute(
    Community_name = Community_name,
    Landrace_name  = Landrace_name,
    Plants         = individuals,
    Birr           = cost_etb,
    GBP            = round(cost_etb / fx),
    No_farmers     = farmers
  ) %>%
  arrange(Community_name, Landrace_name)

write.csv(select_offers, "select_offers.csv", row.names = FALSE)
print(select_offers)
