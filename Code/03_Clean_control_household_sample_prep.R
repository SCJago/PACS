# ============================================================
# PACS CONTROL HOUSEHOLDS: LANDRACE CLEANING AND ANALYSIS
# ============================================================
#
# FILES NEEDED 
#  (not provided for confidentiality - control hh data without coordinates 
#  provided as "clean_control_hh_data_sharenolatlon.csv") 
#  (Required output for next script "hh_boot_35-37E.rds" is provided as this does not
#  contain identifying information)
#
# ----------------------
#   - clean_control_hh_data.csv     (control household data: T1/T2
#                                     landrace lists, coordinates, years)
#   - surveyed_kebeles_PACS.shp     (+ associated shapefile components;
#                                     boundary of surveyed/treated kebeles)
#   - Rift_valley_boundary.kml      (Rift Valley boundary)
#   - agroecological_zones.tif      (raster of Ethiopia's agroecological
#                                     zones)
#
# WHAT THIS SCRIPT DOES
# ----------------------
# - builds the geographically-restricted (35-37 deg E) household
#   bootstrap pool used for pseudo-control community bootstrapping
# - produces the agroecological zone map of the study area, showing
#   treated communities, used/excluded control households, the Rift
#   Valley boundary, and a zoomed inset around the study area
#
# OUTPUTS GENERATED
# ----------------------
#   hh_boot_35_37E.csv / hh_boot_35_37E.rds   bootstrap-ready control
#       household pool restricted to 35-37 deg E (used by the downstream
#       pooled bootstrap script).
#
#   agroecological_map_with_zoom_inset.png    study area map: agroecological
#       zones, treated communities, used/excluded control households,
#       Rift Valley boundary, and zoomed inset.
#
# ============================================================


# ----------------------------
# 0. Set working directory
# ----------------------------
setwd("~/Papers/PACS/Control_households/Clean_code")


# ----------------------------
# 1. Load packages
# ----------------------------
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)  # required by ne_countries() at scale = "medium"
library(terra)
library(tibble)
library(ggspatial)
library(patchwork)


# ----------------------------
# 2. Load control household data
# ----------------------------
# T1_clean / T2_clean are semicolon-separated landrace lists, split
# into list-columns here for use when building the bootstrap pool below.
control_hh_data <- read_csv("clean_control_hh_data.csv") %>%
  mutate(
    T1_clean = map(strsplit(T1_clean, "; "), ~ if (length(.x) == 1 && is.na(.x)) character(0) else .x),
    T2_clean = map(strsplit(T2_clean, "; "), ~ if (length(.x) == 1 && is.na(.x)) character(0) else .x)
  )

# ----------------------------
# 3. Load spatial layers
# ----------------------------
# `hh` is a subset of `control_hh_data` (hh_id, latitude, longitude,
# T1_Year), used by the bootstrap and mapping sections below.
ethiopia <- rnaturalearth::ne_countries(
  country = "Ethiopia",
  scale = "medium",
  returnclass = "sf"
) %>%
  st_make_valid()

# Helper: read a boundary shapefile/KML and harmonise its CRS to Ethiopia's
load_boundary <- function(path, ...) {
  st_read(path, ...) %>%
    st_transform(st_crs(ethiopia))
}

surveyed_kebeles <- load_boundary("surveyed_kebeles_PACS.shp", quiet = TRUE)

hh <- subset(control_hh_data, select = c(hh_id, latitude, longitude, T1_Year))


# ----------------------------
# 4. Restrict to households between 35 and 37 degrees East
# ----------------------------
# Matches the treated study area's longitude range; this restricted
# pool becomes the sampling pool for the downstream bootstrap analysis.
hh_restricted <- hh %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    longitude >= 35,
    longitude <= 37
  )


# ----------------------------
# 5. Create household-level bootstrap pool
# ----------------------------
# This is the pool of control households from which pseudo-villages
# will be sampled with replacement.
hh_boot <- control_hh_data %>%
  filter(hh_id %in% hh_restricted$hh_id) %>%
  select(
    hh_id,
    T1_Year,
    T1_clean,
    T2_clean,
    n_T1,
    n_T2,
    n_lost_row,
    n_new_row,
  ) %>%
  mutate(
    T1_clean = map(T1_clean, ~ str_trim(.x)),
    T2_clean = map(T2_clean, ~ str_trim(.x))
  )


# ----------------------------
# 6. Save readable and R-native versions of bootstrap pool
# ----------------------------
hh_boot_save <- hh_boot %>%
  mutate(
    T1_clean = map_chr(T1_clean, ~ paste(.x, collapse = ", ")),
    T2_clean = map_chr(T2_clean, ~ paste(.x, collapse = ", ")),
    n_lost_row = map_chr(n_lost_row, ~ paste(.x, collapse = ", ")),
    n_new_row  = map_chr(n_new_row, ~ paste(.x, collapse = ", "))
  )

write.csv(hh_boot_save, "hh_boot_35_37E.csv", row.names = FALSE)
saveRDS(hh_boot, "hh_boot_35_37E.rds")


# ----------------------------
# 7. Agroecological lookup table
# ----------------------------
# Sections 7-14 build the agroecological zone map: agroecological
# zones, controls, treated communities, Rift Valley boundary, and a
# zoomed inset.
agro_lookup <- tribble(
  ~CLASS, ~TEXT,
  101, "Temperate / arid",
  102, "Temperate / semi-arid",
  103, "Temperate / sub-humid",
  104, "Temperate / humid",
  211, "Subtropic - warm / arid",
  212, "Subtropic - warm / semi-arid",
  213, "Subtropic - warm / sub-humid",
  214, "Subtropic - warm / humid",
  221, "Subtropic - cool / arid",
  222, "Subtropic - cool / semi-arid",
  223, "Subtropic - cool / sub-humid",
  224, "Subtropic - cool / humid",
  311, "Tropic - warm / arid",
  312, "Tropic - warm / semi-arid",
  313, "Tropic - warm / sub-humid",
  314, "Tropic - warm / humid",
  321, "Tropic - cool / arid",
  322, "Tropic - cool / semi-arid",
  323, "Tropic - cool / sub-humid",
  324, "Tropic - cool / humid",
  400, "Boreal"
)


# ----------------------------
# 8. Load and prepare agroecological raster
# ----------------------------
agro <- rast("agroecological_zones.tif")
agro <- project(agro, st_crs(ethiopia)$wkt)

ethiopia_vect <- vect(ethiopia)

agro_clip <- agro %>%
  crop(ethiopia_vect) %>%
  mask(ethiopia_vect)

agro_df <- as.data.frame(agro_clip, xy = TRUE, na.rm = TRUE)
names(agro_df)[3] <- "CLASS"

agro_df <- agro_df %>%
  mutate(CLASS = as.integer(CLASS)) %>%
  left_join(agro_lookup, by = "CLASS") %>%
  filter(!is.na(TEXT))


# ----------------------------
# 9. Shared overlays used by BOTH maps below
# ----------------------------

treated_centroids <- surveyed_kebeles %>%
  filter(SELECTED == "Yes") %>%
  st_point_on_surface() %>%
  mutate(map_status = "Treated communities (n = 15)")

controls_all <- hh %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(st_crs(ethiopia)) %>%
  mutate(
    map_status = if_else(
      hh_id %in% hh_boot_save$hh_id,
      "Used controls (n = 29)",
      "Excluded controls (n = 23)"
    )
  )

controls_all_jitter <- st_jitter(controls_all, amount = 0.005)

rift_eth <- load_boundary("Rift_valley_boundary.kml") %>%
  st_make_valid() %>%
  st_intersection(ethiopia) %>%
  mutate(line_feature = "Rift Valley")


# ----------------------------
# 10. Define zoom area
# ----------------------------
zoom_geom <- c(
  st_geometry(controls_all %>% filter(map_status == "Used controls (n = 29)")),
  st_geometry(treated_centroids)
)

zoom_bb <- st_bbox(zoom_geom)
pad <- 0.25

zoom_x <- as.numeric(c(zoom_bb["xmin"] - pad, zoom_bb["xmax"] + pad))
zoom_y <- as.numeric(c(zoom_bb["ymin"] - pad, zoom_bb["ymax"] + pad))

zoom_box <- st_as_sfc(
  st_bbox(
    c(
      xmin = zoom_x[1],
      ymin = zoom_y[1],
      xmax = zoom_x[2],
      ymax = zoom_y[2]
    ),
    crs = st_crs(ethiopia)
  )
) %>%
  st_sf() %>%
  mutate(line_feature = "Study area")


# ----------------------------
# 11. Colour and shape palettes
# ----------------------------
agro_cols <- c(
  "Tropic - cool / arid"       = "#6B7FA6",
  "Tropic - cool / semi-arid"  = "#8FB6C9",
  "Tropic - cool / sub-humid"  = "#9BBF9B",
  "Tropic - cool / humid"      = "#B7D7C2",
  "Tropic - warm / arid"       = "#D8B365",
  "Tropic - warm / semi-arid"  = "#E6D39A",
  "Tropic - warm / sub-humid"  = "#C9895B",
  "Tropic - warm / humid"      = "#BFA2C8"
)

agro_df <- agro_df %>%
  mutate(TEXT = factor(TEXT, levels = names(agro_cols)))

point_cols <- c(
  "Excluded controls (n = 23)" = "grey25",
  "Used controls (n = 29)" = "black",
  "Treated communities (n = 15)" = "green4"
)

point_shapes <- c(
  "Excluded controls (n = 23)" = 16,
  "Used controls (n = 29)" = 16,
  "Treated communities (n = 15)" = 17
)

line_types <- c(
  "Rift Valley" = "dashed",
  "Study area" = "solid"
)


# ----------------------------
# 12. Full map
# ----------------------------
p_full <- ggplot() +
  geom_raster(
    data = agro_df,
    aes(x = x, y = y, fill = TEXT),
    alpha = 0.85
  ) +
  geom_sf(
    data = rift_eth,
    aes(linetype = line_feature),
    colour = "black",
    linewidth = 0.35,
    alpha = 0.35
  ) +
  geom_sf(
    data = zoom_box,
    aes(linetype = line_feature),
    fill = NA,
    colour = "black",
    linewidth = 0.75
  ) +
  geom_sf(
    data = ethiopia,
    fill = NA,
    colour = "black",
    linewidth = 0.45
  ) +
  geom_sf(
    data = controls_all_jitter,
    aes(colour = map_status, shape = map_status),
    size = 1.5,
    alpha = 0.9
  ) +
  geom_sf(
    data = treated_centroids,
    aes(colour = map_status, shape = map_status),
    size = 1.8,
    alpha = 0.95
  ) +
  annotation_scale(
    location = "br",
    style = "ticks",
    width_hint = 0.22,
    text_cex = 0.65,
    line_width = 0.45,
    pad_x = unit(0.25, "cm"),
    pad_y = unit(0.6, "cm")
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_orienteering,
    height = unit(0.9, "cm"),
    width = unit(0.9, "cm"),
    pad_x = unit(0.25, "cm"),
    pad_y = unit(0.25, "cm")
  ) +
  scale_fill_manual(
    values = agro_cols,
    name = "Agroecological zone",
    drop = TRUE
  ) +
  scale_colour_manual(
    values = point_cols,
    name = "Households and communities"
  ) +
  scale_shape_manual(
    values = point_shapes,
    name = "Households and communities"
  ) +
  scale_linetype_manual(
    values = line_types,
    name = "Features"
  ) +
  coord_sf(expand = TRUE) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1),
      order = 1
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        colour = "black",
        linewidth = c(0.35, 0.75)
      )
    ),
    colour = guide_legend(
      order = 3,
      override.aes = list(size = 3, alpha = 1)
    ),
    shape = guide_legend(order = 3)
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.6, "cm"),
    axis.title = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )


# ----------------------------
# 13. Zoomed inset map
# ----------------------------
agro_zoom_df <- agro_df %>%
  filter(
    x >= zoom_x[1],
    x <= zoom_x[2],
    y >= zoom_y[1],
    y <= zoom_y[2]
  ) %>%
  droplevels()

agro_cols_zoom <- agro_cols[levels(droplevels(agro_zoom_df$TEXT))]

p_inset <- ggplot() +
  geom_raster(
    data = agro_zoom_df,
    aes(x = x, y = y, fill = TEXT),
    alpha = 0.85
  ) +
  geom_sf(
    data = rift_eth,
    colour = "black",
    linewidth = 0.3,
    linetype = "dashed",
    alpha = 0.25
  ) +
  geom_sf(
    data = ethiopia,
    fill = NA,
    colour = "black",
    linewidth = 0.4
  ) +
  geom_sf(
    data = controls_all_jitter,
    aes(colour = map_status, shape = map_status),
    size = 1,
    alpha = 0.9
  ) +
  geom_sf(
    data = treated_centroids,
    aes(colour = map_status, shape = map_status),
    size = 1,
    alpha = 0.95
  ) +
  coord_sf(
    xlim = zoom_x,
    ylim = zoom_y,
    expand = FALSE
  ) +
  scale_fill_manual(
    values = agro_cols_zoom,
    name = "Agroecological zone",
    drop = TRUE
  ) +
  scale_colour_manual(
    values = point_cols,
    name = "Households and communities"
  ) +
  scale_shape_manual(
    values = point_shapes,
    name = "Households and communities"
  ) +
  scale_x_continuous(
    breaks = pretty(zoom_x, n = 3),
    labels = function(x) paste0(x, "\u00b0E")
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(zoom_y[1] * 2) / 2,
      ceiling(zoom_y[2] * 2) / 2,
      by = 0.25
    ),
    labels = function(y) paste0(y, "\u00b0N")
  ) +
  annotation_scale(
    location = "tr",
    style = "ticks",
    width_hint = 0.25,
    text_cex = 0.6,
    line_width = 0.4,
    pad_x = unit(0.25, "cm"),
    pad_y = unit(0.15, "cm")
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA)
  )


# ----------------------------
# 14. Combine and save
# ----------------------------
p_combined <- p_full +
  inset_element(
    p_inset,
    left = 0.61,
    bottom = 0.6,
    right = 0.98,
    top = 0.98
  )

p_combined

ggsave(
  "agroecological_map_with_zoom_inset.png",
  p_combined,
  width = 9,
  height = 6,
  dpi = 300
)

