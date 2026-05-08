################################################################################
#  Multi-Panel Maps of Study Sites
# Created: 2026-02-23
# Purpose: Two separate multi-panel figures:
#          Figure 1: Snowy Hydro (overview + per-cluster detail panels)
#          Figure 2: Hunter Mining (overview + per-cluster detail panels)
#
# Layout per figure (matching reference image style):
#   TOP:    Overview with convex hull boundaries, numbered cluster markers,
#           Australia inset (top-left), legend, compass arrow
#   BOTTOM: One panel per spatial cluster showing individual plots as
#           grey/forestgreen triangles with auto-calculated scale bars
#
# ADJUSTABLE PARAMETERS:
#   h_km     - clustering distance threshold in km (lower = more clusters)
#   buffer_m - buffer around convex hulls in meters
# 
# NOTE: MUST RUN ANALYSIS SCRIPT FIRST TO GENERATE CLEAN PLOT DATA FOR CASE STUDIES
################################################################################

library(tidyverse)
library(sf)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(grid)

# Fix potential conflict with recode
recode <- dplyr::recode

# Color scheme for publication
color_scheme <- c(Reference = "forestgreen", Rehab = "grey60")

################################################################################
# LOAD HUNTER DATA
################################################################################

proj_root <- here::here()
derived_dir <- file.path(proj_root, "deriveddata")

# Output
plot_dir <- file.path(proj_root, "output/plots")

cat("── Loading Hunter Mining from cache ──\n")
hunter_dat <- readRDS(file.path(derived_dir, "hunter_dat.rds"))

hunter_plots <- hunter_dat %>%
    group_by(Plot_ID, Site_ID, Treat) %>%
    slice(1) %>%
    ungroup() %>%
    select(Plot_ID, Site_ID, Treat, geometry)

cat("\nHunter plots loaded:", nrow(hunter_plots), "plots\n")

################################################################################
# LOAD SNOWIES DATA
################################################################################

cat("── Loading Snowy Hydro from cache ──\n")
snowies_dat <- readRDS(file.path(derived_dir, "snowies_dat.rds"))

snowies_plots <- snowies_dat %>%
    group_by(Plot_ID, Site_ID, Treat) %>%
    slice(1) %>%
    ungroup() %>%
    select(Plot_ID, Site_ID, Treat, geometry)

cat("\nSnowies plots loaded:", nrow(snowies_plots), "plots\n")

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Spatial clustering via hierarchical clustering (complete linkage)
cluster_plots <- function(plots_sf, h_km = 5) {
    # Pairwise great-circle distances in km
    dist_mat <- st_distance(plots_sf)
    dist_km <- as.dist(as.matrix(dist_mat) / 1000)

    # Hierarchical clustering with complete linkage
    hc <- hclust(dist_km, method = "complete")
    plots_sf$Cluster <- cutree(hc, h = h_km)

    # Reorder clusters north-to-south, west-to-east
    order_df <- plots_sf %>%
        group_by(Cluster) %>%
        summarise(geometry = st_union(geometry)) %>%
        st_centroid() %>%
        cbind(st_coordinates(.)) %>%
        st_drop_geometry() %>%
        arrange(desc(Y), X)

    remap <- setNames(seq_len(nrow(order_df)), order_df$Cluster)
    plots_sf$Cluster <- remap[as.character(plots_sf$Cluster)]
    return(plots_sf)
}

################################################################################
# AUSTRALIA INSET
################################################################################

australia <- ne_countries(country = "australia", returnclass = "sf", scale = "medium")
australia_crop <- st_crop(australia, xmin = 112, ymin = -44, xmax = 154, ymax = -10)

make_inset <- function(centroid_sf) {
    ggplot() +
        geom_sf(
            data = australia_crop, fill = "grey95", color = "black",
            linewidth = 0.5
        ) +
        geom_sf(
            data = centroid_sf, color = "red", shape = 8, size = 3,
            stroke = 1.5
        ) +
        coord_sf(expand = FALSE) +
        theme_void() +
        theme(
            panel.border = element_blank(),
            plot.background = element_rect(fill = NA, color = NA),
            panel.background = element_rect(fill = NA, color = NA),
            plot.margin = margin(1, 1, 1, 1)
        )
}

################################################################################
# MAIN FIGURE BUILDER
################################################################################

build_figure <- function(plots_sf, h_km = 5, buffer_m = 800) {
    tmap_mode("plot")

    # ---- 1. Cluster plots ----
    plots_c <- cluster_plots(plots_sf, h_km = h_km)
    n_clust <- max(plots_c$Cluster)
    cat("  Clusters found:", n_clust, "\n")

    # Print cluster composition
    plots_c %>%
        st_drop_geometry() %>%
        count(Cluster, Site_ID) %>%
        print()

    # ---- 2. Convex hull boundaries with buffer ----
    boundaries <- plots_c %>%
        group_by(Cluster) %>%
        summarise(geometry = st_union(geometry)) %>%
        st_convex_hull() %>%
        st_buffer(dist = buffer_m)

    # ---- 3. Cluster centroids for numbered markers ----
    centroids <- plots_c %>%
        group_by(Cluster) %>%
        summarise(geometry = st_union(geometry)) %>%
        st_centroid() %>%
        mutate(label = as.character(Cluster))

    # ---- 4. Study centroid for Australia inset ----
    study_centroid <- plots_c %>%
        summarise(geometry = st_union(geometry)) %>%
        st_centroid()

    # ---- 5. Overview bounding box (square, padded) ----
    ov_bb <- st_bbox(boundaries)
    ov_cx <- mean(c(ov_bb$xmin, ov_bb$xmax))
    ov_cy <- mean(c(ov_bb$ymin, ov_bb$ymax))
    ov_span <- max(
        ov_bb$xmax - ov_bb$xmin,
        ov_bb$ymax - ov_bb$ymin
    ) * 1.8
    ov_bbox <- st_bbox(c(
        xmin = ov_cx - ov_span / 2, xmax = ov_cx + ov_span / 2,
        ymin = ov_cy - ov_span / 2, ymax = ov_cy + ov_span / 2
    ), crs = st_crs(plots_c))

    # Map_type
    map_type <- "Esri.WorldImagery"
    # Esri.WorldTopoMap
    # Esri.NatGeoWorldMap
    # Esri.WorldStreetMap
    # Esri.WorldImagery


    # ---- 6. Overview tmap ----
    # Create the overview map
    map_ov <- tm_basemap(map_type) +
        # Numbered cluster markers only (no individual plots, no boundaries)
        tm_shape(centroids, bbox = ov_bbox) +
        tm_symbols(
            size = 1.8, shape = 21,
            fill = "white", col = "black",
            lwd = 1.5, fill_alpha = 0.9
        ) +
        tm_text("label", size = 1.0, fontface = "bold") +
        # Manual legend for Treatment
        tm_add_legend(
            type = "symbols",
            labels = c("Reference", "Restoration"),
            fill = unname(color_scheme),
            shape = 24,
            col = "white",
            size = 1,
            title = "Treatment"
        ) +
        # Compass
        tm_compass(
            type = "arrow",
            position = c("left", "bottom"),
            size = 2.5, text.size = 1, text.color = "white"
        ) +
        tm_scalebar(
            position = c("left", "bottom"),
            text.size = 1, text.color = "white"
        ) +
        tm_layout(
            legend.position = c("right", "bottom"),
            legend.title.size = 1,
            legend.text.size = .9,
            legend.bg.color = "white",
            legend.bg.alpha = 0.9,
            legend.frame = TRUE,
            frame = TRUE, frame.lwd = 2,
            inner.margins = 0
        )

    # ---- 7. Per-cluster detail maps ----
    detail_maps <- list()
    for (i in seq_len(n_clust)) {
        cl <- plots_c %>% filter(Cluster == i)

        cb <- st_bbox(cl)
        cx <- mean(c(cb$xmin, cb$xmax))
        cy <- mean(c(cb$ymin, cb$ymax))
        cr <- max(cb$xmax - cb$xmin, cb$ymax - cb$ymin)
        cr <- max(cr, 0.003) * 1.8 # minimum extent + 80% padding

        cl_bbox <- st_bbox(c(
            xmin = cx - cr / 2, xmax = cx + cr / 2,
            ymin = cy - cr / 2, ymax = cy + cr / 2
        ), crs = st_crs(cl))

        detail_maps[[i]] <- tm_basemap(map_type) +
            tm_shape(cl, bbox = cl_bbox) +
            tm_symbols(
                fill = "Treat",
                fill.scale = tm_scale(values = color_scheme),
                fill.legend = tm_legend(show = FALSE),
                size = .6, shape = 24,
                col = "#ffffff", lwd = 0.8, fill_alpha = 0.9
            ) +
            tm_scalebar(
                position = c("left", "bottom"),
                text.size = 0.8, text.color = "white"
            ) +
            tm_layout(
                frame = TRUE, frame.lwd = 2,
                inner.margins = 0
            )
    }

    # ---- 8. Convert to grobs and compose ----
    ov_grob <- tmap_grob(map_ov)
    inset_grob <- ggplotGrob(make_inset(study_centroid))

    # Overview panel with Australia inset overlay
    ov_panel <- ggplot() +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        annotation_custom(ov_grob, 0, 1, 0, 1) +
        annotation_custom(inset_grob,
            xmin = 0.01, xmax = 0.22,
            ymin = 0.72, ymax = 0.99
        ) +
        theme_void()

    # Detail panels with cluster number labels
    det_panels <- list()
    for (i in seq_len(n_clust)) {
        g <- tmap_grob(detail_maps[[i]])
        det_panels[[i]] <- ggplot() +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
            annotation_custom(g, 0, 1, 0, 1) +
            annotate("rect",
                xmin = 0.01, xmax = 0.1,
                ymin = 0.87, ymax = 0.99,
                fill = "white", color = "black", linewidth = 0.5
            ) +
            annotate("text",
                x = 0.055, y = 0.93,
                label = as.character(i),
                size = 7, fontface = "bold"
            ) +
            theme_void()
    }

    # ---- 9. Patchwork layout: overview on top, details on bottom ----
    top <- paste(rep("A", n_clust), collapse = "")
    bot <- paste(LETTERS[2:(n_clust + 1)], collapse = "")
    design <- paste(top, top, bot, sep = "\n")

    fig <- wrap_plots(c(list(ov_panel), det_panels)) +
        plot_layout(design = design) &
        theme(plot.margin = margin(3, 3, 3, 3))

    return(list(fig = fig, n_clust = n_clust))
}

################################################################################
# BUILD FIGURES
################################################################################

# --- Figure 1: Snowies ---
# h_km = 8 groups nearby sites; adjust up/down for fewer/more clusters
cat("\n=== Building Snowies Figure ===\n")
res_snowies <- build_figure(snowies_plots, h_km = 15, buffer_m = 400)

# --- Figure 2: Hunter ---
# h_km = 15 merges nearby mines; adjust up/down for fewer/more clusters
cat("\n=== Building Hunter Figure ===\n")
res_hunter <- build_figure(hunter_plots, h_km = 25, buffer_m = 500)

################################################################################
# SAVE FIGURES
################################################################################

# Snowies
ggsave(file.path(plot_dir, "Figure_Map_Snowies.png"),
    res_snowies$fig,
    width = 2 * res_snowies$n_clust, height = 8, dpi = 300
)
# ggsave(file.path(plot_dir, "Figure_Map_Snowies.pdf"),
#   res_snowies$fig,
#  width = 2 * res_snowies$n_clust, height = 10
# )

# Hunter
ggsave(file.path(plot_dir, "Figure_Map_Hunter.png"),
    res_hunter$fig,
    width = 2 * res_hunter$n_clust, height = 10, dpi = 300
)
# ggsave(file.path(plot_dir, "Figure_Map_Hunter.pdf"),
#   res_hunter$fig,
#  width = 2 * res_hunter$n_clust, height = 10
# )


