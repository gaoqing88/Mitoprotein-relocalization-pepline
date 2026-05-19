#!/usr/bin/env Rscript
# Publication-grade FC_Mito vs FC_Cyto scatter plots.
# Three colorings × two layouts = 6 figures per results folder.
#
# Colorings:
#   * Direction      — colored by Direction (only Stable_Dual)
#   * ChangePattern  — colored by Change_Pattern (only Stable_Dual, 9 colors)
#   * Pattern        — colored by Pattern (ALL proteins, 6 colors)
#
# Layouts:
#   A) Multi-panel (one per tissue), 3x3
#   B) Single combined panel, all tissues overlaid
#
# Output:
#   <results_dir>/scatter_redesign/Figure_4_<Coloring>_<A|B>.pdf
#
# Usage:
#   Rscript 05_scatter_redesign.R --results=results_GeneName
#   Rscript 05_scatter_redesign.R --results=results_proteinID

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", all_args[grep("--file=", all_args)])
script_dir <- if (length(script_path) > 0) dirname(script_path) else getwd()
project_dir <- if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir

results_name <- "results_GeneName"
for (arg in args) {
  if (grepl("^--results=", arg)) results_name <- sub("^--results=", "", arg)
}
results_dir <- file.path(project_dir, results_name)
out_dir <- file.path(results_dir, "figure4_redesign")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n")

# ---------------------------------------------------------------------------
# Tissue order
# ---------------------------------------------------------------------------
# (Loaded from palette.R below)

# ---------------------------------------------------------------------------
# Source shared palette (defines tissue_order, tissue_colors, pattern_*, change_*, ALPHA_*)
# ---------------------------------------------------------------------------
source(file.path(script_dir, "palette.R"))
change_levels  <- names(change_colors)
pattern_levels <- names(pattern_colors)

# ---------------------------------------------------------------------------
# Publication-grade theme
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 18) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 6, face = "bold", hjust = 0.5,
                                margin = margin(b = 12)),
      plot.subtitle = element_text(size = base_size + 2, hjust = 0.5, color = "grey25",
                                   margin = margin(b = 16)),
      axis.title.x = element_text(size = base_size + 4, face = "bold", color = "black",
                                  margin = margin(t = 12)),
      axis.title.y = element_text(size = base_size + 4, face = "bold", color = "black",
                                  margin = margin(r = 12)),
      axis.text  = element_text(size = base_size + 2, color = "black"),
      axis.line  = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(5, "pt"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.1),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(1.2, "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 5, face = "bold", color = "black",
                                margin = margin(4, 4, 10, 4)),
      legend.title = element_text(size = base_size + 4, face = "bold"),
      legend.text  = element_text(size = base_size + 3),
      legend.key   = element_blank(),
      legend.key.size = unit(1.3, "lines"),
      legend.spacing.y = unit(0.4, "cm"),
      legend.margin = margin(4, 12, 4, 12),
      legend.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(16, 16, 16, 16)
    )
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df <- read_excel(file.path(results_dir, "relocalization_analysis.xlsx"), sheet = "All_Proteins")
df <- df %>% filter(Tissue %in% tissue_order) %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
clip_outliers <- function(d, ax_max) {
  d %>% mutate(
    FC_Cyto_plot = pmin(pmax(FC_Cyto, -ax_max), ax_max),
    FC_Mito_plot = pmin(pmax(FC_Mito, -ax_max), ax_max),
    is_outlier   = (FC_Cyto != FC_Cyto_plot) | (FC_Mito != FC_Mito_plot)
  )
}

add_jitter <- function(d, ax_max, seed = 42) {
  set.seed(seed)
  jw <- ax_max * 0.04   # ~4% of axis range — stronger spread
  d %>% mutate(
    FC_Cyto_jit = FC_Cyto_plot + runif(n(), -jw, jw),
    FC_Mito_jit = FC_Mito_plot + runif(n(), -jw, jw)
  )
}

# Compute clean axis range from data (capped to avoid cluttered breaks)
get_ax_lim <- function(d, q = 0.99, min_max = 3, max_max = 5) {
  fc_all <- c(d$FC_Mito, d$FC_Cyto)
  q99 <- quantile(abs(fc_all), q, na.rm = TRUE)
  ax_max <- ceiling(q99 * 1.3)
  ax_max <- max(ax_max, min_max)
  ax_max <- min(ax_max, max_max)
  list(lim = c(-ax_max, ax_max), max = ax_max,
       breaks = seq(-ax_max, ax_max, by = 1))
}

# ---------------------------------------------------------------------------
# Universal plotting function
# ---------------------------------------------------------------------------
# coloring: "Direction", "ChangePattern", "Pattern"
# layout:   "A" (per-tissue facets) or "B" (combined)
make_scatter <- function(data, coloring, layout, ax_lim, ax_max, ax_breaks) {

  # Resolve color column, palette, labels
  if (coloring == "Direction") {
    data <- data %>% filter(Pattern == "Stable_Dual") %>%
      mutate(
        is_sig = Significant == "Yes",
        ColorCat = case_when(
          is_sig & Direction == "Enriched_in_Mito" ~ "Enriched_in_Mito",
          is_sig & Direction == "Enriched_in_Cyto" ~ "Enriched_in_Cyto",
          TRUE ~ "NotSig"
        )
      )
    pal <- direction_colors
    lab <- direction_labels
    legend_title <- "Significant"
    has_gray_bg <- TRUE
  } else if (coloring == "ChangePattern") {
    data <- data %>% filter(Pattern == "Stable_Dual") %>%
      mutate(ColorCat = factor(Change_Pattern, levels = change_levels))
    pal <- change_colors
    lab <- change_labels
    legend_title <- "Change pattern"
    has_gray_bg <- FALSE
  } else if (coloring == "Pattern") {
    data <- data %>% mutate(ColorCat = factor(Pattern, levels = pattern_levels))
    pal <- pattern_colors
    lab <- pattern_labels
    legend_title <- "Pattern"
    has_gray_bg <- FALSE
  } else {
    stop("Unknown coloring: ", coloring)
  }

  # Clip and jitter
  data <- clip_outliers(data, ax_max) %>% add_jitter(ax_max)

  # Per-tissue stats for annotation (layout A)
  if (coloring == "Direction") {
    label_df <- data %>%
      group_by(Tissue) %>%
      summarise(
        n = n(),
        n_mito = sum(ColorCat == "Enriched_in_Mito"),
        n_cyto = sum(ColorCat == "Enriched_in_Cyto"),
        S_pct = round((n_mito + n_cyto) / n * 100, 1),
        .groups = "drop"
      ) %>% mutate(label = paste0("Mito: ", n_mito, "\nCyto: ", n_cyto,
                                   "\nn = ", n, " (S = ", S_pct, "%)"))
    combined_label <- with(label_df %>% summarise(
      n = sum(n), n_mito = sum(n_mito), n_cyto = sum(n_cyto),
      S_pct = round((n_mito + n_cyto) / n * 100, 1)
    ), paste0("S = ", S_pct, "%\nn = ", n, "\nMito: ", n_mito, "\nCyto: ", n_cyto))
  } else {
    label_df <- data %>%
      group_by(Tissue) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(label = paste0("n = ", n))
    combined_label <- paste0("n = ", nrow(data))
  }

  # Plot order: NotSig (gray) first if applicable, then colored on top
  if (has_gray_bg) {
    data_bg <- data %>% filter(ColorCat == "NotSig")
    data_fg <- data %>% filter(ColorCat != "NotSig")
  } else {
    # For ChangePattern and Pattern, draw "Stable" / "Stable_Dual" first (background)
    bg_cat <- if (coloring == "ChangePattern") "Stable" else "Stable_Dual"
    data_bg <- data %>% filter(as.character(ColorCat) == bg_cat)
    data_fg <- data %>% filter(as.character(ColorCat) != bg_cat)
  }

  # Base plot
  p <- ggplot(mapping = aes(x = FC_Cyto_jit, y = FC_Mito_jit)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.55) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.55) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.5)

  # Background (gray for Direction, light category color for others)
  if (has_gray_bg) {
    p <- p +
      geom_point(data = data_bg %>% filter(!is_outlier),
                 color = "grey78", size = 1.5, alpha = ALPHA_BACKGROUND, stroke = 0,
                 show.legend = FALSE) +
      geom_point(data = data_bg %>% filter(is_outlier),
                 color = "grey60", size = 1.8, alpha = ALPHA_BACKGROUND + 0.15, shape = 17, stroke = 0,
                 show.legend = FALSE)
  } else {
    # Use aes(color=...) so the bg category also appears in the legend
    p <- p +
      geom_point(data = data_bg %>% filter(!is_outlier),
                 aes(color = ColorCat), size = 1.4, alpha = ALPHA_BACKGROUND, stroke = 0) +
      geom_point(data = data_bg %>% filter(is_outlier),
                 aes(color = ColorCat), size = 1.7, alpha = ALPHA_BACKGROUND + 0.15, shape = 17, stroke = 0)
  }

  # Foreground (colored points)
  p <- p +
    geom_point(data = data_fg %>% filter(!is_outlier),
               aes(color = ColorCat), size = 2.3, alpha = ALPHA_SECONDARY, stroke = 0) +
    geom_point(data = data_fg %>% filter(is_outlier),
               aes(color = ColorCat), size = 2.8, alpha = ALPHA_SECONDARY + 0.10, shape = 17, stroke = 0)

  # Color scale
  if (has_gray_bg) {
    # Direction: only show two real categories in legend
    p <- p + scale_color_manual(values = pal, labels = lab, name = legend_title,
                                breaks = names(pal),
                                drop = FALSE)
  } else {
    p <- p + scale_color_manual(values = pal, labels = lab, name = legend_title,
                                breaks = names(pal),
                                drop = FALSE)
  }

  p <- p +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    coord_fixed(ratio = 1, xlim = ax_lim, ylim = ax_lim) +
    labs(
      x = expression(bold("FC cytoplasm  (log"[2]*" 20m / 3m)")),
      y = expression(bold("FC mitochondria  (log"[2]*" 20m / 3m)"))
    )

  # Title depends on coloring (different data scope)
  plot_title <- if (coloring == "Pattern") {
    "Protein relocalization across tissues"
  } else {
    "Stable dual-localized proteins"
  }

  # Layout
  if (layout == "A") {
    # Tissue name annotation at top-right of each panel
    tissue_label_df <- tibble::tibble(
      Tissue = factor(tissue_order, levels = tissue_order),
      x = ax_lim[2] - 0.25,
      y = ax_lim[2] - 0.25,
      label = as.character(Tissue)
    )

    p <- p + facet_wrap(~Tissue, ncol = 3) +
      # n / Mito / Cyto stats at top-left
      geom_text(data = label_df,
                aes(x = ax_lim[1] + 0.25, y = ax_lim[2] - 0.25, label = label),
                inherit.aes = FALSE, hjust = 0, vjust = 1,
                size = 6.5, fontface = "bold", lineheight = 1.0) +
      # Tissue name at top-right (inside the panel)
      geom_text(data = tissue_label_df,
                aes(x = x, y = y, label = label),
                inherit.aes = FALSE, hjust = 1, vjust = 1,
                size = 8, fontface = "bold", color = "black") +
      labs(title = plot_title) +
      guides(color = guide_legend(override.aes = list(size = 7, alpha = 1), ncol = 1)) +
      theme_pub(base_size = 18) +
      theme(legend.position = "right",
            legend.box.margin = margin(l = 10),
            strip.text = element_blank(),
            strip.background = element_blank())
  } else {
    p <- p +
      annotate("text", x = ax_lim[1] + 0.25, y = ax_lim[2] - 0.25,
               label = combined_label, hjust = 0, vjust = 1,
               size = 7.5, fontface = "bold", lineheight = 1.0) +
      labs(title = plot_title) +
      guides(color = guide_legend(override.aes = list(size = 7, alpha = 1), ncol = 1)) +
      theme_pub(base_size = 20) +
      theme(legend.position = "right", legend.box.margin = margin(l = 10))
  }

  return(p)
}

# ---------------------------------------------------------------------------
# Inverse layout: facet by Pattern/ChangePattern, color by Tissue
# ---------------------------------------------------------------------------
# Each panel shows one pattern type, points colored by tissue.
# Useful to see "which tissues drive each pattern".
make_scatter_byPattern <- function(data, facet_col, ax_lim, ax_max, ax_breaks) {
  # facet_col: "Pattern" or "Change_Pattern"
  if (facet_col == "Pattern") {
    facet_levels <- pattern_levels
    facet_labels <- pattern_labels
    plot_title   <- "Protein relocalization across tissues"
  } else if (facet_col == "Change_Pattern") {
    data <- data %>% filter(Pattern == "Stable_Dual")
    facet_levels <- change_levels
    facet_labels <- change_labels
    plot_title   <- "Stable dual-localized proteins"
  } else {
    stop("Unknown facet_col: ", facet_col)
  }

  # Drop any pattern not present in the data
  facet_levels <- intersect(facet_levels, unique(as.character(data[[facet_col]])))
  data <- data %>%
    filter(.data[[facet_col]] %in% facet_levels) %>%
    mutate(FacetCat = factor(.data[[facet_col]], levels = facet_levels))

  # Clip and jitter
  data <- clip_outliers(data, ax_max) %>% add_jitter(ax_max)

  # Per-facet stats (n, plus tissue breakdown later if needed)
  label_df <- data %>%
    group_by(FacetCat) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0("n = ", n))

  # Facet label = readable name
  facet_label_df <- tibble::tibble(
    FacetCat = factor(facet_levels, levels = facet_levels),
    x = ax_lim[2] - 0.25,
    y = ax_lim[2] - 0.25,
    label = unname(facet_labels[facet_levels])
  )

  # Panel count to choose ncol
  n_panels <- length(facet_levels)
  ncol_use <- if (n_panels <= 4) n_panels else if (n_panels <= 6) 3 else 3

  p <- ggplot(data, aes(x = FC_Cyto_jit, y = FC_Mito_jit)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.55) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.55) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.5) +
    geom_point(data = . %>% filter(!is_outlier),
               aes(color = Tissue), size = 2.3, alpha = ALPHA_SECONDARY, stroke = 0) +
    geom_point(data = . %>% filter(is_outlier),
               aes(color = Tissue), size = 2.8, alpha = ALPHA_SECONDARY + 0.10, shape = 17, stroke = 0) +
    geom_text(data = label_df,
              aes(x = ax_lim[1] + 0.25, y = ax_lim[2] - 0.25, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1,
              size = 6.5, fontface = "bold", lineheight = 1.0) +
    geom_text(data = facet_label_df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1, vjust = 1,
              size = 8, fontface = "bold", color = "black") +
    scale_color_manual(values = tissue_colors, name = "Tissue",
                       breaks = tissue_order, drop = FALSE) +
    scale_x_continuous(breaks = ax_breaks) +
    scale_y_continuous(breaks = ax_breaks) +
    coord_fixed(ratio = 1, xlim = ax_lim, ylim = ax_lim) +
    labs(
      title = plot_title,
      x = expression(bold("FC cytoplasm  (log"[2]*" 20m / 3m)")),
      y = expression(bold("FC mitochondria  (log"[2]*" 20m / 3m)"))
    ) +
    facet_wrap(~FacetCat, ncol = ncol_use) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1), ncol = 1)) +
    theme_pub(base_size = 18) +
    theme(legend.position = "right",
          legend.box.margin = margin(l = 10),
          strip.text = element_blank(),
          strip.background = element_blank())

  return(p)
}

# ---------------------------------------------------------------------------
# Render all combinations
# ---------------------------------------------------------------------------
# Different colorings need different axis ranges:
#   - Direction & ChangePattern: only Stable_Dual, FC values are tame
#   - Pattern: all proteins, FC values for New/Lost can be huge → cap to ±5

ax_stable <- get_ax_lim(df %>% filter(Pattern == "Stable_Dual"), min_max = 3, max_max = 5)
ax_all    <- get_ax_lim(df, min_max = 4, max_max = 5)

configs <- list(
  list(coloring = "ChangePattern", ax = ax_stable, prefix = "ChangePattern"),
  list(coloring = "Pattern",       ax = ax_all,    prefix = "Pattern")
)

for (cfg in configs) {
  for (lay in c("A", "B")) {
    p <- make_scatter(df, cfg$coloring, lay,
                      ax_lim = cfg$ax$lim, ax_max = cfg$ax$max,
                      ax_breaks = cfg$ax$breaks)

    if (lay == "A") {
      w <- 16; h <- 14   # wider to fit right-side legend
    } else {
      w <- 13; h <- 11
    }
    fname <- sprintf("Figure_4_%s_%s.pdf", cfg$prefix, lay)
    ggsave(file.path(out_dir, fname), p, width = w, height = h,
           dpi = 300, useDingbats = FALSE)
    cat(sprintf("  Wrote %s  (%d x %d in)\n", fname, w, h))
  }
}

# ---------------------------------------------------------------------------
# Inverse layout: facet by Pattern, color by Tissue
# ---------------------------------------------------------------------------
inverse_configs <- list(
  list(facet_col = "Change_Pattern", ax = ax_stable, prefix = "ChangePattern"),
  list(facet_col = "Pattern",        ax = ax_all,    prefix = "Pattern")
)

for (cfg in inverse_configs) {
  p <- make_scatter_byPattern(df, cfg$facet_col,
                              ax_lim = cfg$ax$lim, ax_max = cfg$ax$max,
                              ax_breaks = cfg$ax$breaks)
  fname <- sprintf("Figure_4_%s_byTissue.pdf", cfg$prefix)
  # ChangePattern has 9 panels (3x3), Pattern has 6-7 panels (3x?)
  if (cfg$facet_col == "Change_Pattern") {
    w <- 16; h <- 14
  } else {
    w <- 16; h <- 11
  }
  ggsave(file.path(out_dir, fname), p, width = w, height = h,
         dpi = 300, useDingbats = FALSE)
  cat(sprintf("  Wrote %s  (%d x %d in)\n", fname, w, h))
}

cat("\nDone. All 6 figures saved to:", out_dir, "\n")
