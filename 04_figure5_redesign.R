#!/usr/bin/env Rscript
# Figure 5 redesign — two publication-grade views of the relocalization data.
#
#   1) Figure_5_ridge.pdf         — Ridge plot per tissue showing Reloc_Score distribution
#   2) Figure_5_heatmap.pdf       — Top-N proteins by |Reloc_Score| as a heatmap (gene x tissue)
#
# Output:
#   <results_dir>/figure5_redesign/Figure_5_*.pdf
#
# Usage:
#   Rscript 07_figure5_redesign.R --results=results_GeneName
#   Rscript 07_figure5_redesign.R --results=results_proteinID

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggridges)
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
out_dir <- file.path(results_dir, "figure5_redesign")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n")

# ---------------------------------------------------------------------------
# Tissue order and colors
# ---------------------------------------------------------------------------
# Source shared palette
source(file.path(script_dir, "palette.R"))

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 18) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 6, face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size + 1, hjust = 0.5, color = "grey25",
                                   margin = margin(b = 14)),
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
      legend.title = element_text(size = base_size + 3, face = "bold"),
      legend.text  = element_text(size = base_size + 2),
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

stable_df <- df %>%
  filter(Pattern == "Stable_Dual", Tissue %in% tissue_order) %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))

# ===========================================================================
# Figure 5.2 — Ridge plot of Reloc_Score per tissue
# ===========================================================================
# Use Reloc_Score (point estimate from collapsed averages) for the distribution shape.
ridge_df <- stable_df %>%
  mutate(Tissue = factor(Tissue, levels = rev(tissue_order)))   # reverse for top-to-bottom

ridge_lim_raw <- max(quantile(abs(ridge_df$Reloc_Score), 0.99, na.rm = TRUE) * 1.2, 3)
ridge_lim <- min(ceiling(ridge_lim_raw), 6)   # clean integer

p_ridge <- ggplot(ridge_df, aes(x = Reloc_Score, y = Tissue, fill = Tissue)) +
  geom_density_ridges(scale = 2.0, alpha = ALPHA_PRIMARY, color = "white",
                      linewidth = 0.4, rel_min_height = 0.005,
                      quantile_lines = TRUE, quantiles = 0.5,
                      vline_color = "white", vline_linewidth = 1) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = c(-1.5, 1.5), color = "grey45",
             linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = tissue_colors, guide = "none") +
  scale_x_continuous(limits = c(-ridge_lim, ridge_lim),
                     breaks = seq(-ridge_lim, ridge_lim, by = 1)) +
  labs(
    title = "Reloc score distribution by tissue",
    x = expression(bold("Reloc score  =  FC"["mito"] - "FC"["cyto"])),
    y = NULL
  ) +
  theme_pub(base_size = 18) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 20, face = "bold"))

ggsave(file.path(out_dir, "Figure_5_ridge.pdf"),
       p_ridge, width = 8.5, height = 8, dpi = 300, useDingbats = FALSE)
cat("  Wrote Figure_5_ridge.pdf  (8.5 x 8 in)\n")

# ===========================================================================
# Figure 5.3 — Top-N heatmap of |Reloc_Score| (gene × tissue)
# ===========================================================================
# Top-N most-relocalized genes per tissue. Keep N small enough that cells
# are roughly square in the final figure.
TOP_N_PER_TISSUE <- 5

gene_col <- if ("GeneName" %in% colnames(stable_df)) "GeneName" else "ProteinID"

top_genes <- stable_df %>%
  group_by(Tissue) %>%
  arrange(desc(abs(Reloc_Score))) %>%
  slice_head(n = TOP_N_PER_TISSUE) %>%
  ungroup() %>%
  pull(.data[[gene_col]]) %>%
  unique()

heatmap_df <- stable_df %>%
  filter(.data[[gene_col]] %in% top_genes) %>%
  select(all_of(gene_col), Tissue, Reloc_Score, Significant_ttest)

# Order genes by mean |Reloc_Score| across tissues (clustering substitute)
gene_order <- heatmap_df %>%
  group_by(.data[[gene_col]]) %>%
  summarise(score = mean(abs(Reloc_Score), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(score)) %>%
  pull(.data[[gene_col]])

heatmap_df <- heatmap_df %>%
  mutate(
    Gene = factor(.data[[gene_col]], levels = rev(gene_order)),
    Tissue = factor(Tissue, levels = tissue_order),
    score_capped = pmin(pmax(Reloc_Score, -5), 5),
    sig_mark = ifelse(Significant_ttest == "Yes", "*", "")
  )

n_genes <- length(top_genes)

p_heat <- ggplot(heatmap_df, aes(x = Tissue, y = Gene, fill = score_capped)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sig_mark), size = 9, fontface = "bold",
            color = "white", vjust = 0.7) +
  scale_fill_gradient2(low = heatmap_low, mid = "white", high = heatmap_high,
                       midpoint = 0, limits = c(-5, 5),
                       breaks = c(-4, -2, 0, 2, 4),
                       name = "Reloc score",
                       guide = guide_colorbar(barwidth = 2.0, barheight = 18,
                                              title.position = "top",
                                              title.hjust = 0.5,
                                              ticks.colour = "black",
                                              frame.colour = "black",
                                              frame.linewidth = 0.7)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = sprintf("Top %d most-relocalized proteins per tissue", TOP_N_PER_TISSUE),
    x = NULL, y = NULL
  ) +
  theme_pub(base_size = 18) +
  coord_equal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, face = "bold"),
    axis.text.y = element_text(face = "italic", size = 16),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5,
                              margin = margin(b = 14))
  )

# Cell width ~ 1.0 inch each; canvas auto-scales
heatmap_width  <- 13
heatmap_height <- 14
ggsave(file.path(out_dir, "Figure_5_heatmap.pdf"),
       p_heat, width = heatmap_width, height = heatmap_height,
       dpi = 300, useDingbats = FALSE)
cat(sprintf("  Wrote Figure_5_heatmap.pdf  (%d genes x %d tissues, %d x %d in)\n",
            n_genes, length(tissue_order), heatmap_width, heatmap_height))

cat("\nDone. Figure 5 redesigns saved to:", out_dir, "\n")
