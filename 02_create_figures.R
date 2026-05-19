#!/usr/bin/env Rscript
# Create publication figures for protein relocalization analysis.
# Reads <results_dir>/relocalization_analysis.xlsx and the four data/*_merged.txt files,
# writes six PDF figures into <results_dir>/.
#
# Usage:
#   Rscript 02_create_figures.R                        # defaults to results_GeneName
#   Rscript 02_create_figures.R --results=results_proteinID

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(readxl)
  library(ggrepel)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", all_args[grep("--file=", all_args)])
script_dir <- if (length(script_path) > 0) dirname(script_path) else getwd()
project_dir <- if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir
data_dir <- file.path(project_dir, "data")

# Parse --results argument
results_name <- "results_GeneName"  # default
for (arg in args) {
  if (grepl("^--results=", arg)) {
    results_name <- sub("^--results=", "", arg)
  }
}
results_dir <- file.path(project_dir, results_name)
cat("Using results directory:", results_dir, "\n")

# ---------------------------------------------------------------------------
# Shared theme and palettes
# ---------------------------------------------------------------------------
theme_publication <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold", color = "black"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.key = element_blank(),
      panel.grid.major = element_line(color = "grey86", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 12, 12, 12)
    )
}

colors_direction <- c(
  "Enriched_in_Mito" = "#E64B35", "Enriched_in_Cyto" = "#4DBBD5",
  "Moderate" = "#7A7A7A", "Stable" = "#C9C9C9"
)

# Source shared palette
source(file.path(script_dir, "palette.R"))
colors_pattern <- pattern_colors
colors_overlap <- overlap_colors
colors_change  <- change_colors

# ---------------------------------------------------------------------------
# Load analysis table
# ---------------------------------------------------------------------------
df <- read_excel(file.path(results_dir, "relocalization_analysis.xlsx"), sheet = "All_Proteins")
tissue_counts <- df %>% count(Tissue, name = "Count") %>% arrange(desc(Count))
all_tissues <- c("Stomach", "Intestine", "Spleen", "Muscle", "Heart", "Brain", "Liver", "Lung", "Kidney")
major_tissues <- tissue_counts %>% filter(Count >= 20) %>% pull(Tissue)
top_tissues <- major_tissues

# ---------------------------------------------------------------------------
# Figure 2: overall relocalization pattern per tissue
# ---------------------------------------------------------------------------
pattern_dist <- df %>%
  group_by(Tissue, Pattern) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Tissue) %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100) %>%
  ungroup() %>%
  mutate(Tissue = factor(Tissue, levels = all_tissues))

pattern_totals <- pattern_dist %>% distinct(Tissue, Total)

missing_pattern_tissues <- tibble(
  Tissue = factor(setdiff(all_tissues, unique(as.character(pattern_dist$Tissue))), levels = all_tissues),
  y = 50,
  label = "No analyzable\nrecords"
)

p1 <- ggplot(pattern_dist, aes(x = Tissue, y = Percentage, fill = Pattern)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.35, width = 0.78,
           alpha = ALPHA_PRIMARY) +
  geom_text(
    aes(label = ifelse(Percentage >= 3, paste0(round(Percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 3.4, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(values = colors_pattern, labels = pattern_labels, name = "Pattern") +
  geom_text(data = pattern_totals, aes(x = Tissue, y = 103, label = paste0("n=", Total)),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  geom_text(data = missing_pattern_tissues, aes(x = Tissue, y = y, label = label),
            inherit.aes = FALSE, size = 3, color = "grey35", fontface = "bold") +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Figure 2. Overall Relocalization Pattern Across Tissues",
    subtitle = paste0("Key: ", results_name, " | Analyzable non-Other records from 3-month to 20-month"),
    x = "Tissue", y = "Percentage of proteins (%)"
  ) +
  theme_publication(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

ggsave(file.path(results_dir, "Figure_2_pattern_distribution.pdf"), p1, width = 11, height = 7, dpi = 300)

# ---------------------------------------------------------------------------
# Figure 1: mito/cyto proteome overlap at each age (uses raw merged files)
# ---------------------------------------------------------------------------
read_presence <- function(path, value_name) {
  # Collapse to ProteinID + Tissue level for overlap analysis.
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    pivot_longer(cols = -c(ProteinID, GeneName), names_to = "Sample", values_to = value_name) %>%
    mutate(
      Tissue = sub("_Rep[0-9]+$", "", Sample),
      Tissue = case_when(
        grepl("^liver", Tissue, ignore.case = TRUE) ~ "Liver",
        grepl("^lung", Tissue, ignore.case = TRUE) ~ "Lung",
        TRUE ~ paste0(toupper(substr(Tissue, 1, 1)), substr(Tissue, 2, nchar(Tissue)))
      ),
      Present = !is.na(.data[[value_name]]) & .data[[value_name]] > 0
    ) %>%
    group_by(ProteinID, Tissue) %>%
    summarise(Present = any(Present), .groups = "drop")
}

overlap_for_age <- function(age_label, mito_file, cyto_file) {
  mito <- read_presence(file.path(data_dir, mito_file), "Mito_Value") %>% rename(Mito_Present = Present)
  cyto <- read_presence(file.path(data_dir, cyto_file), "Cyto_Value") %>% rename(Cyto_Present = Present)
  full_join(mito, cyto, by = c("ProteinID", "Tissue")) %>%
    mutate(
      Mito_Present = replace_na(Mito_Present, FALSE),
      Cyto_Present = replace_na(Cyto_Present, FALSE),
      Category = case_when(
        Mito_Present & Cyto_Present ~ "Shared",
        Mito_Present & !Cyto_Present ~ "Mito_only",
        !Mito_Present & Cyto_Present ~ "Cyto_only",
        TRUE ~ NA_character_
      ),
      Age = age_label
    ) %>%
    filter(!is.na(Category)) %>%
    count(Age, Tissue, Category, name = "Count")
}

overlap_df <- bind_rows(
  overlap_for_age("3-month young", "3month_mito_merged.txt", "3month_cyto_merged.txt"),
  overlap_for_age("20-month old", "20month_mito_merged.txt", "20month_cyto_merged.txt")
) %>%
  mutate(
    Age = factor(Age, levels = c("3-month young", "20-month old")),
    Tissue = factor(Tissue, levels = all_tissues),
    Category = factor(Category, levels = c("Mito_only", "Shared", "Cyto_only"))
  )

overlap_pct <- overlap_df %>%
  group_by(Age, Tissue) %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100) %>%
  ungroup()

overlap_totals <- overlap_pct %>% distinct(Age, Tissue, Total)

p2 <- ggplot(overlap_pct, aes(x = Tissue, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25, width = 0.78,
           alpha = ALPHA_PRIMARY) +
  geom_text(
    aes(label = ifelse(Percentage >= 8, paste0(round(Percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 3, fontface = "bold", color = "white"
  ) +
  geom_text(data = overlap_totals, aes(x = Tissue, y = 103, label = paste0("n=", Total)),
            inherit.aes = FALSE, size = 2.8, fontface = "bold") +
  facet_wrap(~Age, ncol = 1) +
  scale_fill_manual(values = colors_overlap, labels = c("Mito only", "Shared", "Cyto only"),
                    name = "Protein ID class") +
  labs(
    title = "Figure 1. Mitochondrial and Cytoplasmic Proteome Overlap",
    subtitle = "ProteinID presence in each tissue (overlap by ProteinID + Tissue)",
    x = "Tissue", y = "Percentage of ProteinIDs (%)"
  ) +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0.01))) +
  theme_publication(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

ggsave(file.path(results_dir, "Figure_1_mito_cyto_overlap.pdf"), p2, width = 11, height = 9, dpi = 300)

# ---------------------------------------------------------------------------
# Figure 3: change-pattern distribution within Stable_Dual
# ---------------------------------------------------------------------------
stable_df <- df %>% filter(Pattern == "Stable_Dual", Tissue %in% major_tissues)

stable_change <- stable_df %>%
  count(Tissue, Change_Pattern, name = "Count") %>%
  group_by(Tissue) %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100) %>%
  ungroup() %>%
  mutate(
    Tissue = factor(Tissue, levels = major_tissues),
    Change_Pattern = factor(Change_Pattern, levels = c(
      "Cyto_to_Mito", "Mito_to_Cyto",
      "Cyto_up_and_Mito_up", "Cyto_down_and_Mito_down",
      "Mito_up", "Mito_down", "Cyto_up", "Cyto_down", "Stable"
    ))
  )

stable_totals <- stable_change %>% distinct(Tissue, Total)

p3 <- ggplot(stable_change, aes(x = Tissue, y = Percentage, fill = Change_Pattern)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.35, width = 0.78,
           alpha = ALPHA_PRIMARY) +
  geom_text(
    aes(label = ifelse(Percentage >= 6, paste0(round(Percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 3.3, fontface = "bold", color = "white"
  ) +
  geom_text(data = stable_totals, aes(x = Tissue, y = 103, label = paste0("n=", Total)),
            inherit.aes = FALSE, size = 3.3, fontface = "bold") +
  scale_fill_manual(values = colors_change, labels = change_labels, name = "Change pattern") +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Figure 3. Change Patterns Within Stable Dual-Localized Proteins",
    subtitle = paste0("Key: ", results_name, " | Stable_Dual only"),
    x = "Tissue", y = "Percentage of Stable_Dual proteins (%)"
  ) +
  theme_publication(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

ggsave(file.path(results_dir, "Figure_3_stable_dual_change_pattern.pdf"), p3, width = 10.5, height = 7, dpi = 300)

# ---------------------------------------------------------------------------
# Figure 4: Stable_Dual scatter of FC_Mito vs FC_Cyto per tissue
# ---------------------------------------------------------------------------
plot_list <- list()
for (tissue in top_tissues) {
  tissue_df <- stable_df %>% filter(Tissue == tissue)
  n_mito <- sum(tissue_df$Direction == "Enriched_in_Mito")
  n_cyto <- sum(tissue_df$Direction == "Enriched_in_Cyto")
  p <- ggplot(tissue_df, aes(x = FC_Cyto, y = FC_Mito, color = Direction)) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(aes(alpha = Direction), size = 1.9, stroke = 0) +
    scale_alpha_manual(values = c("Enriched_in_Mito" = 0.8, "Enriched_in_Cyto" = 0.8,
                                  "Moderate" = 0.35, "Stable" = 0.22), guide = "none") +
    scale_color_manual(values = colors_direction, guide = "none") +
    coord_fixed(ratio = 1, xlim = c(-6, 6), ylim = c(-6, 6), clip = "off") +
    annotate("text", x = -5.1, y = 5.35,
             label = paste0("Mito: ", n_mito, "\nCyto: ", n_cyto),
             hjust = 0, size = 3.3, fontface = "bold") +
    labs(title = tissue, x = "FC cytoplasm", y = "FC mitochondria") +
    theme_publication(base_size = 11) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  plot_list[[tissue]] <- p
}

p4 <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "Figure 4. Stable Dual-Localized Proteins: FC_Mito vs FC_Cyto",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(file.path(results_dir, "Figure_4_stable_dual_scatter.pdf"), p4, width = 14, height = 9.5, dpi = 300)

# ---------------------------------------------------------------------------
# Figure 5: Stable_Dual volcano (Reloc_Score vs -log10 P)
# ---------------------------------------------------------------------------
volcano_list <- list()
for (tissue in top_tissues) {
  tissue_df <- stable_df %>%
    filter(Tissue == tissue) %>%
    mutate(
      NegLog10P = pmin(-log10(pmax(P_value, 1e-16)), 6),
      Label = case_when(
        Direction == "Enriched_in_Mito" & rank(-Reloc_Score, ties.method = "first") <= 4 ~ GeneName,
        Direction == "Enriched_in_Cyto" & rank(Reloc_Score, ties.method = "first") <= 4 ~ GeneName,
        TRUE ~ ""
      )
    )
  p <- ggplot(tissue_df, aes(x = Reloc_Score, y = NegLog10P, color = Direction)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.55) +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey40", linewidth = 0.55) +
    geom_point(aes(alpha = Direction), size = 1.9, stroke = 0) +
    scale_alpha_manual(values = c("Enriched_in_Mito" = 0.8, "Enriched_in_Cyto" = 0.8,
                                  "Moderate" = 0.35, "Stable" = 0.22), guide = "none") +
    scale_color_manual(values = colors_direction, guide = "none") +
    coord_cartesian(ylim = c(0, 6), clip = "off") +
    labs(title = tissue, x = "Relocalization score", y = "-log10(P-value), capped at 6") +
    theme_publication(base_size = 11) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  if (sum(tissue_df$Label != "") > 0) {
    p <- p + geom_text_repel(
      aes(label = Label), size = 2.4, max.overlaps = 8,
      box.padding = 0.25, point.padding = 0.15,
      segment.size = 0.25, segment.color = "grey30",
      min.segment.length = 0, show.legend = FALSE
    )
  }
  volcano_list[[tissue]] <- p
}

p5 <- wrap_plots(volcano_list, ncol = 3) +
  plot_annotation(
    title = "Figure 5. Stable Dual-Localized Proteins: Relocalization Volcano",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(file.path(results_dir, "Figure_5_stable_dual_volcano.pdf"), p5, width = 14, height = 9.5, dpi = 300)

# ---------------------------------------------------------------------------
# Figure 6: non-Stable_Dual event summary (new / lost / complete transfer)
# ---------------------------------------------------------------------------
event_df <- df %>%
  filter(Pattern != "Stable_Dual") %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_counts$Tissue),
    Pattern = factor(Pattern, levels = c(
      "New_in_Mito", "New_in_Cyto", "Lost_from_Mito", "Lost_from_Cyto",
      "Complete_Transfer_M2C", "Complete_Transfer_C2M"
    ))
  ) %>%
  count(Tissue, Pattern, name = "Count") %>%
  group_by(Tissue) %>%
  mutate(Total = sum(Count), Percentage = Count / Total * 100) %>%
  ungroup()

event_totals <- event_df %>% distinct(Tissue, Total)

p6 <- ggplot(event_df, aes(x = Tissue, y = Percentage, fill = Pattern)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.35, width = 0.78,
           alpha = ALPHA_PRIMARY) +
  geom_text(
    aes(label = ifelse(Percentage >= 8, paste0(round(Percentage, 1), "%"), "")),
    position = position_stack(vjust = 0.5), size = 3.3, fontface = "bold", color = "white"
  ) +
  geom_text(data = event_totals, aes(x = Tissue, y = 103, label = paste0("n=", Total)),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  scale_fill_manual(values = colors_pattern, labels = pattern_labels, name = "Event") +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Figure 6. New, Lost, and Complete Transfer Events",
    subtitle = paste0("Key: ", results_name, " | Non-Stable_Dual patterns from 3-month to 20-month"),
    x = "Tissue", y = "Percentage of non-Stable_Dual records (%)"
  ) +
  theme_publication(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

ggsave(file.path(results_dir, "Figure_6_new_lost_transfer_events.pdf"), p6, width = 10.5, height = 7, dpi = 300)

cat("Wrote 6 figures to", results_dir, "\n")
