# ===========================================================================
# Shared color palette for all figures.
# Source this file from any figure script:  source("scripts/palette.R")
# ===========================================================================

# ---------------------------------------------------------------------------
# Tissues — 9 colors
# ---------------------------------------------------------------------------
tissue_order <- c("Stomach", "Kidney", "Brain",
                  "Intestine", "Lung", "Spleen",
                  "Muscle", "Liver", "Heart")

tissue_colors <- c(
  "Brain"     = "#E76F51",
  "Heart"     = "#2A9D8F",
  "Intestine" = "#F4A261",
  "Kidney"    = "#E9C46A",
  "Liver"     = "#264653",
  "Lung"      = "#A78BFA",
  "Muscle"    = "#EF476F",
  "Spleen"    = "#06A0C0",
  "Stomach"   = "#8D5524"
)

# ---------------------------------------------------------------------------
# Pattern — 7 categories (Nature-style npg palette)
# ---------------------------------------------------------------------------
pattern_colors <- c(
  "Stable_Dual"             = "#00A087",
  "New_in_Cyto"             = "#4DBBD5",
  "New_in_Mito"             = "#3C5488",
  "Lost_from_Cyto"          = "#E64B35",
  "Lost_from_Mito"          = "#F39B7F",
  "Complete_Transfer_M2C"   = "#B09C85",
  "Complete_Transfer_C2M"   = "#8491B4"
)
pattern_labels <- c(
  "Stable_Dual"             = "Stable dual",
  "New_in_Cyto"             = "New in cyto",
  "New_in_Mito"             = "New in mito",
  "Lost_from_Cyto"          = "Lost from cyto",
  "Lost_from_Mito"          = "Lost from mito",
  "Complete_Transfer_M2C"   = "Mito to cyto",
  "Complete_Transfer_C2M"   = "Cyto to mito"
)

# ---------------------------------------------------------------------------
# Change_Pattern — 9 categories (Stable_Dual subset)
# ---------------------------------------------------------------------------
change_colors <- c(
  "Cyto_to_Mito"             = "#E64B35",
  "Mito_to_Cyto"             = "#4DBBD5",
  "Cyto_up_and_Mito_up"      = "#00A087",
  "Cyto_down_and_Mito_down"  = "#8491B4",
  "Mito_up"                  = "#F39B7F",
  "Mito_down"                = "#B09C85",
  "Cyto_up"                  = "#91D1C2",
  "Cyto_down"                = "#3C5488",
  "Stable"                   = "#C9C9C9"
)
change_labels <- c(
  "Cyto_to_Mito"             = "Cyto to Mito",
  "Mito_to_Cyto"             = "Mito to Cyto",
  "Cyto_up_and_Mito_up"      = "Cyto up and Mito up",
  "Cyto_down_and_Mito_down"  = "Cyto down and Mito down",
  "Mito_up"                  = "Mito up",
  "Mito_down"                = "Mito down",
  "Cyto_up"                  = "Cyto up",
  "Cyto_down"                = "Cyto down",
  "Stable"                   = "Stable"
)

# ---------------------------------------------------------------------------
# Overlap (Figure 1) — uses Pattern palette colors for consistency
#   Mito_only matches Lost_from_Cyto color (red) — stays in mito
#   Cyto_only matches New_in_Cyto color (cyan)   — only in cyto
#   Shared    matches Stable_Dual color (teal)
# ---------------------------------------------------------------------------
overlap_colors <- c(
  "Mito_only" = "#E64B35",
  "Shared"    = "#00A087",
  "Cyto_only" = "#4DBBD5"
)

# ---------------------------------------------------------------------------
# Heatmap — divergent (Reloc score)
# ---------------------------------------------------------------------------
heatmap_low  <- "#118AB2"   # blue (negative score, up in cyto)
heatmap_high <- "#D62828"   # red  (positive score, up in mito)

# ---------------------------------------------------------------------------
# Transparency constants — applied uniformly across figures
# ---------------------------------------------------------------------------
ALPHA_PRIMARY    <- 0.85   # bars, ridge fills, top-layer points
ALPHA_SECONDARY  <- 0.65   # mid-layer scatter points
ALPHA_BACKGROUND <- 0.30   # faded background points / shapes
