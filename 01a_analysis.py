#!/usr/bin/env python3
"""
Build the relocalization analysis table — GeneName version.

Key strategy: Collapse protein-group duplicates by GeneName (max), then merge
mito/cyto and 3m/20m on GeneName + Tissue.

Output: results_GeneName/relocalization_analysis.xlsx
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
PSEUDOCOUNT = 1.0
RELOC_THRESHOLD = 1.5
P_VALUE_THRESHOLD = 0.05
EXCLUDE_OTHER_PATTERN = True

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
script_dir = Path(__file__).parent
project_dir = script_dir.parent
data_dir = project_dir / "data"
results_dir = project_dir / "results_GeneName"
results_dir.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Step 1: Load and average replicates, normalize tissue names
# ---------------------------------------------------------------------------
mito_3m = pd.read_csv(data_dir / "3month_mito_merged.txt", sep="\t")
mito_20m = pd.read_csv(data_dir / "20month_mito_merged.txt", sep="\t")
cyto_3m = pd.read_csv(data_dir / "3month_cyto_merged.txt", sep="\t")
cyto_20m = pd.read_csv(data_dir / "20month_cyto_merged.txt", sep="\t")


def average_replicates(df):
    """Collapse *_Rep1, *_Rep2, ... columns into a single per-tissue mean column."""
    protein_cols = ["ProteinID", "GeneName"]
    tissue_cols = [c for c in df.columns if c not in protein_cols]
    tissues = sorted({c.rsplit("_Rep", 1)[0] for c in tissue_cols})

    result = df[protein_cols].copy()
    for tissue in tissues:
        rep_cols = [c for c in tissue_cols if c.startswith(tissue + "_Rep")]
        if rep_cols:
            result[tissue] = df[rep_cols].mean(axis=1)
    return result


def normalize_tissue_names(df):
    """Merge liver/liverM/liverC -> Liver, lung/lungM/lungC -> Lung, capitalize others."""
    protein_cols = ["ProteinID", "GeneName"]
    rename_map = {}
    for col in df.columns:
        if col in protein_cols:
            continue
        new_name = col
        if col.startswith(("liverM", "liverC", "liver")):
            new_name = col.replace("liverM", "Liver").replace("liverC", "Liver").replace("liver", "Liver")
        elif col.startswith(("lungM", "lungC", "lung")):
            new_name = col.replace("lungM", "Lung").replace("lungC", "Lung").replace("lung", "Lung")
        elif col[0].islower():
            new_name = col.capitalize()
        rename_map[col] = new_name

    df_renamed = df.rename(columns=rename_map)
    tissue_cols = [c for c in df_renamed.columns if c not in protein_cols]

    result = df_renamed[protein_cols].copy()
    for tissue in sorted(set(tissue_cols)):
        matching = [c for c in tissue_cols if c == tissue]
        result[tissue] = df_renamed[matching].mean(axis=1) if len(matching) > 1 else df_renamed[tissue]
    return result


mito_3m_norm = normalize_tissue_names(average_replicates(mito_3m))
mito_20m_norm = normalize_tissue_names(average_replicates(mito_20m))
cyto_3m_norm = normalize_tissue_names(average_replicates(cyto_3m))
cyto_20m_norm = normalize_tissue_names(average_replicates(cyto_20m))

# ---------------------------------------------------------------------------
# Step 2: Collapse protein-group duplicates by GeneName, merge on GeneName + Tissue
# ---------------------------------------------------------------------------


def collapse_duplicate_genes(df):
    """Collapse rows sharing a GeneName by taking the max per tissue cell."""
    protein_cols = ["ProteinID", "GeneName"]
    tissue_cols = [c for c in df.columns if c not in protein_cols]

    df = df.copy()
    df["GeneName"] = df["GeneName"].replace({"nan": None, "": None})
    df = df.dropna(subset=["GeneName"])

    values = df.groupby("GeneName", as_index=False)[tissue_cols].max()
    ids = (df[["GeneName", "ProteinID"]]
           .dropna(subset=["ProteinID"])
           .drop_duplicates(subset=["GeneName"], keep="first"))
    return values.merge(ids, on="GeneName", how="left")[protein_cols + tissue_cols]


mito_3m_norm = collapse_duplicate_genes(mito_3m_norm)
mito_20m_norm = collapse_duplicate_genes(mito_20m_norm)
cyto_3m_norm = collapse_duplicate_genes(cyto_3m_norm)
cyto_20m_norm = collapse_duplicate_genes(cyto_20m_norm)


def melt_data(df, value_name):
    return df.melt(id_vars=["ProteinID", "GeneName"], var_name="Tissue", value_name=value_name)


m3 = melt_data(mito_3m_norm, "Mito_3month").rename(columns={"ProteinID": "PID_m3"})[["GeneName", "Tissue", "PID_m3", "Mito_3month"]]
m20 = melt_data(mito_20m_norm, "Mito_20month").rename(columns={"ProteinID": "PID_m20"})[["GeneName", "Tissue", "PID_m20", "Mito_20month"]]
c3 = melt_data(cyto_3m_norm, "Cyto_3month").rename(columns={"ProteinID": "PID_c3"})[["GeneName", "Tissue", "PID_c3", "Cyto_3month"]]
c20 = melt_data(cyto_20m_norm, "Cyto_20month").rename(columns={"ProteinID": "PID_c20"})[["GeneName", "Tissue", "PID_c20", "Cyto_20month"]]

df = m3.merge(m20, on=["GeneName", "Tissue"], how="outer")
df = df.merge(c3, on=["GeneName", "Tissue"], how="outer")
df = df.merge(c20, on=["GeneName", "Tissue"], how="outer")

df["ProteinID"] = df["PID_m3"].fillna(df["PID_m20"]).fillna(df["PID_c3"]).fillna(df["PID_c20"])
df = df.drop(columns=["PID_m3", "PID_m20", "PID_c3", "PID_c20"])

value_cols = ["Mito_3month", "Mito_20month", "Cyto_3month", "Cyto_20month"]
df[value_cols] = df[value_cols].fillna(0)
df = df[(df[value_cols] > 0).any(axis=1)]

# ---------------------------------------------------------------------------
# Step 3: Classify patterns and compute relocalization metrics
# ---------------------------------------------------------------------------
def classify_pattern(row):
    m3_, m20_ = row["Mito_3month"] > 0, row["Mito_20month"] > 0
    c3_, c20_ = row["Cyto_3month"] > 0, row["Cyto_20month"] > 0

    if m3_ and m20_ and c3_ and c20_:
        return "Stable_Dual"
    if m3_ and m20_ and c3_ and not c20_:
        return "Lost_from_Cyto"
    if m3_ and m20_ and not c3_ and c20_:
        return "New_in_Cyto"
    if m3_ and not m20_ and c3_ and c20_:
        return "Lost_from_Mito"
    if not m3_ and m20_ and c3_ and c20_:
        return "New_in_Mito"
    if m3_ and not m20_ and not c3_ and c20_:
        return "Complete_Transfer_M2C"
    if not m3_ and m20_ and c3_ and not c20_:
        return "Complete_Transfer_C2M"
    return "Other"


df["Pattern"] = df.apply(classify_pattern, axis=1)

if EXCLUDE_OTHER_PATTERN:
    df = df[df["Pattern"] != "Other"].copy()

df["FC_Mito"] = np.log2((df["Mito_20month"] + PSEUDOCOUNT) / (df["Mito_3month"] + PSEUDOCOUNT))
df["FC_Cyto"] = np.log2((df["Cyto_20month"] + PSEUDOCOUNT) / (df["Cyto_3month"] + PSEUDOCOUNT))
df["Reloc_Score"] = df["FC_Mito"] - df["FC_Cyto"]
df["P_value"] = 2 * (1 - stats.norm.cdf(np.abs(df["Reloc_Score"])))


def classify_direction(score):
    if score > RELOC_THRESHOLD:
        return "Enriched_in_Mito"
    if score < -RELOC_THRESHOLD:
        return "Enriched_in_Cyto"
    if abs(score) > 0.5:
        return "Moderate"
    return "Stable"


df["Direction"] = df["Reloc_Score"].apply(classify_direction)
df["Significant"] = np.where(
    (np.abs(df["Reloc_Score"]) > RELOC_THRESHOLD) & (df["P_value"] < P_VALUE_THRESHOLD),
    "Yes", "No",
)

# ---------------------------------------------------------------------------
# Step 4: Enrichment columns
# ---------------------------------------------------------------------------
df["Mito_Ratio_3m"] = df["Mito_3month"] / (df["Mito_3month"] + df["Cyto_3month"] + 0.001)
df["Mito_Ratio_20m"] = df["Mito_20month"] / (df["Mito_20month"] + df["Cyto_20month"] + 0.001)
df["Total_3month"] = df["Mito_3month"] + df["Cyto_3month"]
df["Total_20month"] = df["Mito_20month"] + df["Cyto_20month"]
df["FC_Total"] = np.log2((df["Total_20month"] + PSEUDOCOUNT) / (df["Total_3month"] + PSEUDOCOUNT))


def classify_change(fc):
    if fc > 0.5:
        return "Increase"
    if fc < -0.5:
        return "Decrease"
    return "Stable"


df["Mito_Change"] = df["FC_Mito"].apply(classify_change)
df["Cyto_Change"] = df["FC_Cyto"].apply(classify_change)

CHANGE_PATTERNS = {
    ("Increase", "Increase"): "Cyto_up_and_Mito_up",
    ("Increase", "Decrease"): "Cyto_to_Mito",
    ("Increase", "Stable"): "Mito_up",
    ("Decrease", "Increase"): "Mito_to_Cyto",
    ("Decrease", "Decrease"): "Cyto_down_and_Mito_down",
    ("Decrease", "Stable"): "Mito_down",
    ("Stable", "Increase"): "Cyto_up",
    ("Stable", "Decrease"): "Cyto_down",
    ("Stable", "Stable"): "Stable",
}
df["Change_Pattern"] = df.apply(
    lambda r: CHANGE_PATTERNS.get((r["Mito_Change"], r["Cyto_Change"]), "Other"), axis=1
)


def get_transfer_pattern(row):
    r3, r20 = row["Mito_Ratio_3m"], row["Mito_Ratio_20m"]
    if r3 < 0.4 and r20 > 0.6:
        return "Cyto_to_Mito"
    if r3 > 0.6 and r20 < 0.4:
        return "Mito_to_Cyto"
    if r3 > 0.6 and r20 > 0.6:
        return "Stable_in_Mito"
    if r3 < 0.4 and r20 < 0.4:
        return "Stable_in_Cyto"
    return "Balanced"


df["Transfer_Pattern"] = df.apply(get_transfer_pattern, axis=1)


def get_reloc_intensity(score):
    s = abs(score)
    if s > 3.0:
        return "Very_Strong"
    if s > 2.0:
        return "Strong"
    if s > 1.5:
        return "Moderate"
    if s > 0.5:
        return "Weak"
    return "Very_Weak"


df["Reloc_Intensity"] = df["Reloc_Score"].apply(get_reloc_intensity)


def get_expression_level(total):
    if total > 100000:
        return "Very_High"
    if total > 50000:
        return "High"
    if total > 10000:
        return "Medium"
    if total > 1000:
        return "Low"
    return "Very_Low"


df["Expression_Level_3m"] = df["Total_3month"].apply(get_expression_level)
df["Expression_Level_20m"] = df["Total_20month"].apply(get_expression_level)

df["Fold_Change_Mito"] = (df["Mito_20month"] + PSEUDOCOUNT) / (df["Mito_3month"] + PSEUDOCOUNT)
df["Fold_Change_Cyto"] = (df["Cyto_20month"] + PSEUDOCOUNT) / (df["Cyto_3month"] + PSEUDOCOUNT)
df["Fold_Change_Total"] = (df["Total_20month"] + PSEUDOCOUNT) / (df["Total_3month"] + PSEUDOCOUNT)


def get_significance_level(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


df["Significance_Level"] = df["P_value"].apply(get_significance_level)

# ---------------------------------------------------------------------------
# Step 5: Replicate-based t-test for relocalization significance
# ---------------------------------------------------------------------------
def _load_replicates_long(path):
    """Load merged file → long format (ProteinID, GeneName, Tissue, Rep, Value)."""
    raw = pd.read_csv(path, sep="\t")
    rep_cols = [c for c in raw.columns if "_Rep" in c]
    long = raw.melt(id_vars=["ProteinID", "GeneName"],
                    value_vars=rep_cols, var_name="Sample", value_name="Value")
    long["Tissue"] = long["Sample"].str.rsplit("_Rep", n=1).str[0]
    long["Rep"] = long["Sample"].str.rsplit("_Rep", n=1).str[1].astype(int)
    long = long.drop(columns="Sample")
    def _norm(t):
        if t.startswith(("liverM", "liverC", "liver")): return "Liver"
        if t.startswith(("lungM", "lungC", "lung")): return "Lung"
        return t[0].upper() + t[1:]
    long["Tissue"] = long["Tissue"].apply(_norm)
    return long

def _collapse_reps(long_df, key_col):
    return long_df.groupby([key_col, "Tissue", "Rep"], as_index=False)["Value"].max()

KEY_COL = "GeneName"
print("  Computing replicate-level t-test...")
_m3  = _collapse_reps(_load_replicates_long(data_dir / "3month_mito_merged.txt"), KEY_COL).rename(columns={"Value": "Mito_3m"})
_m20 = _collapse_reps(_load_replicates_long(data_dir / "20month_mito_merged.txt"), KEY_COL).rename(columns={"Value": "Mito_20m"})
_c3  = _collapse_reps(_load_replicates_long(data_dir / "3month_cyto_merged.txt"), KEY_COL).rename(columns={"Value": "Cyto_3m"})
_c20 = _collapse_reps(_load_replicates_long(data_dir / "20month_cyto_merged.txt"), KEY_COL).rename(columns={"Value": "Cyto_20m"})

_rep = (_m3.merge(_m20, on=[KEY_COL, "Tissue", "Rep"], how="outer")
           .merge(_c3,  on=[KEY_COL, "Tissue", "Rep"], how="outer")
           .merge(_c20, on=[KEY_COL, "Tissue", "Rep"], how="outer"))
_rep = _rep.dropna(subset=["Mito_3m", "Mito_20m", "Cyto_3m", "Cyto_20m"])
_rep["Delta_rep"] = (np.log2((_rep["Mito_20m"] + PSEUDOCOUNT) / (_rep["Mito_3m"] + PSEUDOCOUNT))
                   - np.log2((_rep["Cyto_20m"] + PSEUDOCOUNT) / (_rep["Cyto_3m"] + PSEUDOCOUNT)))

def _ttest(g):
    x = g["Delta_rep"].values
    n = len(x)
    if n < 2:
        return pd.Series({"n_rep": n, "Delta_mean": np.nan, "Delta_sd": np.nan, "t_stat": np.nan, "P_ttest": np.nan})
    if np.allclose(x, x[0]):
        return pd.Series({"n_rep": n, "Delta_mean": x.mean(), "Delta_sd": 0.0, "t_stat": np.nan, "P_ttest": 1.0})
    t, p = stats.ttest_1samp(x, 0.0)
    return pd.Series({"n_rep": n, "Delta_mean": x.mean(), "Delta_sd": x.std(ddof=1), "t_stat": t, "P_ttest": p})

_tt = _rep.groupby([KEY_COL, "Tissue"]).apply(_ttest, include_groups=False).reset_index()

# FDR per tissue
fdr_pt = np.full(len(_tt), np.nan)
for _, idx in _tt.groupby("Tissue").groups.items():
    sub = _tt.loc[idx]; valid = sub["P_ttest"].notna()
    if valid.sum() > 0:
        _, p_adj, _, _ = multipletests(sub.loc[valid, "P_ttest"].values, method="fdr_bh")
        fdr_pt[idx[valid].astype(int)] = p_adj
_tt["FDR_BH"] = fdr_pt
_tt["NegLog10_P_ttest"] = -np.log10(_tt["P_ttest"].clip(lower=1e-300))
_tt["NegLog10_FDR"] = -np.log10(_tt["FDR_BH"].clip(lower=1e-300))
_tt["Significant_ttest"] = np.where((_tt["Delta_mean"].abs() > RELOC_THRESHOLD) & (_tt["P_ttest"] < P_VALUE_THRESHOLD), "Yes", "No")
_tt["Significant_FDR"] = np.where((_tt["Delta_mean"].abs() > RELOC_THRESHOLD) & (_tt["FDR_BH"] < P_VALUE_THRESHOLD), "Yes", "No")

# Merge t-test results into main df
ttest_cols = ["n_rep", "Delta_mean", "Delta_sd", "t_stat", "P_ttest", "FDR_BH",
              "NegLog10_P_ttest", "NegLog10_FDR", "Significant_ttest", "Significant_FDR"]
df = df.merge(_tt[[KEY_COL, "Tissue"] + ttest_cols], on=[KEY_COL, "Tissue"], how="left")
print(f"  T-test: {_tt['P_ttest'].notna().sum()} records with p-value, {(_tt['Significant_ttest']=='Yes').sum()} significant")

# ---------------------------------------------------------------------------
# Step 6: Reorder columns and write Excel workbook
# ---------------------------------------------------------------------------
df = df[[
    "ProteinID", "GeneName", "Tissue",
    "Pattern",
    "Mito_3month", "Mito_20month", "Cyto_3month", "Cyto_20month",
    "Total_3month", "Total_20month",
    "FC_Mito", "FC_Cyto", "FC_Total", "Reloc_Score",
    "Fold_Change_Mito", "Fold_Change_Cyto", "Fold_Change_Total",
    "P_value", "Significant", "Significance_Level",
    "n_rep", "Delta_mean", "Delta_sd", "t_stat", "P_ttest", "FDR_BH",
    "NegLog10_P_ttest", "NegLog10_FDR", "Significant_ttest", "Significant_FDR",
    "Direction", "Reloc_Intensity",
    "Transfer_Pattern", "Mito_Ratio_3m", "Mito_Ratio_20m",
    "Mito_Change", "Cyto_Change", "Change_Pattern",
    "Expression_Level_3m", "Expression_Level_20m",
]]

output_file = results_dir / "relocalization_analysis.xlsx"
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    df.to_excel(writer, sheet_name="All_Proteins", index=False)

    stable = df[df["Pattern"] == "Stable_Dual"].sort_values("Reloc_Score", key=abs, ascending=False)
    stable.to_excel(writer, sheet_name="Stable_Dual_Localized", index=False)

    newly = df[df["Pattern"].isin(["New_in_Cyto", "New_in_Mito"])].sort_values(
        ["Pattern", "Reloc_Score"], ascending=[True, False]
    )
    newly.to_excel(writer, sheet_name="Newly_Appeared", index=False)

    lost = df[df["Pattern"].isin(["Lost_from_Cyto", "Lost_from_Mito"])].sort_values(
        ["Pattern", "Reloc_Score"], ascending=[True, False]
    )
    lost.to_excel(writer, sheet_name="Lost_Proteins", index=False)

    transfer = df[df["Pattern"].isin(["Complete_Transfer_M2C", "Complete_Transfer_C2M"])].sort_values(
        ["Pattern", "Reloc_Score"], ascending=[True, False]
    )
    transfer.to_excel(writer, sheet_name="Complete_Transfer", index=False)

    df.groupby(["Tissue", "Direction"]).size().unstack(fill_value=0).to_excel(writer, sheet_name="Tissue_Summary")
    df.groupby(["Tissue", "Pattern"]).size().unstack(fill_value=0).to_excel(writer, sheet_name="Pattern_Stats")
    df.groupby(["Tissue", "Transfer_Pattern"]).size().unstack(fill_value=0).to_excel(writer, sheet_name="Transfer_Stats")
    df.groupby(["Tissue", "Change_Pattern"]).size().unstack(fill_value=0).to_excel(writer, sheet_name="Change_Stats")

print(f"[GeneName key] Wrote {output_file.name}: {len(df)} records, {df['GeneName'].nunique()} genes, {df['Tissue'].nunique()} tissues")
