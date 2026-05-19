#!/usr/bin/env python3
"""Merge per-tissue C-COMPASS Excel files into four tab-delimited inputs."""

from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
DATASET_DIR = DATA_DIR / "dataset"

SHEET_NAME = "数据矩阵"

# Non-sample columns to drop before collecting replicate intensities
EXCLUDE_COLUMNS = {
    "Accession",
    "Gene Name",
    "Symbol",
    "Description",
    "Unique peptides",
    "MitoCarta3.0_MitoPathways",
}

MERGE_JOBS = [
    {"folder": "3-month mitochondria", "output": "3month_mito_merged.txt"},
    {"folder": "20-month mitochondria ", "output": "20month_mito_merged.txt"},
    {"folder": "3-month cytoplsam", "output": "3month_cyto_merged.txt"},
    {"folder": "20-month cytoplasm", "output": "20month_cyto_merged.txt"},
]


def tissue_name_from_file(path: Path) -> str:
    """Extract tissue token (e.g. 'Brain') from a file name like 'Brain_mitochondria_3month.xlsx'."""
    name = path.stem
    for marker in ("_mitochondria_", "_cytoplasm_"):
        if marker in name:
            return name.split(marker, 1)[0]
    raise ValueError(f"Cannot infer tissue name from file name: {path.name}")


def read_tissue_file(path: Path) -> pd.DataFrame:
    """Read one tissue Excel file and rename sample columns to <tissue>_RepN."""
    df = pd.read_excel(path, sheet_name=SHEET_NAME)

    if "Accession" not in df.columns:
        raise ValueError(f"{path} missing required column: Accession")

    gene_col = "Gene Name" if "Gene Name" in df.columns else "Symbol"
    if gene_col not in df.columns:
        raise ValueError(f"{path} missing required column: Gene Name or Symbol")

    value_cols = [col for col in df.columns if col not in EXCLUDE_COLUMNS]
    if not value_cols:
        raise ValueError(f"{path} has no replicate intensity columns")

    tissue = tissue_name_from_file(path)
    keep = df[["Accession", gene_col, *value_cols]].copy()
    keep = keep.rename(columns={"Accession": "ProteinID", gene_col: "GeneName"})
    keep = keep.rename(
        columns={col: f"{tissue}_Rep{i}" for i, col in enumerate(value_cols, start=1)}
    )

    rep_cols = [col for col in keep.columns if col not in {"ProteinID", "GeneName"}]
    keep[rep_cols] = keep[rep_cols].apply(pd.to_numeric, errors="coerce")
    keep["ProteinID"] = keep["ProteinID"].astype(str)
    keep["GeneName"] = keep["GeneName"].astype(str)

    return keep


def merge_folder(folder: Path) -> pd.DataFrame:
    """Outer-merge all tissue files in one condition folder on ProteinID + GeneName."""
    files = sorted(p for p in folder.glob("*.xlsx") if not p.name.startswith("~$"))
    if not files:
        raise FileNotFoundError(f"No .xlsx files found in {folder}")

    merged = None
    for path in files:
        tissue_df = read_tissue_file(path)
        merged = tissue_df if merged is None else merged.merge(
            tissue_df, on=["ProteinID", "GeneName"], how="outer"
        )
    return merged.sort_values(["ProteinID", "GeneName"]).reset_index(drop=True)


def main() -> None:
    for job in MERGE_JOBS:
        folder = DATASET_DIR / job["folder"]
        output = DATA_DIR / job["output"]
        merged = merge_folder(folder)
        merged.to_csv(output, sep="\t", index=False, na_rep="")
        print(f"Wrote {output.name}: {len(merged)} rows, {len(merged.columns) - 2} replicate columns")


if __name__ == "__main__":
    main()
