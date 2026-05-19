#!/usr/bin/env python3
"""
Run the complete relocalization analysis pipeline.

Pipeline order:
  00  Merge raw Excel files → data/*.txt
  01a Analysis table + t-test (GeneName key) → results_GeneName/
  01b Analysis table + t-test (ProteinID key) → results_proteinID/
  02  Basic figures (Figure 1-6) for both versions
  03  Figure 4 redesign (publication-grade scatter) for both versions
  04  Figure 5 redesign (ridge + heatmap) for both versions

Usage:
    python run_all.py
"""

import subprocess
import sys
from pathlib import Path

script_dir = Path(__file__).parent

steps = [
    # Step 0: Merge raw data
    ("00 Merge dataset",
     [sys.executable, str(script_dir / "00_merge_dataset.py")]),

    # Step 1a/1b: Build analysis tables (includes t-test)
    ("01a Analysis + t-test (GeneName key)",
     [sys.executable, str(script_dir / "01a_analysis_GeneName.py")]),
    ("01b Analysis + t-test (ProteinID key)",
     [sys.executable, str(script_dir / "01b_analysis_ProteinID.py")]),

    # Step 2: Basic figures (Figure 1-6)
    ("02 Figures (GeneName)",
     ["Rscript", str(script_dir / "02_create_figures.R"), "--results=results_GeneName"]),
    ("02 Figures (ProteinID)",
     ["Rscript", str(script_dir / "02_create_figures.R"), "--results=results_proteinID"]),

    # Step 3: Figure 4 redesign (publication scatter)
    ("03 Figure 4 redesign (GeneName)",
     ["Rscript", str(script_dir / "03_figure4_redesign.R"), "--results=results_GeneName"]),
    ("03 Figure 4 redesign (ProteinID)",
     ["Rscript", str(script_dir / "03_figure4_redesign.R"), "--results=results_proteinID"]),

    # Step 4: Figure 5 redesign (ridge + heatmap)
    ("04 Figure 5 redesign (GeneName)",
     ["Rscript", str(script_dir / "04_figure5_redesign.R"), "--results=results_GeneName"]),
    ("04 Figure 5 redesign (ProteinID)",
     ["Rscript", str(script_dir / "04_figure5_redesign.R"), "--results=results_proteinID"]),
]

failed = []
for label, cmd in steps:
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    result = subprocess.run(cmd)
    if result.returncode not in (0, 1):
        print(f"  *** FAILED (exit code {result.returncode}) ***")
        failed.append(label)

print(f"\n{'='*60}")
if failed:
    print(f"  COMPLETED WITH ERRORS ({len(failed)} failed):")
    for f in failed:
        print(f"    - {f}")
else:
    print("  ALL STEPS COMPLETED SUCCESSFULLY")
print(f"{'='*60}")
print("\nOutputs:")
print("  results_GeneName/   — GeneName key analysis + figures")
print("  results_proteinID/  — ProteinID key analysis + figures")
