#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch Mitochondrial Morphology Analyzer
Processes Nellie features_organelles.csv files to calculate cell-level metrics.

Author: [Your Name]
License: MIT
"""

import os
import argparse
from typing import Tuple, List
from pathlib import Path

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import filedialog, messagebox


# ---------------- FILE PICKER ----------------
def pick_paths_if_needed(args) -> Tuple[str, str]:
    """Open GUI dialogs to select input and output directories."""
    root = tk.Tk()
    root.withdraw()

    in_dir = args.in_dir or filedialog.askdirectory(
        title="Select folder containing Nellie CSV files"
    )
    if not in_dir:
        messagebox.showerror("Missing folder", "No input folder selected.")
        raise SystemExit(1)

    out_dir = args.out_dir or filedialog.askdirectory(
        title="Select output folder"
    )
    if not out_dir:
        messagebox.showerror("Missing folder", "No output folder selected.")
        raise SystemExit(1)

    root.update()
    root.destroy()
    return in_dir, out_dir


# ---------------- FIND CSV FILES ----------------
def find_csv_files(in_dir: str, recursive: bool = False) -> List[str]:
    """
    Find all CSV files in the input directory.
    
    Args:
        in_dir: Input directory path
        recursive: If True, search subdirectories recursively
    
    Returns:
        List of CSV file paths
    """
    if recursive:
        csv_files = list(Path(in_dir).rglob("*.csv"))
    else:
        csv_files = list(Path(in_dir).glob("*.csv"))
    
    return [str(f) for f in csv_files]


# ---------------- ANALYSIS ----------------
def analyze_cell(csv_path: str) -> pd.DataFrame:
    """
    Analyze mitochondrial morphology from Nellie output CSV.
    
    Args:
        csv_path: Path to features_organelles.csv from Nellie
    
    Returns:
        DataFrame with cell-level metrics
    
    Raises:
        ValueError: If required columns are missing
    """
    df = pd.read_csv(csv_path)

    # ---- EXACT Nellie column names (no guessing) ----
    AREA_COL = "organelle_area_raw"
    MAJOR_COL = "organelle_axis_length_maj_raw"
    MINOR_COL = "organelle_axis_length_min_raw"
    EXTENT_COL = "organelle_extent_raw"
    SOLIDITY_COL = "organelle_solidity_raw"

    required = [AREA_COL, MAJOR_COL, MINOR_COL, EXTENT_COL, SOLIDITY_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}\n"
            f"Available columns are:\n{list(df.columns)}"
        )

    # Convert to numeric
    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=[AREA_COL, MAJOR_COL, MINOR_COL])

    # ---- per-mitochondrion arrays ----
    area = df[AREA_COL].values
    length_maj = df[MAJOR_COL].values
    length_min = df[MINOR_COL].values
    solidity = df[SOLIDITY_COL].values
    extent = df[EXTENT_COL].values

    n_mito = len(area)

    # ---- size metrics ----
    total_area = np.sum(area)
    mean_area = np.mean(area)
    median_area = np.median(area)
    sd_area = np.std(area, ddof=1) if n_mito > 1 else np.nan

    # ---- length metrics ----
    mean_length = np.mean(length_maj)
    max_length = np.max(length_maj)
    sd_length = np.std(length_maj, ddof=1) if n_mito > 1 else np.nan

    # ---- shape metrics ----
    mean_solidity = np.mean(solidity)
    mean_extent = np.mean(extent)

    aspect_ratio = length_maj / length_min
    mean_aspect_ratio = np.mean(aspect_ratio)

    # ---- derived / composite metrics ----
    fragmentation_index = n_mito / total_area if total_area > 0 else np.nan
    network_dominance_index = max_length / mean_length if mean_length > 0 else np.nan
    shape_heterogeneity_index = sd_length / mean_length if mean_length > 0 else np.nan
    tubularity_score = mean_aspect_ratio * (1.0 - mean_solidity)

    out = pd.DataFrame([{
        "filename": os.path.basename(csv_path),
        "n_mitochondria": n_mito,
        "total_mito_area_raw": total_area,
        "mean_mito_area_raw": mean_area,
        "median_mito_area_raw": median_area,
        "sd_mito_area_raw": sd_area,
        "mean_mito_length_raw": mean_length,
        "max_mito_length_raw": max_length,
        "sd_mito_length_raw": sd_length,
        "mean_solidity": mean_solidity,
        "mean_extent": mean_extent,
        "mean_aspect_ratio": mean_aspect_ratio,
        "fragmentation_index": fragmentation_index,
        "network_dominance_index": network_dominance_index,
        "shape_heterogeneity_index": shape_heterogeneity_index,
        "tubularity_score": tubularity_score
    }])

    return out


# ---------------- BATCH PROCESSING ----------------
def batch_process(in_dir: str, out_dir: str, recursive: bool = False):
    """
    Process all CSV files in input directory.
    
    Args:
        in_dir: Input directory containing Nellie CSV files
        out_dir: Output directory for results
        recursive: If True, search subdirectories recursively
    """
    csv_files = find_csv_files(in_dir, recursive)
    
    if not csv_files:
        print(f"No CSV files found in {in_dir}")
        return
    
    print(f"Found {len(csv_files)} CSV file(s) to process")
    
    all_results = []
    failed_files = []
    
    for i, csv_path in enumerate(csv_files, 1):
        print(f"[{i}/{len(csv_files)}] Processing: {os.path.basename(csv_path)}")
        
        try:
            result = analyze_cell(csv_path)
            all_results.append(result)
            
            # Save individual file results
            base = os.path.splitext(os.path.basename(csv_path))[0]
            out_csv = os.path.join(out_dir, f"{base}_cell_metrics.csv")
            result.to_csv(out_csv, index=False)
            print(f"  ✓ Saved: {out_csv}")
            
        except Exception as e:
            print(f"  ✗ Error processing {csv_path}: {e}")
            failed_files.append((csv_path, str(e)))
    
    # Save combined results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined_path = os.path.join(out_dir, "combined_cell_metrics.csv")
        combined.to_csv(combined_path, index=False)
        print(f"\n✓ Combined results saved: {combined_path}")
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"  Successful: {len(all_results)}/{len(csv_files)}")
    if failed_files:
        print(f"  Failed: {len(failed_files)}")
        print("\nFailed files:")
        for path, error in failed_files:
            print(f"  - {os.path.basename(path)}: {error}")
    print(f"{'='*60}")


# ---------------- CLI ----------------
def parse_args():
    """Parse command line arguments."""
    ap = argparse.ArgumentParser(
        description="Batch process Nellie mitochondria CSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # GUI mode (interactive)
  python batch_analyze.py
  
  # Command line mode
  python batch_analyze.py --in /path/to/csvs --out /path/to/output
  
  # Recursive search in subdirectories
  python batch_analyze.py --in /path/to/csvs --out /path/to/output --recursive
        """
    )
    ap.add_argument("--in", dest="in_dir", 
                    help="Input folder containing CSV files")
    ap.add_argument("--out", dest="out_dir", 
                    help="Output folder for results")
    ap.add_argument("--recursive", action="store_true",
                    help="Search for CSV files recursively in subdirectories")
    return ap.parse_args()


def main():
    """Main execution function."""
    args = parse_args()

    if args.in_dir and args.out_dir:
        in_dir, out_dir = args.in_dir, args.out_dir
    else:
        in_dir, out_dir = pick_paths_if_needed(args)

    os.makedirs(out_dir, exist_ok=True)

    batch_process(in_dir, out_dir, recursive=args.recursive)


if __name__ == "__main__":
    main()
