#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mitochondrial Form Factor Calculator for Nellie Output
-------------------------------------------------------
This script processes Nellie organelle feature CSV files and calculates
mitochondrial form factor and related morphological metrics.

Form Factor = Perimeter² / (4π × Area)
- Values close to 1: circular/round mitochondria
- Values > 1: elongated mitochondria
"""

import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path
from typing import Tuple
import warnings

warnings.filterwarnings('ignore')

# Try to import tkinter for GUI support
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox

    HAS_GUI = True
except ImportError:
    HAS_GUI = False
    print("Note: tkinter not available. GUI folder selection disabled.")
    print("Please use --in and --out command line arguments instead.\n")


# ============================================================================
# DIRECTORY PICKER
# ============================================================================
def pick_paths_if_needed(args) -> Tuple[str, str]:
    """
    Open GUI dialogs to select input and output directories if not provided via CLI.

    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments

    Returns:
    --------
    in_dir : str
        Input directory path
    out_dir : str
        Output directory path
    """
    if not HAS_GUI:
        print("\nError: GUI not available (tkinter not installed).")
        print("Please provide directories using command line arguments:")
        print("  python calculate_mito_form_factor.py --in /path/to/input --out /path/to/output\n")
        raise SystemExit(1)

    root = tk.Tk()
    root.withdraw()

    in_dir = args.in_dir or filedialog.askdirectory(
        title="Select folder containing Nellie CSV files"
    )
    if not in_dir:
        messagebox.showerror("Missing folder", "No input folder selected.")
        raise SystemExit(1)

    out_dir = args.out_dir or filedialog.askdirectory(
        title="Select output folder for results"
    )
    if not out_dir:
        messagebox.showerror("Missing folder", "No output folder selected.")
        raise SystemExit(1)

    root.update()
    root.destroy()
    return in_dir, out_dir


# ============================================================================
# FORM FACTOR CALCULATIONS
# ============================================================================
def calculate_ellipse_perimeter(major_axis, minor_axis):
    """
    Calculate approximate perimeter of an ellipse using Ramanujan's approximation.

    Parameters:
    -----------
    major_axis : float or array
        Length of major axis (a)
    minor_axis : float or array
        Length of minor axis (b)

    Returns:
    --------
    perimeter : float or array
        Approximate perimeter
    """
    a = major_axis / 2  # Convert axis length to semi-axis
    b = minor_axis / 2

    # Ramanujan's approximation: P ≈ π(a + b)[1 + 3h/(10 + √(4-3h))]
    # where h = ((a-b)/(a+b))²
    h = ((a - b) / (a + b)) ** 2
    perimeter = np.pi * (a + b) * (1 + (3 * h) / (10 + np.sqrt(4 - 3 * h)))

    return perimeter


def calculate_form_factor(area, perimeter):
    """
    Calculate form factor (circularity index).

    Form Factor = Perimeter² / (4π × Area)
    - Value = 1: perfect circle
    - Value > 1: elongated shape

    Parameters:
    -----------
    area : float or array
        Area of the object
    perimeter : float or array
        Perimeter of the object

    Returns:
    --------
    form_factor : float or array
        Form factor value
    """
    form_factor = (perimeter ** 2) / (4 * np.pi * area)
    return form_factor


def process_nellie_file(filepath):
    """
    Process a single Nellie organelle features CSV file.

    Parameters:
    -----------
    filepath : str or Path
        Path to the Nellie CSV file

    Returns:
    --------
    df : DataFrame
        DataFrame with added form factor columns
    """
    # Read the CSV
    df = pd.read_csv(filepath)

    # Check if required columns exist
    required_cols = ['organelle_area_raw', 'organelle_axis_length_maj_raw',
                     'organelle_axis_length_min_raw']

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Warning: {filepath.name} missing columns: {missing_cols}")
        return None

    # Extract relevant measurements
    area = df['organelle_area_raw']
    major_axis = df['organelle_axis_length_maj_raw']
    minor_axis = df['organelle_axis_length_min_raw']

    # Calculate perimeter (approximation from ellipse)
    perimeter = calculate_ellipse_perimeter(major_axis, minor_axis)

    # Calculate form factor
    form_factor = calculate_form_factor(area, perimeter)

    # Add new columns to dataframe
    df['perimeter_approx'] = perimeter
    df['form_factor'] = form_factor

    # Add aspect ratio if not already present
    if 'aspect_ratio_calc' not in df.columns:
        df['aspect_ratio_calc'] = major_axis / minor_axis

    # Calculate circularity (inverse of form factor, normalized to 0-1)
    df['circularity'] = 1 / form_factor

    # Add filename for tracking
    df['source_file'] = filepath.name

    return df


def summarize_morphology(df, groupby_cols=None):
    """
    Calculate summary statistics for morphological parameters.

    Parameters:
    -----------
    df : DataFrame
        Processed dataframe with form factor
    groupby_cols : list, optional
        Columns to group by (e.g., ['source_file', 't'])

    Returns:
    --------
    summary : DataFrame
        Summary statistics
    """
    metrics = ['organelle_area_raw', 'form_factor', 'aspect_ratio_calc',
               'circularity', 'perimeter_approx']

    if groupby_cols:
        summary = df.groupby(groupby_cols)[metrics].agg([
            'count', 'mean', 'std', 'median', 'min', 'max'
        ])
    else:
        summary = df[metrics].agg([
            'count', 'mean', 'std', 'median', 'min', 'max'
        ])

    return summary


def batch_process_folder(input_folder, output_folder=None, pattern="*features_organelles.csv"):
    """
    Process all Nellie CSV files in a folder.

    Parameters:
    -----------
    input_folder : str or Path
        Folder containing Nellie CSV files
    output_folder : str or Path, optional
        Folder to save results (default: same as input_folder)
    pattern : str
        File pattern to match (default: "*features_organelles.csv")

    Returns:
    --------
    all_data : DataFrame
        Combined dataframe with all processed data
    summary : DataFrame
        Summary statistics per file
    """
    input_path = Path(input_folder)

    if output_folder is None:
        output_path = input_path
    else:
        output_path = Path(output_folder)
        output_path.mkdir(parents=True, exist_ok=True)

    # Find all matching files
    csv_files = list(input_path.glob(pattern))

    if len(csv_files) == 0:
        print(f"No files matching '{pattern}' found in {input_folder}")
        return None, None

    print(f"Found {len(csv_files)} files to process")

    # Process each file
    all_dataframes = []

    for i, filepath in enumerate(csv_files, 1):
        print(f"Processing {i}/{len(csv_files)}: {filepath.name}")

        df = process_nellie_file(filepath)

        if df is not None:
            all_dataframes.append(df)

    # Combine all dataframes
    if len(all_dataframes) == 0:
        print("No valid data processed")
        return None, None

    all_data = pd.concat(all_dataframes, ignore_index=True)

    # Calculate summary statistics per file
    summary = summarize_morphology(all_data, groupby_cols=['source_file'])

    # Save combined results
    output_file = output_path / "all_mitochondria_form_factor.csv"
    all_data.to_csv(output_file, index=False)
    print(f"\n✓ Saved combined data to: {output_file}")

    # Save summary statistics
    summary_file = output_path / "mitochondria_morphology_summary.csv"
    summary.to_csv(summary_file)
    print(f"✓ Saved summary statistics to: {summary_file}")

    # Create a simplified summary (mean values only)
    simple_summary = all_data.groupby('source_file').agg({
        'organelle_area_raw': ['count', 'mean', 'std'],
        'form_factor': ['mean', 'std', 'median'],
        'aspect_ratio_calc': ['mean', 'std'],
        'circularity': ['mean', 'std']
    })
    simple_summary.columns = ['_'.join(col).strip() for col in simple_summary.columns.values]

    simple_summary_file = output_path / "mitochondria_summary_simple.csv"
    simple_summary.to_csv(simple_summary_file)
    print(f"✓ Saved simple summary to: {simple_summary_file}")

    return all_data, summary


def print_interpretation_guide():
    """Print interpretation guide for form factor values."""
    guide = """
    ═══════════════════════════════════════════════════════════════
    FORM FACTOR INTERPRETATION GUIDE
    ═══════════════════════════════════════════════════════════════

    Form Factor = Perimeter² / (4π × Area)

    Values:
    -------
    • 1.0        = Perfect circle (round mitochondria)
    • 1.0 - 2.0  = Slightly elongated
    • 2.0 - 4.0  = Moderately elongated (tubular)
    • > 4.0      = Highly elongated (networked)

    Aspect Ratio (Major axis / Minor axis):
    ----------------------------------------
    • 1.0        = Perfect circle
    • 1.0 - 2.0  = Oval/slightly elongated
    • 2.0 - 3.0  = Moderately elongated
    • > 3.0      = Highly elongated (rod-like)

    Circularity (1 / Form Factor):
    -------------------------------
    • 1.0        = Perfect circle
    • 0.5 - 1.0  = Slightly elongated
    • 0.25 - 0.5 = Moderately elongated
    • < 0.25     = Highly elongated

    ═══════════════════════════════════════════════════════════════
    """
    print(guide)


def parse_args():
    """Parse command line arguments."""
    ap = argparse.ArgumentParser(
        description="Calculate mitochondrial form factor from Nellie CSV files"
    )
    ap.add_argument("--in", dest="in_dir",
                    help="Input folder containing Nellie *features_organelles.csv files")
    ap.add_argument("--out", dest="out_dir",
                    help="Output folder for results")
    ap.add_argument("--pattern", default="*features_organelles.csv",
                    help="File pattern to match (default: *features_organelles.csv)")
    return ap.parse_args()


def main():
    """Main function with CLI and GUI support."""
    args = parse_args()

    # Get directories from CLI or GUI picker
    if args.in_dir and args.out_dir:
        input_folder = args.in_dir
        output_folder = args.out_dir
    else:
        input_folder, output_folder = pick_paths_if_needed(args)

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    print("=" * 70)
    print("MITOCHONDRIAL FORM FACTOR CALCULATOR")
    print("=" * 70)
    print(f"\nInput folder: {input_folder}")
    print(f"Output folder: {output_folder}")
    print(f"File pattern: {args.pattern}")

    # Process all files
    all_data, summary = batch_process_folder(
        input_folder=input_folder,
        output_folder=output_folder,
        pattern=args.pattern
    )

    if all_data is not None:
        print(f"\n{'=' * 70}")
        print(f"PROCESSING COMPLETE")
        print(f"{'=' * 70}")
        print(f"Total mitochondria analyzed: {len(all_data)}")
        print(f"Total files processed: {all_data['source_file'].nunique()}")

        print(f"\n{'OVERALL STATISTICS':^70}")
        print(f"{'-' * 70}")
        print(f"Mean form factor: {all_data['form_factor'].mean():.3f} ± {all_data['form_factor'].std():.3f}")
        print(f"Median form factor: {all_data['form_factor'].median():.3f}")
        print(
            f"Mean aspect ratio: {all_data['aspect_ratio_calc'].mean():.3f} ± {all_data['aspect_ratio_calc'].std():.3f}")
        print(f"Mean area: {all_data['organelle_area_raw'].mean():.2f} ± {all_data['organelle_area_raw'].std():.2f}")

        print_interpretation_guide()


if __name__ == "__main__":
    main()
