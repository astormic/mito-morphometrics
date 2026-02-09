#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Form Factor
==========================================
Statistical comparison of mitochondrial form factor and morphology metrics
between two or three experimental replicates.

Author: Amir Rahmani
License: MIT
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox

    HAS_GUI = True
except ImportError:
    HAS_GUI = False
    print("Note: tkinter not available. GUI file selection disabled.")
    print("Please use --rep1, --rep2, and --out command line arguments instead.\n")


METRIC_LABELS = {
    'organelle_area_raw_count': 'Number of Mitochondria per Cell',
    'organelle_area_raw_mean': 'Mean Mitochondrial Area',
    'organelle_area_raw_std': 'SD of Mitochondrial Area',
    'form_factor_mean': 'Mean Form Factor',
    'form_factor_std': 'SD of Form Factor',
    'form_factor_median': 'Median Form Factor',
    'aspect_ratio_calc_mean': 'Mean Aspect Ratio',
    'aspect_ratio_calc_std': 'SD of Aspect Ratio',
    'circularity_mean': 'Mean Circularity',
    'circularity_std': 'SD of Circularity'
}


def pick_files_if_needed(args):
    """
    Open GUI dialogs to select replicate files and output location if not provided via CLI.

    Parameters:
    -----------
    args : argparse.Namespace
        Command line arguments

    Returns:
    --------
    rep_paths : list[str]
        Paths to replicate CSV files
    out_path : str
        Output prefix path
    """
    if not HAS_GUI:
        print("\nError: GUI not available (tkinter not installed).")
        print("Please provide files using command line arguments:")
        print("  python compare_form_factor_replicates.py --rep1 rep1.csv --rep2 rep2.csv --out results/comparison\n")
        raise SystemExit(1)

    root = tk.Tk()
    root.withdraw()

    rep_paths = []

    rep1_path = args.rep1 or filedialog.askopenfilename(
        title="Select Replicate 1 CSV file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not rep1_path:
        messagebox.showerror("Missing file", "No replicate 1 file selected.")
        raise SystemExit(1)
    rep_paths.append(rep1_path)

    rep2_path = args.rep2 or filedialog.askopenfilename(
        title="Select Replicate 2 CSV file",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not rep2_path:
        messagebox.showerror("Missing file", "No replicate 2 file selected.")
        raise SystemExit(1)
    rep_paths.append(rep2_path)

    rep3_path = args.rep3
    if not rep3_path:
        add_third = messagebox.askyesno(
            "Third replicate",
            "Do you want to select a third replicate file?"
        )
        if add_third:
            rep3_path = filedialog.askopenfilename(
                title="Select Replicate 3 CSV file",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
            )
            if not rep3_path:
                messagebox.showerror("Missing file", "No replicate 3 file selected.")
                raise SystemExit(1)
    if rep3_path:
        rep_paths.append(rep3_path)

    if args.out:
        out_path = args.out
    else:
        out_dir = filedialog.askdirectory(
            title="Select output folder for results"
        )
        if not out_dir:
            messagebox.showerror("Missing folder", "No output folder selected.")
            raise SystemExit(1)

        out_path = str(Path(out_dir) / "form_factor_comparison")

    root.update()
    root.destroy()
    return rep_paths, out_path


def load_replicates(rep_paths: list[str]) -> tuple:
    """Load and label replicate data."""
    replicates = []
    for idx, path in enumerate(rep_paths, start=1):
        rep = pd.read_csv(path)
        rep['replicate'] = f'Replicate {idx}'
        replicates.append(rep)

    combined = pd.concat(replicates, ignore_index=True)

    return replicates, combined


def perform_statistical_tests(replicates: list[pd.DataFrame],
                              metrics: list) -> pd.DataFrame:
    """
    Perform statistical tests on each metric.

    Uses Shapiro-Wilk test to assess normality, then:
    - Welch's t-test or Mann-Whitney U for two groups
    - One-way ANOVA or Kruskal-Wallis for three groups

    Args:
        replicates: List of replicate DataFrames
        metrics: List of metric column names to test

    Returns:
        DataFrame with statistical results
    """
    results = []

    for metric in metrics:
        data_by_rep = [rep[metric].dropna() for rep in replicates]

        means = [data.mean() for data in data_by_rep]
        sems = [data.sem() for data in data_by_rep]
        ns = [len(data) for data in data_by_rep]

        normality_p = [
            stats.shapiro(data)[1] if len(data) >= 3 else None
            for data in data_by_rep
        ]
        all_normal = all(p is not None and p > 0.05 for p in normality_p)

        if any(n < 2 for n in ns):
            test_used = "Insufficient data"
            p_value = np.nan
        elif len(replicates) == 2:
            data1, data2 = data_by_rep
            if all_normal:
                _, p_value = stats.ttest_ind(data1, data2, equal_var=False)
                test_used = "Welch's t-test"
            else:
                _, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
                test_used = "Mann-Whitney U"
        else:
            if all_normal:
                _, p_value = stats.f_oneway(*data_by_rep)
                test_used = "One-way ANOVA"
            else:
                _, p_value = stats.kruskal(*data_by_rep)
                test_used = "Kruskal-Wallis"

        if pd.isna(p_value):
            sig = 'na'
        elif p_value < 0.001:
            sig = '***'
        elif p_value < 0.01:
            sig = '**'
        elif p_value < 0.05:
            sig = '*'
        else:
            sig = 'ns'

        row = {
            'Metric': METRIC_LABELS.get(metric, metric),
            'Test': test_used,
            'p_value': p_value,
            'Significance': sig
        }
        for idx, (mean, sem, n) in enumerate(zip(means, sems, ns), start=1):
            row[f'Rep{idx}_mean'] = mean
            row[f'Rep{idx}_SEM'] = sem
            row[f'Rep{idx}_n'] = n
        results.append(row)

    return pd.DataFrame(results)


def create_comparison_plots(replicates: list[pd.DataFrame],
                            metrics: list,
                            metric_labels: dict, output_prefix: str):
    """Create boxplot visualizations for all and key metrics."""

    n_metrics = len(metrics)
    n_cols = 3
    n_rows = int(np.ceil(n_metrics / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4 * n_rows))
    axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes

    n_reps = len(replicates)
    colors = list(plt.get_cmap("Set2").colors[:n_reps])

    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        data_by_rep = [rep[metric].dropna() for rep in replicates]

        x_pos = np.arange(1, n_reps + 1)
        bp = ax.boxplot(
            data_by_rep,
            positions=x_pos,
            widths=0.6,
            patch_artist=True,
            showfliers=False
        )
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
            patch.set_edgecolor('black')

        np.random.seed(42)
        for rep_idx, data in enumerate(data_by_rep):
            x_jitter = np.random.normal(rep_idx + 1, 0.06, size=len(data))
            ax.scatter(x_jitter, data, alpha=0.4, s=25, color='black', zorder=3)

        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'Rep {idx + 1}\n(n={len(data)})'
                            for idx, data in enumerate(data_by_rep)])
        ax.set_ylabel(metric_labels[metric], fontsize=11, fontweight='bold')
        ax.set_title(metric_labels[metric], fontsize=12, fontweight='bold', pad=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--')

    for idx in range(n_metrics, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_all_metrics_boxplots.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_all_metrics_boxplots.svg', format='svg', bbox_inches='tight')
    plt.close()

    key_metrics = [
        'organelle_area_raw_count',
        'organelle_area_raw_mean',
        'form_factor_mean',
        'form_factor_median',
        'aspect_ratio_calc_mean',
        'circularity_mean'
    ]

    fig2, axes2 = plt.subplots(2, 3, figsize=(16, 10))
    axes2 = axes2.flatten()

    for idx, metric in enumerate(key_metrics):
        ax = axes2[idx]

        data_by_rep = [rep[metric].dropna() for rep in replicates]

        x_pos = np.arange(1, n_reps + 1)
        bp = ax.boxplot(
            data_by_rep,
            positions=x_pos,
            widths=0.6,
            patch_artist=True,
            showfliers=False
        )
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
            patch.set_edgecolor('black')

        np.random.seed(42)
        for rep_idx, data in enumerate(data_by_rep):
            x_jitter = np.random.normal(rep_idx + 1, 0.06, size=len(data))
            ax.scatter(x_jitter, data, alpha=0.5, s=40, color='black',
                       zorder=3, edgecolors='black', linewidth=0.5)

        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'Replicate {idx + 1}\n(n={len(data)})'
                            for idx, data in enumerate(data_by_rep)],
                           fontsize=11, fontweight='bold')
        ax.set_ylabel(metric_labels[metric], fontsize=12, fontweight='bold')
        ax.set_title(metric_labels[metric], fontsize=13, fontweight='bold', pad=15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_key_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_key_metrics.svg', format='svg', bbox_inches='tight')
    plt.close()

    fig3, axes3 = plt.subplots(2, 3, figsize=(16, 10))
    axes3 = axes3.flatten()

    for idx, metric in enumerate(key_metrics):
        ax = axes3[idx]

        data_by_rep = [rep[metric].dropna() for rep in replicates]

        x_pos = np.arange(1, n_reps + 1)
        bp = ax.boxplot(
            data_by_rep,
            positions=x_pos,
            widths=0.6,
            patch_artist=True,
            showfliers=False
        )
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
            patch.set_edgecolor('black')

        np.random.seed(42)
        for rep_idx, data in enumerate(data_by_rep):
            x_jitter = np.random.normal(rep_idx + 1, 0.06, size=len(data))
            ax.scatter(x_jitter, data, alpha=0.4, s=25, color='black', zorder=3)

        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'Replicate {idx + 1}\n(n={len(data)})'
                            for idx, data in enumerate(data_by_rep)],
                           fontsize=11, fontweight='bold')
        ax.set_ylabel(metric_labels[metric], fontsize=12, fontweight='bold')
        ax.set_title(metric_labels[metric], fontsize=13, fontweight='bold', pad=15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_key_metrics_boxplots.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_key_metrics_boxplots.svg', format='svg', bbox_inches='tight')
    plt.close()


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Compare mitochondrial form factor between two or three replicates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (with GUI):
  python compare_form_factor_replicates.py

Example (with CLI):
  python compare_form_factor_replicates.py \\
    --rep1 Rep1_mitochondria_summary_simple.csv \\
    --rep2 Rep2_mitochondria_summary_simple.csv \\
    --rep3 Rep3_mitochondria_summary_simple.csv \\
    --out results/form_factor_comparison
        """
    )
    parser.add_argument("--rep1",
                        help="Path to replicate 1 CSV file (optional, will use GUI if not provided)")
    parser.add_argument("--rep2",
                        help="Path to replicate 2 CSV file (optional, will use GUI if not provided)")
    parser.add_argument("--rep3",
                        help="Path to replicate 3 CSV file (optional)")
    parser.add_argument("--out",
                        help="Output prefix (optional, will use GUI if not provided)")

    args = parser.parse_args()

    # Get file paths from CLI or GUI picker
    if args.rep1 and args.rep2 and args.out:
        rep_paths = [args.rep1, args.rep2]
        if args.rep3:
            rep_paths.append(args.rep3)
        out_path = args.out
    else:
        rep_paths, out_path = pick_files_if_needed(args)

    out_dir = Path(out_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("MITOCHONDRIAL FORM FACTOR REPLICATE COMPARISON")
    print("=" * 100)

    for idx, rep_path in enumerate(rep_paths, start=1):
        print(f"\nReplicate {idx}: {Path(rep_path).name}")
    print(f"Output prefix: {out_path}\n")

    print("Loading replicates...")
    replicates, combined = load_replicates(rep_paths)

    for idx, rep in enumerate(replicates, start=1):
        print(f"Replicate {idx}: n={len(rep)} cells")
    print(f"Total: n={len(combined)} cells\n")

    metrics = [
        'organelle_area_raw_count',
        'organelle_area_raw_mean',
        'organelle_area_raw_std',
        'form_factor_mean',
        'form_factor_std',
        'form_factor_median',
        'aspect_ratio_calc_mean',
        'aspect_ratio_calc_std',
        'circularity_mean',
        'circularity_std'
    ]

    print("Performing statistical tests...")
    results_df = perform_statistical_tests(replicates, metrics)

    results_path = f"{out_path}_statistics.csv"
    results_df.to_csv(results_path, index=False)
    print(f"Statistical results saved: {results_path}\n")

    print("=" * 100)
    print("STATISTICAL COMPARISON SUMMARY")
    print("=" * 100)
    for _, row in results_df.iterrows():
        print(f"\n{row['Metric']}:")
        for idx in range(len(replicates)):
            mean = row.get(f"Rep{idx + 1}_mean")
            sem = row.get(f"Rep{idx + 1}_SEM")
            n_val = row.get(f"Rep{idx + 1}_n")
            mean_text = "nan" if pd.isna(mean) else f"{mean:.3f}"
            sem_text = "nan" if pd.isna(sem) else f"{sem:.3f}"
            n_text = "0" if pd.isna(n_val) else f"{int(n_val)}"
            print(f"  Replicate {idx + 1}: {mean_text} +/- {sem_text} (n={n_text})")
        p_text = "n/a" if pd.isna(row['p_value']) else f"{row['p_value']:.4f}"
        print(f"  {row['Test']}: p = {p_text} {row['Significance']}")

    print("\n" + "=" * 100)
    print("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
    print("=" * 100 + "\n")

    print("Creating visualizations...")
    create_comparison_plots(replicates, metrics, METRIC_LABELS, out_path)
    print(f"Plots saved: {out_path}_all_metrics_boxplots.[png/svg]")
    print(f"Plots saved: {out_path}_key_metrics_boxplots.[png/svg]")

    print("\n" + "=" * 100)
    print("Analysis complete!")
    print("=" * 100)


if __name__ == "__main__":
    main()
