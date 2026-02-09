#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistical Comparison of Mitochondrial Morphology Replicates/Species
Compares two or three sets of cell-level metrics with appropriate statistical tests.

Author: Amir Rahmani
License: MIT
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path


METRIC_LABELS = {
    'n_mitochondria': 'Number of Mitochondria',
    'total_mito_area_raw': 'Total Mitochondrial Area',
    'mean_mito_area_raw': 'Mean Mitochondrial Area',
    'median_mito_area_raw': 'Median Mitochondrial Area',
    'mean_mito_length_raw': 'Mean Mitochondrial Length',
    'max_mito_length_raw': 'Max Mitochondrial Length',
    'mean_solidity': 'Mean Solidity',
    'mean_extent': 'Mean Extent',
    'mean_aspect_ratio': 'Mean Aspect Ratio',
    'fragmentation_index': 'Fragmentation Index',
    'network_dominance_index': 'Network Dominance Index',
    'shape_heterogeneity_index': 'Shape Heterogeneity Index',
    'tubularity_score': 'Tubularity Score'
}


def load_replicates(rep_paths: list[str]) -> tuple:
    """Load and label replicate data from a list of CSV paths."""
    replicates = []
    for idx, path in enumerate(rep_paths, start=1):
        rep = pd.read_csv(path)
        rep['replicate'] = f'Replicate {idx}'
        replicates.append(rep)

    combined = pd.concat(replicates, ignore_index=True)

    return replicates, combined


def perform_statistical_tests(replicates: list[pd.DataFrame], metrics: list) -> pd.DataFrame:
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


def create_comparison_plots(replicates: list[pd.DataFrame], results_df: pd.DataFrame,
                            metrics: list, metric_labels: dict, output_prefix: str):
    """Create bar plot visualizations with statistical annotations."""
    
    n_metrics = len(metrics)
    n_cols = 3
    n_rows = int(np.ceil(n_metrics / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4*n_rows))
    axes = axes.flatten()
    
    n_reps = len(replicates)
    colors = sns.color_palette("Set2", n_reps)
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        
        data_by_rep = [rep[metric].dropna() for rep in replicates]

        means = [data.mean() for data in data_by_rep]
        sems = [data.sem() for data in data_by_rep]
        
        row = results_df[results_df['Metric'] == metric_labels[metric]].iloc[0]
        p_val = row['p_value']
        sig = row['Significance']
        
        # Create bar plot
        x_pos = np.arange(n_reps)
        bars = ax.bar(x_pos, means, yerr=sems, capsize=5,
                      color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        
        # Add individual data points
        np.random.seed(42)
        for rep_idx, data in enumerate(data_by_rep):
            x_jitter = np.random.normal(rep_idx, 0.04, size=len(data))
            ax.scatter(x_jitter, data, alpha=0.4, s=30, color='black', zorder=3)
        
        # Add significance bracket
        y_max = max(mean + sem for mean, sem in zip(means, sems))
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        bracket_height = y_max + 0.1 * y_range
        
        if n_reps == 2 and sig not in ['ns', 'na']:
            ax.plot([0, 0, 1, 1],
                    [bracket_height, bracket_height + 0.02*y_range,
                     bracket_height + 0.02*y_range, bracket_height],
                    'k-', linewidth=1.5)
            ax.text(0.5, bracket_height + 0.04*y_range, sig,
                    ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        # Add p-value
        if pd.isna(p_val):
            p_text = 'p = n/a'
        elif p_val < 0.001:
            p_text = 'p < 0.001'
        else:
            p_text = f'p = {p_val:.3f}'
        
        ax.text(0.5, -0.15, p_text, ha='center', va='top', 
                transform=ax.transAxes, fontsize=9, style='italic')
        
        # Formatting
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
    plt.savefig(f'{output_prefix}_all_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_all_metrics.svg', format='svg', bbox_inches='tight')
    plt.close()
    
    key_metrics = [
        'n_mitochondria',
        'mean_mito_area_raw',
        'mean_mito_length_raw',
        'fragmentation_index',
        'tubularity_score',
        'mean_aspect_ratio'
    ]
    
    fig2, axes2 = plt.subplots(2, 3, figsize=(16, 10))
    axes2 = axes2.flatten()
    
    for idx, metric in enumerate(key_metrics):
        ax = axes2[idx]
        
        data_by_rep = [rep[metric].dropna() for rep in replicates]

        means = [data.mean() for data in data_by_rep]
        sems = [data.sem() for data in data_by_rep]
        
        row = results_df[results_df['Metric'] == metric_labels[metric]].iloc[0]
        p_val = row['p_value']
        sig = row['Significance']
        
        x_pos = np.arange(n_reps)
        bars = ax.bar(x_pos, means, yerr=sems, capsize=8,
                      color=colors, alpha=0.8, edgecolor='black', linewidth=2)
        
        np.random.seed(42)
        for rep_idx, data in enumerate(data_by_rep):
            x_jitter = np.random.normal(rep_idx, 0.05, size=len(data))
            ax.scatter(x_jitter, data, alpha=0.5, s=50, color=colors[rep_idx],
                       zorder=3, edgecolors='black', linewidth=0.5)

        y_max = max(mean + sem for mean, sem in zip(means, sems))
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        bracket_height = y_max + 0.12 * y_range
        
        if n_reps == 2 and sig not in ['ns', 'na']:
            ax.plot([0, 0, 1, 1],
                    [bracket_height, bracket_height + 0.03*y_range,
                     bracket_height + 0.03*y_range, bracket_height],
                    'k-', linewidth=2)
            ax.text(0.5, bracket_height + 0.05*y_range, sig,
                    ha='center', va='bottom', fontsize=14, fontweight='bold')

        if pd.isna(p_val):
            p_text = 'p = n/a'
        elif p_val < 0.001:
            p_text = 'p < 0.001'
        else:
            p_text = f'p = {p_val:.3f}'
        
        ax.text(0.5, -0.18, p_text, ha='center', va='top', 
                transform=ax.transAxes, fontsize=11, style='italic', fontweight='bold')
        
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


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Compare two or three replicates of mitochondrial morphology data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python compare_replicates.py \\
    --rep1 replicate1_combined.csv \\
    --rep2 replicate2_combined.csv \\
    --rep3 replicate3_combined.csv \\
    --out results/comparison
        """
    )
    parser.add_argument("--rep1", required=True, 
                       help="Path to replicate 1 CSV file")
    parser.add_argument("--rep2", required=True, 
                       help="Path to replicate 2 CSV file")
    parser.add_argument("--rep3",
                       help="Optional path to replicate 3 CSV file")
    parser.add_argument("--out", required=True, 
                       help="Output prefix (e.g., results/comparison)")
    
    args = parser.parse_args()
    
    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading replicates...")
    rep_paths = [args.rep1, args.rep2]
    if args.rep3:
        rep_paths.append(args.rep3)
    replicates, combined = load_replicates(rep_paths)

    for idx, rep in enumerate(replicates, start=1):
        print(f"Replicate {idx}: n={len(rep)} cells")
    print(f"Total: n={len(combined)} cells\n")
    
    metrics = [
        'n_mitochondria',
        'total_mito_area_raw',
        'mean_mito_area_raw',
        'median_mito_area_raw',
        'mean_mito_length_raw',
        'max_mito_length_raw',
        'mean_solidity',
        'mean_extent',
        'mean_aspect_ratio',
        'fragmentation_index',
        'network_dominance_index',
        'shape_heterogeneity_index',
        'tubularity_score'
    ]
    
    print("Performing statistical tests...")
    results_df = perform_statistical_tests(replicates, metrics)
    
    results_path = f"{args.out}_statistics.csv"
    results_df.to_csv(results_path, index=False)
    print(f"Statistical results saved: {results_path}\n")
    
    print("="*100)
    print("STATISTICAL COMPARISON SUMMARY")
    print("="*100)
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
    
    print("\n" + "="*100)
    print("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
    print("="*100 + "\n")
    
    # Create plots
    print("Creating visualizations...")
    create_comparison_plots(replicates, results_df, metrics, METRIC_LABELS, args.out)
    print(f"Plots saved: {args.out}_all_metrics.[png/svg]")
    print(f"Plots saved: {args.out}_key_metrics.[png/svg]")
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
