#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistical Comparison of Mitochondrial Morphology Replicates
Compares two sets of cell-level metrics with appropriate statistical tests.

Author: [Your Name]
License: MIT
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path


def load_replicates(rep1_path: str, rep2_path: str) -> tuple:
    """Load and label replicate data."""
    rep1 = pd.read_csv(rep1_path)
    rep2 = pd.read_csv(rep2_path)
    
    rep1['replicate'] = 'Replicate 1'
    rep2['replicate'] = 'Replicate 2'
    
    combined = pd.concat([rep1, rep2], ignore_index=True)
    
    return rep1, rep2, combined


def perform_statistical_tests(rep1: pd.DataFrame, rep2: pd.DataFrame, 
                              metrics: list) -> pd.DataFrame:
    """
    Perform statistical tests on each metric.
    
    Uses Shapiro-Wilk test to assess normality, then:
    - Welch's t-test if both groups are normally distributed
    - Mann-Whitney U test otherwise
    
    Args:
        rep1: Replicate 1 DataFrame
        rep2: Replicate 2 DataFrame
        metrics: List of metric column names to test
    
    Returns:
        DataFrame with statistical results
    """
    results = []
    
    metric_labels = {
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
    
    for metric in metrics:
        data1 = rep1[metric].dropna()
        data2 = rep2[metric].dropna()
        
        # Calculate statistics
        mean1, sem1 = data1.mean(), data1.sem()
        mean2, sem2 = data2.mean(), data2.sem()
        
        # Normality tests (Shapiro-Wilk)
        _, p_norm1 = stats.shapiro(data1) if len(data1) >= 3 else (None, None)
        _, p_norm2 = stats.shapiro(data2) if len(data2) >= 3 else (None, None)
        
        # Choose appropriate test
        if p_norm1 and p_norm2 and p_norm1 > 0.05 and p_norm2 > 0.05:
            stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
            test_used = "Welch's t-test"
        else:
            stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            test_used = "Mann-Whitney U"
        
        # Significance stars
        if p_value < 0.001:
            sig = '***'
        elif p_value < 0.01:
            sig = '**'
        elif p_value < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        
        results.append({
            'Metric': metric_labels.get(metric, metric),
            'Rep1_mean': mean1,
            'Rep1_SEM': sem1,
            'Rep1_n': len(data1),
            'Rep2_mean': mean2,
            'Rep2_SEM': sem2,
            'Rep2_n': len(data2),
            'Test': test_used,
            'p_value': p_value,
            'Significance': sig
        })
    
    return pd.DataFrame(results)


def create_comparison_plots(rep1: pd.DataFrame, rep2: pd.DataFrame, 
                           results_df: pd.DataFrame, metrics: list,
                           metric_labels: dict, output_prefix: str):
    """Create bar plot visualizations with statistical annotations."""
    
    # Comprehensive plot with all metrics
    n_metrics = len(metrics)
    n_cols = 3
    n_rows = int(np.ceil(n_metrics / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4*n_rows))
    axes = axes.flatten()
    
    colors = ['#3498db', '#e74c3c']  # Blue for Rep1, Red for Rep2
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        
        data1 = rep1[metric].dropna()
        data2 = rep2[metric].dropna()
        
        means = [data1.mean(), data2.mean()]
        sems = [data1.sem(), data2.sem()]
        
        row = results_df[results_df['Metric'] == metric_labels[metric]].iloc[0]
        p_val = row['p_value']
        sig = row['Significance']
        
        # Create bar plot
        x_pos = [0, 1]
        bars = ax.bar(x_pos, means, yerr=sems, capsize=5, 
                       color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        
        # Add individual data points
        np.random.seed(42)
        x1_jitter = np.random.normal(0, 0.04, size=len(data1))
        x2_jitter = np.random.normal(1, 0.04, size=len(data2))
        
        ax.scatter(x1_jitter, data1, alpha=0.4, s=30, color='black', zorder=3)
        ax.scatter(x2_jitter, data2, alpha=0.4, s=30, color='black', zorder=3)
        
        # Add significance bracket
        y_max = max(means[0] + sems[0], means[1] + sems[1])
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        bracket_height = y_max + 0.1 * y_range
        
        if sig != 'ns':
            ax.plot([0, 0, 1, 1], 
                    [bracket_height, bracket_height + 0.02*y_range, 
                     bracket_height + 0.02*y_range, bracket_height], 
                    'k-', linewidth=1.5)
            ax.text(0.5, bracket_height + 0.04*y_range, sig, 
                    ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        # Add p-value
        if p_val < 0.001:
            p_text = 'p < 0.001'
        else:
            p_text = f'p = {p_val:.3f}'
        
        ax.text(0.5, -0.15, p_text, ha='center', va='top', 
                transform=ax.transAxes, fontsize=9, style='italic')
        
        # Formatting
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'Rep 1\n(n={len(data1)})', f'Rep 2\n(n={len(data2)})'])
        ax.set_ylabel(metric_labels[metric], fontsize=11, fontweight='bold')
        ax.set_title(metric_labels[metric], fontsize=12, fontweight='bold', pad=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Remove empty subplots
    for idx in range(n_metrics, len(axes)):
        fig.delaxes(axes[idx])
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_all_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_all_metrics.svg', format='svg', bbox_inches='tight')
    plt.close()
    
    # Key metrics plot
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
        
        data1 = rep1[metric].dropna()
        data2 = rep2[metric].dropna()
        
        means = [data1.mean(), data2.mean()]
        sems = [data1.sem(), data2.sem()]
        
        row = results_df[results_df['Metric'] == metric_labels[metric]].iloc[0]
        p_val = row['p_value']
        sig = row['Significance']
        
        x_pos = [0, 1]
        bars = ax.bar(x_pos, means, yerr=sems, capsize=8, 
                       color=colors, alpha=0.8, edgecolor='black', linewidth=2)
        
        # Individual points
        np.random.seed(42)
        x1_jitter = np.random.normal(0, 0.05, size=len(data1))
        x2_jitter = np.random.normal(1, 0.05, size=len(data2))
        
        ax.scatter(x1_jitter, data1, alpha=0.5, s=50, color='darkblue', 
                   zorder=3, edgecolors='black', linewidth=0.5)
        ax.scatter(x2_jitter, data2, alpha=0.5, s=50, color='darkred', 
                   zorder=3, edgecolors='black', linewidth=0.5)
        
        # Significance
        y_max = max(means[0] + sems[0], means[1] + sems[1])
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        bracket_height = y_max + 0.12 * y_range
        
        if sig != 'ns':
            ax.plot([0, 0, 1, 1], 
                    [bracket_height, bracket_height + 0.03*y_range, 
                     bracket_height + 0.03*y_range, bracket_height], 
                    'k-', linewidth=2)
            ax.text(0.5, bracket_height + 0.05*y_range, sig, 
                    ha='center', va='bottom', fontsize=14, fontweight='bold')
        
        if p_val < 0.001:
            p_text = 'p < 0.001'
        else:
            p_text = f'p = {p_val:.3f}'
        
        ax.text(0.5, -0.18, p_text, ha='center', va='top', 
                transform=ax.transAxes, fontsize=11, style='italic', fontweight='bold')
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'Replicate 1\n(n={len(data1)})', 
                            f'Replicate 2\n(n={len(data2)})'], 
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
        description="Compare two replicates of mitochondrial morphology data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python compare_replicates.py \\
    --rep1 replicate1_combined.csv \\
    --rep2 replicate2_combined.csv \\
    --out results/comparison
        """
    )
    parser.add_argument("--rep1", required=True, 
                       help="Path to replicate 1 CSV file")
    parser.add_argument("--rep2", required=True, 
                       help="Path to replicate 2 CSV file")
    parser.add_argument("--out", required=True, 
                       help="Output prefix (e.g., results/comparison)")
    
    args = parser.parse_args()
    
    # Create output directory if needed
    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("Loading replicates...")
    rep1, rep2, combined = load_replicates(args.rep1, args.rep2)
    
    print(f"Replicate 1: n={len(rep1)} cells")
    print(f"Replicate 2: n={len(rep2)} cells")
    print(f"Total: n={len(combined)} cells\n")
    
    # Define metrics
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
    
    metric_labels = {
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
    
    # Perform statistical tests
    print("Performing statistical tests...")
    results_df = perform_statistical_tests(rep1, rep2, metrics)
    
    # Save results
    results_path = f"{args.out}_statistics.csv"
    results_df.to_csv(results_path, index=False)
    print(f"✓ Statistical results saved: {results_path}\n")
    
    # Print summary
    print("="*100)
    print("STATISTICAL COMPARISON SUMMARY")
    print("="*100)
    for _, row in results_df.iterrows():
        print(f"\n{row['Metric']}:")
        print(f"  Replicate 1: {row['Rep1_mean']:.3f} ± {row['Rep1_SEM']:.3f} (n={row['Rep1_n']})")
        print(f"  Replicate 2: {row['Rep2_mean']:.3f} ± {row['Rep2_SEM']:.3f} (n={row['Rep2_n']})")
        print(f"  {row['Test']}: p = {row['p_value']:.4f} {row['Significance']}")
    
    print("\n" + "="*100)
    print("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
    print("="*100 + "\n")
    
    # Create plots
    print("Creating visualizations...")
    create_comparison_plots(rep1, rep2, results_df, metrics, metric_labels, args.out)
    print(f"✓ Plots saved: {args.out}_all_metrics.[png/svg]")
    print(f"✓ Plots saved: {args.out}_key_metrics.[png/svg]")
    
    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()
