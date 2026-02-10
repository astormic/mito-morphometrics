#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Mitochondrial Morphology Between Species
==========================================
Shows all individual cells color-coded by replicate to visualise both
overall distribution and replicate consistency.

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

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox

    HAS_GUI = True
except ImportError:
    HAS_GUI = False
    print("Note: tkinter not available. GUI file selection disabled.")
    print("Please use command line arguments instead.\n")

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def pick_files_if_needed(args):
    """
    Open GUI dialogs to select species replicate files and output location if not provided via CLI.
    """
    if not HAS_GUI:
        print("\nError: GUI not available (tkinter not installed).")
        print("Please provide files using command line arguments.\n")
        raise SystemExit(1)

    root = tk.Tk()
    root.withdraw()

    # Get species 1 files
    species1_paths = []
    if args.species1:
        species1_paths = args.species1
    else:
        messagebox.showinfo(
            "Species 1 Files",
            "Select all replicate CSV files for Species 1 (e.g., L_Fin)\n\n"
            "You can select multiple files at once."
        )
        files = filedialog.askopenfilenames(
            title="Select Species 1 Replicate CSV files",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not files:
            messagebox.showerror("Missing files", "No Species 1 files selected.")
            raise SystemExit(1)
        species1_paths = list(files)

    # Get species 2 files
    species2_paths = []
    if args.species2:
        species2_paths = args.species2
    else:
        messagebox.showinfo(
            "Species 2 Files",
            "Select all replicate CSV files for Species 2 (e.g., H_Fin)\n\n"
            "You can select multiple files at once."
        )
        files = filedialog.askopenfilenames(
            title="Select Species 2 Replicate CSV files",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not files:
            messagebox.showerror("Missing files", "No Species 2 files selected.")
            raise SystemExit(1)
        species2_paths = list(files)

    # Get output path
    if args.out:
        out_path = args.out
    else:
        out_dir = filedialog.askdirectory(
            title="Select output folder for results"
        )
        if not out_dir:
            messagebox.showerror("Missing folder", "No output folder selected.")
            raise SystemExit(1)
        out_path = str(Path(out_dir) / "species_comparison_cell_level")

    root.update()
    root.destroy()
    return species1_paths, species2_paths, out_path


def load_and_label_data(files, species_name, replicate_prefix='Rep'):
    """Load data and add species and replicate labels"""
    all_data = []
    for i, file in enumerate(files, start=1):
        df = pd.read_csv(file)
        df['Species'] = species_name
        df['Replicate'] = f'{replicate_prefix}{i}'
        df['Replicate_ID'] = f"{species_name}_{replicate_prefix}{i}"
        all_data.append(df)
    return pd.concat(all_data, ignore_index=True)


def calculate_replicate_summaries(all_data, parameters):
    """Calculate replicate-level summaries for statistical testing"""
    replicate_summaries = []
    for param in parameters:
        for species in all_data['Species'].unique():
            species_data = all_data[all_data['Species'] == species]
            for rep in species_data['Replicate'].unique():
                rep_data = species_data[species_data['Replicate'] == rep]
                replicate_summaries.append({
                    'Parameter': param,
                    'Species': species,
                    'Replicate': rep,
                    'Replicate_ID': f"{species}_{rep}",
                    'Median': rep_data[param].median(),
                    'Mean': rep_data[param].mean(),
                    'n_cells': len(rep_data)
                })
    return pd.DataFrame(replicate_summaries)


def perform_statistical_tests(replicate_df, parameters):
    """
    Perform statistical tests at the replicate level.
    Uses Mann-Whitney U test (appropriate for small n).
    """
    results = []

    species_list = sorted(replicate_df['Species'].unique())
    if len(species_list) != 2:
        raise ValueError(f"Expected 2 species, got {len(species_list)}")

    species1, species2 = species_list

    for param in parameters:
        param_data = replicate_df[replicate_df['Parameter'] == param]
        species1_values = param_data[param_data['Species'] == species1]['Median'].values
        species2_values = param_data[param_data['Species'] == species2]['Median'].values

        # Mann-Whitney U test (better for small samples)
        u_stat, p_value = stats.mannwhitneyu(species1_values, species2_values,
                                             alternative='two-sided')

        # Effect size (median difference)
        species1_median = np.median(species1_values)
        species2_median = np.median(species2_values)
        percent_diff = ((species2_median - species1_median) / species1_median) * 100

        # Significance level
        if p_value < 0.001:
            sig = '***'
        elif p_value < 0.01:
            sig = '**'
        elif p_value < 0.05:
            sig = '*'
        else:
            sig = 'ns'

        results.append({
            'Parameter': param,
            f'{species1}_median': species1_median,
            f'{species1}_n_replicates': len(species1_values),
            f'{species2}_median': species2_median,
            f'{species2}_n_replicates': len(species2_values),
            'Percent_Difference': percent_diff,
            'Mann_Whitney_U': u_stat,
            'p_value': p_value,
            'Significance': sig
        })

    return pd.DataFrame(results)


def create_cell_level_boxplots(all_data, parameters, param_labels, output_prefix):
    """
    Create boxplots showing all individual cells, color-coded by replicate.
    This shows both the overall distribution AND replicate consistency.
    """

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    species_list = sorted(all_data['Species'].unique())
    color_maps = {
        species_list[0]: ['#1f77b4', '#9467bd', '#17becf', '#8c564b', '#e377c2'],  # Blue, Purple, Cyan, Brown, Pink
        species_list[1]: ['#ff7f0e', '#d62728', '#ff9896', '#bcbd22', '#ffbb78']
        # Orange, Red, Light red, Yellow-green, Light orange
    }

    for idx, param in enumerate(parameters):
        ax = axes[idx]

        # Create boxplots
        bp_data = [all_data[all_data['Species'] == sp][param].values
                   for sp in species_list]
        bp = ax.boxplot(bp_data,
                        positions=[1, 2],
                        widths=0.5,
                        patch_artist=True,
                        tick_labels=species_list,
                        boxprops=dict(facecolor='lightgray', alpha=0.3),
                        medianprops=dict(color='black', linewidth=2.5),
                        showfliers=False)

        # Overlay individual cells, color-coded by replicate
        for species_idx, species in enumerate(species_list):
            species_data = all_data[all_data['Species'] == species]
            replicates = sorted(species_data['Replicate'].unique())

            for rep_idx, rep in enumerate(replicates):
                rep_data = species_data[species_data['Replicate'] == rep][param].values

                # Add jitter for visibility
                np.random.seed(42 + rep_idx)
                x_jitter = np.random.normal(species_idx + 1, 0.08, size=len(rep_data))

                # Color for this replicate
                color = color_maps[species][rep_idx]

                ax.scatter(x_jitter, rep_data,
                           s=30, alpha=0.6, color=color,
                           edgecolors='black', linewidths=0.5, zorder=5,
                           label=f'{species} {rep}' if idx == 0 else '')

        ax.set_ylabel(param_labels[param], fontsize=12, fontweight='bold')
        ax.set_xlabel('Species', fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add sample sizes
        species_counts = {sp: len(all_data[all_data['Species'] == sp])
                          for sp in species_list}
        info_text = '\n'.join([f'{sp}: n={species_counts[sp]} cells'
                               for sp in species_list])
        ax.text(0.02, 0.98, info_text,
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.4))

    # Create legend in the first subplot
    axes[0].legend(loc='upper right', fontsize=8, framealpha=0.9,
                   title='Replicate', title_fontsize=9)

    plt.suptitle('Mitochondrial Morphology: Cell-Level Comparison\n' +
                 'Each dot = one cell (color-coded by replicate)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_cell_level.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_cell_level.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_cell_level.[png/svg]")


def create_split_violin_plots(all_data, parameters, param_labels, output_prefix):
    """
    Create violin plots showing distribution with individual replicates overlaid.
    """

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    species_list = sorted(all_data['Species'].unique())

    # Color palettes - using distinct, saturated colors
    color_maps = {
        species_list[0]: ['#1f77b4', '#9467bd', '#17becf', '#8c564b', '#e377c2'],  # Blue, Purple, Cyan, Brown, Pink
        species_list[1]: ['#ff7f0e', '#d62728', '#ff9896', '#bcbd22', '#ffbb78']
        # Orange, Red, Light red, Yellow-green, Light orange
    }

    for idx, param in enumerate(parameters):
        ax = axes[idx]

        # Violin plots for overall distribution
        parts = ax.violinplot([all_data[all_data['Species'] == sp][param].values
                               for sp in species_list],
                              positions=[1, 2],
                              widths=0.7,
                              showmeans=False,
                              showmedians=True)

        # Color the violins
        for pc, species in zip(parts['bodies'], species_list):
            if species == species_list[0]:
                pc.set_facecolor('lightblue')
            else:
                pc.set_facecolor('lightsalmon')
            pc.set_alpha(0.3)

        # Overlay individual replicate data
        for species_idx, species in enumerate(species_list):
            species_data = all_data[all_data['Species'] == species]
            replicates = sorted(species_data['Replicate'].unique())

            for rep_idx, rep in enumerate(replicates):
                rep_data = species_data[species_data['Replicate'] == rep][param].values

                # Smaller jitter for violin plots
                np.random.seed(42 + rep_idx)
                x_jitter = np.random.normal(species_idx + 1, 0.05, size=len(rep_data))

                color = color_maps[species][rep_idx]

                ax.scatter(x_jitter, rep_data,
                           s=25, alpha=0.5, color=color,
                           edgecolors='black', linewidths=0.3, zorder=5)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(species_list, fontsize=12, fontweight='bold')
        ax.set_ylabel(param_labels[param], fontsize=12, fontweight='bold')
        ax.set_title(param_labels[param], fontsize=13, fontweight='bold', pad=10)
        ax.grid(axis='y', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('Mitochondrial Morphology: Distribution View\n' +
                 'Violin = overall distribution; dots = individual cells by replicate',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_violin.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_violin.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_violin.[png/svg]")


def create_replicate_breakdown(all_data, parameters, param_labels, output_prefix):
    """
    Create separate boxplots for each replicate to clearly show replicate-level variation.
    """

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    axes = axes.flatten()

    species_list = sorted(all_data['Species'].unique())

    for idx, param in enumerate(parameters):
        ax = axes[idx]

        # Get all unique replicate IDs across both species
        all_rep_ids = []
        positions = []
        pos = 1

        for species in species_list:
            species_data = all_data[all_data['Species'] == species]
            replicates = sorted(species_data['Replicate'].unique())

            for rep in replicates:
                all_rep_ids.append((species, rep))
                positions.append(pos)
                pos += 1
            pos += 0.5  # Add gap between species

        # Create boxplots for each replicate
        bp_data = []
        labels = []
        colors = []

        color_map = {species_list[0]: '#4A90E2', species_list[1]: '#E88D3A'}

        for species, rep in all_rep_ids:
            rep_data = all_data[(all_data['Species'] == species) &
                                (all_data['Replicate'] == rep)][param].values
            bp_data.append(rep_data)
            labels.append(f'{species}\n{rep}\n(n={len(rep_data)})')
            colors.append(color_map[species])

        bp = ax.boxplot(bp_data,
                        positions=positions,
                        widths=0.6,
                        patch_artist=True,
                        showfliers=False)

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
            patch.set_edgecolor('black')

        # Add individual points
        for i, ((species, rep), pos) in enumerate(zip(all_rep_ids, positions)):
            rep_data = all_data[(all_data['Species'] == species) &
                                (all_data['Replicate'] == rep)][param].values

            np.random.seed(42)
            x_jitter = np.random.normal(pos, 0.08, size=len(rep_data))

            ax.scatter(x_jitter, rep_data,
                       s=30, alpha=0.4, color='black', zorder=5)

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(param_labels[param], fontsize=12, fontweight='bold')
        ax.set_title(param_labels[param], fontsize=13, fontweight='bold', pad=10)
        ax.grid(axis='y', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add vertical line to separate species
        if len(species_list) == 2:
            n_reps_species1 = len(all_data[all_data['Species'] == species_list[0]]['Replicate'].unique())
            sep_pos = positions[n_reps_species1 - 1] + 0.75
            ax.axvline(sep_pos, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)

    plt.suptitle('Mitochondrial Morphology: Replicate Breakdown\n' +
                 'Each box = one replicate',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_replicate_breakdown.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_replicate_breakdown.svg', format='svg', bbox_inches='tight')
    plt.close()

    print(f"✓ Saved: {output_prefix}_replicate_breakdown.[png/svg]")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Compare mitochondrial morphology showing all cells color-coded by replicate",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example (with GUI):
  python compare_species_cell_level.py

Example (with CLI):
  python compare_species_cell_level.py \\
    --species1 L_Fin_rep1.csv L_Fin_rep2.csv L_Fin_rep3.csv \\
    --species2 H_Fin_rep1.csv H_Fin_rep2.csv \\
    --name1 L_Fin --name2 H_Fin \\
    --out results/species_comparison_cells
        """
    )
    parser.add_argument("--species1", nargs='+',
                        help="Paths to species 1 replicate CSV files")
    parser.add_argument("--species2", nargs='+',
                        help="Paths to species 2 replicate CSV files")
    parser.add_argument("--name1", default="Species1",
                        help="Name for species 1 (default: Species1)")
    parser.add_argument("--name2", default="Species2",
                        help="Name for species 2 (default: Species2)")
    parser.add_argument("--out",
                        help="Output prefix (optional, will use GUI if not provided)")

    args = parser.parse_args()

    # Get file paths from CLI or GUI picker
    if args.species1 and args.species2 and args.out:
        species1_paths = args.species1
        species2_paths = args.species2
        out_path = args.out
    else:
        species1_paths, species2_paths, out_path = pick_files_if_needed(args)

    # Create output directory
    out_dir = Path(out_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("MITOCHONDRIAL MORPHOLOGY SPECIES COMPARISON (CELL-LEVEL VISUALIZATION)")
    print("=" * 100)

    print(f"\n{args.name1} replicates ({len(species1_paths)} files):")
    for idx, path in enumerate(species1_paths, start=1):
        print(f"  Rep{idx}: {Path(path).name}")

    print(f"\n{args.name2} replicates ({len(species2_paths)} files):")
    for idx, path in enumerate(species2_paths, start=1):
        print(f"  Rep{idx}: {Path(path).name}")

    print(f"\nOutput prefix: {out_path}\n")

    # Load data
    print("Loading data...")
    species1_data = load_and_label_data(species1_paths, args.name1)
    species2_data = load_and_label_data(species2_paths, args.name2)
    all_data = pd.concat([species1_data, species2_data], ignore_index=True)

    print(f"{args.name1}: {len(species1_data)} cells across {len(species1_paths)} replicates")
    for i in range(len(species1_paths)):
        n = len(species1_data[species1_data['Replicate'] == f'Rep{i + 1}'])
        print(f"  Rep{i + 1}: {n} cells")

    print(f"\n{args.name2}: {len(species2_data)} cells across {len(species2_paths)} replicates")
    for i in range(len(species2_paths)):
        n = len(species2_data[species2_data['Replicate'] == f'Rep{i + 1}'])
        print(f"  Rep{i + 1}: {n} cells")

    # Parameters to analyze
    parameters = ['organelle_area_raw_mean', 'form_factor_mean',
                  'aspect_ratio_calc_mean', 'circularity_mean']

    param_labels = {
        'organelle_area_raw_mean': 'Mitochondrial Area (μm²)',
        'form_factor_mean': 'Form Factor',
        'aspect_ratio_calc_mean': 'Aspect Ratio',
        'circularity_mean': 'Circularity'
    }

    print("\nCalculating replicate-level summaries for statistical testing...")
    replicate_df = calculate_replicate_summaries(all_data, parameters)
    print("Performing statistical tests (replicate-level)...")
    stats_df = perform_statistical_tests(replicate_df, parameters)
    stats_path = f"{out_path}_statistics.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"Statistical results saved: {stats_path}")

    replicate_path = f"{out_path}_replicate_summaries.csv"
    replicate_df.to_csv(replicate_path, index=False)
    print(f"Replicate summaries saved: {replicate_path}")

    print("\n" + "=" * 100)
    print("STATISTICAL TEST RESULTS (Replicate-level analysis)")
    print("=" * 100)

    for _, row in stats_df.iterrows():
        print(f"\n{row['Parameter']}:")
        print(f"  {args.name1}: median = {row[f'{args.name1}_median']:.3f} "
              f"(n_replicates={int(row[f'{args.name1}_n_replicates'])})")
        print(f"  {args.name2}: median = {row[f'{args.name2}_median']:.3f} "
              f"(n_replicates={int(row[f'{args.name2}_n_replicates'])})")
        print(f"  Difference: {row['Percent_Difference']:+.1f}%")
        print(f"  Mann-Whitney U = {row['Mann_Whitney_U']:.2f}, "
              f"p = {row['p_value']:.4f} {row['Significance']}")

    print("\n" + "=" * 100)

    print("\nCreating cell-level visualizations...")
    create_cell_level_boxplots(all_data, parameters, param_labels, out_path)
    create_split_violin_plots(all_data, parameters, param_labels, out_path)
    create_replicate_breakdown(all_data, parameters, param_labels, out_path)

    print("\n" + "=" * 100)
    print("ANALYSIS COMPLETE!")
    print("=" * 100)
    print("\nGenerated files:")
    print(f"  1. {out_path}_statistics.csv - Statistical test results")
    print(f"  2. {out_path}_replicate_summaries.csv - Replicate-level data")
    print(f"  3. {out_path}_cell_level.[png/svg] - All cells color-coded by replicate")
    print(f"  4. {out_path}_violin.[png/svg] - Distribution view with replicates")
    print(f"  5. {out_path}_replicate_breakdown.[png/svg] - Each replicate shown separately")
    print("\nNOTE:")
    print("  These visualizations show ALL individual cells (not replicate summaries).")
    print("  Colors indicate which replicate each cell came from.")
    print("  Statistics are still calculated at the replicate level (proper approach).")
    print("=" * 100)


if __name__ == "__main__":
    main()
