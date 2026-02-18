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
plt.rcParams["figure.dpi"] = 300


def p_to_sig(p: float) -> str:
    """Convert p-value to significance stars."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def add_sig_bracket(ax, x1: float, x2: float, y: float, h: float, text: str, fontsize: int = 12):
    """Draw a significance bracket between x1 and x2 at height y."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c="black", clip_on=False)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=fontsize)


def finite_minmax(arr: np.ndarray):
    """Return (min,max) of finite values, or (None,None) if none."""
    a = np.asarray(arr, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return None, None
    return float(np.min(a)), float(np.max(a))


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
            "You can select multiple files at once.",
        )
        files = filedialog.askopenfilenames(
            title="Select Species 1 Replicate CSV files",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
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
            "You can select multiple files at once.",
        )
        files = filedialog.askopenfilenames(
            title="Select Species 2 Replicate CSV files",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not files:
            messagebox.showerror("Missing files", "No Species 2 files selected.")
            raise SystemExit(1)
        species2_paths = list(files)

    # Get output path
    if args.out:
        out_path = args.out
    else:
        out_dir = filedialog.askdirectory(title="Select output folder for results")
        if not out_dir:
            messagebox.showerror("Missing folder", "No output folder selected.")
            raise SystemExit(1)
        out_path = str(Path(out_dir) / "species_comparison_cell_level")

    root.update()
    root.destroy()
    return species1_paths, species2_paths, out_path


def load_and_label_data(files, species_name, replicate_prefix="Rep"):
    """Load data and add species and replicate labels"""
    all_data = []
    for i, file in enumerate(files, start=1):
        df = pd.read_csv(file)
        df["Species"] = species_name
        df["Replicate"] = f"{replicate_prefix}{i}"
        df["Replicate_ID"] = f"{species_name}_{replicate_prefix}{i}"
        all_data.append(df)
    return pd.concat(all_data, ignore_index=True)


def calculate_replicate_summaries(all_data, parameters):
    """Calculate replicate-level summaries for statistical testing"""
    replicate_summaries = []
    for param in parameters:
        for species in all_data["Species"].unique():
            species_data = all_data[all_data["Species"] == species]
            for rep in species_data["Replicate"].unique():
                rep_data = species_data[species_data["Replicate"] == rep]
                replicate_summaries.append(
                    {
                        "Parameter": param,
                        "Species": species,
                        "Replicate": rep,
                        "Replicate_ID": f"{species}_{rep}",
                        "Median": rep_data[param].median(),
                        "Mean": rep_data[param].mean(),
                        "n_cells": len(rep_data),
                    }
                )
    return pd.DataFrame(replicate_summaries)


def perform_statistical_tests(replicate_df, parameters):
    """
    Perform statistical tests at the replicate level.
    Uses Mann-Whitney U test (appropriate for small n).
    """
    results = []

    species_list = sorted(replicate_df["Species"].unique())
    if len(species_list) != 2:
        raise ValueError(f"Expected 2 species, got {len(species_list)}")

    species1, species2 = species_list

    for param in parameters:
        param_data = replicate_df[replicate_df["Parameter"] == param]
        species1_values = param_data[param_data["Species"] == species1]["Median"].values
        species2_values = param_data[param_data["Species"] == species2]["Median"].values

        u_stat, p_value = stats.mannwhitneyu(species1_values, species2_values, alternative="two-sided")

        species1_median = np.median(species1_values)
        species2_median = np.median(species2_values)
        percent_diff = ((species2_median - species1_median) / species1_median) * 100

        sig = p_to_sig(float(p_value))

        results.append(
            {
                "Parameter": param,
                f"{species1}_median": species1_median,
                f"{species1}_n_replicates": len(species1_values),
                f"{species2}_median": species2_median,
                f"{species2}_n_replicates": len(species2_values),
                "Percent_Difference": percent_diff,
                "Mann_Whitney_U": u_stat,
                "p_value": p_value,
                "Significance": sig,
            }
        )

    return pd.DataFrame(results)


def create_cell_level_boxplots(all_data, parameters, param_labels, output_prefix, species_display, stats_df):
    """
    Create boxplots showing all individual cells, color-coded by replicate.
    Adds replicate-level p-value + significance bracket.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    species_list = sorted(all_data["Species"].unique())
    color_maps = {
        species_list[0]: ["#1f77b4", "#9467bd", "#17becf", "#8c564b", "#e377c2"],
        species_list[1]: ["#ff7f0e", "#d62728", "#ff9896", "#bcbd22", "#ffbb78"],
    }

    for idx, param in enumerate(parameters):
        ax = axes[idx]

        # boxplots
        bp_data = [all_data[all_data["Species"] == sp][param].values for sp in species_list]
        ax.boxplot(
            bp_data,
            positions=[1, 2],
            widths=0.5,
            patch_artist=True,
            tick_labels=[species_display[sp] for sp in species_list],
            boxprops=dict(facecolor="lightgray", alpha=0.3),
            medianprops=dict(color="black", linewidth=2.5),
            showfliers=False,
        )

        # overlay points (per replicate)
        for species_idx, species in enumerate(species_list):
            species_data = all_data[all_data["Species"] == species]
            replicates = sorted(species_data["Replicate"].unique())

            for rep_idx, rep in enumerate(replicates):
                rep_data = species_data[species_data["Replicate"] == rep][param].values
                np.random.seed(42 + rep_idx)
                x_jitter = np.random.normal(species_idx + 1, 0.08, size=len(rep_data))
                color = color_maps[species][rep_idx]

                ax.scatter(
                    x_jitter,
                    rep_data,
                    s=30,
                    alpha=0.6,
                    color=color,
                    edgecolors="black",
                    linewidths=0.5,
                    zorder=5,
                    label=f"{species_display[species]} {rep}" if idx == 0 else "",
                )

        ax.set_ylabel(param_labels[param], fontsize=12, fontweight="bold")
        ax.set_xlabel("Species", fontsize=12, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        species_counts = {sp: len(all_data[all_data["Species"] == sp]) for sp in species_list}
        info_text = "\n".join([f"{species_display[sp]}: n={species_counts[sp]} cells" for sp in species_list])
        ax.text(
            0.02,
            0.98,
            info_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.4),
        )

        try:
            row = stats_df.loc[stats_df["Parameter"] == param].iloc[0]
            pval = float(row["p_value"])
            sig = str(row.get("Significance", p_to_sig(pval)))

            vals = all_data[param].to_numpy(dtype=float)
            y_min, y_max = finite_minmax(vals)

            if y_min is not None and y_max is not None:
                yrng = y_max - y_min
                if not np.isfinite(yrng) or yrng == 0:
                    yrng = 1.0

                y = y_max + 0.06 * yrng
                h = 0.02 * yrng
                label = f"{sig}\n" + (f"p={pval:.3g}" if pval >= 0.001 else "p<0.001")

                add_sig_bracket(ax, 1, 2, y, h, label, fontsize=11)

                # extend ylim safely
                cur_bottom, cur_top = ax.get_ylim()
                new_top = max(cur_top, y + 0.12 * yrng)
                if np.isfinite(new_top):
                    ax.set_ylim(cur_bottom, new_top)
        except Exception:
            # If anything goes wrong, don't crash—just skip annotations for this panel.
            pass

    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.9, title="Replicate", title_fontsize=9)

    plt.suptitle(
        "Mitochondrial Morphology: Cell-Level Comparison\nEach dot = one cell (color-coded by replicate)",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_cell_level.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_prefix}_cell_level.svg", format="svg", bbox_inches="tight")
    plt.close()
    print(f"✓ Saved: {output_prefix}_cell_level.[png/svg]")


def create_split_violin_plots(all_data, parameters, param_labels, output_prefix, species_display, stats_df):
    """
    Create violin plots showing distribution with individual replicates overlaid.
    Adds replicate-level p-value + significance bracket.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    species_list = sorted(all_data["Species"].unique())
    color_maps = {
        species_list[0]: ["#1f77b4", "#9467bd", "#17becf", "#8c564b", "#e377c2"],
        species_list[1]: ["#ff7f0e", "#d62728", "#ff9896", "#bcbd22", "#ffbb78"],
    }

    for idx, param in enumerate(parameters):
        ax = axes[idx]

        parts = ax.violinplot(
            [all_data[all_data["Species"] == sp][param].values for sp in species_list],
            positions=[1, 2],
            widths=0.7,
            showmeans=False,
            showmedians=True,
        )

        for pc, species in zip(parts["bodies"], species_list):
            pc.set_facecolor("lightblue" if species == species_list[0] else "lightsalmon")
            pc.set_alpha(0.3)

        for species_idx, species in enumerate(species_list):
            species_data = all_data[all_data["Species"] == species]
            replicates = sorted(species_data["Replicate"].unique())

            for rep_idx, rep in enumerate(replicates):
                rep_data = species_data[species_data["Replicate"] == rep][param].values
                np.random.seed(42 + rep_idx)
                x_jitter = np.random.normal(species_idx + 1, 0.05, size=len(rep_data))
                color = color_maps[species][rep_idx]

                ax.scatter(
                    x_jitter,
                    rep_data,
                    s=25,
                    alpha=0.5,
                    color=color,
                    edgecolors="black",
                    linewidths=0.3,
                    zorder=5,
                )

        ax.set_xticks([1, 2])
        ax.set_xticklabels([species_display[sp] for sp in species_list], fontsize=12, fontweight="bold")
        ax.set_ylabel(param_labels[param], fontsize=12, fontweight="bold")
        ax.set_title(param_labels[param], fontsize=13, fontweight="bold", pad=10)
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # --- robust significance annotation ---
        try:
            row = stats_df.loc[stats_df["Parameter"] == param].iloc[0]
            pval = float(row["p_value"])
            sig = str(row.get("Significance", p_to_sig(pval)))

            vals = all_data[param].to_numpy(dtype=float)
            y_min, y_max = finite_minmax(vals)

            if y_min is not None and y_max is not None:
                yrng = y_max - y_min
                if not np.isfinite(yrng) or yrng == 0:
                    yrng = 1.0

                y = y_max + 0.06 * yrng
                h = 0.02 * yrng
                label = f"{sig}\n" + (f"p={pval:.3g}" if pval >= 0.001 else "p<0.001")
                add_sig_bracket(ax, 1, 2, y, h, label, fontsize=11)

                cur_bottom, cur_top = ax.get_ylim()
                new_top = max(cur_top, y + 0.12 * yrng)
                if np.isfinite(new_top):
                    ax.set_ylim(cur_bottom, new_top)
        except Exception:
            pass

    plt.suptitle(
        "Mitochondrial Morphology: Distribution View\nViolin = overall distribution; dots = individual cells by replicate",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_violin.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_prefix}_violin.svg", format="svg", bbox_inches="tight")
    plt.close()
    print(f"✓ Saved: {output_prefix}_violin.[png/svg]")


def create_replicate_breakdown(all_data, parameters, param_labels, output_prefix, species_display):
    """
    Create separate boxplots for each replicate to clearly show replicate-level variation.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    axes = axes.flatten()

    species_list = sorted(all_data["Species"].unique())

    for idx, param in enumerate(parameters):
        ax = axes[idx]

        all_rep_ids = []
        positions = []
        pos = 1

        for species in species_list:
            species_data = all_data[all_data["Species"] == species]
            replicates = sorted(species_data["Replicate"].unique())

            for rep in replicates:
                all_rep_ids.append((species, rep))
                positions.append(pos)
                pos += 1
            pos += 0.5

        bp_data = []
        labels = []
        colors = []
        color_map = {species_list[0]: "#4A90E2", species_list[1]: "#E88D3A"}

        for species, rep in all_rep_ids:
            rep_data = all_data[(all_data["Species"] == species) & (all_data["Replicate"] == rep)][param].values
            bp_data.append(rep_data)
            labels.append(f"{species_display[species]}\n{rep}\n(n={len(rep_data)})")
            colors.append(color_map[species])

        bp = ax.boxplot(bp_data, positions=positions, widths=0.6, patch_artist=True, showfliers=False)

        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
            patch.set_edgecolor("black")

        for (species, rep), p in zip(all_rep_ids, positions):
            rep_data = all_data[(all_data["Species"] == species) & (all_data["Replicate"] == rep)][param].values
            np.random.seed(42)
            x_jitter = np.random.normal(p, 0.08, size=len(rep_data))
            ax.scatter(x_jitter, rep_data, s=30, alpha=0.4, color="black", zorder=5)

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(param_labels[param], fontsize=12, fontweight="bold")
        ax.set_title(param_labels[param], fontsize=13, fontweight="bold", pad=10)
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if len(species_list) == 2:
            n_reps_species1 = len(all_data[all_data["Species"] == species_list[0]]["Replicate"].unique())
            sep_pos = positions[n_reps_species1 - 1] + 0.75
            ax.axvline(sep_pos, color="gray", linestyle="--", alpha=0.5, linewidth=1.5)

    plt.suptitle("Mitochondrial Morphology: Replicate Breakdown\nEach box = one replicate",
                 fontsize=14, fontweight="bold", y=0.995)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_replicate_breakdown.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_prefix}_replicate_breakdown.svg", format="svg", bbox_inches="tight")
    plt.close()
    print(f"✓ Saved: {output_prefix}_replicate_breakdown.[png/svg]")


def main():
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
        """,
    )
    parser.add_argument("--species1", nargs="+", help="Paths to species 1 replicate CSV files")
    parser.add_argument("--species2", nargs="+", help="Paths to species 2 replicate CSV files")
    parser.add_argument("--name1", default="Species1", help="Name for species 1 (default: Species1)")
    parser.add_argument("--name2", default="Species2", help="Name for species 2 (default: Species2)")
    parser.add_argument("--out", help="Output prefix (optional, will use GUI if not provided)")
    args = parser.parse_args()

    if args.species1 and args.species2 and args.out:
        species1_paths = args.species1
        species2_paths = args.species2
        out_path = args.out
    else:
        species1_paths, species2_paths, out_path = pick_files_if_needed(args)

    # Species 1 -> H. antarcticus; Species 2 -> L. pholis (italic in plots)
    species_display = {
        args.name1: r"$\it{H.\ antarcticus}$",
        args.name2: r"$\it{L.\ pholis}$",
    }

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

    print("Loading data...")
    species1_data = load_and_label_data(species1_paths, args.name1)
    species2_data = load_and_label_data(species2_paths, args.name2)
    all_data = pd.concat([species1_data, species2_data], ignore_index=True)

    print(f"{args.name1}: {len(species1_data)} cells across {len(species1_paths)} replicates")
    for i in range(len(species1_paths)):
        n = len(species1_data[species1_data["Replicate"] == f"Rep{i + 1}"])
        print(f"  Rep{i + 1}: {n} cells")

    print(f"\n{args.name2}: {len(species2_data)} cells across {len(species2_paths)} replicates")
    for i in range(len(species2_paths)):
        n = len(species2_data[species2_data["Replicate"] == f"Rep{i + 1}"])
        print(f"  Rep{i + 1}: {n} cells")

    parameters = ["organelle_area_raw_mean", "form_factor_mean", "aspect_ratio_calc_mean", "circularity_mean"]

    param_labels = {
        "organelle_area_raw_mean": "Mitochondrial Area (μm²)",
        "form_factor_mean": "Form Factor",
        "aspect_ratio_calc_mean": "Aspect Ratio",
        "circularity_mean": "Circularity",
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

    print("\nCreating cell-level visualizations...")
    create_cell_level_boxplots(all_data, parameters, param_labels, out_path, species_display, stats_df)
    create_split_violin_plots(all_data, parameters, param_labels, out_path, species_display, stats_df)
    create_replicate_breakdown(all_data, parameters, param_labels, out_path, species_display)

    print("\n" + "=" * 100)
    print("ANALYSIS COMPLETE!")
    print("=" * 100)


if __name__ == "__main__":
    main()
