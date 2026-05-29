"""
Mitochondrial mass comparison

Species:
  - Harpagifer antarcticus
  - Lipophrys pholis        

Metric: total_mito_area_raw (µm²) per cell

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

HARP_FILES = {
    "Rep 1":  "20250805_Skin_488_SYTO24_561_Mito_Orange_single_plane.csv",
    "Rep 2":  "20250915_Skin_488_SYTO24_561_Mito_Orange_zstack.csv",
    "Rep 3":  "20251212_Skin_488_SYTO23_561_Mito_Orange_z_project.csv",
}
LPHOLIS_FILE = "Lpholis_cell_metrics.csv"
COL = "total_mito_area_raw"

# Load Harpagifer replicates
harp_dfs = {}
for label, path in HARP_FILES.items():
    df = pd.read_csv(path, encoding="utf-8-sig")
    harp_dfs[label] = df[COL].dropna().values

# Pool all Harpagifer cells for stats / box
harp_all = np.concatenate(list(harp_dfs.values()))

# Load L. pholis
lph_df = pd.read_csv(LPHOLIS_FILE, encoding="utf-8-sig")
lph_all = lph_df[COL].dropna().values

# Descriptive statistics
def describe(arr, label=""):
    n = len(arr)
    return {
        "n":      n,
        "mean":   np.mean(arr),
        "sem":    stats.sem(arr),
        "median": np.median(arr),
        "sd":     np.std(arr, ddof=1),
        "q25":    np.percentile(arr, 25),
        "q75":    np.percentile(arr, 75),
    }

print("=" * 60)
print("  MITOCHONDRIAL MASS — total_mito_area_raw (µm²) per cell")
print("=" * 60)

for label, arr in {**harp_dfs, "H. antarcticus (all)": harp_all,
                   "L. pholis": lph_all}.items():
    s = describe(arr)
    sw_p = stats.shapiro(arr).pvalue if len(arr) <= 5000 else float("nan")
    print(f"\n{label}")
    print(f"  n          = {s['n']}")
    print(f"  Mean ± SEM = {s['mean']:.2f} ± {s['sem']:.2f} µm²")
    print(f"  Median     = {s['median']:.2f} µm²")
    print(f"  SD         = {s['sd']:.2f} µm²")
    print(f"  IQR        = [{s['q25']:.2f}, {s['q75']:.2f}]")
    print(f"  Shapiro-Wilk p = {sw_p:.4f}")

# Statistical test (Harpagifer pooled vs L. pholis)
sw_h = stats.shapiro(harp_all).pvalue
sw_l = stats.shapiro(lph_all).pvalue
both_normal = sw_h > 0.05 and sw_l > 0.05
lev_p = stats.levene(harp_all, lph_all).pvalue
equal_var = lev_p > 0.05

if both_normal:
    stat, p_val = stats.ttest_ind(harp_all, lph_all, equal_var=equal_var)
    test_name = "Student's t-test" if equal_var else "Welch's t-test"
    stat_label = "t"
else:
    stat, p_val = stats.mannwhitneyu(harp_all, lph_all, alternative="two-sided")
    test_name = "Mann-Whitney U"
    stat_label = "U"

pooled_sd = np.sqrt((np.std(harp_all, ddof=1)**2 + np.std(lph_all, ddof=1)**2) / 2)
cohens_d  = (np.mean(harp_all) - np.mean(lph_all)) / pooled_sd

print(f"\n{'─'*60}")
print(f"Levene's test p = {lev_p:.4f}  (equal variances: {equal_var})")
print(f"\n{test_name}  ({stat_label} = {stat:.4f},  p = {p_val:.4f})")
print(f"Cohen's d = {cohens_d:.3f}")
print("=" * 60)

REP_COLORS  = ["#E07B54", "#C94277", "#7B4FA6"]   
BOX_COLORS  = {"H. antarcticus": "#4C72B0", "L. pholis": "#DDAA77"}

fig, ax = plt.subplots(figsize=(7, 6))

plot_data  = [harp_all, lph_all]
positions  = [1, 2]
box_labels = ["H. antarcticus", "L. pholis"]

bp = ax.boxplot(
    plot_data,
    positions=positions,
    widths=0.45,
    patch_artist=True,
    notch=False,
    showfliers=False,          
    medianprops=dict(color="black", linewidth=2),
    whiskerprops=dict(linewidth=1.2),
    capprops=dict(linewidth=1.2),
    boxprops=dict(linewidth=1.2),
)

for patch, label in zip(bp["boxes"], box_labels):
    patch.set_facecolor(BOX_COLORS[label])
    patch.set_alpha(0.55)

rng = np.random.default_rng(42)
for (rep_label, arr), color in zip(harp_dfs.items(), REP_COLORS):
    jitter = rng.uniform(-0.18, 0.18, size=len(arr))
    ax.scatter(1 + jitter, arr,
               color=color, edgecolors="white", linewidths=0.4,
               s=30, alpha=0.85, zorder=4)

jitter = rng.uniform(-0.18, 0.18, size=len(lph_all))
ax.scatter(2 + jitter, lph_all,
           color=BOX_COLORS["L. pholis"], edgecolors="white", linewidths=0.4,
           s=30, alpha=0.75, zorder=4)

y_top  = max(np.max(harp_all), np.max(lph_all))
y_line = y_top * 1.07
y_text = y_top * 1.11

ax.plot([1, 1, 2, 2],
        [y_line * 0.97, y_line, y_line, y_line * 0.97],
        color="black", linewidth=1.2)

def p_to_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

ax.text(1.5, y_text, p_to_stars(p_val),
        ha="center", va="bottom", fontsize=13)

annot = (
    f"{test_name}\n"
    f"p = {p_val:.4f}\n"
    f"Cohen's d = {cohens_d:.2f}"
)
ax.text(0.97, 0.97, annot, transform=ax.transAxes,
        fontsize=8, va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="grey", alpha=0.85))

for pos, arr in zip(positions, plot_data):
    ax.text(pos, -0.035 * y_top, f"n = {len(arr)}",
            ha="center", va="top", fontsize=9)

rep_patches = [
    mpatches.Patch(color=c, label=lbl, alpha=0.85)
    for lbl, c in zip(harp_dfs.keys(), REP_COLORS)
]
ax.legend(handles=rep_patches, title="H. antarcticus replicates",
          fontsize=8, title_fontsize=8.5,
          loc="upper left", framealpha=0.85)

ax.set_xlim(0.4, 2.6)
ax.set_xticks(positions)
ax.set_xticklabels([f"$\\it{{{sp}}}$" for sp in box_labels], fontsize=11)
ax.set_ylabel("Total mitochondrial area per cell (µm²)", fontsize=11)
ax.set_title("Mitochondrial mass — skin cells", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
out_path = "Skin_comparison_total_mass_.svg"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"\nFigure saved → {out_path}")
plt.show()
