"""
Mitochondrial comparison for cells cultured from the fin tissue
Metric reported: total_mito_area_raw (µm²) per cell.

Replicates
  H. antarcticus:
    Rep 1 — 20251127_Harpagifer_Fin_combined_cell_metrics_zstack.csv
    Rep 2 — 20251127_Harpagifer_Fin_combined_cell_metrics.csv
    Rep 3 — 20260128_Harpagifer_Fin_combined_cell_metrics.csv

  L. pholis:
    Rep 1 — 20251106_Lpholis_Fin_combined_cell_metrics.csv
    Rep 2 — 20251107_Lpholis_Fin_combined_cell_metrics.csv

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

COL = "total_mito_area_raw"

HARP_REPS = {
    "Rep 1": "20251127_Harpagifer_Fin_combined_cell_metrics_zstack.csv",
    "Rep 2": "20251127_Harpagifer_Fin_combined_cell_metrics.csv",
    "Rep 3": "20260128_Harpagifer_Fin_combined_cell_metrics.csv",
}

LPHOLIS_REPS = {
    "Rep 1": "20251106_Lpholis_Fin_combined_cell_metrics.csv",
    "Rep 2": "20251107_Lpholis_Fin_combined_cell_metrics.csv",
}

def load_reps(rep_dict):
    return {label: pd.read_csv(f, encoding="utf-8-sig")[COL].dropna().values
            for label, f in rep_dict.items()}

harp_reps    = load_reps(HARP_REPS)
lpholis_reps = load_reps(LPHOLIS_REPS)

harp_all    = np.concatenate(list(harp_reps.values()))
lpholis_all = np.concatenate(list(lpholis_reps.values()))

# Descriptive statistics
def describe(arr):
    return dict(n=len(arr), mean=np.mean(arr), sem=stats.sem(arr),
                median=np.median(arr), sd=np.std(arr, ddof=1),
                q25=np.percentile(arr, 25), q75=np.percentile(arr, 75))

print("=" * 60)
print("  MITOCHONDRIAL MASS — FIN CELLS  total_mito_area_raw (µm²)")
print("=" * 60)
for label, arr in {**{f"H. antarcticus {k}": v for k, v in harp_reps.items()},
                   "H. antarcticus (pooled)": harp_all,
                   **{f"L. pholis {k}": v for k, v in lpholis_reps.items()},
                   "L. pholis (pooled)": lpholis_all}.items():
    s  = describe(arr)
    sw = stats.shapiro(arr).pvalue if len(arr) >= 3 else float("nan")
    print(f"\n{label}  (n={s['n']})")
    print(f"  Mean ± SEM = {s['mean']:.2f} ± {s['sem']:.2f}  |  Median = {s['median']:.2f}")
    print(f"  SD = {s['sd']:.2f}  |  IQR = [{s['q25']:.2f}, {s['q75']:.2f}]")
    print(f"  Shapiro-Wilk p = {sw:.4f}")

# Statistical test
sw_h = stats.shapiro(harp_all).pvalue
sw_l = stats.shapiro(lpholis_all).pvalue
both_normal = sw_h > 0.05 and sw_l > 0.05
lev_p       = stats.levene(harp_all, lpholis_all).pvalue
equal_var   = lev_p > 0.05

if both_normal:
    stat, p_val = stats.ttest_ind(harp_all, lpholis_all, equal_var=equal_var)
    test_name   = "Student's t-test" if equal_var else "Welch's t-test"
    stat_label  = "t"
else:
    stat, p_val = stats.mannwhitneyu(harp_all, lpholis_all, alternative="two-sided")
    test_name   = "Mann-Whitney U"
    stat_label  = "U"

pooled_sd = np.sqrt((np.std(harp_all, ddof=1)**2 + np.std(lpholis_all, ddof=1)**2) / 2)
cohens_d  = (np.mean(harp_all) - np.mean(lpholis_all)) / pooled_sd

print(f"\n{'─'*60}")
print(f"Levene p = {lev_p:.4f}  |  {test_name}: {stat_label} = {stat:.4f}, p = {p_val:.4f}")
print(f"Cohen's d = {cohens_d:.3f}")
print("=" * 60)

HARP_REP_COLORS    = ["#1B4F8A", "#2E86AB", "#A8DADC"]
LPHOLIS_REP_COLORS = ["#C94A12", "#E07B54", "#F4A261"]

BOX_COLOR_HARP    = "#4C72B0"
BOX_COLOR_LPHOLIS = "#DD8452"

fig, ax = plt.subplots(figsize=(7, 6))

positions  = [1, 2]
plot_data  = [harp_all, lpholis_all]
box_colors = [BOX_COLOR_HARP, BOX_COLOR_LPHOLIS]
sp_labels  = ["H. antarcticus", "L. pholis"]

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
for patch, color in zip(bp["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.45)

rng = np.random.default_rng(42)
for (rep_label, arr), color in zip(harp_reps.items(), HARP_REP_COLORS):
    jitter = rng.uniform(-0.18, 0.18, size=len(arr))
    ax.scatter(1 + jitter, arr, color=color, edgecolors="white",
               linewidths=0.4, s=32, alpha=0.9, zorder=4)

for (rep_label, arr), color in zip(lpholis_reps.items(), LPHOLIS_REP_COLORS):
    jitter = rng.uniform(-0.18, 0.18, size=len(arr))
    ax.scatter(2 + jitter, arr, color=color, edgecolors="white",
               linewidths=0.4, s=32, alpha=0.9, zorder=4)

y_lo, y_hi = ax.get_ylim()
for pos, arr in zip(positions, plot_data):
    ax.text(pos, y_lo - 0.02 * (y_hi - y_lo),
            f"n = {len(arr)}", ha="center", va="top", fontsize=9)

y_top  = max(np.max(harp_all), np.max(lpholis_all))
y_line = y_top * 1.07
y_text = y_top * 1.115

ax.plot([1, 1, 2, 2], [y_line * 0.97, y_line, y_line, y_line * 0.97],
        color="black", linewidth=1.2)

def p_to_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

ax.text(1.5, y_text, p_to_stars(p_val), ha="center", va="bottom", fontsize=13)

annot = f"{test_name}\np = {p_val:.4f}\nCohen's d = {cohens_d:.2f}"
ax.text(0.97, 0.97, annot, transform=ax.transAxes, fontsize=8,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="grey", alpha=0.85))

harp_patches = [mpatches.Patch(color=c, label=lbl, alpha=0.9)
                for lbl, c in zip(harp_reps.keys(), HARP_REP_COLORS)]
lph_patches  = [mpatches.Patch(color=c, label=lbl, alpha=0.9)
                for lbl, c in zip(lpholis_reps.keys(), LPHOLIS_REP_COLORS)]

leg1 = ax.legend(handles=harp_patches, title="H. antarcticus reps",
                 fontsize=8, title_fontsize=8.5,
                 loc="upper left", framealpha=0.85)
ax.add_artist(leg1)
ax.legend(handles=lph_patches, title="L. pholis reps",
          fontsize=8, title_fontsize=8.5,
          loc="upper center", framealpha=0.85)

ax.set_xlim(0.4, 2.6)
ax.set_xticks(positions)
ax.set_xticklabels([f"$\\it{{{sp}}}$" for sp in sp_labels], fontsize=11)
ax.set_ylabel("Total mitochondrial area per cell (µm²)", fontsize=11)
ax.set_title("Mitochondrial mass — fin cells", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
out_path = "mito_fin_boxplot_2RepLpholis.svg"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"\nFigure saved → {out_path}")
plt.show()
