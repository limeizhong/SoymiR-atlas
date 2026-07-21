#!/usr/bin/env python3
"""
Plot: Confidence star-level distribution by miRNA annotation category.
Horizontal stacked bar chart showing proportion of 3-7★ miRNAs within each category.

Input:
  - 05_Figures/input/plotting_data/Figure2/confidence_star_distribution_input.xlsx
    with annotation_category + Confidence columns

Output:
  - 05_Figures/results/intermediate_figures/misc/status_star_distribution.png

Usage:
  cd SoymiR-atlas
  python3 05_Figures/scripts/Figure2/plot_status_star_distribution.py
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from pathlib import Path

plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 22,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 22,
    'legend.title_fontsize': 22,
})

# ── Paths ──
SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parents[1]
S7 = MODULE_DIR / "input/plotting_data/Figure2/confidence_star_distribution_input.xlsx"
OUTPUT = MODULE_DIR / "results/intermediate_figures/misc/status_star_distribution.png"

# ── Load ──
df = pd.read_excel(S7)

# ── Cross-tab ──
ct = pd.crosstab(df['annotation_category'], df['Confidence'])

# Y-axis order: bottom = reference-matched → top = novel (most → least known)
y_order = [
    'unannotated_opposite_arm_variant',
    'unannotated_opposite_arm_product',
    'novel_family_new_locus_variant',
    'novel_family_new_locus',
    'known_family_new_locus_variant',
    'known_family_new_locus',
    'reference_locus_variant',
    'reference_matched',
]
ct = ct.loc[y_order]
totals = ct.sum(axis=1)
ct_pct = ct.div(totals, axis=0) * 100
star_cols = [3, 4, 5, 6, 7]
category_labels = {
    'reference_matched': 'Reference matched',
    'reference_locus_variant': 'Reference locus variant',
    'known_family_new_locus': 'Known-family new locus',
    'known_family_new_locus_variant': 'Known-family new locus variant',
    'novel_family_new_locus': 'Novel-family new locus',
    'novel_family_new_locus_variant': 'Novel-family new locus variant',
    'unannotated_opposite_arm_product': 'Unannotated opposite-arm',
    'unannotated_opposite_arm_variant': 'Unannotated opposite-arm variant',
}
category_colors = {
    'reference_matched': '#4D4D4D',
    'reference_locus_variant': '#0072B2',
    'known_family_new_locus': '#56B4E9',
    'known_family_new_locus_variant': '#009E73',
    'novel_family_new_locus': '#E69F00',
    'novel_family_new_locus_variant': '#F0E442',
    'unannotated_opposite_arm_product': '#D55E00',
    'unannotated_opposite_arm_variant': '#CC79A7',
}

# ── Colors: single-hue blue gradient (light → dark) ──
colors = ['#deebf7', '#9ecae1', '#6baed6', '#3182bd', '#08519c']

# ── Plot ──
fig, ax = plt.subplots(figsize=(13.2, 7.2))

left = np.zeros(len(ct_pct))
for i, star in enumerate(star_cols):
    vals = ct_pct[star].values
    bars = ax.barh(ct_pct.index, vals, left=left, color=colors[i],
                   edgecolor='white', linewidth=0.8, label=f'{star}★')
    for j, (v, l) in enumerate(zip(vals, left)):
        if v >= 4:
            label_color = colors[4] if star in (3, 4) else 'white'
            ax.text(l + v / 2, j, f'{v:.0f}%', ha='center', va='center',
                    fontsize=15, fontweight='normal', color=label_color)
    left += vals

# Total N annotations
for j, total in enumerate(totals):
    ax.text(1.025, j, f'n = {int(total)}', transform=ax.get_yaxis_transform(),
            va='center', ha='left', fontsize=18, color='black', clip_on=False)

# ── Styling ──
ax.set_xlabel('Proportion of miRNAs')
ax.set_xlim(0, 100)
ax.set_yticks(np.arange(len(ct_pct.index)))
ax.set_yticklabels([""] * len(ct_pct.index))
ax.tick_params(axis='y', length=0)
for j, category in enumerate(ct_pct.index):
    ax.add_patch(
        Rectangle(
            (-0.064, j - 0.31),
            0.040,
            0.62,
            transform=ax.get_yaxis_transform(),
            facecolor=category_colors[category],
            edgecolor='none',
            clip_on=False,
        )
    )

legend = ax.legend(loc='lower center', ncol=5,
                   frameon=False,
                   bbox_to_anchor=(0.5, 0.985))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_bounds(0, 100)
ax.xaxis.grid(False)
ax.set_axisbelow(True)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])

plt.tight_layout()
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT, dpi=600, bbox_inches='tight', facecolor='white')
print(f"Saved: {OUTPUT}")
plt.close()
