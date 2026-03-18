"""
Figure Generation: Manhattan + QQ Plots
========================================
Generates publication-ready Manhattan and QQ plots for
primary lymphedema, secondary lymphedema, and lipedema GWAS.

Output:
  figure2_manhattan.png/.svg       — Three-panel Manhattan (Figure 2)
  suppfig1_qq_plots.png/.svg       — Three-panel QQ (Supplementary Figure 1)

Usage: Run locally or in AoU Workbench after GWAS completion
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os
import subprocess

import os

# ── Set this to your workspace results bucket ─────────────────────────────────
# In the AoU Workbench: SECURE_BUCKET = os.environ["WORKSPACE_BUCKET"]
SECURE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "")

GW_SIG     = 5e-8
SUGGEST    = 1e-5
CHR_COLORS = ['#2C5F8A', '#85A9C5']
HIGHLIGHT  = '#E74C3C'
SUGGEST_C  = '#F39C12'

FS_TITLE  = 16
FS_AXIS   = 14
FS_TICK   = 12
FS_LEGEND = 11
FS_ANNOT  = 10
FS_PANEL  = 22


def prep_df(df):
    """Add CHR_NUM, POS as int, and log10p columns."""
    df = df.copy()
    df['CHR_NUM'] = df['CHR'].str.replace('chr', '').str.replace('X', '23')
    df['CHR_NUM'] = pd.to_numeric(df['CHR_NUM'], errors='coerce')
    df = df.dropna(subset=['CHR_NUM', 'POS'])
    df['CHR_NUM'] = df['CHR_NUM'].astype(int)
    df['POS']     = df['POS'].astype(int)
    df['log10p']  = -np.log10(df['p_value'])
    return df


def compute_x_positions(df):
    """Compute cumulative genomic x-positions for Manhattan plot."""
    offsets, mids = {}, {}
    cum = 0
    for chrom in sorted(df['CHR_NUM'].unique()):
        offsets[chrom] = cum
        mids[chrom]    = cum + df[df['CHR_NUM'] == chrom]['POS'].max() / 2
        cum += df[df['CHR_NUM'] == chrom]['POS'].max() + 5_000_000
    df = df.copy()
    df['x_pos'] = df['POS'] + df['CHR_NUM'].map(offsets)
    return df, mids


def compute_lgc(df):
    """Compute genomic inflation factor λGC."""
    chi2 = stats.chi2.ppf(1 - df['p_value'].dropna(), df=1)
    return np.median(chi2) / 0.456


def load_gwas_file(fname, bucket=SECURE_BUCKET, subdir="GWAS_by_chr"):
    """Load GWAS results, downloading from GCS if not local."""
    if not os.path.exists(fname):
        print(f"Downloading {fname}...")
        subprocess.run(f"gsutil cp {bucket}/{subdir}/{fname} .", shell=True)
    df = pd.read_csv(fname, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df['p_value'] = pd.to_numeric(df['p_value'], errors='coerce')
    df = df[df['p_value'] > 0].copy()
    return prep_df(df)


# ============================================================================
# FIGURE 2: THREE-PANEL MANHATTAN
# ============================================================================

def make_figure2(df_primary, df_secondary, df_lipedema,
                 outname="figure2_manhattan"):
    """
    Three-panel Manhattan plot (A=Primary, B=Secondary, C=Lipedema).
    Full width, no QQ inset. QQ plots are in Supplementary Figure 1.
    """
    fig, axes = plt.subplots(3, 1, figsize=(22, 21))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(hspace=0.35)

    panel_data = [
        (df_primary,   "Primary Lymphedema",   "A"),
        (df_secondary, "Secondary Lymphedema", "B"),
        (df_lipedema,  "Lipedema",             "C"),
    ]

    for ax, (df, title, panel_label) in zip(axes, panel_data):
        df, chr_mids = compute_x_positions(df)
        lgc    = compute_lgc(df)
        n_gw   = (df['p_value'] < GW_SIG).sum()
        n_sugg = (df['p_value'] < SUGGEST).sum() - n_gw
        ax.set_facecolor('white')

        # Scatter
        for i, chrom in enumerate(sorted(df['CHR_NUM'].unique())):
            subset = df[df['CHR_NUM'] == chrom]
            color  = CHR_COLORS[i % 2]
            base   = subset[subset['p_value'] >= SUGGEST]
            ax.scatter(base['x_pos'], base['log10p'],
                       c=color, s=5, alpha=0.6, linewidths=0, rasterized=True)
            sugg = subset[(subset['p_value'] < SUGGEST) & (subset['p_value'] >= GW_SIG)]
            ax.scatter(sugg['x_pos'], sugg['log10p'],
                       c=SUGGEST_C, s=22, alpha=0.9, linewidths=0, zorder=4)
            gw = subset[subset['p_value'] < GW_SIG]
            ax.scatter(gw['x_pos'], gw['log10p'],
                       c=HIGHLIGHT, s=50, alpha=1.0,
                       linewidths=0.5, edgecolors='white', zorder=5)

        # Threshold lines
        ax.axhline(-np.log10(GW_SIG),  color=HIGHLIGHT, lw=1.3, ls='--', alpha=0.8)
        ax.axhline(-np.log10(SUGGEST), color=SUGGEST_C, lw=1.0, ls=':', alpha=0.8)

        # Label GW significant hits
        gw_hits = df[df['p_value'] < GW_SIG].nsmallest(10, 'p_value').reset_index(drop=True)
        gw_offsets = [
            {'dy': 1.8, 'dx':  1_000_000, 'ha': 'left'},
            {'dy': 1.2, 'dx': -1_000_000, 'ha': 'right'},
            {'dy': 2.5, 'dx':  1_000_000, 'ha': 'left'},
            {'dy': 2.1, 'dx': -1_000_000, 'ha': 'right'},
        ]
        for i, (_, row) in enumerate(gw_hits.iterrows()):
            off = gw_offsets[i] if i < len(gw_offsets) else \
                  {'dy': 1.0 + i * 0.7, 'dx': 0, 'ha': 'center'}
            lbl = f"chr{row['CHR'].replace('chr','')}:{int(row['POS']):,}"
            ax.annotate(lbl,
                xy=(row['x_pos'], row['log10p']),
                xytext=(row['x_pos'] + off['dx'], row['log10p'] + off['dy']),
                fontsize=FS_ANNOT, color='#C0392B', fontweight='bold',
                ha=off['ha'], va='bottom',
                arrowprops=dict(arrowstyle='->', color='#C0392B',
                                lw=1.0, connectionstyle='arc3,rad=0.1'),
                zorder=6)

        # Label top 3 suggestive hits not overlapping GW labels
        sugg_hits = df[(df['p_value'] < SUGGEST) & (df['p_value'] >= GW_SIG)
                       ].nsmallest(5, 'p_value').reset_index(drop=True)
        labeled_x = [row['x_pos'] for _, row in gw_hits.iterrows()]
        count = 0
        for _, row in sugg_hits.iterrows():
            if count >= 3:
                break
            if any(abs(row['x_pos'] - lx) < 5e7 for lx in labeled_x):
                continue
            lbl = f"chr{row['CHR'].replace('chr','')}:{int(row['POS']):,}"
            ax.annotate(lbl,
                xy=(row['x_pos'], row['log10p']),
                xytext=(row['x_pos'], row['log10p'] + 1.0),
                fontsize=FS_ANNOT - 1, color='#E67E22', fontweight='bold',
                ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color='#E67E22', lw=0.8),
                zorder=5)
            labeled_x.append(row['x_pos'])
            count += 1

        # Chromosome labels
        for chrom, mid in chr_mids.items():
            lbl = 'X' if chrom == 23 else str(chrom)
            ax.text(mid, -0.65, lbl, ha='center', va='top',
                    fontsize=FS_TICK - 1, color='#555555')

        # Axes
        y_max = df['log10p'].max()
        ax.set_xlim(0, df['x_pos'].max() + 5_000_000)
        ax.set_ylim(-0.9, y_max + 3.5)
        ax.set_xlabel("Chromosome", fontsize=FS_AXIS, color='#2C3E50', labelpad=12)
        ax.set_ylabel("\u2212log\u2081\u2080(p)", fontsize=FS_AXIS, color='#2C3E50')
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=FS_TICK)
        ax.set_title(title, fontsize=FS_TITLE, fontweight='bold',
                     color='#1F3864', pad=10)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_color('#CCCCCC')
        ax.spines['left'].set_color('#CCCCCC')

        # λGC annotation
        ax.text(0.99, 0.97, f"\u03BBGC = {lgc:.3f}",
                transform=ax.transAxes, fontsize=FS_LEGEND,
                color='#2C3E50', ha='right', va='top', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#CCCCCC', alpha=0.9))

        # Legend
        handles = [
            mpatches.Patch(color=CHR_COLORS[0], label='Chromosome (odd)'),
            mpatches.Patch(color=CHR_COLORS[1], label='Chromosome (even)'),
            mpatches.Patch(color=SUGGEST_C,
                label=f'Suggestive (p\u00A0<\u00A01\u00D710\u207B\u2075, n\u00A0=\u00A0{n_sugg:,})'),
            mpatches.Patch(color=HIGHLIGHT,
                label=f'Genome-wide significant (p\u00A0<\u00A05\u00D710\u207B\u2078, n\u00A0=\u00A0{n_gw})'),
        ]
        ax.legend(handles=handles, loc='upper left',
                  fontsize=FS_LEGEND, framealpha=0.9, edgecolor='#CCCCCC')

        # Panel label
        ax.text(-0.02, 1.05, panel_label, transform=ax.transAxes,
                fontsize=FS_PANEL, fontweight='bold', color='#1F3864', va='top')

    for ext in ['png', 'svg']:
        plt.savefig(f"{outname}.{ext}",
                    dpi=300 if ext == 'png' else None,
                    bbox_inches='tight', facecolor='white', format=ext)
    subprocess.run(
        f"gsutil cp {outname}.png {outname}.svg {SECURE_BUCKET}/figures/",
        shell=True
    )
    print(f"✓ Figure 2 saved: {outname}.png/.svg")
    plt.show()


# ============================================================================
# SUPPLEMENTARY FIGURE 1: THREE-PANEL QQ
# ============================================================================

def make_suppfig1(df_primary, df_secondary, df_lipedema,
                  outname="suppfig1_qq_plots"):
    """
    Three-panel QQ plot (A=Primary, B=Secondary, C=Lipedema).
    λGC annotated on each panel. Legend on panel A only.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(wspace=0.35)

    qq_data = [
        (df_primary,   "Primary Lymphedema",   "A"),
        (df_secondary, "Secondary Lymphedema", "B"),
        (df_lipedema,  "Lipedema",             "C"),
    ]

    for ax, (df, title, panel_label) in zip(axes, qq_data):
        ax.set_facecolor('white')
        lgc = compute_lgc(df)
        n   = len(df)

        observed = np.sort(df['log10p'].values)[::-1]
        expected = -np.log10(np.arange(1, n + 1) / n)
        ci_upper = -np.log10(
            stats.beta.ppf(0.025, np.arange(1, n + 1), np.arange(n, 0, -1))
        )
        ci_lower = -np.log10(
            stats.beta.ppf(0.975, np.arange(1, n + 1), np.arange(n, 0, -1))
        )

        ax.fill_between(expected, ci_lower, ci_upper,
                        alpha=0.2, color='#85A9C5',
                        label='95% CI' if panel_label == 'A' else None)
        ax.scatter(expected, observed,
                   c=np.where(observed >= -np.log10(GW_SIG),  HIGHLIGHT,
                     np.where(observed >= -np.log10(SUGGEST), SUGGEST_C,
                              '#2C5F8A')),
                   s=4, alpha=0.7, linewidths=0, rasterized=True)
        ax.plot([0, expected.max()], [0, expected.max()],
                color='#E74C3C', lw=1.3, ls='--', zorder=5,
                label='Expected (null)' if panel_label == 'A' else None)

        ax.set_xlabel("Expected \u2212log\u2081\u2080(p)",
                      fontsize=FS_AXIS, color='#2C3E50')
        ax.set_ylabel("Observed \u2212log\u2081\u2080(p)",
                      fontsize=FS_AXIS, color='#2C3E50')
        ax.tick_params(labelsize=FS_TICK)
        ax.set_title(title, fontsize=FS_TITLE, fontweight='bold',
                     color='#1F3864', pad=10)
        ax.text(0.05, 0.95, f"\u03BBGC = {lgc:.3f}",
                transform=ax.transAxes, fontsize=FS_LEGEND,
                fontweight='bold', color='#2C3E50',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#CCCCCC', alpha=0.9))
        if panel_label == 'A':
            ax.legend(fontsize=FS_LEGEND, framealpha=0.9,
                      edgecolor='#CCCCCC', loc='upper left')
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_color('#CCCCCC')
        ax.spines['left'].set_color('#CCCCCC')
        ax.text(-0.08, 1.05, panel_label, transform=ax.transAxes,
                fontsize=FS_PANEL, fontweight='bold', color='#1F3864', va='top')

    for ext in ['png', 'svg']:
        plt.savefig(f"{outname}.{ext}",
                    dpi=300 if ext == 'png' else None,
                    bbox_inches='tight', facecolor='white', format=ext)
    subprocess.run(
        f"gsutil cp {outname}.png {outname}.svg {SECURE_BUCKET}/figures/",
        shell=True
    )
    print(f"✓ Supplementary Figure 1 saved: {outname}.png/.svg")
    plt.show()


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    # Load data
    df_primary   = load_gwas_file("gwas_PRIMARY_genome_wide_full.tsv")
    df_secondary = load_gwas_file("gwas_SECONDARY_genome_wide_clean.tsv")
    df_lipedema  = load_gwas_file("gwas_LIPEDEMA_genome_wide.tsv",
                                  subdir="")  # root of bucket

    print(f"Primary:   {len(df_primary):,} variants")
    print(f"Secondary: {len(df_secondary):,} variants")
    print(f"Lipedema:  {len(df_lipedema):,} variants")

    make_figure2(df_primary, df_secondary, df_lipedema)
    make_suppfig1(df_primary, df_secondary, df_lipedema)
