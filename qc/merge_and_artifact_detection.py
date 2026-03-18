"""
QC: Merge Per-Chromosome TSVs + Artifact Detection
====================================================
Merges per-chromosome GWAS result TSVs into genome-wide files
and screens for segmental duplication artifacts using beta
uniformity analysis.

Artifact signature (chr10 secondary lymphedema example):
  - Span > 1 Mb across suggestive variants
  - Beta CV < 0.05 (nearly identical effect sizes)
  - All betas same directionality

Usage: Run locally or in AoU Workbench after per-chromosome GWAS completes
"""

import pandas as pd
import numpy as np
import subprocess
import os

import os

# ── Set this to your workspace results bucket ─────────────────────────────────
# In the AoU Workbench: SECURE_BUCKET = os.environ["WORKSPACE_BUCKET"]
SECURE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "")
CHROMS = [str(i) for i in range(1, 23)] + ['X']

SUGGESTIVE_SIG  = 1e-5
GENOME_WIDE_SIG = 5e-8

# ============================================================================
# MERGE PER-CHROMOSOME FILES
# ============================================================================

def merge_per_chr(condition, bucket=SECURE_BUCKET, subdir="GWAS_by_chr"):
    """
    Download and merge per-chromosome GWAS TSVs for a given condition.
    Returns merged DataFrame with schema-validated rows only.
    """
    outfile = f"gwas_{condition}_genome_wide_full.tsv"
    if os.path.exists(outfile):
        print(f"✓ {outfile} already exists — loading")
        df = pd.read_csv(outfile, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        df['p_value'] = pd.to_numeric(df['p_value'], errors='coerce')
        return df[df['p_value'] > 0].copy()

    print(f"\n=== Merging {condition} ===")
    dfs = []
    header_cols = None

    for chrom in CHROMS:
        fname = f"gwas_{condition}_chr{chrom}.tsv"
        gcs_path = f"{bucket}/{subdir}/{fname}"

        ret = subprocess.run(
            f"gsutil -q cp {gcs_path} /tmp/{fname}",
            shell=True, capture_output=True, text=True
        )
        if ret.returncode != 0:
            print(f"  ✗ chr{chrom}: not found — skipping")
            continue

        try:
            df = pd.read_csv(f"/tmp/{fname}", sep="\t")
            df.columns = [c.strip() for c in df.columns]

            if header_cols is None:
                header_cols = df.columns.tolist()
                print(f"  Schema: {header_cols}")
            elif df.columns.tolist() != header_cols:
                print(f"  ⚠ chr{chrom}: schema mismatch — skipping")
                continue

            n_sugg = (df['p_value'] < SUGGESTIVE_SIG).sum() \
                if 'p_value' in df.columns else 0
            n_gw   = (df['p_value'] < GENOME_WIDE_SIG).sum() \
                if 'p_value' in df.columns else 0
            print(f"  ✓ chr{chrom}: {len(df):,} variants "
                  f"(suggestive={n_sugg}, GW sig={n_gw})")
            dfs.append(df)

        except Exception as e:
            print(f"  ✗ chr{chrom}: error — {e}")
            continue

    if not dfs:
        raise RuntimeError(f"No files loaded for {condition}")

    merged = pd.concat(dfs, ignore_index=True)
    merged = merged[merged['p_value'].notna()]
    merged['p_value'] = pd.to_numeric(merged['p_value'], errors='coerce')
    merged = merged[merged['p_value'] > 0]

    print(f"\n  Total merged: {len(merged):,}")
    print(f"  GW sig: {(merged['p_value'] < GENOME_WIDE_SIG).sum()}")
    print(f"  Suggestive: {(merged['p_value'] < SUGGESTIVE_SIG).sum()}")

    merged.to_csv(outfile, sep="\t", index=False)
    subprocess.run(
        f"gsutil cp {outfile} {bucket}/{subdir}/",
        shell=True
    )
    print(f"  ✓ Saved: {outfile}")
    return merged


# ============================================================================
# ARTIFACT DETECTION
# ============================================================================

def screen_for_artifacts(df, condition, cv_threshold=0.05, span_threshold=1_000_000,
                          n_threshold=10):
    """
    Screen suggestive variants for segmental duplication artifacts.

    Artifact signature:
      - CV of beta < cv_threshold (uniform effect sizes)
      - Genomic span > span_threshold bp (too wide for real signal)
      - n > n_threshold variants (enough to compute reliable CV)

    Real LD blocks also show low CV but are confined to small spans (<100kb).
    The combined span + CV criterion distinguishes artifacts from real signals.

    Returns dict of flagged regions with artifact statistics.
    """
    print(f"\n=== Artifact screen: {condition} ===")
    print(f"  Criteria: CV < {cv_threshold} AND span > {span_threshold/1e6:.1f}Mb AND n > {n_threshold}")

    sugg = df[df['p_value'] < SUGGESTIVE_SIG].copy()
    flagged = {}

    for chrom in sorted(sugg['CHR'].unique()):
        sub = sugg[sugg['CHR'] == chrom]
        if len(sub) < 2:
            continue

        span = sub['POS'].max() - sub['POS'].min()
        mean_beta = sub['beta'].mean()
        cv = sub['beta'].std() / abs(mean_beta) if mean_beta != 0 else 0
        all_same_dir = (sub['beta'] > 0).all() or (sub['beta'] < 0).all()

        is_artifact = (cv < cv_threshold and
                       span > span_threshold and
                       len(sub) > n_threshold)

        status = "⚠ ARTIFACT?" if is_artifact else "OK"
        print(f"  {chrom:6s}: n={len(sub):3d}  "
              f"span={span/1e6:6.2f}Mb  "
              f"beta={mean_beta:.3f}±{sub['beta'].std():.3f}  "
              f"CV={cv:.3f}  "
              f"same_dir={all_same_dir}  {status}")

        if is_artifact:
            flagged[chrom] = {
                'n': len(sub),
                'span_mb': span / 1e6,
                'pos_min': sub['POS'].min(),
                'pos_max': sub['POS'].max(),
                'beta_mean': mean_beta,
                'beta_std': sub['beta'].std(),
                'cv': cv,
                'all_same_dir': all_same_dir
            }

    return flagged


def exclude_artifacts(df, artifact_regions, buffer=10_000):
    """
    Remove variants in flagged artifact regions from DataFrame.

    artifact_regions: dict of {chrom: {'pos_min': int, 'pos_max': int}}
    buffer: bp to add on each side of region boundary
    """
    mask = pd.Series(False, index=df.index)

    for chrom, region in artifact_regions.items():
        chrom_mask = (
            (df['CHR'] == chrom) &
            (df['POS'] >= region['pos_min'] - buffer) &
            (df['POS'] <= region['pos_max'] + buffer)
        )
        n_removed = chrom_mask.sum()
        print(f"  Removing {n_removed:,} variants from {chrom} artifact region "
              f"({region['pos_min']:,}–{region['pos_max']:,})")
        mask |= chrom_mask

    df_clean = df[~mask].copy()
    print(f"  Total removed: {mask.sum():,}  Remaining: {len(df_clean):,}")
    return df_clean


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    from scipy import stats

    # Merge all three conditions
    df_primary   = merge_per_chr("PRIMARY")
    df_secondary = merge_per_chr("SECONDARY")

    # Screen for artifacts
    primary_artifacts   = screen_for_artifacts(df_primary,   "PRIMARY")
    secondary_artifacts = screen_for_artifacts(df_secondary, "SECONDARY")

    # Apply artifact exclusions
    if secondary_artifacts:
        print("\nExcluding secondary lymphedema artifacts...")
        df_secondary_clean = exclude_artifacts(df_secondary, secondary_artifacts)

        # Recompute lambda after exclusion
        chi2 = stats.chi2.ppf(1 - df_secondary_clean['p_value'].dropna(), df=1)
        lgc  = np.median(chi2) / 0.456
        print(f"  λGC after exclusion: {lgc:.3f}")

        df_secondary_clean.to_csv(
            "gwas_SECONDARY_genome_wide_clean.tsv", sep="\t", index=False
        )
        subprocess.run(
            f"gsutil cp gwas_SECONDARY_genome_wide_clean.tsv "
            f"{SECURE_BUCKET}/GWAS_by_chr/",
            shell=True
        )
        print("  ✓ Clean secondary file saved")
    else:
        print("\nNo artifacts detected in secondary lymphedema")
        df_secondary_clean = df_secondary

    # Summary
    print("\n=== FINAL SUMMARY ===")
    for label, df in [("Primary",   df_primary),
                      ("Secondary", df_secondary_clean)]:
        chi2 = stats.chi2.ppf(1 - df['p_value'].dropna(), df=1)
        lgc  = np.median(chi2) / 0.456
        print(f"  {label}: {len(df):,} variants  "
              f"GW sig={(df['p_value'] < GENOME_WIDE_SIG).sum()}  "
              f"Suggestive={(df['p_value'] < SUGGESTIVE_SIG).sum()}  "
              f"λGC={lgc:.3f}")
