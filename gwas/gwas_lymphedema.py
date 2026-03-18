"""
Lymphedema Common Variant GWAS Pipeline
========================================
All of Us Controlled Tier Dataset v8 | Hail v0.2 | GRCh38

Runs parallel GWAS for primary and secondary lymphedema vs matched controls.

Design:
  - Logistic regression, Wald test, additive model
  - Covariates: age, sex_at_birth, PC1-PC4
  - chrX sex-aware dosage with PAR exclusion
  - MAF >= 0.01, call rate >= 0.95, HWE >= 1e-6
  - Per-chromosome TSVs merged to genome-wide file

Cluster: 10 workers x 32 CPUs x 208 GB RAM + 416 GB driver

Usage: Run in All of Us Researcher Workbench Jupyter notebook with Hail cluster
"""

import hail as hl
import pandas as pd
import numpy as np
import os
import ast
import time
from datetime import datetime

# ============================================================================
# CONFIGURATION
# ============================================================================

WORKSPACE_BUCKET  = os.environ.get("WORKSPACE_BUCKET")
WORKSPACE_CDR     = os.environ.get("WORKSPACE_CDR")
GOOGLE_PROJECT    = os.environ.get("GOOGLE_PROJECT")

COHORT_FILE    = "lymphedema_cohort_final.csv"
ANCESTRY_GCS   = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
VDS_PATH       = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
# ── Set this to your workspace results bucket ─────────────────────────────────
# In the AoU Workbench: RESULTS_BUCKET = os.environ["WORKSPACE_BUCKET"]
# Or use a separate controlled-access results bucket if preferred.
RESULTS_BUCKET = WORKSPACE_BUCKET

# QC thresholds
MAF_THRESHOLD       = 0.01
CALL_RATE_THRESHOLD = 0.95
HWE_THRESHOLD       = 1e-6

REGRESSION_TEST = 'wald'
ALL_COVARIATES  = ['sex_at_birth_num', 'age_at_gwas', 'PC1', 'PC2', 'PC3', 'PC4']

CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX']
START_CHROMOSOME_INDEX = 0

GENOME_WIDE_SIG = 5e-8
SUGGESTIVE_SIG  = 1e-5

# ============================================================================
# INITIALIZE HAIL
# ============================================================================

print("=" * 70)
print("LYMPHEDEMA COMMON VARIANT GWAS")
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# Cluster: 10 workers x 32 CPUs x 208 GB RAM + 416 GB driver
hl.init(
    log='hail_lymphedema_gwas.log',
    min_block_size=128,
    local_tmpdir='/tmp',
    spark_conf={
        'spark.driver.memory':             '200g',
        'spark.driver.maxResultSize':      '64g',
        'spark.executor.memory':           '180g',
        'spark.executor.cores':            '32',
        'spark.executor.instances':        '10',
        'spark.executor.memoryOverhead':   '20g',
        'spark.memory.fraction':           '0.8',
        'spark.memory.storageFraction':    '0.3',
        'spark.default.parallelism':       '640',
        'spark.sql.shuffle.partitions':    '640',
        'spark.dynamicAllocation.enabled': 'false',
        'spark.shuffle.service.enabled':   'false',
        'spark.local.dir':                 '/tmp',
    },
    idempotent=True
)
hl.default_reference('GRCh38')
print("✓ Hail initialized — 10 workers x 32 CPUs x 180 GB RAM")

# ============================================================================
# LOAD COHORT AND ANCESTRY PCs
# ============================================================================

print("\n1. Loading lymphedema cohort...")
if not os.path.exists(COHORT_FILE):
    os.system(f"gsutil -u {GOOGLE_PROJECT} cp {WORKSPACE_BUCKET}/{COHORT_FILE} .")

cohort = pd.read_csv(COHORT_FILE)
print(f"  {len(cohort):,} individuals")
print(cohort.groupby(["condition", "group"]).size().to_string())

cohort["phenotype_primary"]   = (cohort["condition"] == "PRIMARY").astype(int)
cohort["phenotype_secondary"] = (cohort["condition"] == "SECONDARY").astype(int)
cohort["phenotype_any"]       = (cohort["group"] == "Case").astype(int)
cohort["sex_at_birth_num"]    = (
    cohort["sex_at_birth"].str.lower().str.contains("female", na=False)
).astype(int)

print("\n2. Loading ancestry PCs...")
ANCESTRY_LOCAL = "ancestry_preds.tsv"
if not os.path.exists(ANCESTRY_LOCAL):
    os.system(f"gsutil -u {GOOGLE_PROJECT} cp {ANCESTRY_GCS} {ANCESTRY_LOCAL}")

ancestry = pd.read_csv(ANCESTRY_LOCAL, sep="\t")
ancestry.columns = [c.strip().lower() for c in ancestry.columns]
ID_COL = next(
    (c for c in ancestry.columns if 'research' in c or 'person' in c),
    ancestry.columns[0]
)
ancestry = ancestry.rename(columns={ID_COL: "person_id"})
ancestry["person_id"] = ancestry["person_id"].astype(int)

def parse_pca(val):
    try:
        return ast.literal_eval(str(val))[:4]
    except Exception:
        return [np.nan, np.nan, np.nan, np.nan]

pcs = ancestry["pca_features"].apply(parse_pca)
ancestry[["PC1","PC2","PC3","PC4"]] = pd.DataFrame(pcs.tolist(), index=ancestry.index)
ancestry_pcs = ancestry[["person_id","PC1","PC2","PC3","PC4"]].copy()

cohort["person_id"] = cohort["person_id"].astype(int)
cohort = cohort.merge(ancestry_pcs, on="person_id", how="left")
cohort = cohort.dropna(subset=["PC1","PC2","PC3","PC4"])
print(f"  Final GWAS cohort: {len(cohort):,}")

# ============================================================================
# LOAD VDS
# ============================================================================

print("\n3. Loading VDS...")
vds = hl.vds.read_vds(VDS_PATH)
sample_table = hl.Table.from_pandas(
    pd.DataFrame({"s": cohort["person_id"].astype(str).tolist()}), key="s"
)
vds_filtered = hl.vds.filter_samples(vds, sample_table, remove_dead_alleles=True)
print(f"  ✓ VDS filtered to {vds_filtered.variant_data.count_cols():,} samples")

# ============================================================================
# GWAS FUNCTION
# ============================================================================

def run_gwas_chromosome(vds_filtered, pheno_table, chromosome, analysis_label):
    """
    Run Wald logistic regression for one chromosome and analysis.
    analysis_label: 'PRIMARY' or 'SECONDARY'
    """
    chr_start = time.time()
    print(f"\n[{analysis_label}] {chromosome}")

    # Filter + densify
    vds_chr = hl.vds.filter_intervals(
        vds_filtered,
        [hl.parse_locus_interval(chromosome, reference_genome='GRCh38')]
    )
    mt = hl.vds.to_dense_mt(vds_chr)
    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))
    n_var, _ = mt.count()

    # QC
    mt = hl.variant_qc(mt)
    mt_qc = mt.filter_rows(
        (mt.variant_qc.call_rate >= CALL_RATE_THRESHOLD) &
        (mt.variant_qc.AF[1] >= MAF_THRESHOLD) &
        (mt.variant_qc.AF[1] <= 1 - MAF_THRESHOLD) &
        (mt.variant_qc.p_value_hwe >= HWE_THRESHOLD)
    )
    n_qc = mt_qc.count_rows()
    print(f"  QC: {n_qc:,}/{n_var:,} variants")
    if n_qc == 0:
        return None

    # Annotate covariates
    mt_analysis = mt_qc.annotate_cols(**pheno_table[mt_qc.s])

    # chrX sex-aware dosage (PAR filter BEFORE dosage expression)
    if chromosome == 'chrX':
        par1 = hl.parse_locus_interval('chrX:10001-2781479', reference_genome='GRCh38')
        par2 = hl.parse_locus_interval('chrX:155701383-156030895', reference_genome='GRCh38')
        mt_analysis = mt_analysis.filter_rows(
            ~par1.contains(mt_analysis.locus) &
            ~par2.contains(mt_analysis.locus)
        )

    x_dosage = hl.float64(mt_analysis.GT.n_alt_alleles())

    # Logistic regression
    phenotype_col = (mt_analysis.phenotype_primary
                     if analysis_label == 'PRIMARY'
                     else mt_analysis.phenotype_secondary)
    gwas = hl.logistic_regression_rows(
        test=REGRESSION_TEST,
        y=phenotype_col,
        x=x_dosage,
        covariates=[1.0] + [mt_analysis[cov] for cov in ALL_COVARIATES]
    )

    # Annotate and export
    gwas = gwas.annotate(
        CHR=gwas.locus.contig,
        POS=gwas.locus.position,
        REF=gwas.alleles[0],
        ALT=gwas.alleles[1],
        analysis=analysis_label
    )
    output_path = f"{RESULTS_BUCKET}/GWAS_by_chr/gwas_{analysis_label}_{chromosome}.tsv"
    gwas.export(output_path)
    print(f"  ✓ Exported ({(time.time()-chr_start)/60:.1f} min)")
    return output_path

# ============================================================================
# MAIN LOOP — BOTH ANALYSES
# ============================================================================

for analysis_label in ["PRIMARY", "SECONDARY"]:
    print(f"\n{'='*70}")
    print(f"ANALYSIS: {analysis_label} LYMPHEDEMA")
    print(f"{'='*70}")

    # Build phenotype table for this analysis
    pheno_table = hl.Table.from_pandas(
        cohort[["person_id", "phenotype_primary", "phenotype_secondary",
                "sex_at_birth_num", "age_at_gwas",
                "PC1","PC2","PC3","PC4"]].assign(
            s=cohort["person_id"].astype(str)
        ),
        key="s"
    )

    failed = []
    for chrom_idx, chromosome in enumerate(CHROMOSOMES[START_CHROMOSOME_INDEX:],
                                            start=START_CHROMOSOME_INDEX):
        try:
            run_gwas_chromosome(vds_filtered, pheno_table, chromosome, analysis_label)
        except Exception as e:
            print(f"  ✗ ERROR on {chromosome}: {e}")
            failed.append(chromosome)

    # Merge
    merge_cmd = (
        f"gsutil cat {RESULTS_BUCKET}/GWAS_by_chr/gwas_{analysis_label}_chr*.tsv | "
        f"awk 'NR==1 || !/^locus/' > gwas_{analysis_label}_genome_wide_full.tsv"
    )
    os.system(merge_cmd)
    os.system(f"gsutil cp gwas_{analysis_label}_genome_wide_full.tsv "
              f"{RESULTS_BUCKET}/GWAS_by_chr/")
    print(f"\n✓ {analysis_label} genome-wide results merged")
    if failed:
        print(f"  ⚠ Failed chromosomes: {failed}")

print("\n✓ Lymphedema GWAS complete")
