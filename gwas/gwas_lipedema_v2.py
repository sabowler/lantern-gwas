"""
Lipedema Common Variant GWAS Pipeline v2
=========================================
All of Us Controlled Tier Dataset v8 | Hail v0.2 | GRCh38

Fixes from v1:
  - FIXED chrX source mismatch: PAR filter_rows moved BEFORE x_dosage
    definition so all expressions reference the same MT object.
  - FIXED Firth non-convergence at N~110: switched to Wald regression.
    Firth's iterative solver fails wholesale at this case N.
  - FIXED null model explosion: removed sex_at_birth_num covariate.
    Lipedema cohort is female-only cases with female-matched controls,
    making sex a constant that causes Newton iteration to explode.
  - FIXED summary stats KeyError: column name normalisation maps Hail
    export variants to consistent names before display.

Design:
  - Lipedema cases (Tier 1 + Tier 2) vs 1:10 matched controls
  - Wald logistic regression (additive model)
  - Covariates: age, PC1-PC4 (sex removed — female-only cohort)
  - chrX sex-aware dosage with PAR exclusion
  - MAF >= 0.01, call rate >= 0.95, HWE >= 1e-6
  - Per-chromosome TSVs merged to genome-wide file

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

WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET")
WORKSPACE_CDR    = os.environ.get("WORKSPACE_CDR")
GOOGLE_PROJECT   = os.environ.get("GOOGLE_PROJECT")

COHORT_FILE    = "lipedema_final_cohort_v4.csv"
ANCESTRY_GCS   = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
OUTPUT_PREFIX  = "gwas_LIPEDEMA"
VDS_PATH       = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"

# ── Set this to your workspace results bucket ─────────────────────────────────
# In the AoU Workbench, your workspace bucket is available as:
#   os.environ["WORKSPACE_BUCKET"]
# You can also use a separate controlled-access results bucket if preferred.
# Example: RESULTS_FOLDER = f"{WORKSPACE_BUCKET}/gwas_lipedema_results"
RESULTS_FOLDER = f"{WORKSPACE_BUCKET}/gwas_lipedema_results"

# QC thresholds
MAF_THRESHOLD       = 0.01
CALL_RATE_THRESHOLD = 0.95
HWE_THRESHOLD       = 1e-6

# Regression: Wald (Firth fails to converge at N~110)
REGRESSION_TEST = 'wald'

# Covariates: sex_at_birth_num EXCLUDED — female-only cohort, constant column
# causes logistic null model to explode at Newton iteration 1
ALL_COVARIATES = ['age_at_gwas', 'PC1', 'PC2', 'PC3', 'PC4']

CHROMOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX']
START_CHROMOSOME_INDEX = 0   # set > 0 to resume after interruption

GENOME_WIDE_SIG = 5e-8
SUGGESTIVE_SIG  = 1e-5

# ============================================================================
# INITIALIZE HAIL
# ============================================================================

print("=" * 70)
print("LIPEDEMA COMMON VARIANT GWAS — v2")
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# Cluster: 8 workers x 16 CPUs x 104 GB RAM + 104 GB driver
hl.init(
    log='hail_lipedema_gwas.log',
    min_block_size=128,
    local_tmpdir='/tmp',
    spark_conf={
        'spark.driver.memory':             '80g',
        'spark.driver.maxResultSize':      '32g',
        'spark.executor.memory':           '85g',
        'spark.executor.cores':            '16',
        'spark.executor.instances':        '8',
        'spark.executor.memoryOverhead':   '10g',
        'spark.memory.fraction':           '0.8',
        'spark.memory.storageFraction':    '0.3',
        'spark.default.parallelism':       '256',
        'spark.sql.shuffle.partitions':    '256',
        'spark.dynamicAllocation.enabled': 'false',
        'spark.shuffle.service.enabled':   'false',
        'spark.local.dir':                 '/tmp',
    },
    idempotent=True
)
hl.default_reference('GRCh38')
print("✓ Hail initialized — 8 workers x 16 CPUs x 85 GB RAM")

# ============================================================================
# LOAD COHORT
# ============================================================================

print("\n1. Loading lipedema cohort...")

if not os.path.exists(COHORT_FILE):
    os.system(f"gsutil -u {GOOGLE_PROJECT} cp {WORKSPACE_BUCKET}/{COHORT_FILE} .")

cohort = pd.read_csv(COHORT_FILE)
cohort = cohort.rename(columns={"cohort_label": "Group"})
cohort["Group"] = cohort["Group"].str.capitalize()

# Keep Tier 1 + Tier 2 cases and matched controls
ACTIVE_TIERS = [
    "Tier 1 — Verified (code + WHR + BMI)",
    "Tier 2 — Code confirmed, measurements missing",
    "Tier 2 — Code confirmed, anthro does not fit profile",
    "Control",
]
cohort = cohort[
    (cohort["Group"] == "Control") |
    (cohort["confidence_tier"].isin(ACTIVE_TIERS))
].copy()

cohort["phenotype"] = (cohort["Group"] == "Case").astype(int)

# Age
if "current_age" in cohort.columns and cohort["current_age"].notna().any():
    cohort["age_at_gwas"] = cohort["current_age"]
else:
    cohort["date_of_birth"] = pd.to_datetime(cohort["date_of_birth"], utc=True)
    cohort["age_at_gwas"] = (
        (pd.Timestamp.now(tz="UTC") - cohort["date_of_birth"]).dt.days / 365.25
    )

print(f"  Cases:    {cohort['phenotype'].sum()}")
print(f"  Controls: {(cohort['phenotype'] == 0).sum()}")

# ============================================================================
# LOAD ANCESTRY PCs
# ============================================================================

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
ancestry[["PC1", "PC2", "PC3", "PC4"]] = pd.DataFrame(
    pcs.tolist(), index=ancestry.index
)
ancestry_pcs = ancestry[["person_id", "PC1", "PC2", "PC3", "PC4"]].copy()

cohort["person_id"] = cohort["person_id"].astype(int)
cohort = cohort.merge(ancestry_pcs, on="person_id", how="left")
cohort = cohort.dropna(subset=["PC1", "PC2", "PC3", "PC4"])

# Verify sex covariate is constant (confirms correct removal)
sex_unique = cohort["sex_at_birth"].str.lower().str.contains("female", na=False).nunique()
print(f"  sex_at_birth unique values: {sex_unique} (expected 1 for female-only cohort)")
print(f"  Final GWAS cohort: {len(cohort)} individuals")

cohort.to_csv("lipedema_gwas_cohort.csv", index=False)
os.system(f"gsutil cp lipedema_gwas_cohort.csv {WORKSPACE_BUCKET}/")
print("  ✓ GWAS cohort saved")

# ============================================================================
# LOAD VDS AND FILTER TO COHORT
# ============================================================================

print("\n3. Loading VDS and filtering to cohort samples...")

vds = hl.vds.read_vds(VDS_PATH)

sample_ids_list = cohort["person_id"].astype(str).tolist()
sample_table = hl.Table.from_pandas(
    pd.DataFrame({"s": sample_ids_list}), key="s"
)

vds_filtered = hl.vds.filter_samples(
    vds, sample_table, remove_dead_alleles=True
)
print(f"  ✓ VDS filtered to {vds_filtered.variant_data.count_cols():,} samples")

pheno_table = hl.Table.from_pandas(
    cohort[["person_id", "phenotype", "age_at_gwas",
            "PC1", "PC2", "PC3", "PC4"]].assign(
        s=cohort["person_id"].astype(str)
    ),
    key="s"
)

# ============================================================================
# GWAS FUNCTION
# ============================================================================

def run_gwas_chromosome(vds_filtered, pheno_table, chromosome):
    """
    Run Wald logistic regression GWAS for a single chromosome.

    Key implementation notes:
    - LGT → GT conversion required after VDS densification
    - chrX: PAR filter BEFORE x_dosage definition (source mismatch fix)
    - sex covariate excluded (female-only cohort)
    """
    chr_start = time.time()
    print(f"\n{'='*60}")
    print(f"CHROMOSOME: {chromosome}")
    print(f"{'='*60}")

    # 1. Filter and densify
    vds_chr = hl.vds.filter_intervals(
        vds_filtered,
        [hl.parse_locus_interval(chromosome, reference_genome='GRCh38')]
    )
    mt = hl.vds.to_dense_mt(vds_chr)
    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))
    n_var, n_samp = mt.count()
    print(f"  Loaded: {n_var:,} variants, {n_samp:,} samples")

    # 2. Variant QC
    mt = hl.variant_qc(mt)
    mt_qc = mt.filter_rows(
        (mt.variant_qc.call_rate >= CALL_RATE_THRESHOLD) &
        (mt.variant_qc.AF[1] >= MAF_THRESHOLD) &
        (mt.variant_qc.AF[1] <= 1 - MAF_THRESHOLD) &
        (mt.variant_qc.p_value_hwe >= HWE_THRESHOLD)
    )
    n_after_qc = mt_qc.count_rows()
    print(f"  After QC: {n_after_qc:,} variants ({100*n_after_qc/n_var:.1f}% retained)")

    if n_after_qc == 0:
        print(f"  ⚠ No variants passed QC — skipping {chromosome}")
        return None

    # 3. Annotate covariates
    mt_analysis = mt_qc.annotate_cols(**pheno_table[mt_qc.s])

    # 4. Sex-aware dosage for chrX
    # CRITICAL: PAR filter_rows MUST precede x_dosage definition.
    # Moving filter_rows after annotate_cols creates a new MT object;
    # covariate expressions remain bound to the pre-filter MT, causing
    # 'source mismatch' in logistic_regression_rows.
    if chromosome == 'chrX':
        par1 = hl.parse_locus_interval(
            'chrX:10001-2781479', reference_genome='GRCh38'
        )
        par2 = hl.parse_locus_interval(
            'chrX:155701383-156030895', reference_genome='GRCh38'
        )
        mt_analysis = mt_analysis.filter_rows(
            ~par1.contains(mt_analysis.locus) &
            ~par2.contains(mt_analysis.locus)
        )

    x_dosage = hl.float64(mt_analysis.GT.n_alt_alleles())

    # 5. Logistic regression (Wald)
    print(f"  Running {REGRESSION_TEST} regression...")
    gwas = hl.logistic_regression_rows(
        test=REGRESSION_TEST,
        y=mt_analysis.phenotype,
        x=x_dosage,
        covariates=[1.0] + [mt_analysis[cov] for cov in ALL_COVARIATES]
    )

    # 6. Annotate and export
    gwas = gwas.annotate(
        CHR=gwas.locus.contig,
        POS=gwas.locus.position,
        REF=gwas.alleles[0],
        ALT=gwas.alleles[1],
        analysis="LIPEDEMA"
    )
    output_path = f"{RESULTS_FOLDER}/{OUTPUT_PREFIX}_{chromosome}.tsv"
    gwas.export(output_path)

    elapsed = time.time() - chr_start
    print(f"  ✓ Exported: {output_path} ({elapsed/60:.1f} min)")
    return output_path

# ============================================================================
# MAIN LOOP
# ============================================================================

print("\n" + "="*70)
print("STARTING GWAS — ALL CHROMOSOMES")
print("="*70)

overall_start     = time.time()
completed_paths   = []
failed_chromosomes = []

for chrom_idx, chromosome in enumerate(
    CHROMOSOMES[START_CHROMOSOME_INDEX:],
    start=START_CHROMOSOME_INDEX
):
    print(f"\nChromosome {chrom_idx+1}/{len(CHROMOSOMES)}: {chromosome}")
    try:
        path = run_gwas_chromosome(vds_filtered, pheno_table, chromosome)
        if path:
            completed_paths.append(path)
    except Exception as e:
        print(f"  ✗ ERROR on {chromosome}: {e}")
        failed_chromosomes.append(chromosome)
        continue

# ============================================================================
# MERGE PER-CHROMOSOME RESULTS
# ============================================================================

print("\n" + "="*70)
print("MERGING PER-CHROMOSOME RESULTS")
print("="*70)

merge_cmd = (
    f"gsutil cat {RESULTS_FOLDER}/{OUTPUT_PREFIX}_chr*.tsv | "
    f"awk 'NR==1 || !/^locus/' > gwas_LIPEDEMA_genome_wide.tsv"
)
os.system(merge_cmd)
os.system(f"gsutil cp gwas_LIPEDEMA_genome_wide.tsv {RESULTS_FOLDER}/")
print("✓ Genome-wide results uploaded")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

if os.path.exists("gwas_LIPEDEMA_genome_wide.tsv"):
    results = pd.read_csv("gwas_LIPEDEMA_genome_wide.tsv", sep="\t")
    results.columns = [c.strip() for c in results.columns]

    # Normalise column names across Hail versions
    col_map = {}
    for col in results.columns:
        cl = col.lower()
        if cl == 'p_value' or cl.endswith('.p_value'):
            col_map[col] = 'p_value'
        elif cl in ('standard_error', 'se') or cl.endswith('.standard_error'):
            col_map[col] = 'standard_error'
        elif cl == 'beta' or cl.endswith('.beta'):
            col_map[col] = 'beta'
    results = results.rename(columns=col_map)

    total     = len(results)
    converged = results["p_value"].notna().sum() if "p_value" in results.columns else 0
    gw_sig    = (results["p_value"] <= GENOME_WIDE_SIG).sum() if converged > 0 else 0
    sugg      = ((results["p_value"] > GENOME_WIDE_SIG) &
                 (results["p_value"] <= SUGGESTIVE_SIG)).sum() if converged > 0 else 0
    min_p     = results["p_value"].min() if converged > 0 else float('nan')

    print(f"  Total variants:          {total:,}")
    print(f"  Converged:               {converged:,} ({100*converged/total:.1f}%)")
    print(f"  GW significant:          {gw_sig}")
    print(f"  Suggestive (p < 1e-5):   {sugg}")
    if pd.notna(min_p):
        print(f"  Min p-value:             {min_p:.3e}")

    display_cols = [c for c in ["locus", "alleles", "beta", "standard_error", "p_value"]
                    if c in results.columns]
    if gw_sig > 0:
        print(f"\n  GW significant hits:")
        print(results[results["p_value"] <= GENOME_WIDE_SIG]
              [display_cols].sort_values("p_value").to_string(index=False))
    elif sugg > 0:
        print(f"\n  Top 10 suggestive hits:")
        print(results.nsmallest(10, "p_value")[display_cols].to_string(index=False))

if failed_chromosomes:
    print(f"\n  ⚠ Failed chromosomes: {failed_chromosomes}")

total_elapsed = (time.time() - overall_start) / 60
print(f"\nTotal runtime: {total_elapsed:.0f} min")
print("✓ Lipedema GWAS complete")
