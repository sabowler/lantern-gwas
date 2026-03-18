"""
LD Calculation Between Lead Variants
======================================
Computes Pearson r² between two lead variants from the lipedema GWAS
using dosage values across cohort samples.

Lead variants:
  chr19:2,738,482 G>C  (rs111311797, p = 1.5e-8, OR = 4.74)
  chr19:2,740,851 G>A  (rs111567842, p = 4.85e-8, OR = 4.39)

Usage: Run in AoU Workbench with Hail cluster after GWAS completion
"""

import hail as hl
import pandas as pd
import numpy as np
import os
import time

WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET")
GOOGLE_PROJECT   = os.environ.get("GOOGLE_PROJECT")
VDS_PATH = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"

LEAD_1 = 2738482   # rs111311797
LEAD_2 = 2740851   # rs111567842

# ============================================================================
# INITIALIZE HAIL
# ============================================================================

hl.init(quiet=True)
hl.default_reference('GRCh38')

# ============================================================================
# LOAD COHORT
# ============================================================================

print("Loading cohort...")
cohort = pd.read_csv("lipedema_gwas_cohort.csv")
cohort_ids = cohort['person_id'].astype(str).tolist()
print(f"  {len(cohort_ids):,} samples")

# ============================================================================
# LOAD VDS, FILTER SAMPLES, FILTER TO LEAD VARIANTS
# ============================================================================

print("Loading VDS...")
vds = hl.vds.read_vds(VDS_PATH)

cohort_table = hl.Table.parallelize(
    [hl.struct(s=sid) for sid in cohort_ids],
    schema=hl.tstruct(s=hl.tstr),
    key='s'
)

print("Filtering samples...")
vds_filtered = hl.vds.filter_samples(vds, cohort_table, keep=True)

print("Converting to dense MT...")
mt = hl.vds.to_dense_mt(vds_filtered)

print("Filtering to 2 lead variants...")
mt_ld = mt.filter_rows(
    (mt.locus == hl.locus('chr19', LEAD_1, reference_genome='GRCh38')) |
    (mt.locus == hl.locus('chr19', LEAD_2, reference_genome='GRCh38'))
)
mt_ld = mt_ld.annotate_entries(
    GT=hl.vds.lgt_to_gt(mt_ld.LGT, mt_ld.LA)
)
mt_ld = mt_ld.annotate_entries(
    dosage=hl.float64(mt_ld.GT.n_alt_alleles())
)

n = mt_ld.count_rows()
print(f"  Found {n} variants (expected 2)")

# ============================================================================
# EXTRACT DOSAGES AND COMPUTE r²
# ============================================================================

print("Extracting dosages...")
t1 = time.time()

entries = mt_ld.entries()
entries = entries.key_by()
entries = entries.select('s', 'var_key', 'dosage')

# Add variant key
mt_ld_keyed = mt_ld.annotate_rows(
    var_key=hl.str(mt_ld.locus.position)
)
mt_ld_keyed = mt_ld_keyed.annotate_entries(
    var_key=mt_ld_keyed.var_key
)

entries = mt_ld_keyed.entries()
entries = entries.key_by()
entries = entries.select(
    's',
    'var_key',
    dosage=hl.float64(entries.GT.n_alt_alleles())
)

# Group per sample: collect dosages for both variants
sample_dosages = (
    entries
    .group_by('s')
    .aggregate(
        dosages=hl.agg.collect(hl.struct(var=entries.var_key, d=entries.dosage))
    )
)

print(f"  Aggregation defined ({time.time()-t1:.1f}s)")
print("  Collecting to pandas (this is the slow step)...")

df = sample_dosages.to_pandas()
print(f"  ✓ Collected {len(df):,} samples ({time.time()-t1:.1f}s)")

# ============================================================================
# COMPUTE r²
# ============================================================================

print("\nComputing r²...")

def parse_dosages(row, var_key1, var_key2):
    d = {item['var']: item['d'] for item in row}
    return d.get(var_key1, np.nan), d.get(var_key2, np.nan)

var1_key = str(LEAD_1)
var2_key = str(LEAD_2)

parsed = df['dosages'].apply(
    lambda x: parse_dosages(x, var1_key, var2_key)
)
df['d1'] = parsed.apply(lambda x: x[0])
df['d2'] = parsed.apply(lambda x: x[1])

df_complete = df[['d1','d2']].dropna()
print(f"  Samples with both dosages: {len(df_complete):,}")

# r² = Pearson correlation squared
r = np.corrcoef(df_complete['d1'], df_complete['d2'])[0, 1]
r2 = r ** 2

print(f"\n{'='*50}")
print(f"  chr19:{LEAD_1} vs chr19:{LEAD_2}")
print(f"  r  = {r:.4f}")
print(f"  r² = {r2:.4f}")
print(f"{'='*50}")

# Dosage distributions
print(f"\nDosage distribution — chr19:{LEAD_1}:")
print(df_complete['d1'].value_counts().sort_index().to_string())

print(f"\nDosage distribution — chr19:{LEAD_2}:")
print(df_complete['d2'].value_counts().sort_index().to_string())

# Save result
result = pd.DataFrame({
    'var1': [f'chr19:{LEAD_1}'],
    'var2': [f'chr19:{LEAD_2}'],
    'r': [r],
    'r2': [r2],
    'n_samples': [len(df_complete)]
})
result.to_csv("lipedema_ld_result.csv", index=False)
import subprocess
WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "")
subprocess.run(
    f"gsutil cp lipedema_ld_result.csv {WORKSPACE_BUCKET}/",
    shell=True
)
print(f"\n✓ LD result saved")
print(f"  r² = {r2:.4f} — "
      f"{'likely same signal' if r2 > 0.8 else 'possibly independent signals' if r2 < 0.2 else 'moderate LD'}")
