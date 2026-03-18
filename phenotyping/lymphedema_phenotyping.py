"""
Lymphedema Cohort Phenotyping Pipeline
=======================================
Identifies primary and secondary lymphedema cases from the All of Us
Controlled Tier Dataset v8 using OMOP/SNOMED condition codes.

Primary lymphedema:  participants with any hereditary lymphedema concept
Secondary lymphedema: participants with non-hereditary lymphedema concepts

Inclusion criteria:
  - Adult (age >= 18)
  - Self-identified White (8527) or Black/African American (8516)
  - Whole-genome sequencing data available

Exclusion criteria (controls):
  - Any lymphedema diagnosis
  - Documented immunosuppressant exposure

Matching: 1:10, race + sex + age (nearest neighbor by days)

Usage: Run in All of Us Researcher Workbench (Controlled Tier v8)
"""

import pandas as pd
import numpy as np
from datetime import datetime

# ── BigQuery / CDR ────────────────────────────────────────────────────────────
import os
WORKSPACE_CDR   = os.environ["WORKSPACE_CDR"]
WORKSPACE_BUCKET = os.environ["WORKSPACE_BUCKET"]
GOOGLE_PROJECT  = os.environ["GOOGLE_PROJECT"]

# ── Configuration ─────────────────────────────────────────────────────────────
RACE_CONCEPTS = [8527, 8516]   # White, Black/African American
MIN_AGE       = 18
MATCH_RATIO   = 10             # 1:10 case:control

# ── Step 1: Extract participant demographics ──────────────────────────────────
demo_query = f"""
SELECT
    p.person_id,
    p.birth_datetime,
    p.sex_at_birth_concept_id,
    p.race_concept_id,
    p.ethnicity_concept_id,
    DATE_DIFF(CURRENT_DATE(), DATE(p.birth_datetime), YEAR) AS age_at_analysis,
    c_sex.concept_name   AS sex_at_birth,
    c_race.concept_name  AS race,
    c_eth.concept_name   AS ethnicity
FROM `{WORKSPACE_CDR}.person` p
LEFT JOIN `{WORKSPACE_CDR}.concept` c_sex
    ON p.sex_at_birth_concept_id = c_sex.concept_id
LEFT JOIN `{WORKSPACE_CDR}.concept` c_race
    ON p.race_concept_id = c_race.concept_id
LEFT JOIN `{WORKSPACE_CDR}.concept` c_eth
    ON p.ethnicity_concept_id = c_eth.concept_id
WHERE p.race_concept_id IN ({",".join(str(r) for r in RACE_CONCEPTS)})
  AND DATE_DIFF(CURRENT_DATE(), DATE(p.birth_datetime), YEAR) >= {MIN_AGE}
"""

print("Extracting demographics...")
demo_df = pd.read_gbq(demo_query, project_id=GOOGLE_PROJECT,
                      dialect="standard")
print(f"  {len(demo_df):,} eligible participants")

# ── Step 2: Identify lymphedema cases ─────────────────────────────────────────
# Query cb_criteria for hereditary lymphedema concept IDs
hereditary_query = f"""
SELECT DISTINCT person_id
FROM `{WORKSPACE_CDR}.cb_search_all_events`
WHERE is_standard = 1
  AND concept_id IN (
      SELECT concept_id FROM `{WORKSPACE_CDR}.cb_criteria`
      WHERE full_text LIKE '%hereditary lymphedema%'
         OR full_text LIKE '%primary lymphedema%'
         OR concept_name LIKE '%hereditary lymphedema%'
  )
"""

secondary_query = f"""
SELECT DISTINCT person_id
FROM `{WORKSPACE_CDR}.condition_occurrence`
WHERE condition_concept_id IN (
    SELECT concept_id FROM `{WORKSPACE_CDR}.concept`
    WHERE (concept_name LIKE '%lymphedema%'
       OR concept_name LIKE '%lymphoedema%')
      AND vocabulary_id IN ('SNOMED', 'ICD10CM', 'ICD9CM')
      AND standard_concept = 'S'
)
"""

print("Identifying primary lymphedema cases...")
primary_ids = pd.read_gbq(hereditary_query, project_id=GOOGLE_PROJECT,
                           dialect="standard")["person_id"].tolist()
print(f"  Primary candidates: {len(primary_ids):,}")

print("Identifying secondary lymphedema cases...")
all_lymph_ids = pd.read_gbq(secondary_query, project_id=GOOGLE_PROJECT,
                              dialect="standard")["person_id"].tolist()
secondary_ids = [i for i in all_lymph_ids if i not in primary_ids]
print(f"  Secondary candidates: {len(secondary_ids):,}")

# ── Step 3: Check WGS availability ────────────────────────────────────────────
wgs_query = f"""
SELECT DISTINCT person_id
FROM `{WORKSPACE_CDR}.cb_search_all_events`
WHERE is_standard = 1
  AND concept_id = 1333101   -- Has WGS data
"""

wgs_ids = set(pd.read_gbq(wgs_query, project_id=GOOGLE_PROJECT,
                            dialect="standard")["person_id"].tolist())
print(f"\nParticipants with WGS: {len(wgs_ids):,}")

primary_ids   = [i for i in primary_ids   if i in wgs_ids]
secondary_ids = [i for i in secondary_ids if i in wgs_ids]
print(f"Primary with WGS:   {len(primary_ids):,}")
print(f"Secondary with WGS: {len(secondary_ids):,}")

# ── Step 4: Build case dataframes ─────────────────────────────────────────────
demo_df = demo_df[demo_df["person_id"].isin(wgs_ids)].copy()

primary_df = demo_df[demo_df["person_id"].isin(primary_ids)].copy()
primary_df["condition"] = "PRIMARY"

secondary_df = demo_df[demo_df["person_id"].isin(secondary_ids)].copy()
secondary_df["condition"] = "SECONDARY"

cases_df = pd.concat([primary_df, secondary_df], ignore_index=True)
print(f"\nFinal cases: {len(cases_df):,} "
      f"({len(primary_df):,} primary, {len(secondary_df):,} secondary)")

# ── Step 5: Identify eligible controls ────────────────────────────────────────
# Exclude: any lymphedema Dx, any immunosuppressant
immunosuppressant_query = f"""
SELECT DISTINCT person_id
FROM `{WORKSPACE_CDR}.drug_exposure`
WHERE drug_concept_id IN (
    SELECT concept_id FROM `{WORKSPACE_CDR}.concept`
    WHERE concept_class_id = 'Ingredient'
      AND (concept_name LIKE '%tacrolimus%'
        OR concept_name LIKE '%cyclosporine%'
        OR concept_name LIKE '%mycophenolate%'
        OR concept_name LIKE '%azathioprine%'
        OR concept_name LIKE '%sirolimus%'
        OR concept_name LIKE '%prednisone%'
        OR concept_name LIKE '%methotrexate%')
)
"""

print("\nIdentifying ineligible controls...")
immuno_ids = set(pd.read_gbq(immunosuppressant_query,
                              project_id=GOOGLE_PROJECT,
                              dialect="standard")["person_id"].tolist())
exclude_ids = set(all_lymph_ids) | immuno_ids

controls_df = demo_df[
    ~demo_df["person_id"].isin(exclude_ids) &
    ~demo_df["person_id"].isin(cases_df["person_id"])
].copy()
controls_df["condition"] = "CONTROL"
print(f"  Eligible controls: {len(controls_df):,}")

# ── Step 6: 1:10 matching ─────────────────────────────────────────────────────
def match_cases_to_controls(cases, controls, ratio=10):
    """
    Match cases to controls 1:ratio on race + sex, nearest age.
    Controls are removed from pool once matched (sampling without replacement).
    """
    matched_pairs = []
    controls_pool = controls.copy().reset_index(drop=True)
    controls_pool["age_days"] = (
        pd.to_datetime(controls_pool["birth_datetime"]).apply(
            lambda x: (datetime.now() - x.replace(tzinfo=None)).days
        )
    )

    cases = cases.copy()
    cases["age_days"] = (
        pd.to_datetime(cases["birth_datetime"]).apply(
            lambda x: (datetime.now() - x.replace(tzinfo=None)).days
        )
    )

    for _, case_row in cases.iterrows():
        # Filter pool: same race + sex
        pool = controls_pool[
            (controls_pool["race_concept_id"] == case_row["race_concept_id"]) &
            (controls_pool["sex_at_birth_concept_id"] ==
             case_row["sex_at_birth_concept_id"])
        ].copy()

        if len(pool) == 0:
            print(f"  WARNING: No eligible controls for case {case_row['person_id']}")
            continue

        # Sort by age distance, take up to ratio matches
        pool["age_diff"] = abs(pool["age_days"] - case_row["age_days"])
        pool = pool.sort_values("age_diff").head(ratio)

        # Add matched controls to output
        for _, ctrl_row in pool.iterrows():
            matched_pairs.append({
                "case_person_id":    case_row["person_id"],
                "control_person_id": ctrl_row["person_id"],
                "condition":         case_row["condition"],
                "age_diff_days":     ctrl_row["age_diff"]
            })
            # Remove from pool (sampling without replacement)
            controls_pool = controls_pool[
                controls_pool["person_id"] != ctrl_row["person_id"]
            ]

    return pd.DataFrame(matched_pairs)

print("\nRunning 1:10 matching...")
matched = match_cases_to_controls(cases_df, controls_df, ratio=MATCH_RATIO)
print(f"  Matched pairs: {len(matched):,}")

n_unique_controls = matched["control_person_id"].nunique()
print(f"  Unique controls: {n_unique_controls:,}")
print(f"  Unique cases:    {matched['case_person_id'].nunique():,}")

# ── Step 7: Build final cohort file ──────────────────────────────────────────
all_ids = (
    list(matched["case_person_id"].unique()) +
    list(matched["control_person_id"].unique())
)

final_cohort = demo_df[demo_df["person_id"].isin(all_ids)].copy()
final_cohort["group"] = final_cohort["person_id"].apply(
    lambda pid: "Case" if pid in cases_df["person_id"].values else "Control"
)
final_cohort["condition"] = final_cohort["person_id"].apply(
    lambda pid: cases_df[cases_df["person_id"] == pid]["condition"].values[0]
    if pid in cases_df["person_id"].values else "CONTROL"
)

print(f"\nFinal cohort:")
print(final_cohort.groupby(["condition", "group"]).size().to_string())

# Save
out_file = "lymphedema_cohort_final.csv"
final_cohort.to_csv(out_file, index=False)
os.system(f"gsutil cp {out_file} {WORKSPACE_BUCKET}/")
print(f"\n✓ Saved: {WORKSPACE_BUCKET}/{out_file}")
