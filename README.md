# LANTERN Collaboratory — Lymphedema & Lipedema GWAS Pipeline

Genome-wide association study pipeline for primary lymphedema, secondary
lymphedema, and lipedema using the All of Us Research Program Controlled
Tier Dataset v8.

## Citation

Bowler SA, Hampilos NH, Ndhlovu LC, and on behalf of the LANTERN Collaboratory.
*Divergent Genetic Architectures Underlying Primary and Secondary Lymphedema
and Lipedema: A Genome-wide Association Study.* (In preparation)

---

## Study Overview

Three parallel GWAS conducted on All of Us v8 WGS data (GRCh38, Hail v0.2):

| Condition | Cases | Controls | GW Significant | Suggestive |
|-----------|------:|--------:|---------------:|-----------:|
| Primary lymphedema | 93 | 930 | 1 (chr15:36,333,087 — *MEIS2* proximal) | 30 |
| Secondary lymphedema | 4,336 | 43,360 | 0 | 49 |
| Lipedema | 110 | 1,100 | 2 (*SLC39A3*, chr19p13.3) | 325 |

**Total participants: 49,929**

---

## Repository Structure

```
lantern_gwas/
├── phenotyping/
│   ├── lymphedema_phenotyping.py      # Primary/secondary case ID + 1:10 matching
│   ├── lipedema_phenotyping_v4.py     # Lipedema case ID, tiering, anthropometric QC
│   └── lipedema_controls.py          # Hierarchical control matching
├── gwas/
│   ├── gwas_lymphedema.py            # Primary + secondary GWAS (Hail, Wald)
│   └── gwas_lipedema_v2.py           # Lipedema GWAS (Hail, Wald, female-only)
├── qc/
│   └── merge_and_artifact_detection.py  # Per-chr merge + segmental dup screening
├── figures/
│   ├── manhattan_qq.py               # Figure 2 + Supplementary Figure 1
│   └── slc39a3_gene_plot.py          # Figure 3 — SLC39A3 locus plot
├── utils/
│   ├── annotation.py                 # Ensembl VEP variant annotation
│   ├── pathway_enrichment.py         # Enrichr pathway analysis
│   └── ld_calculation.py             # LD r² between lead variants (Hail)
├── .gitignore
└── README.md
```

---

## Requirements

### Data Access

All analyses were performed on the **All of Us Researcher Workbench** using
the Controlled Tier Dataset v8. Access requires:

1. Registration at [researchallofus.org](https://researchallofus.org)
2. Completion of required ethics and data use training
3. Approval for Controlled Tier access
4. A registered workspace on the Researcher Workbench

These scripts cannot be run outside the All of Us environment — individual-level
genomic data is not publicly available and must remain within the Workbench.

### Software

All analyses run within the All of Us Jupyter notebook environment.
The following are available pre-installed or installable within the Workbench:

| Package | Version | Purpose |
|---------|---------|---------|
| Hail | 0.2 | GWAS, VDS operations, variant QC |
| Python | 3.10 | All scripting |
| pandas | >=1.5 | Data manipulation |
| numpy | >=1.23 | Numerical operations |
| matplotlib | >=3.6 | Figure generation |
| scipy | >=1.9 | Statistical tests, QQ plots |
| requests | >=2.28 | Ensembl VEP API calls |

### Installation within AoU Workbench

No installation is required for Hail — it is configured as part of the
Dataproc cluster when you create a Hail-enabled notebook environment.

For any missing Python packages, install within a notebook cell:

```python
import subprocess
subprocess.run(["pip", "install", "requests", "--quiet"])
```

### Cluster Configuration

**Lymphedema pipeline** (larger cohort, N ~48,000):
- 10 worker nodes, 32 CPUs, 208 GB RAM each
- Driver: 32 CPUs, 416 GB RAM

**Lipedema pipeline** (smaller cohort, N ~1,200):
- 8 worker nodes, 16 CPUs, 104 GB RAM each
- Driver: 16 CPUs, 104 GB RAM

---

## Usage

Scripts should be run in the following order within the All of Us Workbench:

### 1. Phenotyping

```python
# Identify lymphedema cases and match controls
phenotyping/lymphedema_phenotyping.py

# Identify lipedema cases (tiered phenotyping)
phenotyping/lipedema_phenotyping_v4.py

# Match lipedema controls (hierarchical)
phenotyping/lipedema_controls.py
```

### 2. GWAS

```python
# Run lymphedema GWAS (primary + secondary in parallel)
gwas/gwas_lymphedema.py

# Run lipedema GWAS
gwas/gwas_lipedema_v2.py
```

### 3. QC

```python
# Merge per-chromosome TSVs and screen for artifacts
qc/merge_and_artifact_detection.py
```

### 4. Annotation and Enrichment

```python
# Annotate suggestive variants via Ensembl VEP
utils/annotation.py

# Pathway enrichment via Enrichr API
utils/pathway_enrichment.py

# LD between lead variants
utils/ld_calculation.py
```

### 5. Figures

```python
# Manhattan plots (Figure 2) + QQ plots (Supplementary Figure 1)
figures/manhattan_qq.py

# SLC39A3 locus plot (Figure 3)
figures/slc39a3_gene_plot.py
```

---

## Key Implementation Notes

**gwas_lipedema_v2.py contains three important fixes relative to earlier versions:**

1. **chrX source mismatch** — PAR region `filter_rows` must occur *before*
   dosage expressions are defined. Moving it after `annotate_cols` creates a
   new MT object; covariate expressions remain bound to the pre-filter MT,
   causing `source mismatch` in `logistic_regression_rows`.

2. **Firth non-convergence** — Hail's Firth solver fails wholesale at N ~110
   cases. Switched to Wald regression, which is standard GWAS practice at
   this sample size.

3. **Null model explosion** — `sex_at_birth_num` is a constant in the
   female-only lipedema cohort, causing Newton iteration to explode. Removed
   from covariates. Only `age` and `PC1-PC4` are used.

---

## Data

Input data paths are available to registered All of Us researchers within the
Controlled Tier Workbench environment. Key inputs:

- All of Us v8 WGS VDS (GRCh38, Hail format)
- All of Us ancestry PCs (`ancestry_preds.tsv`, controlled tier)
- Cohort CSV files generated by the phenotyping scripts above

Output files are stored in your personal workspace bucket on Google Cloud
Storage, accessible only within your registered Workbench workspace.
Do not share workspace bucket IDs or paths publicly.

---

## License

Code in this repository is released under the **MIT License**.

Individual-level data from the All of Us Research Program is not included
and is not publicly available. Access to the underlying data requires
independent registration and approval through the All of Us Researcher
Workbench.

---

## Support

For questions about the analysis pipeline, please open a GitHub Issue or
contact the Ndhlovu Lab:

- GitHub: [NdhlovuLab](https://github.com/NdhlovuLab)
- Institution: Division of Infectious Disease, Department of Medicine,
  Weill Cornell Medicine, New York, NY

For questions about All of Us data access:
- [All of Us Support](https://support.researchallofus.org)
- [Researcher Workbench Documentation](https://support.researchallofus.org/hc/en-us/categories/360001554391)
