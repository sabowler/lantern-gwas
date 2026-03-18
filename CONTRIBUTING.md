# Contributing to lantern-gwas

Thank you for your interest in contributing to this project. This repository
accompanies a peer-reviewed manuscript and is primarily maintained as a
reproducibility resource. We welcome contributions that improve usability,
fix bugs, or extend the pipeline to new conditions or datasets.

---

## Ways to Contribute

### Bug Reports

If you find a bug — incorrect output, a script that fails to run, or
behavior that doesn't match the documentation — please open a
[GitHub Issue](https://github.com/NdhlovuLab/lantern-gwas/issues) with:

- A clear description of the problem
- The script name and relevant section
- The error message or unexpected output
- Your Hail version and All of Us dataset version (if applicable)

### Code Improvements

We welcome pull requests for:

- Bug fixes
- Improved error handling or logging
- Documentation improvements
- Performance optimizations
- Support for new All of Us dataset versions (e.g. v9+)

### Extending the Pipeline

If you have adapted this pipeline for a related condition or dataset and
would like to contribute your extensions, please open an Issue first to
discuss before submitting a pull request. We are particularly interested
in extensions to:

- Additional lymphatic or adipose conditions
- Other large EHR-linked biobanks (UK Biobank, BioMe, etc.)
- Updated Hail versions

---

## Getting Started

1. **Fork** the repository on GitHub
2. **Clone** your fork locally
3. Create a **feature branch**: `git checkout -b fix/your-fix-name`
4. Make your changes
5. **Test** your changes in an All of Us Workbench environment if possible
6. **Commit** with a clear message (see below)
7. **Push** to your fork and open a **Pull Request**

---

## Commit Message Format

Please use clear, descriptive commit messages. We loosely follow the
conventional commits format:

```
fix: correct chrX PAR filter ordering in lipedema pipeline
feat: add support for Hail v0.2.x annotation API
docs: update cluster configuration for AoU v9
refactor: consolidate QC thresholds into shared config module
```

---

## Code Style

- Python 3.10+
- Follow PEP 8 conventions
- Use descriptive variable names — this is scientific code that others
  will need to read and audit
- Include docstrings on all functions
- Comment any Hail-specific workarounds or non-obvious implementation
  choices (see the chrX PAR fix in `gwas_lipedema_v2.py` as an example)
- Do not hardcode bucket paths, person IDs, or any workspace-specific
  identifiers — use `os.environ["WORKSPACE_BUCKET"]` or equivalent

---

## Data and Privacy

**Never commit individual-level data.** This includes:

- CSV files containing person IDs or phenotype data
- GWAS summary statistics with small cell sizes
- Any file that could be used to re-identify participants

The `.gitignore` is configured to exclude common data file types.
If you are unsure whether a file is safe to commit, do not commit it.

All contributors are expected to comply with the
[All of Us Data Use Agreement](https://www.researchallofus.org/data-tools/data-access/)
when working with or extending this pipeline.

---

## Questions

For questions about the scientific methods, please open a GitHub Issue
or contact the Ndhlovu Lab at Weill Cornell Medicine.

For questions about All of Us data access, see:
[All of Us Support](https://support.researchallofus.org)
