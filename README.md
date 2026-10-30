# All of Us orthopedic recovery GWAS analysis code

This repository contains analysis code for the manuscript *Genome-wide association study of maladaptive recovery after elective orthopedic surgery in the All of Us Research Program*.

## Scope

Included code covers cohort construction, phenotype and covariate derivation, genotype quality control, ancestry-specific REGENIE inputs, trans-ancestry meta-analysis, candidate-locus heterogeneity, LD score regression preparation, and published polygenic-score scoring.

This repository intentionally contains no participant-level data, summary statistics, results, figures, formatted manuscript tables, manuscript files, or manuscript/figure/table-generation utilities. CSV and TSV files produced at run time are analytical outputs, not bundled data.

All participant-level operations must run inside an authorized All of Us Researcher Workbench Controlled Tier workspace. Users must supply their own authorized environment variables and source data. Nothing in this repository grants data access.

## Layout

- `src/cohort_roster.py`: nine-procedure cohort construction, eligibility, ancestry, and follow-up roster.
- `src/phenotype_derivation.py`: persistent opioid-use, region-specific pain, complication-exclusion, and sensitivity phenotypes.
- `src/covariate_derivation.py`: demographic, clinical, BMI, EHR-density, ancestry-PC, and sample-QC derivation.
- `notebooks/05_genomic_prep.ipynb`: ACAF extraction/QC, ancestry-specific PLINK2 export, and positive-control preparation.
- `src/gwas_inputs.py`: REGENIE phenotype/covariate input construction and command generation.
- `src/meta_analysis.py` and `src/gwas_meta_driver.py`: fixed/random-effects meta-analysis, calibration, genome-wide hits, and candidate-locus readout.
- `src/candidate_heterogeneity.py`: cross-procedure Cochran Q and I-squared analyses.
- `src/ldsc_prep.py`: HapMap3 annotation and LD score regression input/log parsing.
- `src/prs_scoring.py`: scoring of PGS002767, PGS002763, PGS004926, and PGS003585 within the secure workspace.

## Required environment

Python dependencies are listed in `requirements.txt`. The workflow also requires the All of Us Workbench command-line environment, PLINK2, REGENIE 4.1.2, BGENIX, and LDSC. Set `GOOGLE_PROJECT`, `WORKSPACE_CDR`, and `WORKSPACE_BUCKET` to resources in the authorized workspace. Do not replace these with public or local participant-level files.

## Reproducibility and disclosure controls

The code rounds counts to the nearest 20 and suppresses counts of 20 or fewer at model-visible output surfaces. Participant-level intermediates are written only to `/tmp` or the authorized workspace bucket. Review all outputs under current All of Us dissemination policy before export.

## License

MIT. See `LICENSE`.

