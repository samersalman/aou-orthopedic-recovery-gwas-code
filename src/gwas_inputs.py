"""
s5_part_a.py — Session 5, Part A driver (primary recovery GWAS + trans-ancestry meta)

Builds REGENIE phenotype + covariate files from the frozen Session-3 analytic parquet
(schema verified live 2026-07-20: analytic_pooled.parquet, 12,592 x 44), generates
REGENIE Step-1/Step-2 command lines, computes lambda_GC. AGGREGATE-ONLY at every
model-visible surface: never prints rows; counts round-to-20 (<=20 suppressed);
GWAS reporting is MAC>=20 (REGENIE Step 2 --minMAC 20).

Frozen upstream (do NOT re-derive):
  Phenotypes (already 0/1/NA with strict-control def baked in, S3 DoD MET):
    pheno_primary_composite (PRIMARY), pheno_opioid_only, pheno_pain_only, pheno_allcause
  Covars: age, age2, sex, PC1..PC16, surgical_year, ehr_density, bmi,
          mh_burden = cx_depression|cx_anxiety, pain_burden = cx_pain_region|cx_pain_chronic,
          proc_class one-hot pc_* (pooled only; ref = pc_arthroplasty)
  IDs (S4): FID=0, IID=person_id ; ancestries {eur,afr,amr}; AMR meta-only (RUS)
  Per-cohort scope: pooled parquet is one-row/person (earliest-index cohort), so a per-cohort
    phenotype column = composite restricted to that person's assigned cohort, NA elsewhere.
    Joint REGENIE Step-1 fits all these columns at once (one sample set per ancestry).
"""

import os
import numpy as np
import pandas as pd

ANCESTRIES = ["eur", "afr", "amr"]
PHENO_MAP = {
    "composite": "pheno_primary_composite",
    "opioid": "pheno_opioid_only",
    "pain": "pheno_pain_only",
    "allcause": "pheno_allcause",
}
COHORTS = ["TKA", "THA", "rotator_cuff", "TSA", "cervical_acdf",
           "lumbar_fusion", "lumbar_decompression", "ankle", "ACL"]
PC_COLS = [f"PC{i}" for i in range(1, 17)]
PROC_DUMMIES = ["pc_foot_ankle", "pc_soft_tissue",
                "pc_spine_decompression", "pc_spine_fusion"]   # ref = pc_arthroplasty
BASE_COVARS = (["age", "age2", "sex"] + PC_COLS +
               ["surgical_year", "ehr_density", "bmi", "mh_burden", "pain_burden"])


def round20(n):
    n = int(n)
    if n <= 0:
        return 0
    if n <= 20:
        return "<=20"
    return int(round(n / 20.0) * 20)


def _to01(s):
    if s.dtype == bool:
        return s.astype("float64")
    return pd.to_numeric(s, errors="coerce").astype("float64")


def _encode_sex(s):
    return s.astype(str).map({"Male": 1, "M": 1, "Female": 0, "F": 0}).astype("float64")


def build_inputs(pooled_parquet, out_dir):
    """Write pheno_<anc>.tsv, covar_<anc>.tsv (+_nopc), return aggregate census DataFrame.
    Row-level TSVs -> out_dir on the VM only (never model-visible)."""
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_parquet(pooled_parquet)
    if "qc_pass" in df.columns:
        df = df[df["qc_pass"] == True].copy()

    # derived burden covariates
    df["mh_burden"] = ((_to01(df["cx_depression"]) > 0) | (_to01(df["cx_anxiety"]) > 0)).astype("float64")
    df["pain_burden"] = ((_to01(df["cx_pain_region"]) > 0) | (_to01(df["cx_pain_chronic"]) > 0)).astype("float64")

    # BMI: impute missing to ancestry median + missingness indicator (avoid REGENIE listwise deletion;
    #      BMI is the only covariate with missingness -- ~2,985 rows, EUR ~2,262)
    df["bmi_miss"] = df["bmi"].isna().astype("float64")
    df["bmi"] = df.groupby("ancestry")["bmi"].transform(lambda s: s.fillna(s.median()))

    census = []
    for anc in ANCESTRIES:
        sub = df[df["ancestry"].astype(str).str.upper() == anc.upper()].copy()
        if len(sub) == 0:
            continue
        fid = pd.Series(0, index=sub.index, name="FID")
        iid = sub["person_id"].astype(str).rename("IID")

        # ----- phenotype file -----
        ph = pd.DataFrame({"FID": fid, "IID": iid})
        for short, col in PHENO_MAP.items():
            ph[short] = _to01(sub[col])
        comp = _to01(sub["pheno_primary_composite"])
        for c in COHORTS:
            ph[f"comp_{c}"] = comp.where(sub["cohort"].astype(str) == c, np.nan)
        ph.to_csv(f"{out_dir}/pheno_{anc}.tsv", sep="\t", index=False, na_rep="NA")

        # ----- covariate files (with / without proc_class) -----
        cov = pd.DataFrame({"FID": fid, "IID": iid})
        cov["age"] = sub["age"].astype("float64")
        cov["age2"] = sub["age2"].astype("float64")
        cov["sex"] = _encode_sex(sub["sex"])
        for pc in PC_COLS:
            cov[pc] = sub[pc].astype("float64")
        cov["surgical_year"] = sub["surgical_year"].astype("float64")
        cov["ehr_density"] = sub["ehr_density"].astype("float64")
        cov["bmi"] = sub["bmi"].astype("float64")
        cov["mh_burden"] = sub["mh_burden"]
        cov["pain_burden"] = sub["pain_burden"]
        cov["bmi_miss"] = sub["bmi_miss"]
        cov_nopc = cov.copy()
        for d in PROC_DUMMIES:
            cov[d] = _to01(sub[d]) if d in sub.columns else 0.0
        cov.to_csv(f"{out_dir}/covar_{anc}.tsv", sep="\t", index=False, na_rep="NA")
        cov_nopc.to_csv(f"{out_dir}/covar_{anc}_nopc.tsv", sep="\t", index=False, na_rep="NA")

        # ----- aggregate census (round20) -----
        for short in PHENO_MAP:
            y = ph[short]
            census.append(dict(ancestry=anc, scope="pooled", pheno=short,
                               N=round20(int(y.notna().sum())),
                               cases=round20(int((y == 1).sum())),
                               controls=round20(int((y == 0).sum()))))
        for c in COHORTS:
            y = ph[f"comp_{c}"]
            n = int(y.notna().sum())
            if n == 0:
                continue
            census.append(dict(ancestry=anc, scope=f"cohort-{c}", pheno="composite",
                               N=round20(n),
                               cases=round20(int((y == 1).sum())),
                               controls=round20(int((y == 0).sum()))))
    return pd.DataFrame(census)


# ---- REGENIE command helpers (paths resolve at runtime) ----
def step1_cmd(anc, geno_dir, pheno, covar, out):
    return (f"regenie --step 1 --bed {geno_dir}/step1_pruned_{anc} "
            f"--phenoFile {pheno} --covarFile {covar} --bt --bsize 1000 "
            f"--lowmem --lowmem-prefix /tmp/rg_{anc} --threads 15 --out {out}")


def step2_cmd(anc, geno_dir, pheno, covar, pred, out, pheno_col=None):
    sel = f"--phenoColList {pheno_col} " if pheno_col else ""
    return (f"regenie --step 2 --bgen {geno_dir}/step2_{anc}.bgen "
            f"--sample {geno_dir}/step2_{anc}.sample --ref-first --pred {pred} "
            f"--phenoFile {pheno} --covarFile {covar} {sel}--bt --firth --approx "
            f"--pThresh 0.01 --minMAC 20 --bsize 400 --threads 15 --gz --out {out}")


def lambda_gc(pvals):
    import scipy.stats as st
    p = np.asarray(pvals, float)
    p = p[(p > 0) & (p <= 1)]
    return float(np.median(st.chi2.isf(p, 1)) / st.chi2.ppf(0.5, 1))
