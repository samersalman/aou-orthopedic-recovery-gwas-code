"""
s5_het.py -- cross-cohort heterogeneity at the pre-registered candidate loci (Option-B anchor).

Reads per-cohort EUR composite sumstats restricted to candidate-region variants, picks a lead
variant per candidate gene (max cohorts reporting, then min pooled-EUR p), computes an
inverse-variance fixed-effect meta ACROSS COHORTS + Cochran's Q + I^2 + p_het, and writes an
aggregate per-locus heterogeneity result file. Low I^2 => shared recovery
architecture; high I^2 => procedure-specific. Aggregate-only; MAC>=20 enforced per cohort.
"""
import os
import numpy as np, pandas as pd
import scipy.stats as st
import meta_analysis as M

COHORTS = ["TKA", "THA", "rotator_cuff", "TSA", "cervical_acdf",
           "lumbar_fusion", "lumbar_decompression", "ankle", "ACL"]

perco = {}
for c in COHORTS:
    p = f"ss/cohort/eur_{c}_composite.regenie.gz"
    if os.path.exists(p):
        d = M.read_regenie(p)
        if len(d):
            perco[c] = d.set_index("ID")

pooled = M.read_regenie("ss/pooled_eur_composite.regenie.gz").set_index("ID")

rows = []
for g, (chrom, s, e) in M.CANDIDATES.items():
    win = pooled[(pooled["CHR"] == chrom) & (pooled["POS"].between(s, e))].copy()
    if len(win) == 0:
        continue
    win["ncoh"] = [sum(vid in perco[c].index for c in perco) for vid in win.index]
    src = win[win["ncoh"] >= 4]
    if len(src) == 0:
        src = win[win["ncoh"] >= 1]
    if len(src) == 0:
        continue
    lead_id = src["P"].idxmin()
    ref_a1 = pooled.loc[lead_id, "ALLELE1"]
    betas, ses, labs = [], [], []
    for c in COHORTS:
        if c in perco and lead_id in perco[c].index:
            r = perco[c].loc[lead_id]
            b, se = float(r["BETA"]), float(r["SE"])
            if r["ALLELE1"] != ref_a1:
                b = -b
            if np.isfinite(b) and np.isfinite(se) and se > 0:
                betas.append(b); ses.append(se); labs.append(c)
    if len(betas) < 2:
        continue
    betas, ses = np.array(betas), np.array(ses)
    w = 1.0 / ses ** 2
    bfe = float(np.sum(w * betas) / np.sum(w))
    sefe = float(np.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (betas - bfe) ** 2))
    dfree = len(betas) - 1
    I2 = max(0.0, (Q - dfree) / Q) * 100 if Q > 0 else 0.0
    pQ = float(st.chi2.sf(Q, max(dfree, 1)))
    verdict = "shared" if I2 < 40 else ("moderate" if I2 < 75 else "procedure-specific")
    rows.append(dict(gene=g, lead=lead_id, chr=chrom, n_cohorts=len(betas),
                     beta_meta=round(bfe, 4), se_meta=round(sefe, 4),
                     Q=round(Q, 2), I2=round(I2, 1), p_het=pQ, verdict=verdict))

het = pd.DataFrame(rows)
het.to_csv("results/cross_cohort_heterogeneity.csv", index=False)
print("=== cross-cohort heterogeneity at candidate loci (EUR composite) ===")
hp = het.copy()
hp["p_het"] = hp["p_het"].map(lambda x: f"{x:.2e}")
print(hp.to_string(index=False))

print("HET_DONE")
