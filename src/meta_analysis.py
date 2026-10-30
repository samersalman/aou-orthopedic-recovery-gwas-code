"""
s5_meta.py — Session 5, Part A analysis backend (post-REGENIE): lambda_GC,
trans-ancestry meta-analysis (fixed + random effects), cross-cohort heterogeneity (Q, I^2),
genome-wide hit table, and pre-registered candidate-locus readout.

All model-visible output is AGGREGATE (betas/SEs/p, lambda, counts). REGENIE Step 2 already
enforced --minMAC 20, so every reported variant has MAC>=20; no small per-genotype cells.
Row-level sumstats stay on the VM / workspace bucket; nothing is downloaded off-Workbench.

REGENIE .regenie(.gz) columns (v3/v4, space-sep):
  CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ [INFO] N TEST BETA SE CHISQ LOG10P EXTRA
  -> effect allele = ALLELE1 (ALT / tested); BETA = log-OR per ALLELE1; LOG10P = -log10(p).
"""

import gzip
import numpy as np
import pandas as pd
import scipy.stats as st

GW = 5e-8
CANDIDATES = {  # pre-registered recovery candidates + positive-control OA loci (GRCh38)
    "ASTN2":   (9, 118_000_000, 120_300_000),
    "GABRB2":  (5, 160_700_000, 161_050_000),
    "OPRM1":   (6, 154_010_000, 154_250_000),
    "COMT":    (22, 19_940_000, 19_970_000),
    "SCN11A":  (3, 38_870_000, 38_950_000),
    "EFEMP1":  (2, 55_860_000, 56_010_000),
    "GDF5":    (20, 35_430_000, 35_460_000),   # posctrl OA
    "COL27A1": (9, 114_100_000, 114_320_000),  # posctrl OA
    "COL11A1": (1, 102_870_000, 103_110_000),  # posctrl OA
    "SLC39A8": (4, 102_250_000, 102_460_000),  # posctrl OA
}


def read_regenie(path):
    """Parse a REGENIE sumstats (.regenie or .regenie.gz). Adds P (from LOG10P).
    Keeps essential columns; row-level -> VM memory only."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        df = pd.read_csv(fh, delim_whitespace=True, comment="#")  # C parser (fast)
    df.columns = [c.upper() for c in df.columns]
    ren = {"CHROM": "CHR", "GENPOS": "POS"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    if "LOG10P" in df.columns:
        df["P"] = np.power(10.0, -df["LOG10P"].astype(float))
    df["CHR"] = pd.to_numeric(df["CHR"].astype(str).str.replace("chr", "", case=False), errors="coerce")
    for c in ("POS", "BETA", "SE", "A1FREQ", "P"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["BETA", "SE", "P"])


def lambda_gc(p):
    p = np.asarray(p, float)
    p = p[(p > 0) & (p <= 1)]
    return float(np.median(st.chi2.isf(p, 1)) / st.chi2.ppf(0.5, 1))


def _align(dfs):
    """Align sumstats on shared IDs via hash reindex (fast; avoids slow string merges) and
    orient BETA to the reference (first) file's ALLELE1. Returns (merged, cols)."""
    common = set(dfs[0]["ID"])
    for d in dfs[1:]:
        common &= set(d["ID"])
    base = dfs[0].drop_duplicates("ID").set_index("ID")
    idx = base.index.intersection(list(common))
    base = base.reindex(idx)
    merged = pd.DataFrame({"ID": idx.values, "CHR": base["CHR"].values, "POS": base["POS"].values,
                           "ALLELE0": base["ALLELE0"].values, "ALLELE1": base["ALLELE1"].values})
    ref_a1 = base["ALLELE1"].values
    cols = []
    for i, d in enumerate(dfs):
        di = d.drop_duplicates("ID").set_index("ID").reindex(idx)
        b = di["BETA"].values.astype(float).copy()
        s = di["SE"].values.astype(float)
        f = di["A1FREQ"].values.astype(float).copy()
        swap = di["ALLELE1"].values != ref_a1
        b[swap] = -b[swap]; f[swap] = 1 - f[swap]
        merged[f"B_{i}"], merged[f"S_{i}"], merged[f"F_{i}"] = b, s, f
        cols.append(i)
    return merged, cols


def meta_analyze(dfs):
    """Inverse-variance FIXED + DerSimonian-Laird RANDOM effects meta across >=2 sumstats.
    Returns an aggregate DataFrame (one row/variant) with Q, I^2, p_het."""
    merged, cols = _align(dfs)
    B = np.vstack([merged[f"B_{i}"].values for i in cols]).T
    S = np.vstack([merged[f"S_{i}"].values for i in cols]).T
    W = 1.0 / (S ** 2)
    sw = np.nansum(W, axis=1)
    beta_fe = np.nansum(W * B, axis=1) / sw
    se_fe = np.sqrt(1.0 / sw)
    z_fe = beta_fe / se_fe
    p_fe = 2 * st.norm.sf(np.abs(z_fe))
    k = np.sum(~np.isnan(B), axis=1)
    Q = np.nansum(W * (B - beta_fe[:, None]) ** 2, axis=1)
    dfree = np.maximum(k - 1, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        I2 = np.where(Q > 0, np.maximum(0, (Q - dfree) / Q) * 100, 0.0)
        p_het = st.chi2.sf(Q, np.maximum(dfree, 1))
        sw2 = np.nansum(W ** 2, axis=1)
        tau2 = np.where(sw - sw2 / sw > 0, np.maximum(0, (Q - dfree) / (sw - sw2 / sw)), 0.0)
    Wr = 1.0 / (S ** 2 + tau2[:, None])
    swr = np.nansum(Wr, axis=1)
    beta_re = np.nansum(Wr * B, axis=1) / swr
    se_re = np.sqrt(1.0 / swr)
    p_re = 2 * st.norm.sf(np.abs(beta_re / se_re))
    out = merged[["ID", "CHR", "POS", "ALLELE0", "ALLELE1"]].copy()
    out["n_anc"] = k
    out["beta_fe"], out["se_fe"], out["p_fe"] = beta_fe, se_fe, p_fe
    out["beta_re"], out["se_re"], out["p_re"] = beta_re, se_re, p_re
    out["Q"], out["I2"], out["p_het"] = Q, I2, p_het
    out["a1freq_mean"] = np.nanmean(np.vstack([merged[f"F_{i}"].values for i in cols]).T, axis=1)
    return out


def hit_table(df, pcol="P", thresh=GW):
    """Genome-wide hits (aggregate). Clumps by nearest within 1 Mb keeping min-p lead."""
    d = df[df[pcol] < thresh].sort_values(pcol)
    leads = []
    for _, r in d.iterrows():
        if all(not (r["CHR"] == L["CHR"] and abs(r["POS"] - L["POS"]) < 1e6) for L in leads):
            leads.append(r)
    return pd.DataFrame(leads)


def candidate_readout(df, pcol="P"):
    """Min-p variant within each pre-registered candidate window (report regardless of sig)."""
    rows = []
    for g, (c, s, e) in CANDIDATES.items():
        sub = df[(df["CHR"] == c) & (df["POS"].between(s, e))]
        if len(sub) == 0:
            rows.append(dict(gene=g, chr=c, n_var=0, lead=None, minp=np.nan, beta=np.nan, se=np.nan))
            continue
        lead = sub.loc[sub[pcol].idxmin()]
        rows.append(dict(gene=g, chr=c, n_var=int(len(sub)), lead=lead["ID"],
                         minp=float(lead[pcol]),
                         beta=float(lead.get("BETA", lead.get("beta_fe", np.nan))),
                         se=float(lead.get("SE", lead.get("se_fe", np.nan)))))
    return pd.DataFrame(rows)
