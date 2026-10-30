"""
s5_partB.py -- Session 5 Part B backend (genetic-architecture backbone).

Prepares Part A REGENIE sumstats for LDSC (annotates rsID from the LD-score reference,
which REGENIE chr:pos:a1:a2 IDs lack), parses LDSC h2/rg logs, and assembles the
candidate-gene deep-dive + the "Option C floor" summary. Aggregate-only; sumstats stay
in the bucket. LDSC itself (py2.7 conda env) is invoked from the shell; this module does
the munge-prep + output parsing + reporting.

LDSC needs SNP(rsID), A1(effect), A2, N, and a signed statistic. Our REGENIE output has
CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P. We map (chr,bp)->rsID
from the eur_w_ld_chr *.l2.ldscore.gz files (cols include CHR, SNP, BP), restricting to the
HapMap3 SNP set LDSC uses -- exactly the SNPs needed for h2.
"""
import glob
import gzip
import re
import numpy as np
import pandas as pd
import meta_analysis as M


def build_snp_map(ldscore_dir):
    """chr:pos -> rsID from LDSC per-chr .l2.ldscore.gz (has CHR, SNP, BP)."""
    m = {}
    for f in sorted(glob.glob(f"{ldscore_dir}/*.l2.ldscore.gz")):
        d = pd.read_csv(f, sep=r"\s+", usecols=lambda c: c.upper() in ("CHR", "SNP", "BP"))
        d.columns = [c.upper() for c in d.columns]
        for chrom, bp, snp in zip(d["CHR"], d["BP"], d["SNP"]):
            m[(int(chrom), int(bp))] = snp
    return m


def munge_prep(regenie_path, snp_map, out_path):
    """Annotate a REGENIE sumstats with rsID (restrict to HapMap3 map) -> LDSC-ready TSV.
    Columns: SNP A1 A2 N BETA SE P  (A1 = ALLELE1 = effect allele)."""
    d = M.read_regenie(regenie_path)
    d["rsid"] = [snp_map.get((int(c), int(p))) for c, p in zip(d["CHR"], d["POS"])]
    d = d.dropna(subset=["rsid"])
    out = pd.DataFrame({
        "SNP": d["rsid"], "A1": d["ALLELE1"], "A2": d["ALLELE0"],
        "N": d["N"] if "N" in d.columns else np.nan,
        "BETA": d["BETA"], "SE": d["SE"], "P": d["P"],
    })
    out.to_csv(out_path, sep="\t", index=False)
    return len(out)


def h2_from_log(path):
    """Parse 'Total Observed scale h2: X (SE)' from an LDSC h2 log."""
    txt = open(path).read()
    m = re.search(r"Total Observed scale h2:\s*([\-0-9.]+)\s*\(([\-0-9.]+)\)", txt)
    lam = re.search(r"Lambda GC:\s*([\-0-9.]+)", txt)
    rat = re.search(r"Ratio:\s*([\-0-9.NAna]+)", txt)
    if not m:
        return None
    h2, se = float(m.group(1)), float(m.group(2))
    z = h2 / se if se > 0 else np.nan
    return dict(h2=h2, se=se, z=round(z, 2),
                lambda_gc=float(lam.group(1)) if lam else np.nan,
                ratio=(rat.group(1) if rat else "NA"))


def rg_from_log(path):
    """Parse the rg summary table at the end of an LDSC --rg log."""
    lines = open(path).read().splitlines()
    hdr = next((i for i, l in enumerate(lines) if l.split()[:2] == ["p1", "p2"]), None)
    if hdr is None:
        return pd.DataFrame()
    rows = []
    for l in lines[hdr + 1:]:
        parts = l.split()
        if len(parts) < 6 or parts[0].startswith("Analysis"):
            break
        rows.append(parts)
    cols = lines[hdr].split()
    df = pd.DataFrame(rows, columns=cols[:len(rows[0])]) if rows else pd.DataFrame()
    return df


# Public reference sumstats for rg (public, non-AoU; download to reference/). Verify URLs live.
RG_TRAITS = {
    "OA_allOA":   "OA (GO consortium / Tachmazidou)",
    "chronic_pain": "multisite chronic pain (UKB)",
    "OUD":        "opioid use disorder (PGC/MVP)",
    "depression": "MDD (PGC)",
    "BMI":        "BMI (GIANT/UKB)",
}
