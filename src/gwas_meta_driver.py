"""
s5_runA.py — Part A meta driver. For each phenotype: trans-ancestry inverse-variance
FIXED + DerSimonian-Laird RANDOM meta across EUR/AFR/AMR; meta lambda_GC; genome-wide
hit table (5e-8); and pre-registered candidate-locus readout.

Writes: full meta sumstats -> ss/ (VM, bucket-bound; NOT model-visible / not local-egress);
aggregate tables -> results/ (>=20, MAC>=20 already enforced by REGENIE Step 2).
Prints only aggregate summaries.
"""
import os
import pandas as pd
import meta_analysis as M

ANC = ["eur", "afr", "amr"]
PHENOS = ["composite", "opioid", "pain", "allcause"]
os.makedirs("results", exist_ok=True)
os.makedirs("ss/meta", exist_ok=True)

summary = []
cand_frames = []
for ph in PHENOS:
    dfs = {a: M.read_regenie(f"ss/pooled_{a}_{ph}.regenie.gz") for a in ANC}
    # per-ancestry lambda + candidate readout
    for a in ANC:
        d = dfs[a]
        summary.append(dict(pheno=ph, scope=a, n_var=len(d),
                            lam=round(M.lambda_gc(d["P"]), 4),
                            n_gwsig=int((d["P"] < 5e-8).sum())))
    # trans-ancestry meta (FE + RE + Q/I2); reuse a completed meta file if present (resumable)
    meta_path = f"ss/meta/{ph}_meta.tsv.gz"
    if os.path.exists(meta_path):
        meta = pd.read_csv(meta_path, sep="\t")
    else:
        meta = M.meta_analyze([dfs["eur"], dfs["afr"], dfs["amr"]])
        meta.to_csv(meta_path + ".tmp", sep="\t", index=False, compression="gzip")
        os.replace(meta_path + ".tmp", meta_path)  # atomic -> never leaves a truncated file
    lam_fe = M.lambda_gc(meta["p_fe"])
    summary.append(dict(pheno=ph, scope="meta_fe", n_var=len(meta),
                        lam=round(lam_fe, 4), n_gwsig=int((meta["p_fe"] < 5e-8).sum())))
    summary.append(dict(pheno=ph, scope="meta_re", n_var=len(meta),
                        lam=round(M.lambda_gc(meta["p_re"]), 4),
                        n_gwsig=int((meta["p_re"] < 5e-8).sum())))
    # candidate readout on the META (FE p) for every phenotype
    cr = M.candidate_readout(meta.rename(columns={"p_fe": "P"}), pcol="P")
    cr.insert(0, "pheno", ph)
    cand_frames.append(cr)
    # hit table (meta FE)
    hits = M.hit_table(meta.rename(columns={"p_fe": "P"}), pcol="P")
    if len(hits):
        hits.to_csv(f"results/hits_{ph}_meta.csv", index=False)

summ = pd.DataFrame(summary)
summ.to_csv("results/lambda_and_hits_summary.csv", index=False)
cand = pd.concat(cand_frames, ignore_index=True)
cand.to_csv("results/candidate_loci_meta.csv", index=False)

print("=== per-scope lambda_GC + genome-wide (5e-8) hit counts ===")
print(summ.to_string(index=False))
print("\n=== candidate loci (min-p in window, trans-ancestry META FE) ===")
cc = cand.copy()
cc["minp"] = cc["minp"].map(lambda x: f"{x:.2e}" if pd.notna(x) else "NA")
cc["beta"] = cc["beta"].map(lambda x: f"{x:+.3f}" if pd.notna(x) else "NA")
print(cc[["pheno", "gene", "chr", "n_var", "lead", "minp", "beta"]].to_string(index=False))

print("META_DONE")
