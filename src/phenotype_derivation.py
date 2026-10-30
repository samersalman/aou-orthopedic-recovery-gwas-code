#!/usr/bin/env python3
# ============================================================================
# s3_part_a.py  --  SESSION 3 / PART A
# Composite PRIMARY maladaptive-recovery phenotype + 3 SECONDARY phenotypes, per cohort.
#   opioid arm      : reuse _transfer/outcomes.py machinery (RxNorm ingredient set +
#                     concept_ancestor descendant expansion + drug_strength + CDC-2022 MME).
#                     Windows RE-POINTED to the composite's two discrete post-op windows
#                     91-180 & 181-365 d; persist_opioid = fill in BOTH (>=2 windows).
#                     Brummett new_persistent_use (periop + 90-180) kept as sensitivity flag.
#   pain arm  (NEW) : region-linked pain/symptom dx via condition_source_concept_id (ICD-10-CM).
#                     McDonald durability: >=2 encounters >30 d apart, both in post-op days 90-365.
#   structural (NEW): same-year (index..index+365) competing-explanation purge (revision /
#                     infection / dislocation / periprosthetic fx / pseudarthrosis / malignancy).
#   compose         : primary = (persist_opioid | persist_pain) & ~struct_excl ; 3 secondaries.
# DUCC 3.1 AGGREGATE-ONLY: every printed number via round20()/prev(); rows stay in-VM.
# Reads /tmp/analytic_rosters.parquet ; writes /tmp/phenotypes_<cohort>.parquet (Part A->Part B handoff).
# ============================================================================
import os, re, subprocess, datetime as dt
import numpy as np, pandas as pd

GOOGLE_PROJECT=(os.environ.get("GOOGLE_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT") or os.environ.get("GCP_PROJECT"))
def wb_resolve(name):
    for args in (["wb","resource","resolve","--name",name],["wb","resolve","--name",name]):
        try:
            out=subprocess.run(args,capture_output=True,text=True,timeout=90); v=(out.stdout or "").strip()
            if v: return v
        except Exception: pass
    return None
CDR=os.environ.get("WORKSPACE_CDR") or wb_resolve("C2025Q4R6")
from google.cloud import bigquery
_bq=bigquery.Client(project=GOOGLE_PROJECT)
def q(sql): return _bq.query(sql.replace("{CDR}",CDR)).result().to_dataframe()
def round20(n):
    try: n=int(n)
    except Exception: return n
    if n==0: return 0
    if n<=20: return "<=20 (suppressed)"
    return int(round(n/20.0)*20)
def prev(k,n):
    r=round20(int(k)); pctv=(100.0*int(k)/n) if n else float('nan')
    return f"{r} ({pctv:.1f}%)"
il=lambda xs: ",".join("'"+x+"'" for x in xs)
ANC_ORDER=['EUR','AFR','AMR','EAS','SAS','MID','OTH']
order=['TKA','THA','rotator_cuff','TSA','cervical_acdf','lumbar_fusion','lumbar_decompression','ankle','ACL']
print("CDR:",CDR,"| GOOGLE_PROJECT:",GOOGLE_PROJECT)

# ---------- load frozen analytic roster (S2->S3 handoff) ----------
rost=pd.read_parquet('/tmp/analytic_rosters.parquet')
rost['person_id']=rost['person_id'].astype('int64'); rost['index_date']=pd.to_datetime(rost['index_date'])
if 'fu365' not in rost.columns: rost['fu365']=True
print("roster: analytic rows",round20(len(rost)),"| persons",round20(rost['person_id'].nunique()),
      "| cohorts",sorted(rost['cohort'].dropna().unique()))
person_ids=sorted(rost['person_id'].unique()); inids="("+",".join(str(i) for i in person_ids)+")"  # in-VM only

# ---------- opioid ingredient set + MME factors (VERBATIM from outcomes.py) ----------
opioids=['codeine','hydrocodone','oxycodone','tramadol','morphine','hydromorphone','oxymorphone',
 'fentanyl','methadone','tapentadol','buprenorphine','meperidine','levorphanol','butorphanol',
 'nalbuphine','pentazocine','dihydrocodeine']
FAC={'codeine':0.15,'hydrocodone':1.0,'morphine':1.0,'oxycodone':1.5,'hydromorphone':4.0,
     'oxymorphone':3.0,'tramadol':0.1,'tapentadol':0.4}

# ---------- region-linked pain ICD-10-CM anchors (source codes; expanded per cohort; verify live) ----------
region_pain={
 'TKA':['M25.56','M23','M17'],                                   # knee pain / internal derangement / knee OA
 'THA':['M25.55','M16'],                                         # hip pain / hip OA
 'rotator_cuff':['M25.51','M75'],                                # shoulder pain / rotator-cuff & shoulder lesions
 'TSA':['M25.51','M75'],                                         # shoulder pain / shoulder lesions
 'cervical_acdf':['M54.2','M54.12','M54.13','M50'],              # cervicalgia / cervical radiculopathy / cervical disc
 'lumbar_fusion':['M54.5','M54.4','M54.16','M54.17','M51','M43.1'],       # LBP / sciatica / lumbar radic / disc / spondylolisthesis
 'lumbar_decompression':['M54.5','M54.4','M54.16','M54.17','M51','M48.06'],# LBP / sciatica / lumbar radic / disc / lumbar stenosis
 'ankle':['M25.57','M25.67'],                                    # ankle/foot pain & stiffness
 'ACL':['M25.56','M23'],                                         # knee pain / internal derangement
}
G89=['G89']  # chronic pain (region-agnostic) -- SUPPORTING signal only (sensitivity toggle)

# ---------- structural / infectious competing-explanation set (ICD-10-CM; same-year purge) ----------
struct_prefixes=['T84.0','T84.5','T81.4','M00','M86','M97','M96.0','C']
#   T84.0 mechanical complication of internal joint prosthesis (incl. dislocation T84.02x)
#   T84.5 infection/inflammatory reaction to internal joint prosthesis ; T81.4 post-procedural infection
#   M00 pyogenic arthritis ; M86 osteomyelitis ; M97 periprosthetic fracture around prosthetic joint
#   M96.0 pseudarthrosis after fusion/arthrodesis (spine cohorts) ; C* any malignancy

# ============================================================================
# (1) ONE-PASS pulls for all roster persons (split per cohort in pandas; window per cohort index_date).
# ============================================================================
FILLS=("WITH oi AS (SELECT concept_id, LOWER(concept_name) ing FROM `{CDR}.concept` "
 "WHERE vocabulary_id='RxNorm' AND concept_class_id='Ingredient' AND LOWER(concept_name) IN ("+il(opioids)+")), "
 "od AS (SELECT ca.descendant_concept_id d, oi.ing FROM `{CDR}.concept_ancestor` ca JOIN oi ON ca.ancestor_concept_id=oi.concept_id) "
 "SELECT de.person_id, od.ing AS ingredient, de.drug_exposure_start_date AS sd, "
 "de.quantity, de.days_supply, ds.amount_value, ds.numerator_value "
 "FROM `{CDR}.drug_exposure` de JOIN od ON de.drug_concept_id=od.d "
 "LEFT JOIN `{CDR}.drug_strength` ds ON ds.drug_concept_id=de.drug_concept_id "
 "WHERE de.person_id IN "+inids)
f_all=q(FILLS); f_all['person_id']=f_all['person_id'].astype('int64')
f_all['sd']=pd.to_datetime(f_all['sd'].astype(str),errors='coerce')
for col in ('quantity','days_supply','amount_value','numerator_value'): f_all[col]=pd.to_numeric(f_all[col],errors='coerce')
f_all['strength_mg']=f_all['amount_value'].where(f_all['amount_value'].notna(),f_all['numerator_value'])
f_all['factor']=f_all['ingredient'].map(FAC)
f_all['mme']=f_all['strength_mg']*f_all['quantity']*f_all['factor']
f_all.loc[f_all['ingredient'].isin(['fentanyl','methadone','buprenorphine']),'mme']=np.nan  # non-standard factors: flagged, excluded from oral MME
print("opioid fills pulled (round20 rows):",round20(len(f_all)))

allpx=sorted({p for v in region_pain.values() for p in v} | set(G89))   # include G89 for the sensitivity variant
like_pain=" OR ".join("c.concept_code LIKE '"+p+"%'" for p in allpx)
PAIN=("SELECT co.person_id, c.concept_code AS code, co.condition_start_date AS sd "
 "FROM `{CDR}.condition_occurrence` co JOIN `{CDR}.concept` c ON c.concept_id=co.condition_source_concept_id "
 "WHERE co.person_id IN "+inids+" AND c.vocabulary_id='ICD10CM' AND ("+like_pain+")")
pain_all=q(PAIN); pain_all['person_id']=pain_all['person_id'].astype('int64')
pain_all['sd']=pd.to_datetime(pain_all['sd'].astype(str),errors='coerce')
print("region-pain dx rows pulled (round20):",round20(len(pain_all)))

like_s=" OR ".join("c.concept_code LIKE '"+p+"%'" for p in struct_prefixes)
STRUCT=("SELECT co.person_id, c.concept_code AS code, co.condition_start_date AS sd "
 "FROM `{CDR}.condition_occurrence` co JOIN `{CDR}.concept` c ON c.concept_id=co.condition_source_concept_id "
 "WHERE co.person_id IN "+inids+" AND c.vocabulary_id='ICD10CM' AND ("+like_s+")")
struct_all=q(STRUCT); struct_all['person_id']=struct_all['person_id'].astype('int64')
struct_all['sd']=pd.to_datetime(struct_all['sd'].astype(str),errors='coerce')
print("structural-excl dx rows pulled (round20):",round20(len(struct_all)))

# ============================================================================
# (2) Per-cohort phenotype build.
# ============================================================================
def mcdonald(days):
    d=sorted(set(int(x) for x in days if pd.notna(x)))
    return (len(d)>=2) and ((d[-1]-d[0])>30)   # >=2 encounters >30 d apart

pheno={}
day_audit=[]
for c in order:
    r=rost[rost.cohort==c][['person_id','index_date','ancestry','fu365']].drop_duplicates('person_id').copy()
    if not len(r):
        print("  (skip empty cohort",c,")"); continue
    ids_c=set(r['person_id']); idx=r[['person_id','index_date']]
    g=r.copy()

    # --- opioid arm ---
    f=f_all[f_all['person_id'].isin(ids_c)].merge(idx,on='person_id',how='inner')
    f['day']=(f['sd']-f['index_date']).dt.days
    def win(lo,hi): return (f['day']>=lo)&(f['day']<hi)
    def flag(mask): s=set(f[mask]['person_id']); return g['person_id'].isin(s)
    g['fill_w1']=flag(win(91,181))          # 91-180 d
    g['fill_w2']=flag(win(181,366))         # 181-365 d
    g['persist_opioid']=g['fill_w1']&g['fill_w2']
    g['periop_fill']=flag(win(-7,15)); g['fill_90_180']=flag(win(90,180))
    g['new_persistent_use']=g['periop_fill']&g['fill_90_180']   # Brummett sensitivity
    # descriptive MME (aggregate only): median daily MME over 0-180 among persons with any oral-MME fill
    mme0=f[win(0,180)].groupby('person_id')['mme'].sum(min_count=1)
    g['daily_mme_0_180']=(g['person_id'].map(mme0)/180.0)

    # --- region pain arm (McDonald >=2 codes >30 d apart, days 90-365) ---
    pref=tuple(region_pain[c])
    sub=pain_all[pain_all['person_id'].isin(ids_c)].merge(idx,on='person_id',how='inner')
    sub['day']=(sub['sd']-sub['index_date']).dt.days
    sub=sub[(sub['day']>=90)&(sub['day']<=365)]
    reg=sub[sub['code'].str.startswith(pref)]                 # region-anchored dx only (primary)
    pc=reg
    pp=reg.groupby('person_id')['day'].apply(mcdonald) if len(reg) else pd.Series(dtype=bool)
    g['persist_pain']=g['person_id'].map(pp).fillna(False)
    # sensitivity comparators (reported only; NOT in the frozen handoff columns):
    def mcd3(days):
        d=sorted(set(int(x) for x in days if pd.notna(x))); return (len(d)>=3) and ((d[-1]-d[0])>30)
    pp3=reg.groupby('person_id')['day'].apply(mcd3) if len(reg) else pd.Series(dtype=bool)
    g['persist_pain_ge3']=g['person_id'].map(pp3).fillna(False)   # >=3 codes variant
    rg=sub[sub['code'].str.startswith(pref+tuple(G89))]           # region + G89 chronic-pain variant
    ppg=rg.groupby('person_id')['day'].apply(mcdonald) if len(rg) else pd.Series(dtype=bool)
    g['persist_pain_g89']=g['person_id'].map(ppg).fillna(False)

    # --- structural/infectious same-year exclusion (index..index+365) ---
    sc_=struct_all[struct_all['person_id'].isin(ids_c)].merge(idx,on='person_id',how='inner')
    sc_['day']=(sc_['sd']-sc_['index_date']).dt.days
    sc_=sc_[(sc_['day']>=0)&(sc_['day']<=365)]
    g['struct_excl']=g['person_id'].isin(set(sc_['person_id']))

    # --- compose ---
    g['pheno_opioid_only']=g['persist_opioid']
    g['pheno_pain_only']=g['persist_pain']
    g['pheno_allcause']=g['persist_opioid']|g['persist_pain']
    g['pheno_primary_composite']=(g['persist_opioid']|g['persist_pain'])&(~g['struct_excl'])
    g['adequate_followup']=g['fu365'].fillna(True).astype(bool)
    g['cohort']=c
    keep=['person_id','cohort','index_date','ancestry','persist_opioid','persist_pain','struct_excl',
          'pheno_primary_composite','pheno_opioid_only','pheno_pain_only','pheno_allcause',
          'adequate_followup','new_persistent_use']
    g[keep].to_parquet(f'/tmp/phenotypes_{c}.parquet',index=False)
    pheno[c]=g
    # window sanity (aggregate): day-offset ranges among contributing events
    w1d=f[win(91,181)]['day']; w2d=f[win(181,366)]['day']; pnd=pc['day']
    day_audit.append((c,round20(len(w1d)),round20(len(w2d)),
                      (int(pnd.min()) if len(pnd) else None),(int(pnd.max()) if len(pnd) else None)))

# ============================================================================
# (3) AGGREGATE-ONLY report: case/control/prevalence per cohort x ancestry (round20; <=20 suppressed).
# ============================================================================
print("\n=== PART A -- phenotype prevalence (aggregate-only; round20; <=20 suppressed) ===")
summary_rows=[]
for c in order:
    if c not in pheno: continue
    d=pheno[c]; n=len(d)
    cases=int(d['pheno_primary_composite'].sum())
    controls=int(((~d['pheno_primary_composite'])&(~d['struct_excl'])&d['adequate_followup']).sum())
    print(f"\n{c}: N={round20(n)} | PRIMARY composite cases={prev(cases,n)} | controls(strict)={round20(controls)} | struct_excl={prev(int(d['struct_excl'].sum()),n)}")
    print(f"    opioid_only={prev(int(d['pheno_opioid_only'].sum()),n)} | pain_only={prev(int(d['pheno_pain_only'].sum()),n)} | allcause={prev(int(d['pheno_allcause'].sum()),n)} | Brummett_npu={prev(int(d['new_persistent_use'].sum()),n)}")
    for a in ANC_ORDER:
        da=d[d['ancestry']==a]; na=len(da)
        if na>20:
            print(f"      {a:4s} N={round20(na):>6} composite={prev(int(da['pheno_primary_composite'].sum()),na)}")
        elif na>0:
            print(f"      {a:4s} N=<=20 (suppressed)")
    row={'cohort':c,'N':round20(n),'cases_composite':round20(cases),
         'pct_composite':round(100.0*cases/n,1) if n else None,
         'opioid_only':round20(int(d['pheno_opioid_only'].sum())),
         'pain_only':round20(int(d['pheno_pain_only'].sum())),
         'allcause':round20(int(d['pheno_allcause'].sum())),
         'struct_excl':round20(int(d['struct_excl'].sum())),
         'brummett_npu':round20(int(d['new_persistent_use'].sum()))}
    summary_rows.append(row)

print("\n=== window sanity (aggregate: contributing-event counts + pain day-range 90-365) ===")
for c,w1n,w2n,pmn,pmx in day_audit:
    print(f"  {c:22s} opioid w1(91-180) n={w1n} w2(181-365) n={w2n} | pain day range [{pmn},{pmx}]")

# --- SENSITIVITY VARIANTS (DoD): report prevalence of each vs the frozen primary (aggregate) ---
print("\n=== SENSITIVITY VARIANTS -- prevalence vs frozen primary (aggregate; round20) ===")
print("  (1) opioid rule: >=2-windows [persist_opioid] vs Brummett periop+90-180 [new_persistent_use]")
print("  (2) pain rule:   >=2 codes/30d [persist_pain] vs >=3 codes/30d [persist_pain_ge3]")
print("  (3) pain anchor: region-only [persist_pain] vs region+G89 [persist_pain_g89]")
print("  (4) window:      90-365 (primary) vs 90-180 short-follow-up (noted; recompute toggle)")
print("  (5) exclusion:   strict [pheno_primary_composite] vs all-cause [pheno_allcause]")
for c in order:
    if c not in pheno: continue
    d=pheno[c]; n=len(d)
    print(f"  {c:22s} persist_opioid={prev(int(d['persist_opioid'].sum()),n)}  npu={prev(int(d['new_persistent_use'].sum()),n)}"
          f"  pain>=2={prev(int(d['persist_pain'].sum()),n)}  pain>=3={prev(int(d['persist_pain_ge3'].sum()),n)}"
          f"  pain+G89={prev(int(d['persist_pain_g89'].sum()),n)}")

pd.DataFrame(summary_rows).to_csv('/tmp/s3_part_a_prevalence.csv',index=False)
print("\nliterature sanity: TKA chronic pain ~20%; persistent opioid ortho ~13.7% / spine ~23.6% (band 14-24%).")
print("wrote /tmp/phenotypes_<cohort>.parquet (9) + /tmp/s3_part_a_prevalence.csv")
print("DONE Part A")
