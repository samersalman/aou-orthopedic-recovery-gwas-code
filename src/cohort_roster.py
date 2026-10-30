#!/usr/bin/env python3
# ============================================================================
# s3_roster.py  --  regenerate the frozen S2 analytic roster IN-VM (byte-equivalent
# to s2_part_a.py + s2_part_b.py), using the FROZEN code sets from
# decisions/2026-07-19-cohort-definitions.md, then VALIDATE vs S2 logged aggregates.
# Writes /tmp/base_cohorts.parquet + /tmp/analytic_rosters.parquet (S2->S3 handoff).
# DUCC 3.1 AGGREGATE-ONLY.
# ============================================================================
import os, re, subprocess
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
def sc(s): return s.value_counts(dropna=False).apply(round20)
il=lambda xs: ",".join("'"+x+"'" for x in xs)
order=['TKA','THA','rotator_cuff','TSA','cervical_acdf','lumbar_fusion','lumbar_decompression','ankle','ACL']
ANC_ORDER=['EUR','AFR','AMR','EAS','SAS','MID','OTH']

# --- FROZEN 9-cohort code sets (decisions/2026-07-19-cohort-definitions.md; = s2_part_a COH) ---
COH={
 'TKA':{'cpt':['27447'],'pcs':['0SRC','0SRD']},
 'THA':{'cpt':['27130'],'pcs':['0SR9','0SRB'],'excl_cpt':['27125']},
 'rotator_cuff':{'cpt':['23410','23412','29827'],'pcs':[]},
 'TSA':{'cpt':['23472'],'pcs':['0RRJ','0RRK']},
 'cervical_acdf':{'cpt':['22551','22552','22554','22590','22595','22600','22614'],'pcs':['0RG0','0RG1','0RG2','0RG4']},
 'lumbar_fusion':{'cpt':['22558','22612','22630','22632','22633','22634'],'pcs':['0SG0','0SG1','0SG3']},
 'lumbar_decompression':{'cpt':['63005','63012','63017','63030','63035','63047','63048'],'pcs':['0SB2','0SB4']},
 'ankle':{'cpt':['27870','27702'],'pcs':[]},
 'ACL':{'cpt':['29888'],'pcs':[]},
}
whens=[]
for lab in order:
    spec=COH[lab]
    if spec['cpt']:
        cond="vocabulary_id='CPT4' AND concept_code IN ("+il(spec['cpt'])+")"
        if spec.get('excl_cpt'): cond+=" AND concept_code NOT IN ("+il(spec['excl_cpt'])+")"
        whens.append("WHEN "+cond+" THEN '"+lab+"'")
    if spec['pcs']:
        whens.append("WHEN vocabulary_id='ICD10PCS' AND SUBSTR(concept_code,1,4) IN ("+il(spec['pcs'])+") THEN '"+lab+"'")
all_cpt=sorted({c for v in COH.values() for c in v['cpt']})
all_pcs=sorted({p for v in COH.values() for p in v['pcs']})
SRC_CTE=("WITH cohort_src AS (SELECT concept_id, vocabulary_id, concept_code, CASE "
 +" ".join(whens)+" END AS cohort FROM `{CDR}.concept` "
 "WHERE (vocabulary_id='CPT4' AND concept_code IN ("+il(all_cpt)+")) "
 "OR (vocabulary_id='ICD10PCS' AND SUBSTR(concept_code,1,4) IN ("+il(all_pcs)+"))) ")

# --- base membership + eligibility (>=18, >=365d pre-index EHR, >=2 pre-index visits) ---
ELIG=(SRC_CTE+
 ", po AS (SELECT p.person_id,p.procedure_date,c.cohort FROM `{CDR}.procedure_occurrence` p "
 "JOIN cohort_src c ON p.procedure_source_concept_id=c.concept_id WHERE c.cohort IS NOT NULL), "
 "cohort_long AS (SELECT person_id,cohort,MIN(procedure_date) index_date FROM po GROUP BY person_id,cohort), "
 "ob AS (SELECT person_id, MIN(observation_period_start_date) min_start FROM `{CDR}.observation_period` GROUP BY person_id), "
 "vis AS (SELECT cl.person_id, cl.cohort, COUNT(DISTINCT vo.visit_occurrence_id) n_pre_visits "
 "FROM cohort_long cl JOIN `{CDR}.visit_occurrence` vo ON vo.person_id=cl.person_id "
 "AND vo.visit_start_date < cl.index_date GROUP BY cl.person_id, cl.cohort) "
 "SELECT cl.person_id, cl.cohort, cl.index_date, per.year_of_birth, "
 "EXTRACT(YEAR FROM cl.index_date)-per.year_of_birth AS age_at_index, "
 "ob.min_start, COALESCE(vis.n_pre_visits,0) n_pre_visits "
 "FROM cohort_long cl JOIN `{CDR}.person` per ON per.person_id=cl.person_id "
 "LEFT JOIN ob ON ob.person_id=cl.person_id "
 "LEFT JOIN vis ON vis.person_id=cl.person_id AND vis.cohort=cl.cohort")
raw=q(ELIG); raw['person_id']=raw['person_id'].astype('int64')
for c in ('index_date','min_start'): raw[c]=pd.to_datetime(raw[c].astype(str),errors='coerce')
raw['f_age']=raw['age_at_index']>=18
raw['f_pre']=raw['min_start']<=(raw['index_date']-pd.Timedelta(days=365))
raw['f_visit']=raw['n_pre_visits']>=2
raw['eligible']=raw['f_age']&raw['f_pre']&raw['f_visit']
base=raw[raw['eligible']][['person_id','cohort','index_date']].copy()
base['index_date_pooled']=base.groupby('person_id')['index_date'].transform('min')
base.to_parquet('/tmp/base_cohorts.parquet',index=False)
print("BASE per-cohort N (distinct persons; round20):")
print(sc(base['cohort']).reindex(order).to_string())
print("distinct persons any cohort:",round20(base['person_id'].nunique()),"| person-cohort rows:",round20(len(base)))

# --- ancestry + srWGS membership (ancestry_preds.tsv; parse_pcs verbatim) ---
GEN=wb_resolve("v9-genomics-folder") or "gs://vwb-aou-datasets-controlled/v9"
ANC=GEN.rstrip('/')+"/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
subprocess.run(["gsutil","-u",GOOGLE_PROJECT,"cp",ANC,"/tmp/ancestry.tsv"],check=True)
anc=pd.read_csv('/tmp/ancestry.tsv',sep='\t',dtype=str)
anc2=anc[['research_id','ancestry_pred']].rename(columns={'research_id':'person_id','ancestry_pred':'ancestry'}).drop_duplicates('person_id')
anc2['person_id']=anc2['person_id'].astype('int64'); anc2['ancestry']=anc2['ancestry'].astype(str).str.upper()
anc_ids=set(anc2['person_id'])

# --- funnel: opioid-naive [index-365,index-30] + >=365d follow-up ---
op=['codeine','hydrocodone','oxycodone','tramadol','morphine','hydromorphone','oxymorphone','fentanyl',
 'methadone','tapentadol','buprenorphine','meperidine','levorphanol','butorphanol','nalbuphine',
 'pentazocine','dihydrocodeine','opium']
FLAGS=(SRC_CTE+
 ", po AS (SELECT p.person_id,p.procedure_date,c.cohort FROM `{CDR}.procedure_occurrence` p "
 "JOIN cohort_src c ON p.procedure_source_concept_id=c.concept_id WHERE c.cohort IS NOT NULL), "
 "cl AS (SELECT person_id,cohort,MIN(procedure_date) index_date FROM po GROUP BY person_id,cohort), "
 "oi AS (SELECT concept_id FROM `{CDR}.concept` WHERE vocabulary_id='RxNorm' AND concept_class_id='Ingredient' AND LOWER(concept_name) IN ("+il(op)+")), "
 "od AS (SELECT DISTINCT ca.descendant_concept_id d FROM `{CDR}.concept_ancestor` ca JOIN oi ON ca.ancestor_concept_id=oi.concept_id), "
 "pre AS (SELECT cl.person_id, cl.cohort, COUNT(*) n_pre FROM cl JOIN `{CDR}.drug_exposure` de ON de.person_id=cl.person_id "
 "JOIN od ON de.drug_concept_id=od.d WHERE de.drug_exposure_start_date BETWEEN DATE_SUB(cl.index_date,INTERVAL 365 DAY) "
 "AND DATE_SUB(cl.index_date,INTERVAL 30 DAY) GROUP BY cl.person_id, cl.cohort), "
 "ob AS (SELECT cl.person_id, cl.cohort, cl.index_date, MAX(op2.observation_period_end_date) max_end "
 "FROM cl LEFT JOIN `{CDR}.observation_period` op2 ON op2.person_id=cl.person_id GROUP BY cl.person_id, cl.cohort, cl.index_date) "
 "SELECT cl.person_id, cl.cohort, cl.index_date, COALESCE(pre.n_pre,0) n_pre, ob.max_end "
 "FROM cl LEFT JOIN pre ON pre.person_id=cl.person_id AND pre.cohort=cl.cohort "
 "LEFT JOIN ob ON ob.person_id=cl.person_id AND ob.cohort=cl.cohort")
fl=q(FLAGS); fl['person_id']=fl['person_id'].astype('int64')
for c in ('index_date','max_end'): fl[c]=pd.to_datetime(fl[c].astype(str),errors='coerce')
fl['naive']=fl['n_pre'].fillna(0)==0
fl['fu365']=fl['max_end']>=(fl['index_date']+pd.Timedelta(days=365))

m=base.merge(fl[['person_id','cohort','naive','fu365']],on=['person_id','cohort'],how='left')
for c in ('naive','fu365'): m[c]=m[c].fillna(False)
m['srwgs']=m['person_id'].isin(anc_ids)
m=m.merge(anc2[['person_id','ancestry']],on='person_id',how='left')
m['ancestry']=m['ancestry'].where(m['srwgs'])
m['analytic']=m['srwgs']&m['naive']&m['fu365']
rost=m[m['analytic']][['person_id','index_date','cohort','ancestry','srwgs','naive','fu365']].copy()
rost.to_parquet('/tmp/analytic_rosters.parquet',index=False)

# --- VALIDATE vs S2 logged aggregates ---
an_ids=set(rost['person_id'])
print("\n=== ROSTER VALIDATION vs S2 (round20) ===")
print("analytic person-cohort rows:",round20(len(rost)),"| distinct analytic persons:",round20(len(an_ids)),
      "  [S2: ~14,800 rows / ~12,620 persons]")
pooled_anc=rost.drop_duplicates('person_id').groupby('ancestry').size()
print("analytic distinct persons by ancestry (round20):")
for a in ANC_ORDER: print(f"  {a:4s} {round20(int(pooled_anc.get(a,0)))}")
print("  [S2: EUR 9,760 / AFR 1,420 / AMR 1,200 / EAS 120 / SAS 40 / MID 40]")
print("per-cohort analytic N (round20):")
print(sc(rost['cohort']).reindex(order).to_string())
print("DONE roster -> /tmp/base_cohorts.parquet + /tmp/analytic_rosters.parquet")
