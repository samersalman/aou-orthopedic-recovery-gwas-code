#!/usr/bin/env python3
# ============================================================================
# s3_part_b.py  --  SESSION 3 / PART B
# Covariates + genetic ancestry/PC1-16 + sample QC -> one analytic parquet per cohort
# + pooled analytic parquet. REGENIE-ready phenotype/covariate/PC files.
#   covariates : REUSE _transfer/covariates.py core (age, sex, BMI baseline-window mean,
#                baseline ICD-10-CM conditions -> comorbidity_count) + NEW age2, surgical_year,
#                ehr_density, procedure_class dummy (pooled model).
#   ancestry   : REUSE _transfer/exposure.py parse_pcs 16-PC loader VERBATIM (research_id=person_id).
#   sample QC  : read /tmp/qc_flags.parquet (built live from aux/ QC files) -> drop qc_pass==False;
#                RETAIN related (flag only; REGENIE handles relatedness).
# DUCC 3.1 AGGREGATE-ONLY: Table-1 cells >=20 (round20); rows stay in-VM.
# Reads /tmp/phenotypes_<cohort>.parquet (Part A) + /tmp/qc_flags.parquet ;
# writes /tmp/analytic_<cohort>.parquet (9) + /tmp/analytic_pooled.parquet (S3->S4 handoff).
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
def pct(k,n): return f"{round20(int(k))} ({100.0*int(k)/n:.1f}%)" if n else "NA"
def msd(x):
    x=pd.to_numeric(x,errors='coerce').dropna()
    return f"{x.mean():.1f} +/- {x.std():.1f}" if len(x)>20 else "<=20"
ANC_ORDER=['EUR','AFR','AMR','EAS','SAS','MID','OTH']
order=['TKA','THA','rotator_cuff','TSA','cervical_acdf','lumbar_fusion','lumbar_decompression','ankle','ACL']
PROC_CLASS={'TKA':'arthroplasty','THA':'arthroplasty','TSA':'arthroplasty',
 'lumbar_fusion':'spine_fusion','cervical_acdf':'spine_fusion',
 'lumbar_decompression':'spine_decompression','rotator_cuff':'soft_tissue','ACL':'soft_tissue',
 'ankle':'foot_ankle'}
print("CDR:",CDR,"| GOOGLE_PROJECT:",GOOGLE_PROJECT)

# ============================================================================
# (1) Load Part A phenotypes (all 9; long person-cohort).
# ============================================================================
phe=pd.concat([pd.read_parquet(f'/tmp/phenotypes_{c}.parquet') for c in order
               if os.path.exists(f'/tmp/phenotypes_{c}.parquet')],ignore_index=True)
phe['person_id']=phe['person_id'].astype('int64'); phe['index_date']=pd.to_datetime(phe['index_date'])
ids=sorted(phe['person_id'].unique()); inids="("+",".join(str(i) for i in ids)+")"   # in-VM only
print("phenotypes loaded | person-cohort rows",round20(len(phe)),"| persons",round20(phe['person_id'].nunique()))

# ============================================================================
# (2) Demographics (age at index per cohort, sex) -- covariates.py verbatim core.
# ============================================================================
dem=q("SELECT p.person_id, p.year_of_birth, c.concept_name AS sex FROM `{CDR}.person` p "
 "LEFT JOIN `{CDR}.concept` c ON c.concept_id=p.sex_at_birth_concept_id WHERE p.person_id IN "+inids)
dem['person_id']=dem['person_id'].astype('int64')
cov=phe.merge(dem,on='person_id',how='left')
cov['age']=cov['index_date'].dt.year - pd.to_numeric(cov['year_of_birth'],errors='coerce')
cov['age2']=cov['age']**2
cov['sex']=cov['sex'].where(cov['sex'].isin(['Male','Female']),'Other/Unknown')
cov['surgical_year']=cov['index_date'].dt.year

# ============================================================================
# (3) BMI baseline-window mean [index-365, index+30] (per person-cohort).
# ============================================================================
bmi=q("SELECT m.person_id, m.measurement_date AS md, m.value_as_number AS v FROM `{CDR}.measurement` m "
 "JOIN `{CDR}.concept` c ON c.concept_id=m.measurement_concept_id "
 "WHERE m.person_id IN "+inids+" AND LOWER(c.concept_name) LIKE '%body mass index%' "
 "AND m.value_as_number BETWEEN 10 AND 80")
bmi['person_id']=bmi['person_id'].astype('int64'); bmi['md']=pd.to_datetime(bmi['md'].astype(str),errors='coerce')
bt=cov[['person_id','cohort','index_date']].merge(bmi,on='person_id',how='left')
bt=bt[(bt['md']>=bt['index_date']-pd.Timedelta(days=365))&(bt['md']<=bt['index_date']+pd.Timedelta(days=30))]
bmi_pc=bt.groupby(['person_id','cohort'])['v'].mean().rename('bmi').reset_index()
cov=cov.merge(bmi_pc,on=['person_id','cohort'],how='left')

# ============================================================================
# (4) Baseline conditions via ICD-10-CM source prefixes (condition_start_date <= index) -- covariates.py.
# ============================================================================
grp={'depression':['F32','F33'],'anxiety':['F41'],'sud_oud':['F11','F10','F12','F13','F14','F15','F16','F18','F19'],
 'pain_region':['M25','M54','M75','M51','M50','M23','M16','M17','M48'],'pain_chronic':['G89'],
 'tobacco':['F17','Z87891'],'diabetes':['E10','E11','E13'],'obesity_dx':['E66'],
 'chf':['I50'],'copd':['J44'],'ckd':['N18']}
allpx=sorted({p for v in grp.values() for p in v})
like=" OR ".join("c.concept_code LIKE '"+p+"%'" for p in allpx)
cond=q("SELECT co.person_id, c.concept_code AS code, co.condition_start_date AS sd "
 "FROM `{CDR}.condition_occurrence` co JOIN `{CDR}.concept` c ON c.concept_id=co.condition_source_concept_id "
 "WHERE co.person_id IN "+inids+" AND c.vocabulary_id='ICD10CM' AND ("+like+")")
cond['person_id']=cond['person_id'].astype('int64'); cond['sd']=pd.to_datetime(cond['sd'].astype(str),errors='coerce')
ct=cov[['person_id','cohort','index_date']].merge(cond,on='person_id',how='left')
ct=ct[ct['sd']<=ct['index_date']]
flags=cov[['person_id','cohort']].drop_duplicates().copy()
for name,px in grp.items():
    m=ct['code'].str.startswith(tuple(px))
    s=set(zip(ct[m]['person_id'],ct[m]['cohort']))
    flags['cx_'+name]=[(pid,coh) in s for pid,coh in zip(flags['person_id'],flags['cohort'])]
cov=cov.merge(flags,on=['person_id','cohort'],how='left')
como=['cx_diabetes','cx_chf','cx_copd','cx_ckd','cx_depression','cx_anxiety','cx_sud_oud']
cov['comorbidity_count']=cov[como].sum(axis=1)

# ============================================================================
# (5) ehr_density = distinct pre-index visit-days / pre-index observation-years.
# ============================================================================
vis=q("SELECT vo.person_id, vo.visit_start_date AS vd FROM `{CDR}.visit_occurrence` vo WHERE vo.person_id IN "+inids)
vis['person_id']=vis['person_id'].astype('int64'); vis['vd']=pd.to_datetime(vis['vd'].astype(str),errors='coerce')
vt=cov[['person_id','cohort','index_date']].merge(vis,on='person_id',how='left')
vt=vt[vt['vd']<vt['index_date']]
vday=vt.groupby(['person_id','cohort'])['vd'].nunique().rename('n_pre_vday').reset_index()
obs=q("SELECT person_id, MIN(observation_period_start_date) s FROM `{CDR}.observation_period` "
 "WHERE person_id IN "+inids+" GROUP BY person_id")
obs['person_id']=obs['person_id'].astype('int64'); obs['s']=pd.to_datetime(obs['s'].astype(str),errors='coerce')
cov=cov.merge(vday,on=['person_id','cohort'],how='left').merge(obs,on='person_id',how='left')
cov['pre_years']=((cov['index_date']-cov['s']).dt.days/365.25).clip(lower=0.5)
cov['ehr_density']=(cov['n_pre_vday'].fillna(0)/cov['pre_years']).round(3)

# ============================================================================
# (6) procedure_class dummy (pooled model).
# ============================================================================
cov['procedure_class']=cov['cohort'].map(PROC_CLASS)

# ============================================================================
# (7) Genetic ancestry + PC1-16 (exposure.py parse_pcs loader VERBATIM).
# ============================================================================
GEN=wb_resolve("v9-genomics-folder") or "gs://vwb-aou-datasets-controlled/v9"
ANC=GEN.rstrip('/')+"/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
subprocess.run(["gsutil","-u",GOOGLE_PROJECT,"cp",ANC,"/tmp/ancestry.tsv"],check=True)
anc=pd.read_csv('/tmp/ancestry.tsv',sep='\t',dtype=str)
def parse_pcs(x):
    s=str(x).strip().strip('[]'); vals=[]
    for p in re.split(r'[,\s]+', s.strip()):
        if p=='': continue
        try: vals.append(float(p))
        except Exception: pass
    return (vals+[np.nan]*16)[:16]
pcdf=pd.DataFrame(anc['pca_features'].map(parse_pcs).tolist(),columns=[f'PC{i}' for i in range(1,17)])
anc2=pd.concat([anc[['research_id','ancestry_pred']].reset_index(drop=True),pcdf],axis=1)
anc2=anc2.rename(columns={'research_id':'person_id','ancestry_pred':'ancestry_file'}).drop_duplicates('person_id')
anc2['person_id']=anc2['person_id'].astype('int64'); anc2['ancestry_file']=anc2['ancestry_file'].astype(str).str.upper()
cov=cov.merge(anc2,on='person_id',how='left')
mism=int((cov['ancestry'].notna()&cov['ancestry_file'].notna()&(cov['ancestry']!=cov['ancestry_file'])).sum())
print("ancestry reconcile (roster vs file) mismatches (round20):",round20(mism))
cov['ancestry']=cov['ancestry'].fillna(cov['ancestry_file'])
pc_cols=[f'PC{i}' for i in range(1,17)]
pc_ok=cov[pc_cols].notna().all(axis=1)
print("all-16-PCs present (round20 person-cohort rows):",round20(int(pc_ok.sum())),
      "| PC-missing rows:",round20(int((~pc_ok).sum())))

# ============================================================================
# (8) Sample QC: read /tmp/qc_flags.parquet built live from aux/ (person_id, qc_pass, related).
#     Drop qc_pass==False genomes; RETAIN related (flag only). Default all-pass w/ caveat if absent.
# ============================================================================
if os.path.exists('/tmp/qc_flags.parquet'):
    qcf=pd.read_parquet('/tmp/qc_flags.parquet'); qcf['person_id']=qcf['person_id'].astype('int64')
    cov=cov.merge(qcf[['person_id','qc_pass','related']],on='person_id',how='left')
    cov['qc_pass']=cov['qc_pass'].fillna(True).astype(bool)
    cov['related']=cov['related'].fillna(False).astype(bool)
    print("QC flags applied from /tmp/qc_flags.parquet")
else:
    cov['qc_pass']=True; cov['related']=False
    print("WARNING: /tmp/qc_flags.parquet absent -> qc_pass defaulted TRUE (finalize QC live before freezing).")
n_drop=int((~cov.drop_duplicates('person_id')['qc_pass']).sum())
print("genomes failing QC (dropped; round20 persons):",round20(n_drop),
      "| related retained+flagged (round20 persons):",round20(int(cov.drop_duplicates('person_id')['related'].sum())))

# ============================================================================
# (9) Assemble per-cohort analytic parquet (one row/person/cohort; keep qc_pass).
# ============================================================================
ANALYTIC_COLS=(['person_id','cohort','index_date','ancestry']+pc_cols+
 ['pheno_primary_composite','pheno_opioid_only','pheno_pain_only','pheno_allcause',
  'age','age2','sex','bmi','surgical_year','ehr_density','procedure_class',
  'cx_depression','cx_anxiety','cx_sud_oud','cx_pain_region','cx_pain_chronic','comorbidity_count',
  'qc_pass','related'])
analytic={}
for c in order:
    d=cov[cov.cohort==c].drop_duplicates('person_id').copy()
    d=d[d['qc_pass']]
    d=d[[col for col in ANALYTIC_COLS if col in d.columns]]
    assert d['person_id'].is_unique, f"{c}: not one row/person"
    d.to_parquet(f'/tmp/analytic_{c}.parquet',index=False)
    analytic[c]=d
    print(f"  analytic_{c}: rows(round20)={round20(len(d))} one-row/person={d['person_id'].is_unique}")

# ============================================================================
# (9b) Pooled analytic parquet: concat 9, dedup to EARLIEST index/person (index_date_pooled rule).
# ============================================================================
pooled=(pd.concat([analytic[c] for c in order if c in analytic],ignore_index=True)
        .sort_values('index_date').drop_duplicates('person_id',keep='first'))
for cls in sorted(set(PROC_CLASS.values())):
    pooled['pc_'+cls]=(pooled['procedure_class']==cls).astype(int)
pooled.to_parquet('/tmp/analytic_pooled.parquet',index=False)
print("POOLED analytic: persons(round20)=",round20(len(pooled)),"| one-row/person=",pooled['person_id'].is_unique,
      "| procedure_class one-hot:",[c for c in pooled.columns if c.startswith('pc_')])

# missingness (aggregate)
print("\nmissingness (round20 person-cohort rows): age",round20(int(cov['age'].isna().sum())),
      "| bmi",round20(int(cov['bmi'].isna().sum())),
      "| PC-any-missing",round20(int((~pc_ok).sum())),
      "| sex Other/Unknown",round20(int((cov['sex']=='Other/Unknown').sum())))
print("wrote /tmp/analytic_<cohort>.parquet (9) + /tmp/analytic_pooled.parquet")
print("DONE Part B")
