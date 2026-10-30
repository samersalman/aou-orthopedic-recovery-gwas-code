# s6_rescore.py -- regenerate per-person PRS from genotypes ALREADY IN THE WORKSPACE BUCKET.
# No new CDR/genomic pull, no Hail, no Dataproc: reuses gwas/geno/step2/*.bgen extracted in Session 5
# and the PGS Catalog weight files in reference/pgs/. Row-level scores stay in-VM (/tmp) and are
# persisted to the workspace bucket only (never egressed). Prints aggregate diagnostics only.
import os, sys, gzip, subprocess
import pandas as pd
import numpy as np

WSB = os.environ.get('WORKSPACE_BUCKET')
if not WSB:
    raise RuntimeError('WORKSPACE_BUCKET must be set to the authorized workspace bucket URI')
PGS = {'kneeOA': 'PGS002767', 'hipOA': 'PGS002763', 'backpain': 'PGS004926', 'mdd': 'PGS003585'}
PGDIR = '/tmp/pgs'
GDIR = '/tmp/geno'
for d in (PGDIR, GDIR, '/tmp/score'):
    os.makedirs(d, exist_ok=True)
TAB = chr(9)
NL = chr(10)


def round20(n):
    n = int(n)
    if n <= 0:
        return 0
    if n <= 20:
        return '<=20'
    return int(round(n / 20.0) * 20)


def sh(c):
    return subprocess.run(c, shell=True, capture_output=True, text=True)


def free_gb():
    s = os.statvfs('/tmp')
    return s.f_bavail * s.f_frsize / 1e9


for tag, pid in PGS.items():
    f = PGDIR + '/' + pid + '_hmPOS_GRCh38.txt.gz'
    if not os.path.exists(f):
        sh('gsutil -u $GOOGLE_PROJECT cp ' + WSB + '/reference/pgs/' + pid + '_hmPOS_GRCh38.txt.gz ' + PGDIR + '/')
print('pgs weight files present: ' + str(sum(os.path.exists(PGDIR + '/' + p + '_hmPOS_GRCh38.txt.gz')
                                             for p in PGS.values())) + ' of 4', flush=True)


def build_score(tag, pid):
    out = '/tmp/score/' + tag + '.score'
    if os.path.exists(out):
        return out
    ix = None
    n = 0
    cch = cps = None
    with gzip.open(PGDIR + '/' + pid + '_hmPOS_GRCh38.txt.gz', 'rt') as f, open(out, 'w') as o:
        o.write('ID' + TAB + 'EA' + TAB + 'W' + NL)
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip(NL).split(TAB)
            if ix is None:
                ix = {c: i for i, c in enumerate(p)}
                cch = 'hm_chr' if 'hm_chr' in ix else 'chr_name'
                cps = 'hm_pos' if 'hm_pos' in ix else 'chr_position'
                continue
            ch = p[ix[cch]]
            ps = p[ix[cps]]
            ea = p[ix['effect_allele']]
            w = p[ix['effect_weight']]
            if ch and ps and ea and w and ea not in ('.', 'NA'):
                o.write(ch + ':' + ps + TAB + ea + TAB + w + NL)
                n += 1
    print('  built ' + tag + '.score (' + pid + '): ' + str(n) + ' variants', flush=True)
    return out


for tag, pid in PGS.items():
    build_score(tag, pid)


def do_anc(anc):
    for ext in ('bgen', 'sample'):
        f = GDIR + '/step2_' + anc + '.' + ext
        if not os.path.exists(f):
            sh('gsutil -u $GOOGLE_PROJECT -m cp ' + WSB + '/gwas/geno/step2/step2_' + anc + '.' + ext + ' ' + GDIR + '/')
    pf = '/tmp/' + anc
    if not os.path.exists(pf + '.pgen'):
        print('[' + anc + '] bgen->pgen  free_gb=' + str(round(free_gb())), flush=True)
        r = subprocess.run(['plink2', '--bgen', GDIR + '/step2_' + anc + '.bgen', 'ref-first',
                            '--sample', GDIR + '/step2_' + anc + '.sample',
                            '--set-all-var-ids', '@:#', '--rm-dup', 'exclude-all', '--make-pgen',
                            '--out', pf, '--threads', '8', '--memory', '45000'],
                           capture_output=True, text=True)
        tail = r.stderr.strip().split(NL)[-1] if r.stderr.strip() else ''
        print('  pgen rc ' + str(r.returncode) + '  ' + tail[:160], flush=True)
        if r.returncode == 0:
            for ext in ('bgen', 'sample'):
                try:
                    os.remove(GDIR + '/step2_' + anc + '.' + ext)
                except OSError:
                    pass
    if not os.path.exists(pf + '.pgen'):
        print('[' + anc + '] no pgen -- skipping', flush=True)
        return None
    res = None
    for tag in PGS:
        o = '/tmp/prs_' + anc + '_' + tag
        r = subprocess.run(['plink2', '--pfile', pf, '--score', '/tmp/score/' + tag + '.score',
                            '1', '2', '3', 'header', 'cols=+scoresums,+denom',
                            '--out', o, '--threads', '8', '--memory', '45000'],
                           capture_output=True, text=True)
        sf = o + '.sscore'
        if not os.path.exists(sf):
            print('  ' + anc + '/' + tag + ': FAILED rc' + str(r.returncode) + ': '
                  + r.stderr.strip().split(NL)[-1][:140], flush=True)
            continue
        s = pd.read_csv(sf, sep=TAB)
        idc = [c for c in s.columns if c.lstrip('#') == 'IID'][0]
        avgc = [c for c in s.columns if c.endswith('_AVG')][0]
        dnc = [c for c in s.columns if c.upper().endswith('DENOM')]
        s = s.rename(columns={idc: 'IID', avgc: 'prs_' + tag})
        dn = int(s[dnc[0]].median()) if dnc else -1
        print('  ' + anc + '/' + tag + ': scored_persons=' + str(round20(s.shape[0]))
              + ' median_matched_variants=' + str(dn), flush=True)
        s2 = s[['IID', 'prs_' + tag]]
        res = s2 if res is None else res.merge(s2, on='IID', how='outer')
        try:
            os.remove(sf)
        except OSError:
            pass
    for ext in ('pgen', 'pvar', 'psam'):
        try:
            os.remove(pf + '.' + ext)
        except OSError:
            pass
    if res is not None:
        res['ancestry'] = anc
        print('[' + anc + '] persons=' + str(round20(res.shape[0])), flush=True)
    return res


parts = []
for anc in ('amr', 'afr', 'eur'):
    try:
        r = do_anc(anc)
        if r is not None:
            parts.append(r)
    except Exception as e:
        print('[' + anc + '] EXC ' + type(e).__name__ + ' ' + str(e)[:140], flush=True)

if parts:
    allp = pd.concat(parts, ignore_index=True)
    allp = allp.rename(columns={'IID': 'person_id'})
    allp['person_id'] = pd.to_numeric(allp['person_id'], errors='coerce')
    allp = allp.dropna(subset=['person_id'])
    allp['person_id'] = allp['person_id'].astype('int64')
    allp.to_parquet('/tmp/prs_scores.parquet', index=False)
    sh('gsutil -u $GOOGLE_PROJECT cp /tmp/prs_scores.parquet ' + WSB + '/gwas/prs/prs_scores.parquet')
    print('PRS TOTAL persons=' + str(round20(allp.shape[0])) + ' cols='
          + str([c for c in allp.columns if c != 'person_id']), flush=True)
    print('persisted -> <W>/gwas/prs/prs_scores.parquet', flush=True)
else:
    print('NO PRS PARTS PRODUCED', flush=True)
print('RESCORE DONE', flush=True)
