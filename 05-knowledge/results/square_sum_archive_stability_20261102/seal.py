#!/usr/bin/env python3
"""Seal/check packet bytes and cited source snapshots; not a proof checker."""
import hashlib,json,re,subprocess,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MANIFEST=HERE/'MANIFEST.json'
BASE='66cf347030dc'
INCOMING='d7bbbf3f05'
FINAL_INCOMING='c190e11b77'
DEPENDENCIES=[
 '05-knowledge/results/session137_overnight_proof_frontiers_20261102/discriminant/GLOBAL_PROOF.md',
 '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/TRIANGLE_VARIANCE.md',
 '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/MULTISTAR_STABILITY.md',
 '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/H1_STABILITY.md',
 '05-knowledge/results/session138_stability_and_actual_contexts_20261102/stability/SHARP_NINE_GRAPH_CERTIFICATE.md',
]
def need(b,s):
    if not b:raise RuntimeError(s)
def sha(b):return hashlib.sha256(b).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT,stderr=subprocess.DEVNULL)
def paths():
    out=set(DEPENDENCIES)
    for f in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',f.read_text()):
            if '://' in target or target.startswith('#'):continue
            p=(f.parent/target.split('#')[0]).resolve()
            if p==HERE or HERE in p.parents:continue
            try:out.add(str(p.relative_to(ROOT)))
            except ValueError:raise RuntimeError('nonportable citation '+str(p))
    return sorted(out)
def make():
    pins=[]
    for p in paths():
        for revision in [BASE,INCOMING,FINAL_INCOMING]:
            try:raw=git('show',revision+':'+p);break
            except subprocess.CalledProcessError:continue
        else:raise RuntimeError('missing source '+p)
        pins.append({'path':p,'revision':git('rev-parse',revision).decode().strip(),'blob':git('rev-parse',revision+':'+p).decode().strip(),'sha256':sha(raw),'role':'stability proof dependency' if p in DEPENDENCIES else 'recovery or incoming reference'})
    local={f.name:sha(f.read_bytes()) for f in sorted(HERE.iterdir()) if f.is_file() and f.name!='MANIFEST.json' and f.suffix in {'.md','.py','.cpp','.json'}}
    return {'scope':'Integrity of this packet and explicit cited sources. Does not certify mathematics, all historical files, or transitive dependencies.','files':local,'sources':pins}
def main():
    if sys.argv[1:]==['--write']:
        MANIFEST.write_text(json.dumps(make(),indent=2,sort_keys=True)+'\n')
    else:
        need(not sys.argv[1:] or sys.argv[1:]==['--check'],'usage: seal.py [--write|--check]')
    r=json.loads(MANIFEST.read_text())
    for f,h in r['files'].items():need(sha((HERE/f).read_bytes())==h,'packet hash '+f)
    for row in r['sources']:
        raw=git('show',row['revision']+':'+row['path'])
        need(sha(raw)==row['sha256'],'source hash '+row['path'])
        need(git('rev-parse',row['revision']+':'+row['path']).decode().strip()==row['blob'],'source object '+row['path'])
    print(json.dumps({'status':'passed','packet_files':len(r['files']),'source_snapshots':len(r['sources']),'mathematical_proof_audit':False},sort_keys=True))
if __name__=='__main__':main()
