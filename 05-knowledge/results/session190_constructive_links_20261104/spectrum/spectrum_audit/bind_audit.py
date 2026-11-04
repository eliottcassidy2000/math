#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
import json
from independent_audit import profile, need, relabel, B
D=Path(__file__).resolve().parent
P=D.parent/'spectrum'
if not P.exists():P=D.parent/'session190-spectrum'
PROOF='d1e5debdf607c3eb83f8054c08f94ec1137ef0b4b17aac775973f878429d5479'

def digest(p):return sha256(p.read_bytes()).hexdigest()
def main():
    release=json.loads((P/'RELEASE_FILES.json').read_text())
    for f in release['files']:need(digest(P/f['path'])==f['sha256'],'producer release pin '+f['path'])
    need(digest(P/'PROOF.md')==PROOF,'reviewed proof version')
    own=json.loads((D/'local_audit_results.json').read_text());producer=json.loads((P/'receipt.json').read_text())
    need(relabel(B,(2,0,3,1,5,4,6))==tuple(producer['patch'][-1]['rows']),'original/final seed transport')
    for ours,theirs in zip(own['selected']['records'],producer['patch']):
        need(ours['rows']==theirs['rows'],'physical seed')
        for key in ['s','e','F','q']:need(ours[key]==theirs['profile'][key],'local coordinate '+key)
        need(ours['H']==theirs['profile']['h'],'local Hamilton count')
    larger=[]
    for twin in producer['larger_twins']:
        records=[]
        for key,delta in [('A',0),('B',2)]:
            p=profile(twin[key],twin['root']);need(p['strong'],'larger twin strong')
            need((p['H'],p['s'],p['F'])==(twin['H']+delta,1,twin['F']),'larger twin exact profile')
            records.append(dict(side=key,rows=twin[key],root=twin['root'],profile=p))
        need(twin['necessary_parent_H']==(twin['H']-2)//3 and 3*twin['necessary_parent_F']==twin['F'],'inverse coefficients')
        larger.append(records)
    for base in ['audit','wrapper_audit','census_audit']:
        need((D/(base+'.stdout')).read_bytes()==(D/(base+'.optimized.stdout')).read_bytes(),'normal optimized '+base)
    need((D/'census_audit_results.json').read_bytes()==(D/'census_audit_normal.json').read_bytes(),'normal optimized census data')
    inputs={name:digest(P/name) for name in release['publication_files']}
    artifacts=['REVIEW.md','README.md','independent_audit.py','wrapper_audit.py','census_audit.py','independent_census.cpp','bind_audit.py','local_audit_results.json','wrapper_audit_results.json','census_audit_results.json','census_audit_normal.json','inherited_cover8.json','independent_low9_bank.json','audit.stdout','audit.optimized.stdout','wrapper_audit.stdout','wrapper_audit.optimized.stdout','census_audit.stdout','census_audit.optimized.stdout','small.out','bank8.out','extend8.out','extend9.out']
    out=dict(status='PASS',approved_scope='Explicit all-clock patch, stated sharp finite local bounds, and two fixed-parent exclusions; no global spectrum closure.',reviewed_proof_sha256=PROOF,producer_pins=inputs,audit_pins={name:digest(D/name) for name in artifacts},independent_larger_twin_profiles=larger)
    (D/'audit_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS final producer hash binding, every local coordinate, independent larger twins, and normal/optimized audit agreement.')
if __name__=='__main__':main()
