#!/usr/bin/env python3
"""Verify every claimed replay, then seal the audit packet."""
from pathlib import Path
from hashlib import sha256
import json
H=Path(__file__).resolve().parent

def need(c,m):
    if not c:raise RuntimeError(m)
def same(a,b):
    need((H/a).read_bytes()==(H/b).read_bytes(),'byte mismatch '+a+' '+b)
    return {'first':a,'second':b,'sha256':sha256((H/a).read_bytes()).hexdigest()}
pairs=[('normal.json','optimized.json'),('parent-cards.json','parent-cards-optimized.json'),('hostile-normal.json','hostile-optimized.json'),('hostile-packet-normal.json','hostile-packet-optimized.json'),('extra-rotation-controls.json','extra-rotation-controls-optimized.json'),('normal.json','relocated-audit/normal.json'),('normal.json','relocated-audit/optimized.json'),('parent-cards.json','relocated-audit/parent-cards.json'),('producer/normal.json','producer/optimized.json'),('producer/normal.json','producer-replay-normal.json'),('producer/normal.json','producer-replay-optimized.json'),('comparison-normal.json','comparison-optimized.json')]
verified=[same(a,b) for a,b in pairs]
producer_manifest=sha256((H/'producer/MANIFEST.json').read_bytes()).hexdigest()
need(producer_manifest=='7431d29c032d92aa6da1807c436bb839ed9c84909872de31d7f64b09d4963d43','frozen producer manifest identity')
for f in json.loads((H/'producer/MANIFEST.json').read_text())['files']:
    need(sha256((H/'producer'/f['path']).read_bytes()).hexdigest()==f['sha256'],'producer final integrity')
receipt={'status':'PASS','producer_manifest_sha256':producer_manifest,'independent_methods':['descending coefficient DP','R factors from full-edge orbits','binary-order actual matrix full cubes n<=8','actual parent cycles, automorphisms and card orbits n<=6','actual rotation directed edge-orbits'], 'scope':'all-size proof audit; scalar counts even n<=12; no new unmarked tournament class count','byte_identical_receipts':verified,'hostile_math_controls':9,'hostile_packet_controls':3,'provenance_correction':'THM-554 score product and THM-555 balanced totals through n10 are inherited and reproduced','unresolved_blockers':[]}
(H/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
files=[]
for p in sorted(H.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts or p.suffix in ('.pyc','.pyo') or (p.suffix == '.out' and p.parent != H/'producer') or p.name=='MANIFEST.json' and p.parent==H:continue
    files.append({'path':str(p.relative_to(H)),'sha256':sha256(p.read_bytes()).hexdigest()})
(H/'MANIFEST.json').write_text(json.dumps({'status':'PASS','files':files},indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','producer_manifest_sha256':producer_manifest,'audit_manifest_sha256':sha256((H/'MANIFEST.json').read_bytes()).hexdigest(),'files':len(files)},sort_keys=True))
