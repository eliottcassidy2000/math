#!/usr/bin/env python3
"""Verify exact publication whitelist, pinned bytes, and paired audit receipts."""
import hashlib,json
from pathlib import Path
base=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise RuntimeError(message)

manifest=json.loads((base/'MANIFEST.json').read_text())
actual={str(p.relative_to(base)) for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
require(actual==set(manifest),'publication artifact whitelist mismatch')
for path,pin in manifest.items():
    b=(base/path).read_bytes()
    require(len(b)==pin['bytes'],'byte-count mismatch: '+path)
    require(hashlib.sha256(b).hexdigest()==pin['sha256'],'checksum mismatch: '+path)
for stem in ['magic_fence','theta','ramsey']:
    a=(base/(stem+'.normal.json')).read_bytes();b=(base/(stem+'.optimized.json')).read_bytes()
    require(a==b,'normal/optimized mismatch: '+stem)
    data=json.loads(a);require(data['status']=='PASS','failed receipt: '+stem)
    if stem in ['magic_fence','theta']:
        source=base/('check_'+stem+'.py')
        require(data['checker_sha256']==hashlib.sha256(source.read_bytes()).hexdigest(),'checker identity mismatch: '+stem)
magic=json.loads((base/'magic_fence.normal.json').read_text())
require(magic['checked_conditions']==7160,'magic/fence gate count')
require(magic['magic']['full_permutations']==362880,'full target census')
require(magic['magic']['rational_parameter_pairs']==2025,'rational parameter count')
theta=json.loads((base/'theta.normal.json').read_text())
require(theta['checked_conditions']==15,'theta gate count')
ramsey=json.loads((base/'ramsey.normal.json').read_text())
require(ramsey['closed_component_cases']==496,'closed-component amendment controls')
require(ramsey['paley11']['H']==95095 and ramsey['paley11']['direct_cycles_mod_rotation']==5505,'Paley controls')
require(ramsey['support_graphs']['ordered_pairs']==2208196,'linear-forest pair census')
receipt=json.loads((base/'AUDIT_RECEIPTS.json').read_text())
require(receipt['status']=='PASS','discovery replay status')
require(receipt['ramsey_separate_independent']['publication_checked_conditions']==112989,'final corrected checker gates')
require(receipt['ramsey_separate_independent']['added_closed_component_equalities']==76,'closed-component addition gates')
require(receipt['ramsey_separate_independent']['added_extra_circle_hostiles']==57,'extra-circle hostile gates')
print(json.dumps({'status':'PASS','sealed_files':len(actual)+1,'normal_optimized_equal':True,'scope_corrections_verified':True},sort_keys=True))
