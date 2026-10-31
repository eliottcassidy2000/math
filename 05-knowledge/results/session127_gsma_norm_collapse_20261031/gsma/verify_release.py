#!/usr/bin/env python3
import hashlib,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
checks=0
for line in (root/'SHA256SUMS').read_text().splitlines():
    want,name=line.split('  ',1)
    got=hashlib.sha256((root/name).read_bytes()).hexdigest()
    if got!=want:raise RuntimeError('hash mismatch: '+name)
    checks+=1
expected=(root/'unit_kernel.normal.json').read_bytes()
if expected!=(root/'unit_kernel.optimized.json').read_bytes():raise RuntimeError('normal/O mismatch')
for flags in [[],['-O']]:
    actual=subprocess.check_output([sys.executable,'-B',*flags,str(root/'check_unit_kernel.py')],cwd=root)
    if actual!=expected:raise RuntimeError('executable replay mismatch')
print(json.dumps({'status':'PASS','files_hashed':checks,'normal_and_optimized_replay':True,'exact_checks_each':json.loads(expected)['checks']},sort_keys=True))
