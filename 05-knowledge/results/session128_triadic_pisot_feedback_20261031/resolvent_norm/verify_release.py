#!/usr/bin/env python3
"""Verify file hashes and replay both modes after a local relocation."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

root=Path(__file__).resolve().parent
checks=0
for line in (root/'SHA256SUMS').read_text().splitlines():
    expected,name=line.split('  ',1)
    actual=hashlib.sha256((root/name).read_bytes()).hexdigest()
    if actual != expected: raise RuntimeError('hash mismatch: '+name)
    checks+=1
normal=(root/'normal.json').read_bytes()
optimized=(root/'optimized.json').read_bytes()
if normal != optimized: raise RuntimeError('saved mode mismatch')
with tempfile.TemporaryDirectory(prefix='.relocated-',dir=root) as dest:
    target=Path(dest)/'check_resolvent.py'
    shutil.copyfile(root/'check_resolvent.py',target)
    outputs=[]
    for flags in (['-B'],['-B','-O']):
        out=subprocess.check_output([sys.executable,*flags,str(target)],cwd=dest)
        if out != normal: raise RuntimeError('relocated replay mismatch')
        outputs.append(hashlib.sha256(out).hexdigest())
print(json.dumps({'status':'PASS','hashed_files':checks,'checks':json.loads(normal)['checks'],
                  'replay_modes':['normal','optimized'],'relocated':True,
                  'output_sha256':outputs},indent=2,sort_keys=True))
