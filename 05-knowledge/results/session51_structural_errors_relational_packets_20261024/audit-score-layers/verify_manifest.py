#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent
m=json.loads((p/'MANIFEST.json').read_text())['files']
a={str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in p.rglob('*') if f.is_file() and f.name!='MANIFEST.json' and '__pycache__' not in f.parts}
if a!=m:raise RuntimeError('frozen audit manifest file-set/hash mismatch')
print(json.dumps({'status':'PASS','files':len(m)},sort_keys=True))
