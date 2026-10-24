#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent
m=json.loads((p/'MANIFEST.json').read_text())
expected={r['path']:r['sha256'] for r in m['files']}
actual={str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in p.rglob('*') if f.is_file() and f.name!='MANIFEST.json'}
if expected!=actual:raise RuntimeError('packet file closure or hash mismatch')
print(json.dumps({'status':'PASS','files_verified':len(actual)},sort_keys=True))
