#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
m=json.loads((HERE/'MANIFEST.json').read_text())
expected={r['path'] for r in m['files']}
actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.relative_to(HERE).as_posix() not in {'MANIFEST.json','replay-receipt.json'}}
if actual!=expected:raise ValueError('exact packet membership mismatch')
for r in m['files']:
 p=HERE/r['path']
 if p.stat().st_size!=r['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise ValueError('manifest byte mismatch: '+r['path'])
print(json.dumps({'status':'PASS','files':len(expected)},sort_keys=True))
