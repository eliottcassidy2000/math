#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
expected={r['path'] for r in manifest['files']}
if actual!=expected:raise RuntimeError(f'manifest file-set mismatch: {actual^expected}')
for r in manifest['files']:
    if hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()!=r['sha256']:raise RuntimeError('SHA256 mismatch: '+r['path'])
print(json.dumps({'status':'PASS','files':len(expected)},sort_keys=True))
