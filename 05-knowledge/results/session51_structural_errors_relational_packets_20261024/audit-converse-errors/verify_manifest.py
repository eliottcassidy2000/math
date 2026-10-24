#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
expected=set(manifest['files'])
if actual!=expected:raise RuntimeError('audit manifest file set mismatch: '+str(actual^expected))
for p,sha in manifest['files'].items():
    if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=sha:raise RuntimeError('audit SHA256 mismatch '+p)
print(json.dumps({'status':'PASS','files':len(expected)},sort_keys=True))
