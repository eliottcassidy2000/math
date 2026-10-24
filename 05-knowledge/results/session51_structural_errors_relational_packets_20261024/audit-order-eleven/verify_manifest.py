#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
manifest=json.loads((HERE/'MANIFEST.json').read_text())
expected={e['path'] for e in manifest['files']}
actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.relative_to(HERE).as_posix() not in {'MANIFEST.json','replay-receipt.json'}}
if actual!=expected:raise ValueError('exact audit membership mismatch')
for entry in manifest['files']:
 name=entry['path'];path=HERE/name
 if name.startswith('/') or '..' in Path(name).parts:raise ValueError('unsafe manifest path')
 if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('manifest byte mismatch: '+name)
print(json.dumps({'status':'PASS','files':len(manifest['files'])},sort_keys=True))
