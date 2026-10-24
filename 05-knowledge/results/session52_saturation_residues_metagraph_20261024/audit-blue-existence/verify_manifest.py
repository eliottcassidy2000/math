#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
m=json.loads((root/'MANIFEST.json').read_text())
paths=[x['path'] for x in m['files']]
if len(paths)!=len(set(paths)): raise RuntimeError('duplicate manifest path')
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}-{'MANIFEST.json'}
if set(paths)!=actual: raise RuntimeError('manifest not closed')
for record in m['files']:
    p=Path(record['path'])
    if p.is_absolute() or '..' in p.parts: raise RuntimeError('unsafe path')
    if hashlib.sha256((root/p).read_bytes()).hexdigest()!=record['sha256']: raise RuntimeError('hash mismatch: '+str(p))
print(json.dumps({'status':'PASS','files_verified':len(paths)},sort_keys=True))
